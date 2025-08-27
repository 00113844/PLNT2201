# app.R -- corrected Shiny app for PRO4SAIL (fixes LAI scoping bug)
# Install prosail if needed:
# remotes::install_gitlab("jbferet/prosail")

library(shiny)
library(prosail)
library(ggplot2)

# Silence R CMD check notes for NSE/global data objects when packaging/sharing
if (getRversion() >= "2.15.1") utils::globalVariables(c(
  "SpecSOIL", "SpecPROSPECT_FullRange", "refl", "LAI", "trans"))

ui <- fluidPage(
  titlePanel("Retrieve canopy extinction k from PRO4SAIL (population → LAI → absorptance → k)"),
  sidebarLayout(
    sidebarPanel(
      h4("Population → LAI"),
      numericInput("plants_per_m2", "Plants per m² (population)", value = 200, min = 1),
      numericInput("leaf_area_per_plant", "Leaf area per plant (m² plant⁻¹)", value = 0.02, min = 0.0001, step = 0.005),
      hr(),
      h4("Leaf angle distribution (LID)"),
      selectInput("TypeLidf", "LIDF type", choices = list("Ellipsoidal (Type 2)" = 2, "Verhoef-style (Type 1)" = 1), selected = 2),
      numericInput("LIDFa", "LIDFa (mean leaf angle or param)", value = 30, min = 0, max = 90),
      conditionalPanel(
        condition = "input.TypeLidf == '1'",
        numericInput("LIDFb", "LIDFb (bimodality param, use with Type 1)", value = 0, min = -1, max = 1)
      ),
      hr(),
      h4("PROSPECT (leaf) parameters"),
      numericInput("CHL", "Chlorophyll (µg/cm²)", value = 40, min = 0),
      numericInput("CAR", "Carotenoids (µg/cm²)", value = 8, min = 0),
      numericInput("EWT", "Equivalent water thickness (cm)", value = 0.01, min = 0.0001, step = 0.005),
      numericInput("LMA", "LMA (g cm⁻²)", value = 0.009, min = 0.0001, step = 0.001),
      numericInput("N", "Leaf structure N (PROSPECT)", value = 1.5, min = 1, step = 0.1),
      hr(),
      h4("Geometry & soil"),
      sliderInput("sza", "Solar zenith angle (deg)", min = 0, max = 80, value = 30),
      selectInput("soil", "Soil reflectance", choices = c("Dry" = "Dry_Soil", "Wet" = "Wet_Soil")),
      numericInput("q", "Hot-spot parameter q", value = 0.1, min = 0, step = 0.01),
      actionButton("run", "Run PRO4SAIL")
    ),
    mainPanel(
      # Results moved to bottom; start with solar helper
      hr(),
      h4("Solar position helper (compute SZA)"),
      fluidRow(
        column(4, dateInput("date_local", "Date", value = Sys.Date())),
        column(2, sliderInput("hour_local", "Hour", min = 0, max = 23, value = 12)),
        column(2, sliderInput("minute_local", "Minute", min = 0, max = 59, value = 0)),
        column(4, numericInput("lat", "Latitude (deg, N positive; Perth -31.95)", value = -31.95, step = 0.01))
      ),
      fluidRow(
        column(4, numericInput("lon", "Longitude (deg, E positive; Perth 115.86)", value = 115.86, step = 0.01)),
        column(4, numericInput("tz_offset", "UTC offset hours", value = 8, step = 1)),
        column(4, actionButton("calc_sza", "Compute solar zenith"))
      ),
      helpText("Uses NOAA solar position formula (eqn of time + declination) to derive solar zenith; updates SZA slider used in PRO4SAIL."),
      verbatimTextOutput("solar_pos"),
      hr(),
      h4("Field measurement mode (optional)"),
      fluidRow(
        column(3, numericInput("Io_field", HTML("I<sub>0</sub> (incoming above canopy)"), value = NA, min = 0)),
        column(3, numericInput("Ir_field", HTML("I<sub>r</sub> (reflected above canopy)"), value = NA, min = 0)),
        column(3, numericInput("Ic_field", HTML("I<sub>c</sub> (transmitted below canopy)"), value = NA, min = 0)),
        column(3, numericInput("LAI_field", "LAI (measured)", value = NA, min = 0))
      ),
      verbatimTextOutput("field_calc"),
      fluidRow(
        column(6, plotOutput("specPlot", height = "320px")),
        column(6, plotOutput("kLAIPlot", height = "320px"))
      ),
      hr(),
  h4("Results"),
  verbatimTextOutput("summary"),
  hr(),
  h5("Notes"),
  p("Beer–Lambert canopy form: F = 1 - exp(-k * LAI); here F ~ fcover (fractional cover) or fAPAR (fraction of absorbed PAR)."),
  p("Transmittance T = exp(-k*LAI). We solve k = -ln(T)/LAI with T chosen as (1 - fcover) or (1 - fAPAR)."),
  p("Absorptance spectra: A(λ) = 1 - R(λ) - T(λ) when both reflectance and transmittance are available; no assumption that A = 1 - R."),
  p("k can exceed 1 at high solar zenith angles or strongly planophile canopies; lower values occur in erectophile canopies."),
  p(HTML("Field method: Transmission I<sub>c</sub>; Intercepted fraction F = I<sub>im</sub> = ((I<sub>0</sub> - I<sub>r</sub>) - I<sub>c</sub>)/I<sub>0</sub>. Beer–Lambert with reflection accounted: I<sub>c</sub> = (I<sub>0</sub>-I<sub>r</sub>) exp(-k LAI). Thus k = -ln(I<sub>c</sub>/(I<sub>0</sub>-I<sub>r</sub>))/LAI."))
    )
  )
)

server <- function(input, output, session) {

  # Compute solar position (solar zenith) using NOAA algorithm when button pressed
  observeEvent(input$calc_sza, {
    req(input$date_local, input$hour_local, input$minute_local, input$lat, input$lon, input$tz_offset)
    # Convert to day of year
    d <- as.Date(input$date_local)
    doy <- as.integer(strftime(d, format = "%j"))
    # Fractional hour local
    local_time <- input$hour_local + input$minute_local/60
    # Fractional year gamma (radians) (NOAA approximate)
    gamma <- 2*pi/365 * (doy - 1 + (local_time - 12)/24)
    # Equation of time (minutes)
    eqtime <- 229.18 * (0.000075 + 0.001868*cos(gamma) - 0.032077*sin(gamma) - 0.014615*cos(2*gamma) - 0.040849*sin(2*gamma))
    # Solar declination (radians)
    decl <- 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - 0.006758*cos(2*gamma) + 0.000907*sin(2*gamma) - 0.002697*cos(3*gamma) + 0.00148*sin(3*gamma)
    # Time offset (minutes)
    time_offset <- eqtime + 4*input$lon - 60*input$tz_offset
    # True solar time (minutes)
    tst <- local_time * 60 + time_offset
    # Hour angle (deg)
    ha <- (tst / 4) - 180
    # Convert degrees to radians for lat and hour angle
    lat_rad <- input$lat * pi/180
    ha_rad <- ha * pi/180
    # Solar zenith angle (radians)
    cos_zen <- sin(lat_rad)*sin(decl) + cos(lat_rad)*cos(decl)*cos(ha_rad)
    cos_zen <- max(-1,min(1,cos_zen))
    zen_rad <- acos(cos_zen)
    zen_deg <- zen_rad * 180/pi
    # Solar elevation
    elev_deg <- 90 - zen_deg
    # Update slider
    updateSliderInput(session, "sza", value = round(zen_deg,2))
    output$solar_pos <- renderPrint({
      cat("Solar position (NOAA algorithm):\n")
      cat(sprintf("  Day of year = %d\n", doy))
      cat(sprintf("  Equation of time (min) = %.2f\n", eqtime))
      cat(sprintf("  Declination (deg) = %.2f\n", decl*180/pi))
      cat(sprintf("  Hour angle (deg) = %.2f\n", ha))
      cat(sprintf("  Solar zenith (deg) = %.2f\n", zen_deg))
      cat(sprintf("  Solar elevation (deg) = %.2f\n", elev_deg))
      cat("(Slider 'Solar zenith angle' updated.)\n")
    })
  })
  
  # reactive LAI - explicit name to avoid confusion with other 'LAI' symbols
  LAI_val_reactive <- reactive({
    req(input$plants_per_m2, input$leaf_area_per_plant)
    lai_val <- input$plants_per_m2 * input$leaf_area_per_plant
    return(lai_val)
  })
  
  # PROSPECT input builder
  prospect_input <- reactive({
    data.frame(Cab = input$CHL,
               Car = input$CAR,
               Ant = 0.0,
               Cw  = input$EWT,
               Cm  = input$LMA,
               N   = input$N)
  })
  
  # choose soil spectrum from package data safely
  soil_spectrum <- reactive({
    # safe check: prosail package provides SpecSOIL, SpecPROSPECT etc.
    if (exists("SpecSOIL", where = asNamespace("prosail"), inherits = FALSE)) {
      if (input$soil == "Dry_Soil") {
        return(SpecSOIL$Dry_Soil)
      } else {
        return(SpecSOIL$Wet_Soil)
      }
    } else {
      # fallback: return NULL and PRO4SAIL will use its default soil
      return(NULL)
    }
  })
  
  # Run PRO4SAIL on button press; return a well-formed list always
  prosail_result <- eventReactive(input$run, {
    # ensure prosail is loaded
    if (!"prosail" %in% loadedNamespaces()) {
      showNotification("prosail package not loaded. Install and restart app.", type = "error", duration = 6)
      return(NULL)
    }
    
    LAI_val <- LAI_val_reactive()      # explicit variable name
    Input_PROSPECT <- prospect_input()
    soil_spec <- soil_spectrum()
    
    # Call PRO4SAIL: use explicit argument names as in prosail vignette
    # Note: the exact argument names depend on prosail version; these names follow the vignette.
    Ref <- tryCatch({
      PRO4SAIL(
        Input_PROSPECT = Input_PROSPECT,
        TypeLidf = as.integer(input$TypeLidf),
        LIDFa = as.numeric(input$LIDFa),
        LIDFb = ifelse(is.null(input$LIDFb), 0, as.numeric(input$LIDFb)),
        lai = as.numeric(LAI_val),
        q = as.numeric(input$q),
        tts = as.numeric(input$sza),
        tto = 0,
        psi = 0,
        rsoil = soil_spec
      )
    }, error = function(e) {
      showNotification(paste("PRO4SAIL call failed:", e$message), type = "error", duration = 8)
      return(NULL)
    })
    
    if (is.null(Ref)) return(NULL)
    
    # gather wavelengths from package if available (safe fallback)
    if (exists("SpecPROSPECT_FullRange", where = asNamespace("prosail"), inherits = FALSE)) {
      wvl <- SpecPROSPECT_FullRange$lambda
    } else if (!is.null(Ref$wavelengths)) {
      wvl <- Ref$wavelengths
    } else {
      # fallback vector (not ideal but avoid crash)
      wvl <- seq(400, 2500, by = 1)
    }
    par_mask <- (wvl >= 400 & wvl <= 700)
    
    # compute I_r if rsdt present
    I_r <- NA_real_
    if (!is.null(Ref$rsdt)) {
      rsdt_vec <- Ref$rsdt
      if (length(rsdt_vec) == length(wvl)) {
        I_r <- mean(rsdt_vec[par_mask], na.rm = TRUE)
      } else {
        I_r <- mean(rsdt_vec, na.rm = TRUE)
      }
    }
    
    # absorbed direct / hemispherical if available
    I_abs_dir <- ifelse(!is.null(Ref$abs_dir), as.numeric(Ref$abs_dir), NA_real_)
    # Hemispherical (all-direction) absorbed fraction (preferred for Beer-Lambert F)
    Abs_hem <- NA_real_
    if (!is.null(Ref$abs_hem)) {
      # Some versions may return a scalar; if vector, average over PAR
      if (length(Ref$abs_hem) == 1) {
        Abs_hem <- as.numeric(Ref$abs_hem)
      } else if (exists("par_mask")) {
        # if spectral length matches wavelengths mask
        if (length(Ref$abs_hem) == length(wvl)) {
          Abs_hem <- mean(Ref$abs_hem[par_mask], na.rm = TRUE)
        } else {
          Abs_hem <- mean(Ref$abs_hem, na.rm = TRUE)
        }
      }
    }
    
    # Compute fAPAR if function available (prosail >= versions providing Compute_fAPAR)
    fAPAR <- NA_real_
    if ("Compute_fAPAR" %in% ls("package:prosail")) {
      fAPAR <- tryCatch({
        Compute_fAPAR(abs_dir = Ref$abs_dir,
                      abs_hem = Ref$abs_hem,
                      tts = as.numeric(input$sza))
      }, error = function(e) NA_real_)
      if (!is.null(fAPAR) && length(fAPAR) > 1) fAPAR <- fAPAR[1]
    }

  # beam transmittance: single source policy -> always use fAPAR when available (legacy)
    beam_trans <- NA_real_
    source_used <- "fAPAR"
    if (!is.na(fAPAR)) {
      beam_trans <- 1 - fAPAR
    } else if (!is.null(Ref$fcover)) { # fallback only if fAPAR unavailable
      beam_trans <- as.numeric(1 - Ref$fcover)
      source_used <- "fcover(fallback)"
    }
    if (!is.na(beam_trans)) beam_trans <- max(0, min(1, beam_trans))
    
    # alternative I_t (from I_r and I_abs_dir if both present)
    I_t_alt <- NA_real_
    if (!is.na(I_r) && !is.na(I_abs_dir)) {
      I_t_alt <- max(0, 1 - I_r - I_abs_dir)
    }
    
    # compute k if LAI > 0
    k_from_beam <- NA_real_
    k_from_alt <- NA_real_
    k_from_abs_hem <- NA_real_
    if (!is.na(Abs_hem) && is.finite(Abs_hem) && LAI_val > 0) {
      # Beer-Lambert: F = Abs_hem = 1 - exp(-k * LAI) ⇒ k = -ln(1 - F)/LAI
      if (Abs_hem >= 1) {
        k_from_abs_hem <- Inf
      } else if (Abs_hem <= 0) {
        k_from_abs_hem <- 0
      } else {
        k_from_abs_hem <- -log(1 - Abs_hem)/LAI_val
      }
    }
    if (!is.na(beam_trans) && LAI_val > 0) {
      eps <- 1e-8
      if (beam_trans > eps) k_from_beam <- -log(beam_trans) / LAI_val else k_from_beam <- Inf
    }
    if (!is.na(I_t_alt) && is.finite(I_t_alt) && LAI_val > 0) {
      if (I_t_alt > 0) k_from_alt <- -log(I_t_alt) / LAI_val else k_from_alt <- Inf
    }
    
    # If fcover missing, estimate from k_from_alt (Beer-Lambert: I_t = exp(-k*LAI); fcover ≈ 1 - I_t)
    fcover_est_k_alt <- NA_real_
    beam_trans_est_k_alt <- NA_real_
  if (is.null(Ref$fcover) && !is.na(k_from_alt) && is.finite(k_from_alt) && LAI_val > 0) {
      beam_trans_est_k_alt <- exp(-k_from_alt * LAI_val)
      beam_trans_est_k_alt <- max(0, min(1, beam_trans_est_k_alt))
      fcover_est_k_alt <- 1 - beam_trans_est_k_alt
    }
    
  list(Ref = Ref,
         wavelengths = wvl,
         par_mask = par_mask,
         I_r = I_r,
         I_abs_dir = I_abs_dir,
         Abs_hem = Abs_hem,
         beam_trans = beam_trans,
         I_t_alt = I_t_alt,
         k_from_beam = k_from_beam,
         k_from_alt = k_from_alt,
         k_from_abs_hem = k_from_abs_hem,
     fAPAR = fAPAR,
         beam_source = source_used,
         fcover_est_k_alt = fcover_est_k_alt,
         beam_trans_est_k_alt = beam_trans_est_k_alt,
         LAI_val = LAI_val)
  }, ignoreNULL = FALSE)
  
  output$summary <- renderPrint({
    res <- prosail_result()
    if (is.null(res)) {
      cat("No result yet or PRO4SAIL failed. Click Run PRO4SAIL and check notifications.\n")
      return()
    }
    
    cat("Derived inputs & results:\n")
    cat(sprintf("  LAI (plants/m² × leaf area per plant) = %.4f\n", res$LAI_val))
    cat("\nPRO4SAIL scalar outputs (if available):\n")
    if (!is.null(res$Ref$fcover)) cat(sprintf("  fcover (fraction green cover) = %.4f\n", as.numeric(res$Ref$fcover)))
    if (!is.na(res$fAPAR)) {
      cat(sprintf("  fAPAR (Compute_fAPAR; sole source for k) = %.4f\n", res$fAPAR))
    } else {
      cat("  fAPAR not available; USING fcover as temporary fallback for k.\n")
    }
    if (!is.na(res$beam_trans)) cat(sprintf("  beam_trans (1 - source) = %.4f  [source=%s]\n", res$beam_trans, res$beam_source))
    if (is.na(res$fAPAR) && is.null(res$Ref$fcover) && !is.na(res$fcover_est_k_alt)) {
      cat(sprintf("  fcover (estimated from k_alt & LAI) = %.4f\n", res$fcover_est_k_alt))
      cat(sprintf("  beam_trans (estimated) = %.4f\n", res$beam_trans_est_k_alt))
    }
    if (!is.na(res$I_abs_dir)) cat(sprintf("  abs_dir (absorptance for direct flux) = %.4f\n", res$I_abs_dir)) else cat("  abs_dir not returned\n")
    if (!is.na(res$I_r)) cat(sprintf("  I_r (PAR-averaged rsdt) = %.4f\n", res$I_r)) else cat("  I_r not computed (rsdt missing)\n")
    
    cat("\nEstimated extinction coefficients k (Beer-Lambert):\n")
    if (!is.na(res$k_from_abs_hem)) {
      cat(sprintf("  k (from abs_hem; primary) = %s\n", ifelse(is.infinite(res$k_from_abs_hem), "Inf", sprintf("%.4f", res$k_from_abs_hem))))
    } else {
      cat("  k (abs_hem) not available\n")
    }
    if (!is.na(res$k_from_beam)) {
      cat(sprintf("  k (from fAPAR/fcover beam trans) = %.4f\n", res$k_from_beam))
    }
    if (!is.na(res$k_from_alt)) {
      cat(sprintf("  k (from I_r/I_abs cross-check) = %.4f\n", res$k_from_alt))
    }
    
  cat("\nInterpretation: k is the extinction coefficient for direct beam under current solar zenith and leaf angle distribution.\n")
  cat("k may exceed 1 at high SZA or planophile canopies (G/µ_s > 1).\n")
  # Use is.null instead of is.na because Ref$fcover may be NULL (is.na(NULL) -> logical(0))
  if (is.null(res$Ref$fcover) && !is.na(res$fcover_est_k_alt)) {
      cat("fcover shown above is an inferred value assuming Beer-Lambert and fcover ≈ 1 - exp(-k*LAI).\n")
    }
  })

  # Field measurement calculations
  output$field_calc <- renderPrint({
    Io <- input$Io_field; Ir <- input$Ir_field; Ic <- input$Ic_field; LAI_f <- input$LAI_field
    if (all(is.na(c(Io,Ir,Ic,LAI_f)))) {
      cat("(Enter field measurements above to compute intercepted fraction and k_field)\n"); return()
    }
    if (any(is.na(c(Io,Ir,Ic))) || is.na(LAI_f) || LAI_f <= 0) {
      cat("Provide Io, Ir, Ic (>=0) and LAI_field > 0.\n"); return()
    }
    if (Io <= 0) { cat("Io must be > 0.\n"); return() }
    avail <- Io - Ir
    if (avail <= 0) { cat("Io - Ir must be > 0 (otherwise reflected >= incoming).\n"); return() }
    F_im <- ((Io - Ir) - Ic) / Io
    F_im <- max(0, min(1, F_im))
    k_field <- NA_real_
    if (Ic > 0) {
      ratio <- Ic / avail
      ratio <- max(1e-12, min(1, ratio))
      k_field <- -log(ratio) / LAI_f
    } else {
      k_field <- Inf
    }
    cat("Field-derived metrics (with reflection correction):\n")
    cat(sprintf("  Intercepted fraction F_im = %.4f\n", F_im))
    cat(sprintf("  k_field = %s\n", ifelse(is.infinite(k_field), "Inf (complete interception)", sprintf("%.4f", k_field))))
    if (F_im < 0.01) cat("  (Very low interception; check sensor alignment.)\n")
    if (F_im > 0.99) cat("  (Near-total interception; k may appear very large.)\n")
  })
  
  output$specPlot <- renderPlot({
    res <- prosail_result(); req(res)
    R <- NULL; T <- NULL; label_mode <- ""
    if (!is.null(res$Ref$rddt)) { R <- res$Ref$rddt; label_mode <- "bi-hemispherical"; if (!is.null(res$Ref$tddt)) T <- res$Ref$tddt }
    else if (!is.null(res$Ref$rsdt)) { R <- res$Ref$rsdt; label_mode <- "directional-hemispherical"; if (!is.null(res$Ref$tsdt)) T <- res$Ref$tsdt }
    if (is.null(R)) { plot.new(); text(0.5,0.5,"No spectral data", cex=1.2); return() }
    wv <- if (length(R) == length(res$wavelengths)) res$wavelengths else seq_len(length(R))
    R[!is.finite(R)] <- NA
    have_T <- !is.null(T) && length(T) == length(R)
    if (have_T) T[!is.finite(T)] <- NA else T <- rep(NA_real_, length(R))
    A <- if (have_T) { tmp <- 1 - R - T; tmp[tmp<0] <- 0; tmp[tmp>1] <- 1; tmp } else rep(NA_real_, length(R))
    df <- data.frame(wvl = wv, Reflectance = R, Transmittance = T, Absorptance = A)
    keep <- names(df)[colSums(!is.na(df)) > 1]
    df_long <- tidyr::pivot_longer(df[, keep, drop = FALSE], cols = -wvl, names_to = "Quantity", values_to = "Value")
    df_long <- df_long[is.finite(df_long$Value), ]
    title_txt <- paste("Canopy spectra —", label_mode)
    ggplot(df_long, aes(x = wvl, y = Value, colour = Quantity)) +
      geom_line(linewidth = 0.7, na.rm = TRUE) +
      geom_vline(xintercept = c(400,700), linetype = "dashed", alpha = 0.35) +
      labs(x = "Wavelength (nm)", y = "Fraction", title = title_txt, colour = "") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  })
  
  output$kLAIPlot <- renderPlot({
    res <- prosail_result(); req(res)
    k_primary <- res$k_from_abs_hem
    k_secondary <- res$k_from_beam
    ks <- c(k_primary, k_secondary, res$k_from_alt)
    names(ks) <- c("k_abs_hem", "k_beam", "k_alt")
    ks <- ks[!is.na(ks) & is.finite(ks)]
    if (length(ks) == 0) { plot.new(); text(0.5,0.5,"k not available yet (run PRO4SAIL).", cex=1.2); return() }
    LAI_vals <- seq(0, 8, length.out = 201)
    build_curve <- function(k){
      if (is.infinite(k) || k > 50) return(ifelse(LAI_vals==0,1,0))
      exp(-k * LAI_vals)
    }
    df_list <- lapply(names(ks), function(nm){ data.frame(LAI = LAI_vals, trans = build_curve(ks[[nm]]), source = nm) })
    df_plot <- do.call(rbind, df_list)
    source_labs <- c(k_abs_hem = "abs_hem", k_beam = "fAPAR/fcover", k_alt = "I_r/I_abs")
    cols <- c(k_abs_hem = "#7570b3", k_beam = "#1b9e77", k_alt = "#d95f02")
    ggplot(df_plot, aes(x = LAI, y = trans, colour = source)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = cols, labels = source_labs[names(cols) %in% unique(df_plot$source)]) +
      labs(x = "LAI", y = "Transmitted fraction I/I0", colour = "k source",
           title = "Transmitted fraction vs LAI for available k estimates") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
}

shinyApp(ui, server)
