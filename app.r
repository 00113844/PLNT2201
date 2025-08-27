# app.R -- corrected Shiny app for PRO4SAIL (fixes LAI scoping bug)
# Install prosail if needed:
# remotes::install_gitlab("jbferet/prosail")

library(shiny)
library(prosail)
library(ggplot2)

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
      h4("Results"),
      verbatimTextOutput("summary"),
      fluidRow(
        column(6, plotOutput("specPlot", height = "320px")),
        column(6, plotOutput("kLAIPlot", height = "320px"))
      ),
      hr(),
      h5("Notes"),
      p("k is computed as k = -ln(I_t) / LAI using beam transmittance estimate from PRO4SAIL (beam_trans = 1 - fcover)."),
      p("If spectral reflectance rsdt is available, the app also computes I_r (reflectance over PAR) and recomputes I_t = 1 - I_r - I_abs for cross-check.")
    )
  )
)

server <- function(input, output, session) {
  
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
    
    # absorbed direct if available
    I_abs_dir <- ifelse(!is.null(Ref$abs_dir), as.numeric(Ref$abs_dir), NA_real_)
    
    # beam transmittance estimate from fcover (vignette: fcover = 1 - beam transmittance)
    beam_trans <- NA_real_
    if (!is.null(Ref$fcover)) {
      beam_trans <- as.numeric(1 - Ref$fcover)
      # clip to [0,1]
      beam_trans <- max(0, min(1, beam_trans))
    }
    
    # alternative I_t (from I_r and I_abs_dir if both present)
    I_t_alt <- NA_real_
    if (!is.na(I_r) && !is.na(I_abs_dir)) {
      I_t_alt <- max(0, 1 - I_r - I_abs_dir)
    }
    
    # compute k if LAI > 0
    k_from_beam <- NA_real_
    k_from_alt <- NA_real_
    if (!is.na(beam_trans) && LAI_val > 0) {
      if (beam_trans > 0) k_from_beam <- -log(beam_trans) / LAI_val else k_from_beam <- Inf
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
         beam_trans = beam_trans,
         I_t_alt = I_t_alt,
         k_from_beam = k_from_beam,
         k_from_alt = k_from_alt,
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
    if (!is.null(res$Ref$fcover)) {
      cat(sprintf("  fcover (fraction green cover) = %.4f\n", as.numeric(res$Ref$fcover)))
      cat(sprintf("  beam_trans (1 - fcover) = %.4f\n", res$beam_trans))
    } else {
      cat("  fcover not returned by PRO4SAIL\n")
      if (!is.na(res$fcover_est_k_alt)) {
        cat(sprintf("  fcover (estimated from k_alt & LAI) = %.4f\n", res$fcover_est_k_alt))
        cat(sprintf("  beam_trans (estimated) = %.4f\n", res$beam_trans_est_k_alt))
      }
    }
    if (!is.na(res$I_abs_dir)) cat(sprintf("  abs_dir (absorptance for direct flux) = %.4f\n", res$I_abs_dir)) else cat("  abs_dir not returned\n")
    if (!is.na(res$I_r)) cat(sprintf("  I_r (PAR-averaged rsdt) = %.4f\n", res$I_r)) else cat("  I_r not computed (rsdt missing)\n")
    
    cat("\nEstimated extinction coefficients k (Beer-Lambert):\n")
    if (!is.na(res$k_from_beam)) {
      cat(sprintf("  k (from fcover → beam_trans) = %.4f\n", res$k_from_beam))
    } else {
      cat("  k (from fcover) not available\n")
    }
    if (!is.na(res$k_from_alt)) {
      cat(sprintf("  k (from I_r/I_abs cross-check) = %.4f\n", res$k_from_alt))
    } else {
      cat("  k (from I_r/I_abs) not available\n")
    }
    
    cat("\nInterpretation: k is the extinction coefficient for the direct beam under the chosen geometry and LAD.\n")
  # Use is.null instead of is.na because Ref$fcover may be NULL (is.na(NULL) -> logical(0))
  if (is.null(res$Ref$fcover) && !is.na(res$fcover_est_k_alt)) {
      cat("fcover shown above is an inferred value assuming Beer-Lambert and fcover ≈ 1 - exp(-k*LAI).\n")
    }
  })
  
  output$specPlot <- renderPlot({
    res <- prosail_result()
    req(res)  # require result present
    # choose spectral reflectance to plot
    spec_vec <- NULL
    lab <- ""
    if (!is.null(res$Ref$rddt)) {
      spec_vec <- res$Ref$rddt; lab <- "rddt (bi-hemispherical)"
    } else if (!is.null(res$Ref$rsdt)) {
      spec_vec <- res$Ref$rsdt; lab <- "rsdt (directional-hemispherical)"
    }
    
    if (is.null(spec_vec)) {
      plot.new(); text(0.5,0.5,"No spectral reflectance returned to plot", cex = 1.2)
      return()
    }
    
    # ensure wavelengths length match, else make index
    if (length(spec_vec) == length(res$wavelengths)) {
      df <- data.frame(wvl = res$wavelengths, refl = spec_vec)
    } else {
      df <- data.frame(wvl = seq_len(length(spec_vec)), refl = spec_vec)
    }
    ggplot(df, aes(x = wvl, y = refl)) +
      geom_line() +
      geom_vline(xintercept = c(400,700), linetype = "dashed", alpha = 0.4) +
      labs(x = "Wavelength (nm)", y = "Reflectance", title = paste("Simulated spectral reflectance —", lab)) +
      theme_minimal()
  })
  
  output$kLAIPlot <- renderPlot({
    res <- prosail_result()
    req(res)
    k_use <- if (!is.na(res$k_from_beam)) res$k_from_beam else res$k_from_alt
    if (is.na(k_use)) {
      plot.new(); text(0.5,0.5,"k not available yet (run PRO4SAIL to compute).", cex = 1.2); return()
    }
    LAI_vals <- seq(0, 8, length.out = 201)
    trans <- exp(-k_use * LAI_vals)
    df2 <- data.frame(LAI = LAI_vals, trans = trans)
    ggplot(df2, aes(x = LAI, y = trans)) +
      geom_line() +
      labs(x = "LAI", y = "Transmitted fraction I/I0 (direct beam approx)",
           title = sprintf("Transmitted fraction vs LAI (k = %.3f)", k_use)) +
      theme_minimal()
  })
}

shinyApp(ui, server)
