# ---- Packages ----

library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- App data ----
# Self-contained bundle built by scripts/jobs/build_shiny_app_data.R; re-run that
# script when posteriors or data change. The app has no other runtime data dependency.

app_data         <- readRDS("shiny_data.RDS")
pred_curves      <- app_data$pred_curves
observations     <- app_data$observations
coral_cat_mid    <- app_data$coral_cat_mid     # named num vector: category label -> midpoint (proportion 0-1)
coral_cat_breaks <- app_data$coral_cat_breaks  # bin edges (proportion 0-1), same order as coral_cat_mid
uncertainty_comparison <- app_data$uncertainty_comparison  # RHIS-inclusive descriptive method comparison (reef level)
uncertainty_cots_log1p <- app_data$uncertainty_cots_log1p  # COTS-only companion, mean/SD computed on log1p(count) per unit
uncertainty_pooled            <- app_data$uncertainty_pooled            # one row per method x quantity, pooled across ALL reefs
uncertainty_pooled_cots_log1p <- app_data$uncertainty_pooled_cots_log1p # pooled COTS companion, mean/SD on log1p(count) per unit
uncertainty_band_comparison   <- app_data$uncertainty_band_comparison   # grouped by reef-ID band (e.g. "21-XXX")
uncertainty_band_cots_log1p   <- app_data$uncertainty_band_cots_log1p   # band-grouped COTS companion, mean/SD on log1p(count) per unit
cots_log1p_breaks      <- app_data$cots_log1p_breaks       # original-count axis breaks, log1p-transformed at plot time

# Harmonise display labels in the bundle once, so every tab/plot/table agrees.
clean_label <- function(x) {
  x <- gsub("RS Surface", "ReefScan Surface", x, fixed = TRUE)
  x <- gsub("RS Deep", "ReefScan Deep", x, fixed = TRUE)
  gsub(" -- ", " \u2013 ", x, fixed = TRUE)
}
pred_curves  <- pred_curves  |>
  mutate(across(c(from, to, xlab, ylab), clean_label),
         # COTS rows are NA here; treated as not a percent-scale axis.
         pct_axes = coalesce(pct_axes, FALSE))
observations <- observations |> mutate(across(c(from, to), clean_label))

n_reefs_total <- n_distinct(uncertainty_comparison$reef_label)
data_built    <- format(as.Date(app_data$built_at), "%d %B %Y")

# ---- Branding / app metadata ----
# Set `contact` and `citation` to show them on the About tab and footer. CSIRO logo
# goes in www/ as csiro-logo.svg or .png (picked up automatically).
app_meta <- list(
  title    = "Rosetta-GBR",
  subtitle = "Method Calibration Explorer (CCIP-BR-04)",
  version  = "1.0",
  year     = "2026",
  contact_name = "Dan Gladish",
  contact  = "dan.gladish@csiro.au",
  citation = NULL    # TBD; set to a string to show it on the About tab
)

csiro <- list(
  midnight = "#00313C",
  midday   = "#00A9CE",
  deep     = "#00778F",   # darker teal-blue: accessible contrast with white text
  tint     = "#E6F6FA"
)
col_curve <- csiro$midnight
col_band  <- csiro$midday
col_pt    <- "#B2182B"
method_cols <- c(
  "Manta tow"        = csiro$midnight,
  "ReefScan Surface" = csiro$midday,
  "ReefScan Deep"    = "#E69F00",
  "RHIS"             = "#7B5EA7"
)

theme_csiro <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = csiro$tint, colour = "grey70"),
      strip.text       = element_text(face = "bold", colour = csiro$midnight),
      plot.title       = element_text(face = "bold", colour = csiro$midnight),
      panel.grid.minor = element_blank()
    )
}

csiro_theme <- bs_theme(
  version = 5,
  primary = csiro$deep,
  "navbar-bg" = csiro$midnight,
  "link-color" = csiro$deep
) |>
  bs_add_rules("
    .navbar-brand { white-space: normal; }
    .brand-title  { font-weight: 600; font-size: 1.15rem; display: block; line-height: 1.2; }
    .brand-sub    { font-size: 0.75rem; opacity: 0.8; display: block; font-weight: 400; }
    .brand-logo   { height: 38px; margin-right: 12px; vertical-align: middle; }
    .brand-wrap   { display: inline-flex; align-items: center; }
    .navbar       { border-bottom: 4px solid #00A9CE; }
    .app-footer   { padding: 12px 24px; font-size: 0.8rem; color: #555;
                    border-top: 1px solid #ddd; background: #f8f9fa; }
    .app-footer a { color: #00778F; }
    .about-body h5 { color: #00313C; margin-top: 1rem; }
  ")

logo_file <- sort(list.files("www", pattern = "^csiro-logo\\.(svg|png)$"), decreasing = TRUE)  # prefer svg
logo_tag  <- if (length(logo_file) > 0) {
  tags$img(src = logo_file[1], class = "brand-logo", alt = "CSIRO logo")
}

family_choices <- c("Coral cover", "COTS cross-method", "eDNA")
level_choices  <- c("Tow", "Cull-site", "Reef")
calib_families <- family_choices

# Labels for the family selector (left); internal values (right) drive the logic.
family_labels <- c(
  "Coral cover"                                                  = "Coral cover",
  "COTS cross-method (Manta / ReefScan Surface / ReefScan Deep)" = "COTS cross-method",
  "COTS vs. eDNA (Manta / ReefScan Surface vs. eDNA)"            = "eDNA"
)

# Shared caveat wording (all three calibration tabs)
msg_edna <- function(plural = TRUE) paste0(
  "Only ~17\u201324 reefs have both measurements, so ",
  if (plural) "these curves are" else "this curve is",
  " indicative only."
)
msg_deep_cull <- "ReefScan Deep cull-site coral models have no reef effect (too few reefs) and are exploratory."
msg_prob <- paste(
  "Probability curves: the prediction interval covers one new 0/1 tow, so it spans nearly",
  "0\u2013100%; the credible interval is more informative."
)
msg_pi_missing <- "No prediction interval for some curves; the credible interval is shown for those."

# Optional detail lives in tooltips on the controls, not in the caveat box.
tip_box <- "Boxplots (with raw points) of observations grouped by Manta tow's 11 ordinal coral-cover categories."
tip_se  <- paste(
  "Analytic standard error of each cull-site/reef observation (not the model's interval).",
  "COTS errors are Poisson-based, so they are lower bounds."
)

# ---- UI ----

about_ui <- nav_panel(
  title = "About",
  icon  = icon("circle-info"),
  div(
    class = "about-body container-fluid", style = "max-width: 1000px;",
    layout_columns(
      col_widths = c(12),
      card(
        card_header("Purpose"),
        p("Reef Joint Field Management Program (RJFMP) voyages in 2025\u201326 measured coral cover and",
          "crown-of-thorns starfish (COTS) with several methods. This tool translates a reading from",
          "one method into the value another would be expected to report, with uncertainty, to support",
          "the COTS Control Program. Methods: manta tow, ReefScan Surface (towed camera), ReefScan Deep,",
          "eDNA and, descriptively only, the Reef Health and Impact Survey (RHIS).")
      ),
      card(
        card_header("How to use"),
        tags$ul(
          tags$li(tags$b("Predicting to / from / from and to"), " \u2013 pick a target, a source, or one specific pairing."),
          tags$li(tags$b("Uncertainty comparison"), " \u2013 spread of readings across methods (includes RHIS)."),
          tags$li("Move the \u201cObserved value\u201d slider to read the expected value (red point) and its interval: ",
                  "a ", tags$b("prediction interval"), " (where a single new observation is expected to fall) or ",
                  "a ", tags$b("credible interval"), " (uncertainty in the average relationship).")
        )
      ),
      card(
        card_header("Methods in brief"),
        p("Bayesian mixed models, fitted for each direction of each method pair at three scales: manta-tow",
          "segments (tow), cull sites and reefs. Coral cover uses Beta regression, COTS counts negative",
          "binomial with an effort offset, and presence/absence binomial; predictions assume a typical",
          "sampling effort for the predicted method. Manta coral cover is ordinal (converted to category",
          "midpoints), and the \u201ccombined\u201d manta COTS predictor adds the feeding-scar count to the COTS count (scars for FMPCRS 2025-02 imputed from the present/common category). The table's",
          "\u201cslope\u201d is on the model's link scale.")
      ),
      card(
        card_header("Data and scope"),
        p(sprintf(
          "Analyses use %d reefs from five voyages (FMPCRS 2025-02 to 2026-02); three of the eight",
          n_reefs_total),
          "voyages (2024-05, 2025-01, 2026-03) were excluded after a data-quality sensitivity analysis.",
          "Models involving ReefScan Deep were the least stable under it."),
        p(sprintf("Model outputs bundled on %s.", data_built))
      ),
      card(
        card_header("Limitations"),
        tags$ul(
          tags$li("Small samples: reef level often has only 13\u201332 reefs, fewer for pairings with ReefScan Deep or eDNA."),
          tags$li("Curves describe these data only: don't extrapolate beyond the plotted range, or use one spatial level in place of another."),
          tags$li("COTS counts are sparse and overdispersed; observation standard errors are Poisson-based lower bounds.")
        )
      ),
      if (!is.null(app_meta$citation) || !is.null(app_meta$contact)) {
        card(
          card_header("Citation and contact"),
          if (!is.null(app_meta$citation)) p(app_meta$citation),
          if (!is.null(app_meta$contact)) p("Contact: ", app_meta$contact_name, " (", tags$a(href = paste0("mailto:", app_meta$contact), app_meta$contact), ")")
        )
      }
    )
  )
)

footer_ui <- div(
  class = "app-footer",
  HTML(sprintf(
    "&copy; %s Commonwealth Scientific and Industrial Research Organisation (CSIRO) &middot; v%s",
    app_meta$year, app_meta$version
  )),
  tags$br(),
  "Provided for information; results are statistical estimates with uncertainty and are not a substitute for expert judgement. ",
  if (!is.null(app_meta$contact)) tags$a(href = paste0("mailto:", app_meta$contact), "Contact")
)

ui <- page_navbar(
  title = tags$span(
    class = "brand-wrap", logo_tag,
    tags$span(
      tags$span(class = "brand-title", app_meta$title),
      tags$span(class = "brand-sub", app_meta$subtitle)
    )
  ),
  window_title = paste("CSIRO |", app_meta$title),
  lang   = "en",
  theme  = csiro_theme,
  navbar_options = navbar_options(bg = csiro$midnight, theme = "dark"),
  footer = footer_ui,

  nav_panel(
    title = "Predicting to",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,

        selectInput(
          "family", "Comparison type:",
          choices = family_labels, selected = "Coral cover"
        ),

        conditionalPanel(
          condition = "input.family != 'eDNA'",
          input_switch("show_all_levels", "Show all spatial levels together", value = FALSE),
          conditionalPanel(
            condition = "!input.show_all_levels",
            selectInput("level", "Spatial level:", choices = level_choices, selected = "Reef")
          )
        ),

        uiOutput("to_ui"),

        conditionalPanel(
          condition = "input.family == 'Coral cover'",
          tooltip(
            input_switch("box_manta", "Show Manta tow coral-category boxplots", value = FALSE),
            tip_box
          )
        ),

        uiOutput("slider_ui"),

        radioButtons(
          "interval_type", "Uncertainty band shown:",
          choices = c(
            "Prediction interval (range for a new observation)" = "pi",
            "Credible interval (uncertainty on the mean curve)"  = "ci"
          ),
          selected = "pi"
        ),

        tooltip(
          input_switch("show_obs_se", "Show \u00b1SE error bars (cull-site/reef level only)", value = FALSE),
          tip_se
        ),

        hr(),
        tags$small(textOutput("caveat"))
      ),

      card(
        full_screen = TRUE,
        card_header("Calibration curves \u2013 every source method that can predict the chosen target"),
        uiOutput("plot_ui")
      ),
      uiOutput("secondary_plot_ui"),
      card(
        card_header("Prediction at the chosen value"),
        tableOutput("pred_table")
      )
    )
  ),

  nav_panel(
    title = "Predicting from",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,

        selectInput(
          "pf_family", "Comparison type:",
          choices = family_labels, selected = "Coral cover"
        ),

        conditionalPanel(
          condition = "input.pf_family != 'eDNA'",
          input_switch("pf_show_all_levels", "Show all spatial levels together", value = FALSE),
          conditionalPanel(
            condition = "!input.pf_show_all_levels",
            selectInput("pf_level", "Spatial level:", choices = level_choices, selected = "Reef")
          )
        ),

        uiOutput("pf_from_ui"),

        conditionalPanel(
          condition = "input.pf_family == 'Coral cover' && input.pf_from == 'Manta tow'",
          tooltip(
            input_switch("pf_box_manta", "Show Manta tow coral-category boxplots", value = FALSE),
            tip_box
          )
        ),

        uiOutput("pf_slider_ui"),

        radioButtons(
          "pf_interval_type", "Uncertainty band shown:",
          choices = c(
            "Prediction interval (range for a new observation)" = "pi",
            "Credible interval (uncertainty on the mean curve)"  = "ci"
          ),
          selected = "pi"
        ),

        tooltip(
          input_switch("pf_show_obs_se", "Show \u00b1SE error bars (cull-site/reef level only)", value = FALSE),
          tip_se
        ),

        hr(),
        tags$small("Each panel has its own scale; units differ across targets."),
        tags$br(),
        tags$small(textOutput("pf_caveat"))
      ),

      card(
        full_screen = TRUE,
        card_header("Calibration curves \u2013 the chosen source method predicting each target"),
        uiOutput("pf_plot_ui")
      ),
      uiOutput("pf_secondary_plot_ui"),
      card(
        card_header("Prediction at the chosen value"),
        tableOutput("pf_pred_table")
      )
    )
  ),

  nav_panel(
    title = "Predicting from and to",
    layout_sidebar(
      sidebar = sidebar(
        width = 360,

        selectInput(
          "ft_family", "Comparison type:",
          choices = family_labels, selected = "Coral cover"
        ),

        conditionalPanel(
          condition = "input.ft_family != 'eDNA'",
          input_switch("ft_show_all_levels", "Show all spatial levels together", value = FALSE),
          conditionalPanel(
            condition = "!input.ft_show_all_levels",
            selectInput("ft_level", "Spatial level:", choices = level_choices, selected = "Reef")
          )
        ),

        uiOutput("ft_from_ui"),
        uiOutput("ft_to_ui"),

        conditionalPanel(
          condition = "input.ft_family == 'Coral cover' && input.ft_from == 'Manta tow'",
          tooltip(
            input_switch("ft_box_manta", "Show Manta tow coral-category boxplots", value = FALSE),
            tip_box
          )
        ),

        uiOutput("ft_slider_ui"),

        radioButtons(
          "ft_interval_type", "Uncertainty band shown:",
          choices = c(
            "Prediction interval (range for a new observation)" = "pi",
            "Credible interval (uncertainty on the mean curve)"  = "ci"
          ),
          selected = "pi"
        ),

        tooltip(
          input_switch("ft_show_obs_se", "Show \u00b1SE error bars (cull-site/reef level only)", value = FALSE),
          tip_se
        ),

        hr(),
        tags$small(textOutput("ft_caveat"))
      ),

      card(
        full_screen = TRUE,
        card_header("Calibration curve \u2013 chosen source and target"),
        uiOutput("ft_plot_ui")
      ),
      uiOutput("ft_secondary_plot_ui"),
      card(
        card_header("Prediction at the chosen value"),
        tableOutput("ft_pred_table")
      )
    )
  ),

  nav_panel(
    title = "Uncertainty comparison",
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        radioButtons("uc_quantity", "Quantity:",
                     choices = c("Hard coral cover (%)", "COTS count (per tow or survey)"),
                     selected = "Hard coral cover (%)"),
        conditionalPanel(
          condition = "input.uc_quantity == 'COTS count (per tow or survey)'",
          input_switch(
            "uc_cots_log1p",
            "Show COTS counts on log1p scale (compresses high-count outliers)",
            value = FALSE
          )
        ),
        radioButtons(
          "uc_grouping", "Grouping:",
          choices = c(
            "Per reef"                               = "reef",
            "Per reef-ID band (e.g. 10-XXX, 21-XXX)" = "band",
            "Pooled across all reefs"                = "all"
          ),
          selected = "reef"
        ),
        conditionalPanel(
          condition = "input.uc_grouping != 'all'",
          radioButtons(
            "uc_arrange", "Arrange by:",
            choices = c(
              "Reef / band, with tools side by side"  = "unit",
              "Tool, with reefs / bands listed within" = "tool"
            ),
            selected = "unit"
          )
        ),
        conditionalPanel(
          condition = "input.uc_grouping == 'reef'",
          textInput("uc_reef_filter", "Filter reefs (substring match, blank = all):", value = "")
        ),
        conditionalPanel(
          condition = "input.uc_grouping == 'band'",
          textInput("uc_band_filter", "Filter bands (substring match, blank = all):", value = "")
        ),
        uiOutput("uc_point_desc"),
        hr(),
        tags$small(textOutput("uc_caveat"))
      ),

      card(
        full_screen = TRUE,
        card_header(textOutput("uc_card_header")),
        uiOutput("uc_plot_ui")
      )
    )
  ),

  nav_spacer(),
  about_ui
)

# ---- Server ----

# Wrap a plot's (long) subtitle so it doesn't run off the panel
wrap_subtitle <- function(p, width = 80) {
  s <- tryCatch(ggplot2::get_labs(p)$subtitle, error = function(e) p$labels$subtitle)
  if (is.character(s) && length(s) == 1) {
    p <- p + ggplot2::labs(subtitle = paste(strwrap(s, width), collapse = "\n"))
  }
  p
}

server <- function(input, output, session) {

  ## ---- Convenience reactives for the view toggles ----
  multi_level_active <- reactive({
    isTRUE(input$show_all_levels) && isTRUE(input$family %in% c("Coral cover", "COTS cross-method"))
  })

  # Primary boxplot: Manta tow's raw value is ordinal, so facets whose source ("from")
  # method is Manta tow box observations by category instead of scattering them.
  # The split is done per facet inside pred_plot's renderPlot.
  box_toggle_active <- reactive({
    isTRUE(input$box_manta) && isTRUE(input$family == "Coral cover")
  })

  # Secondary boxplot: when Manta tow is the predicted ("to") method, bin its response
  # into the 11 ordinal categories and box the source method's raw values within each.
  secondary_applicable <- reactive({
    isTRUE(input$box_manta) && isTRUE(input$family == "Coral cover") &&
      isTRUE(input$to == "Manta tow")
  })

  ## ---- "Predict to" choices: depend on family (and level, unless "show all levels" is on) ----
  output$to_ui <- renderUI({
    req(input$family)
    d <- pred_curves |> filter(family_app == input$family)
    if (input$family != "eDNA" && !multi_level_active()) {
      req(input$level)
      d <- d |> filter(level == input$level)
    }
    choices <- d |> pull(to) |> unique() |> sort()
    selectInput("to", "Predict to:", choices = choices)
  })

  ## ---- Curves/observations for the current selection (all levels kept when "show all levels" is on) ----
  current_curves <- reactive({
    req(input$family, input$to)
    d <- pred_curves |> filter(family_app == input$family, to == input$to)
    if (input$family != "eDNA" && !multi_level_active()) {
      req(input$level)
      d <- d |> filter(level == input$level)
    }
    d
  })

  current_obs <- reactive({
    req(input$family, input$to)
    d <- observations |> filter(family_app == input$family, to == input$to)
    if (input$family != "eDNA" && !multi_level_active()) {
      req(input$level)
      d <- d |> filter(level == input$level)
    }
    d
  })

  ## ---- Slider range follows the selection's observed x-range (union across "from" facets/levels) ----
  output$slider_ui <- renderUI({
    d <- current_curves()
    req(nrow(d) > 0)
    rng  <- range(d$x, na.rm = TRUE)
    step <- signif(diff(rng) / 200, 3)
    if (!is.finite(step) || step <= 0) step <- 0.001
    sliderInput(
      "obs_value", "Observed value:",
      min = signif(rng[1], 4), max = signif(rng[2], 4),
      value = signif(rng[1], 4), step = step
    )
  })

  ## ---- Plot height grows to fit 3 stacked levels when that view is on ----
  output$plot_ui <- renderUI({
    plotOutput("pred_plot", height = if (multi_level_active()) "820px" else "500px")
  })

  ## ---- Interpolate mean + interval at the slider value, per (level, from) ----
  interp_df <- reactive({
    d <- current_curves()
    req(nrow(d) > 0, !is.null(input$obs_value))
    use_pi <- input$interval_type == "pi"

    d <- d |> mutate(level_key = if_else(is.na(level), "single", level))
    combos <- d |> distinct(level, level_key, from)

    map_dfr(seq_len(nrow(combos)), function(i) {
      dd <- d |> filter(from == combos$from[i], level_key == combos$level_key[i]) |> arrange(x)
      pi_avail <- !all(is.na(dd$pi_lo95))
      lo_col <- if (use_pi && pi_avail) dd$pi_lo95 else dd$lo95
      hi_col <- if (use_pi && pi_avail) dd$pi_hi95 else dd$hi95
      tibble(
        level = combos$level[i], from = combos$from[i], x = input$obs_value,
        n = dd$n[1], beta_mean = dd$beta_mean[1], beta_lo95 = dd$beta_lo95[1], beta_hi95 = dd$beta_hi95[1],
        is_probability = isTRUE(dd$is_probability[1]),
        pred = approx(dd$x, dd$mean, xout = input$obs_value, rule = 2)$y,
        lo   = approx(dd$x, lo_col,  xout = input$obs_value, rule = 2)$y,
        hi   = approx(dd$x, hi_col,  xout = input$obs_value, rule = 2)$y,
        pi_available = pi_avail,
        used_interval = if (use_pi && pi_avail) "Prediction" else "Credible"
      )
    })
  })

  ## ---- Helper: bin coral cover into the 11 ordinal categories (returns category midpoint) ----
  bin_to_category_mid <- function(x) {
    cat <- cut(x, breaks = coral_cat_breaks, labels = names(coral_cat_mid), include.lowest = TRUE)
    unname(coral_cat_mid[as.character(cat)])
  }

  ## ---- Primary plot: facets by "from" (or level x "from" in all-levels view) ----
  ## Coral-cover facets with a Manta tow source can show observations as boxplots by ordinal
  ## category. x is free per facet, y is shared. pct_axes rows scale x and y by 100;
  ## is_probability rows scale only y by 100.
  output$pred_plot <- renderPlot({
    d <- current_curves()
    o <- current_obs()
    req(nrow(d) > 0)

    multi_level <- multi_level_active()
    box_on      <- box_toggle_active()

    if (multi_level) {
      d <- d |> mutate(level = factor(level, levels = level_choices))
      o <- o |> mutate(level = factor(level, levels = level_choices))
    }

    is_pct  <- isTRUE(unique(d$pct_axes)[1])
    is_prob <- isTRUE(unique(d$is_probability)[1])
    mult_x  <- if (is_pct) 100 else 1
    mult_y  <- if (is_pct || is_prob) 100 else 1

    d2 <- d |> mutate(x = x * mult_x, mean = mean * mult_y, lo95 = lo95 * mult_y, hi95 = hi95 * mult_y,
                       pi_lo95 = pi_lo95 * mult_y, pi_hi95 = pi_hi95 * mult_y)
    o2 <- o |> mutate(x = x * mult_x, y = y * mult_y,
                       x_se = if ("x_se" %in% names(o)) x_se * mult_x else NA_real_,
                       y_se = if ("y_se" %in% names(o)) y_se * mult_y else NA_real_)
    interp <- interp_df() |> mutate(x = x * mult_x, pred = pred * mult_y, lo = lo * mult_y, hi = hi * mult_y)
    if (multi_level) interp <- interp |> mutate(level = factor(level, levels = level_choices))

    has_pi_any <- any(!is.na(d2$pi_lo95))
    use_pi     <- isTRUE(input$interval_type == "pi") && has_pi_any
    band_lo <- if (use_pi) d2$pi_lo95 else d2$lo95
    band_hi <- if (use_pi) d2$pi_hi95 else d2$hi95
    d2 <- d2 |> mutate(band_lo = band_lo, band_hi = band_hi)

    p <- ggplot(d2, aes(x = x))

    if (is_pct) {
      p <- p + geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4)
    }
    if (has_pi_any) {
      # Reference line for the interval not shown as the ribbon
      if (use_pi) {
        p <- p + geom_line(aes(y = lo95), colour = "grey45", linetype = "dashed", na.rm = TRUE) +
          geom_line(aes(y = hi95), colour = "grey45", linetype = "dashed", na.rm = TRUE)
      } else {
        p <- p +
          geom_line(aes(y = pi_lo95), colour = "grey45", linetype = "dashed", na.rm = TRUE) +
          geom_line(aes(y = pi_hi95), colour = "grey45", linetype = "dashed", na.rm = TRUE)
      }
    }

    p <- p + geom_ribbon(aes(ymin = band_lo, ymax = band_hi), fill = col_band, alpha = 0.25)

    o_manta_box <- box_on && any(o$from == "Manta tow")
    o_rest      <- if (o_manta_box) o2 |> filter(from != "Manta tow") else o2

    if (o_manta_box) {
      o_box <- o |>
        filter(from == "Manta tow") |>
        mutate(cat_mid = bin_to_category_mid(x) * mult_x, y = y * mult_y) |>
        filter(!is.na(cat_mid))
      if (multi_level) o_box <- o_box |> mutate(level = factor(level, levels = level_choices))
      p <- p +
        geom_boxplot(
          data = o_box, aes(x = cat_mid, y = y, group = cat_mid),
          width = 1, fill = "grey90", colour = "grey35", alpha = 0.8, outlier.shape = NA
        ) +
        geom_jitter(
          data = o_box, aes(x = cat_mid, y = y),
          width = 0.35, height = 0, alpha = 0.3, size = 0.8, shape = 16, colour = "grey20"
        )
    }

    if (nrow(o_rest) > 0) {
      if (isTRUE(input$show_obs_se) && any(!is.na(o_rest$x_se) | !is.na(o_rest$y_se))) {
        p <- p +
          geom_errorbar(data = o_rest, aes(x = x, ymin = y - y_se, ymax = y + y_se), width = 0, colour = "grey55", alpha = 0.5, na.rm = TRUE) +
          geom_errorbar(data = o_rest, aes(y = y, xmin = x - x_se, xmax = x + x_se), orientation = "y", width = 0, colour = "grey55", alpha = 0.5, na.rm = TRUE)
      }
      p <- p + geom_point(data = o_rest, aes(x = x, y = y), shape = 1, colour = "grey30", alpha = 0.5)
    }

    p <- p +
      geom_line(aes(y = mean), colour = col_curve, linewidth = 1) +
      geom_vline(xintercept = unique(interp$x), linetype = "dashed", colour = "grey40") +
      geom_point(data = interp, aes(x = x, y = pred), size = 3, colour = col_pt) +
      geom_errorbar(data = interp, aes(x = x, ymin = lo, ymax = hi), width = 0, colour = col_pt)

    p <- p + (if (multi_level) {
      facet_wrap(vars(level, from), ncol = 2, scales = "free_x")
    } else {
      facet_wrap(~ from, scales = "free_x", ncol = 2)
    })

    # x: each facet is a different source measuring the same quantity, so name it
    # explicitly. eDNA sources differ in units, so use each row's own xlab.
    xlab_display <- if (input$family == "Coral cover") {
      "Source method hard coral cover (%)"
    } else if (input$family == "COTS cross-method") {
      "Source method COTS density (count / effort)"
    } else if (n_distinct(d$xlab) == 1) {
      unique(d$xlab)[1]
    } else {
      "Source method value (units differ by panel; see panel titles)"
    }

    # y: the fixed "to" target; append the unit suffix.
    ylab_base <- if (n_distinct(d$ylab) == 1) unique(d$ylab)[1] else input$to
    ylab_display <- paste0(ylab_base, if (is_pct || is_prob) " (%)" else "")

    p <- p +
      labs(
        x = xlab_display,
        y = ylab_display,
        title = paste("Predicting:", input$to),
        subtitle = if (multi_level) {
          "All spatial levels"
        } else if (input$family != "eDNA") {
          paste("Level:", input$level)
        } else NULL
      ) +
      theme_csiro(13)

    if (is_pct) {
      p <- p +
        scale_x_continuous(labels = function(v) paste0(v, "%")) +
        scale_y_continuous(labels = function(v) paste0(v, "%"))
    } else if (is_prob) {
      p <- p + scale_y_continuous(labels = function(v) paste0(v, "%"))
    }

    p
  })

  ## ---- Secondary plot: Manta tow is the "to" method and boxplot toggle is on ----
  ## Descriptive: bins Manta's response into 11 ordinal categories and boxes the source
  ## method's raw values, faceted by "from". Dashed line = current slider value.
  output$secondary_plot_ui <- renderUI({
    if (!secondary_applicable()) return(NULL)
    card(
      full_screen = TRUE,
      card_header(
        "Manta tow as the PREDICTED method \u2014 empirical distribution of the source method(s) by Manta's coral-cover category"
      ),
      plotOutput("secondary_plot", height = "380px")
    )
  })

  output$secondary_plot <- renderPlot({
    req(secondary_applicable())
    o <- current_obs()
    req(nrow(o) > 0)

    multi_level <- multi_level_active()
    if (multi_level) o <- o |> mutate(level = factor(level, levels = level_choices))

    # Axes swapped vs the primary plot: Manta (o$y) is categorised, source (o$x) is boxed.
    o_box <- o |>
      mutate(cat_mid = bin_to_category_mid(y) * 100, from_pct = x * 100) |>
      filter(!is.na(cat_mid))
    req(nrow(o_box) > 0)

    p <- ggplot(o_box, aes(x = cat_mid, y = from_pct)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4) +
      geom_boxplot(aes(group = cat_mid), width = 1, fill = "grey90", colour = "grey35",
                   alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.35, height = 0, alpha = 0.3, size = 0.8, shape = 16, colour = "grey20") +
      geom_hline(yintercept = input$obs_value * 100, linetype = "dashed", colour = "grey40") +
      labs(
        x = "Manta tow coral cover category (midpoint, %)",
        y = "Source method hard coral cover (%; see facet label)",
        subtitle = paste0(
          "Dashed line = current slider value (", round(input$obs_value * 100, 1), "%). ",
          if (multi_level) "All spatial levels; " else paste0("Level: ", input$level, "; "),
          "boxplots are empirical (raw data), not the fitted model curve."
        )
      ) +
      scale_x_continuous(labels = function(v) paste0(v, "%")) +
      scale_y_continuous(labels = function(v) paste0(v, "%")) +
      theme_csiro(12)

    if (multi_level) {
      p <- p + facet_wrap(vars(level, from), ncol = 3)
    } else {
      p <- p + facet_wrap(~ from, ncol = 3)
    }

    p
  })

  ## ---- "Predicting from" tab: fixed source, facets over targets ----
  ## Free x and y scales per facet, since target units and ranges differ.

  pf_multi_level_active <- reactive({
    isTRUE(input$pf_show_all_levels) && isTRUE(input$pf_family %in% c("Coral cover", "COTS cross-method"))
  })

  # Main-plot boxplot: only when the fixed source is Manta tow (ordinal x). Coral cover only,
  # so pct_axes/is_probability are uniform across facets.
  pf_box_toggle_active <- reactive({
    isTRUE(input$pf_box_manta) && isTRUE(input$pf_family == "Coral cover") && isTRUE(input$pf_from == "Manta tow")
  })

  # Secondary plot: when the source is RS Surface/Deep and Manta tow is a faceted target,
  # Manta's ordinal scale gets its own boxplot-by-category view.
  pf_secondary_applicable <- reactive({
    isTRUE(input$pf_box_manta) && isTRUE(input$pf_family == "Coral cover") &&
      isTRUE(input$pf_from %in% c("ReefScan Surface", "ReefScan Deep"))
  })

  output$pf_from_ui <- renderUI({
    req(input$pf_family)
    d <- pred_curves |> filter(family_app == input$pf_family)
    if (input$pf_family != "eDNA" && !pf_multi_level_active()) {
      req(input$pf_level)
      d <- d |> filter(level == input$pf_level)
    }
    choices <- d |> pull(from) |> unique() |> sort()
    selectInput("pf_from", "Predict from:", choices = choices)
  })

  pf_current_curves <- reactive({
    req(input$pf_family, input$pf_from)
    d <- pred_curves |> filter(family_app == input$pf_family, from == input$pf_from)
    if (input$pf_family != "eDNA" && !pf_multi_level_active()) {
      req(input$pf_level)
      d <- d |> filter(level == input$pf_level)
    }
    d
  })

  pf_current_obs <- reactive({
    req(input$pf_family, input$pf_from)
    d <- observations |> filter(family_app == input$pf_family, from == input$pf_from)
    if (input$pf_family != "eDNA" && !pf_multi_level_active()) {
      req(input$pf_level)
      d <- d |> filter(level == input$pf_level)
    }
    d
  })

  ## pct_axes is uniform for the whole selection (set per family, not per model row).
  pf_slider_range <- reactive({
    d <- pf_current_curves()
    req(nrow(d) > 0)
    range(d$x, na.rm = TRUE)
  })

  output$pf_slider_ui <- renderUI({
    rng  <- pf_slider_range()
    step <- signif(diff(rng) / 200, 3)
    if (!is.finite(step) || step <= 0) step <- 0.001
    sliderInput(
      "pf_obs_value", "Observed value:",
      min = signif(rng[1], 4), max = signif(rng[2], 4),
      value = signif(rng[1], 4), step = step
    )
  })

  output$pf_plot_ui <- renderUI({
    plotOutput("pf_pred_plot", height = if (pf_multi_level_active()) "820px" else "500px")
  })

  pf_interp_df <- reactive({
    d <- pf_current_curves()
    req(nrow(d) > 0, !is.null(input$pf_obs_value))
    use_pi <- input$pf_interval_type == "pi"

    d <- d |> mutate(level_key = if_else(is.na(level), "single", level))
    combos <- d |> distinct(level, level_key, to, pct_axes, is_probability)

    map_dfr(seq_len(nrow(combos)), function(i) {
      dd <- d |> filter(to == combos$to[i], level_key == combos$level_key[i]) |> arrange(x)
      pi_avail <- !all(is.na(dd$pi_lo95))
      lo_col <- if (use_pi && pi_avail) dd$pi_lo95 else dd$lo95
      hi_col <- if (use_pi && pi_avail) dd$pi_hi95 else dd$hi95
      tibble(
        level = combos$level[i], to = combos$to[i], x = input$pf_obs_value,
        n = dd$n[1], beta_mean = dd$beta_mean[1], beta_lo95 = dd$beta_lo95[1], beta_hi95 = dd$beta_hi95[1],
        pct_axes = combos$pct_axes[i], is_probability = combos$is_probability[i],
        pred = approx(dd$x, dd$mean, xout = input$pf_obs_value, rule = 2)$y,
        lo   = approx(dd$x, lo_col,  xout = input$pf_obs_value, rule = 2)$y,
        hi   = approx(dd$x, hi_col,  xout = input$pf_obs_value, rule = 2)$y,
        pi_available = pi_avail,
        used_interval = if (use_pi && pi_avail) "Prediction" else "Credible"
      )
    })
  })

  ## ---- Helper: unit tag per facet (target units can differ when the source is fixed) ----
  pf_unit_tag <- function(is_probability, pct_axes) {
    dplyr::case_when(
      is_probability ~ " [probability, %]",
      pct_axes       ~ " [%]",
      TRUE           ~ " [count/rate]"
    )
  }

  output$pf_pred_plot <- renderPlot({
    d <- pf_current_curves()
    o <- pf_current_obs()
    req(nrow(d) > 0)

    multi_level <- pf_multi_level_active()
    box_on      <- pf_box_toggle_active()

    if (multi_level) {
      d <- d |> mutate(level = factor(level, levels = level_choices))
      o <- o |> mutate(level = factor(level, levels = level_choices))
    }

    # pct_axes is one flag for the shared x-scale; is_probability varies per "to" facet,
    # so y gets a row-wise multiplier.
    is_pct_family <- isTRUE(unique(d$pct_axes)[1])
    mult_x <- if (is_pct_family) 100 else 1

    d2 <- d |>
      mutate(
        mult_y = if_else(pct_axes | is_probability, 100, 1),
        facet_lab = paste0(to, pf_unit_tag(is_probability, pct_axes)),
        x = x * mult_x, mean = mean * mult_y, lo95 = lo95 * mult_y, hi95 = hi95 * mult_y,
        pi_lo95 = pi_lo95 * mult_y, pi_hi95 = pi_hi95 * mult_y
      )

    flag_lookup <- d |> distinct(level, to, is_probability)
    o2 <- o |>
      left_join(flag_lookup, by = intersect(c("level", "to"), names(o))) |>
      mutate(
        mult_y = if_else(is_pct_family | is_probability, 100, 1),
        facet_lab = paste0(to, pf_unit_tag(is_probability, is_pct_family)),
        x = x * mult_x, y = y * mult_y,
        x_se = if ("x_se" %in% names(o)) x_se * mult_x else NA_real_,
        y_se = if ("y_se" %in% names(o)) y_se * mult_y else NA_real_
      )

    interp <- pf_interp_df() |>
      mutate(
        mult_y = if_else(pct_axes | is_probability, 100, 1),
        facet_lab = paste0(to, pf_unit_tag(is_probability, pct_axes)),
        x = x * mult_x, pred = pred * mult_y, lo = lo * mult_y, hi = hi * mult_y
      )
    if (multi_level) interp <- interp |> mutate(level = factor(level, levels = level_choices))

    has_pi_any <- any(!is.na(d2$pi_lo95))
    use_pi     <- isTRUE(input$pf_interval_type == "pi") && has_pi_any
    d2 <- d2 |> mutate(band_lo = if (use_pi) pi_lo95 else lo95, band_hi = if (use_pi) pi_hi95 else hi95)

    p <- ggplot(d2, aes(x = x))

    if (is_pct_family) {
      p <- p + geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4)
    }
    if (has_pi_any) {
      if (use_pi) {
        p <- p + geom_line(aes(y = lo95), colour = "grey45", linetype = "dashed", na.rm = TRUE) +
          geom_line(aes(y = hi95), colour = "grey45", linetype = "dashed", na.rm = TRUE)
      } else {
        p <- p +
          geom_line(aes(y = pi_lo95), colour = "grey45", linetype = "dashed", na.rm = TRUE) +
          geom_line(aes(y = pi_hi95), colour = "grey45", linetype = "dashed", na.rm = TRUE)
      }
    }

    p <- p + geom_ribbon(aes(ymin = band_lo, ymax = band_hi), fill = col_band, alpha = 0.25)

    if (box_on) {
      o_box <- o |>
        mutate(cat_mid = bin_to_category_mid(x) * mult_x, y = y * 100) |>
        filter(!is.na(cat_mid)) |>
        left_join(flag_lookup, by = intersect(c("level", "to"), names(o))) |>
        mutate(facet_lab = paste0(to, pf_unit_tag(is_probability, is_pct_family)))
      if (multi_level) o_box <- o_box |> mutate(level = factor(level, levels = level_choices))
      p <- p +
        geom_boxplot(
          data = o_box, aes(x = cat_mid, y = y, group = cat_mid),
          width = 1, fill = "grey90", colour = "grey35", alpha = 0.8, outlier.shape = NA
        ) +
        geom_jitter(
          data = o_box, aes(x = cat_mid, y = y),
          width = 0.35, height = 0, alpha = 0.3, size = 0.8, shape = 16, colour = "grey20"
        )
    } else {
      if (isTRUE(input$pf_show_obs_se) && any(!is.na(o2$x_se) | !is.na(o2$y_se))) {
        p <- p +
          geom_errorbar(data = o2, aes(x = x, ymin = y - y_se, ymax = y + y_se), width = 0, colour = "grey55", alpha = 0.5, na.rm = TRUE) +
          geom_errorbar(data = o2, aes(y = y, xmin = x - x_se, xmax = x + x_se), orientation = "y", width = 0, colour = "grey55", alpha = 0.5, na.rm = TRUE)
      }
      p <- p + geom_point(data = o2, aes(x = x, y = y), shape = 1, colour = "grey30", alpha = 0.5)
    }

    p <- p +
      geom_line(aes(y = mean), colour = col_curve, linewidth = 1) +
      geom_vline(xintercept = input$pf_obs_value * mult_x, linetype = "dashed", colour = "grey40") +
      geom_point(data = interp, aes(x = x, y = pred), size = 3, colour = col_pt) +
      geom_errorbar(data = interp, aes(x = x, ymin = lo, ymax = hi), width = 0, colour = col_pt)

    p <- p + (if (multi_level) {
      facet_wrap(vars(level, facet_lab), ncol = 2, scales = "free")
    } else {
      facet_wrap(~ facet_lab, scales = "free", ncol = 2)
    })

    # x: the fixed source; xlab is one explicit string.
    xlab_display <- paste0(unique(d$xlab)[1], if (is_pct_family) " (%)" else "")

    # y: the varying target; name the possible quantities and point to the per-facet unit tag.
    ylab_display <- if (input$pf_family == "Coral cover") {
      "Predicted hard coral cover (%)"
    } else if (input$pf_family == "COTS cross-method") {
      "Predicted value (units in panel titles)"
    } else {
      "Predicted value (units in panel titles)"
    }

    p <- p +
      labs(
        x = xlab_display,
        y = ylab_display,
        title = paste("Predicting from:", input$pf_from),
        subtitle = if (multi_level) {
          "All spatial levels; each facet has its own independent x/y scale"
        } else if (input$pf_family != "eDNA") {
          paste0("Level: ", input$pf_level, "; each facet has its own independent x/y scale")
        } else {
          "Each facet has its own independent x/y scale"
        }
      ) +
      theme_csiro(13)

    if (is_pct_family) {
      p <- p + scale_x_continuous(labels = function(v) paste0(v, "%"))
    }

    p
  })

  ## ---- Secondary plot: RS Surface/Deep source, Manta tow target, boxplot toggle on (descriptive) ----
  output$pf_secondary_plot_ui <- renderUI({
    if (!pf_secondary_applicable()) return(NULL)
    card(
      full_screen = TRUE,
      card_header(paste0(
        "Manta tow as one of the PREDICTED targets \u2014 empirical distribution of ",
        input$pf_from, " by Manta's coral-cover category"
      )),
      plotOutput("pf_secondary_plot", height = "380px")
    )
  })

  output$pf_secondary_plot <- renderPlot({
    req(pf_secondary_applicable())
    o <- pf_current_obs() |> filter(to == "Manta tow")
    req(nrow(o) > 0)

    multi_level <- pf_multi_level_active()
    if (multi_level) o <- o |> mutate(level = factor(level, levels = level_choices))

    # Axes swapped vs the primary plot: Manta (o$y) is categorised, source (o$x) is boxed.
    o_box <- o |>
      mutate(cat_mid = bin_to_category_mid(y) * 100, from_pct = x * 100) |>
      filter(!is.na(cat_mid))
    req(nrow(o_box) > 0)

    p <- ggplot(o_box, aes(x = cat_mid, y = from_pct)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4) +
      geom_boxplot(aes(group = cat_mid), width = 1, fill = "grey90", colour = "grey35",
                   alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.35, height = 0, alpha = 0.3, size = 0.8, shape = 16, colour = "grey20") +
      geom_hline(yintercept = input$pf_obs_value * 100, linetype = "dashed", colour = "grey40") +
      labs(
        x = "Manta tow coral cover category (midpoint, %)",
        y = paste(input$pf_from, "hard coral cover (%)"),
        subtitle = paste0(
          "Dashed line = current slider value (", round(input$pf_obs_value * 100, 1), "%). ",
          if (multi_level) "All spatial levels; " else paste0("Level: ", input$pf_level, "; "),
          "boxplots are empirical (raw data), not the fitted model curve."
        )
      ) +
      scale_x_continuous(labels = function(v) paste0(v, "%")) +
      scale_y_continuous(labels = function(v) paste0(v, "%")) +
      theme_csiro(12)

    if (multi_level) p <- p + facet_wrap(~ level, ncol = 3)

    p
  })

  ## ---- Table ----
  output$pf_pred_table <- renderTable({
    pf_interp_df() |>
      transmute(
        `Level`              = if_else(is.na(level), "\u2014", level),
        `Predicts`           = to,
        `Predicted value`    = if_else(is_probability, paste0(signif(pred * 100, 3), "%"), as.character(signif(pred, 3))),
        `Lower`              = if_else(is_probability, paste0(signif(lo * 100, 3), "%"), as.character(signif(lo, 3))),
        `Upper`              = if_else(is_probability, paste0(signif(hi * 100, 3), "%"), as.character(signif(hi, 3))),
        `Interval shown`     = used_interval,
        `n (model fit)`      = n,
        `Slope (beta, model scale)` = paste0(
          signif(beta_mean, 2), " [", signif(beta_lo95, 2), ", ", signif(beta_hi95, 2), "]"
        )
      )
  }, striped = TRUE, bordered = TRUE)

  ## ---- Caveats (Predicting from tab) ----
  output$pf_caveat <- renderText({
    msgs <- character(0)
    d_cur <- tryCatch(pf_current_curves(), error = function(e) NULL)

    if (isTRUE(input$pf_family == "eDNA")) {
      msgs <- c(msgs, msg_edna(TRUE))
    }

    if (isTRUE(input$pf_family == "Coral cover") &&
        (pf_multi_level_active() || isTRUE(input$pf_level == "Cull-site")) &&
        !is.null(d_cur) && any(grepl("Deep", c(d_cur$from, d_cur$to)))) {
      msgs <- c(msgs, msg_deep_cull)
    }

    if (!is.null(d_cur) && nrow(d_cur) > 0 && any(d_cur$is_probability)) {
      msgs <- c(msgs, msg_prob)
    }

    if (!is.null(d_cur) && nrow(d_cur) > 0 && "pi_note" %in% names(d_cur) && any(!is.na(d_cur$pi_note))) {
      msgs <- c(msgs, unique(na.omit(d_cur$pi_note)))
    }

    interp <- tryCatch(pf_interp_df(), error = function(e) NULL)
    if (!is.null(interp) && input$pf_interval_type == "pi" && any(!interp$pi_available)) {
      msgs <- c(msgs, msg_pi_missing)
    }

    paste(msgs, collapse = " ")
  })

  ## ---- "Predicting from and to" tab: both ends fixed, one curve (per level if shown together) ----

  ft_multi_level_active <- reactive({
    isTRUE(input$ft_show_all_levels) && isTRUE(input$ft_family %in% c("Coral cover", "COTS cross-method"))
  })

  ft_box_toggle_active <- reactive({
    isTRUE(input$ft_box_manta) && isTRUE(input$ft_family == "Coral cover") && isTRUE(input$ft_from == "Manta tow")
  })

  ft_secondary_applicable <- reactive({
    isTRUE(input$ft_box_manta) && isTRUE(input$ft_family == "Coral cover") &&
      isTRUE(input$ft_from %in% c("ReefScan Surface", "ReefScan Deep")) &&
      isTRUE(input$ft_to == "Manta tow")
  })

  output$ft_from_ui <- renderUI({
    req(input$ft_family)
    d <- pred_curves |> filter(family_app == input$ft_family)
    if (input$ft_family != "eDNA" && !ft_multi_level_active()) {
      req(input$ft_level)
      d <- d |> filter(level == input$ft_level)
    }
    choices <- d |> pull(from) |> unique() |> sort()
    selectInput("ft_from", "Source method (from):", choices = choices)
  })

  output$ft_to_ui <- renderUI({
    req(input$ft_family, input$ft_from)
    d <- pred_curves |> filter(family_app == input$ft_family, from == input$ft_from)
    if (input$ft_family != "eDNA" && !ft_multi_level_active()) {
      req(input$ft_level)
      d <- d |> filter(level == input$ft_level)
    }
    choices <- d |> pull(to) |> unique() |> sort()
    selectInput("ft_to", "Target method (to):", choices = choices)
  })

  ft_current_curves <- reactive({
    req(input$ft_family, input$ft_from, input$ft_to)
    d <- pred_curves |> filter(family_app == input$ft_family, from == input$ft_from, to == input$ft_to)
    if (input$ft_family != "eDNA" && !ft_multi_level_active()) {
      req(input$ft_level)
      d <- d |> filter(level == input$ft_level)
    }
    d
  })

  ft_current_obs <- reactive({
    req(input$ft_family, input$ft_from, input$ft_to)
    d <- observations |> filter(family_app == input$ft_family, from == input$ft_from, to == input$ft_to)
    if (input$ft_family != "eDNA" && !ft_multi_level_active()) {
      req(input$ft_level)
      d <- d |> filter(level == input$ft_level)
    }
    d
  })

  ft_slider_range <- reactive({
    d <- ft_current_curves()
    req(nrow(d) > 0)
    range(d$x, na.rm = TRUE)
  })

  output$ft_slider_ui <- renderUI({
    rng  <- ft_slider_range()
    step <- signif(diff(rng) / 200, 3)
    if (!is.finite(step) || step <= 0) step <- 0.001
    sliderInput(
      "ft_obs_value", "Observed value:",
      min = signif(rng[1], 4), max = signif(rng[2], 4),
      value = signif(rng[1], 4), step = step
    )
  })

  output$ft_plot_ui <- renderUI({
    plotOutput("ft_pred_plot", height = if (ft_multi_level_active()) "620px" else "420px")
  })

  ft_interp_df <- reactive({
    d <- ft_current_curves()
    req(nrow(d) > 0, !is.null(input$ft_obs_value))
    use_pi <- input$ft_interval_type == "pi"

    d <- d |> mutate(level_key = if_else(is.na(level), "single", level))
    combos <- d |> distinct(level, level_key, pct_axes, is_probability)

    map_dfr(seq_len(nrow(combos)), function(i) {
      dd <- d |> filter(level_key == combos$level_key[i]) |> arrange(x)
      pi_avail <- !all(is.na(dd$pi_lo95))
      lo_col <- if (use_pi && pi_avail) dd$pi_lo95 else dd$lo95
      hi_col <- if (use_pi && pi_avail) dd$pi_hi95 else dd$hi95
      tibble(
        level = combos$level[i], x = input$ft_obs_value,
        n = dd$n[1], beta_mean = dd$beta_mean[1], beta_lo95 = dd$beta_lo95[1], beta_hi95 = dd$beta_hi95[1],
        pct_axes = combos$pct_axes[i], is_probability = combos$is_probability[i],
        pred = approx(dd$x, dd$mean, xout = input$ft_obs_value, rule = 2)$y,
        lo   = approx(dd$x, lo_col,  xout = input$ft_obs_value, rule = 2)$y,
        hi   = approx(dd$x, hi_col,  xout = input$ft_obs_value, rule = 2)$y,
        pi_available = pi_avail,
        used_interval = if (use_pi && pi_avail) "Prediction" else "Credible"
      )
    })
  })

  output$ft_pred_plot <- renderPlot({
    d <- ft_current_curves()
    o <- ft_current_obs()
    req(nrow(d) > 0)

    multi_level <- ft_multi_level_active()
    box_on      <- ft_box_toggle_active()

    if (multi_level) {
      d <- d |> mutate(level = factor(level, levels = level_choices))
      o <- o |> mutate(level = factor(level, levels = level_choices))
    }

    is_pct  <- isTRUE(unique(d$pct_axes)[1])
    is_prob <- isTRUE(unique(d$is_probability)[1])
    mult_x  <- if (is_pct) 100 else 1
    mult_y  <- if (is_pct || is_prob) 100 else 1

    d2 <- d |> mutate(x = x * mult_x, mean = mean * mult_y, lo95 = lo95 * mult_y, hi95 = hi95 * mult_y,
                       pi_lo95 = pi_lo95 * mult_y, pi_hi95 = pi_hi95 * mult_y)
    o2 <- o |> mutate(x = x * mult_x, y = y * mult_y,
                       x_se = if ("x_se" %in% names(o)) x_se * mult_x else NA_real_,
                       y_se = if ("y_se" %in% names(o)) y_se * mult_y else NA_real_)
    interp <- ft_interp_df() |> mutate(x = x * mult_x, pred = pred * mult_y, lo = lo * mult_y, hi = hi * mult_y)
    if (multi_level) interp <- interp |> mutate(level = factor(level, levels = level_choices))

    has_pi_any <- any(!is.na(d2$pi_lo95))
    use_pi     <- isTRUE(input$ft_interval_type == "pi") && has_pi_any
    d2 <- d2 |> mutate(band_lo = if (use_pi) pi_lo95 else lo95, band_hi = if (use_pi) pi_hi95 else hi95)

    p <- ggplot(d2, aes(x = x))

    if (is_pct) {
      p <- p + geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4)
    }
    if (has_pi_any) {
      if (use_pi) {
        p <- p + geom_line(aes(y = lo95), colour = "grey45", linetype = "dashed", na.rm = TRUE) +
          geom_line(aes(y = hi95), colour = "grey45", linetype = "dashed", na.rm = TRUE)
      } else {
        p <- p +
          geom_line(aes(y = pi_lo95), colour = "grey45", linetype = "dashed", na.rm = TRUE) +
          geom_line(aes(y = pi_hi95), colour = "grey45", linetype = "dashed", na.rm = TRUE)
      }
    }

    p <- p + geom_ribbon(aes(ymin = band_lo, ymax = band_hi), fill = col_band, alpha = 0.25)

    if (box_on) {
      o_box <- o |>
        mutate(cat_mid = bin_to_category_mid(x) * mult_x, y = y * mult_y) |>
        filter(!is.na(cat_mid))
      if (multi_level) o_box <- o_box |> mutate(level = factor(level, levels = level_choices))
      p <- p +
        geom_boxplot(
          data = o_box, aes(x = cat_mid, y = y, group = cat_mid),
          width = 1, fill = "grey90", colour = "grey35", alpha = 0.8, outlier.shape = NA
        ) +
        geom_jitter(
          data = o_box, aes(x = cat_mid, y = y),
          width = 0.35, height = 0, alpha = 0.3, size = 0.8, shape = 16, colour = "grey20"
        )
    } else {
      if (isTRUE(input$ft_show_obs_se) && any(!is.na(o2$x_se) | !is.na(o2$y_se))) {
        p <- p +
          geom_errorbar(data = o2, aes(x = x, ymin = y - y_se, ymax = y + y_se), width = 0, colour = "grey55", alpha = 0.5, na.rm = TRUE) +
          geom_errorbar(data = o2, aes(y = y, xmin = x - x_se, xmax = x + x_se), orientation = "y", width = 0, colour = "grey55", alpha = 0.5, na.rm = TRUE)
      }
      p <- p + geom_point(data = o2, aes(x = x, y = y), shape = 1, colour = "grey30", alpha = 0.5)
    }

    p <- p +
      geom_line(aes(y = mean), colour = col_curve, linewidth = 1) +
      geom_vline(xintercept = unique(interp$x), linetype = "dashed", colour = "grey40") +
      geom_point(data = interp, aes(x = x, y = pred), size = 3, colour = col_pt) +
      geom_errorbar(data = interp, aes(x = x, ymin = lo, ymax = hi), width = 0, colour = col_pt)

    if (multi_level) p <- p + facet_wrap(~ level, ncol = 3, scales = "free_x")

    # One from/to pairing is shown, so xlab/ylab are single strings.
    xlab_display <- paste0(unique(d$xlab)[1], if (is_pct) " (%)" else "")
    ylab_display <- paste0(unique(d$ylab)[1], if (is_pct || is_prob) " (%)" else "")

    p <- p +
      labs(
        x = xlab_display,
        y = ylab_display,
        title = paste0(input$ft_from, " \u2192 ", input$ft_to),
        subtitle = if (multi_level) {
          "All spatial levels"
        } else if (input$ft_family != "eDNA") {
          paste("Level:", input$ft_level)
        } else NULL
      ) +
      theme_csiro(13)

    if (is_pct) {
      p <- p +
        scale_x_continuous(labels = function(v) paste0(v, "%")) +
        scale_y_continuous(labels = function(v) paste0(v, "%"))
    } else if (is_prob) {
      p <- p + scale_y_continuous(labels = function(v) paste0(v, "%"))
    }

    p
  })

  ## ---- Secondary plot: RS Surface/Deep source, Manta tow target, boxplot toggle on (descriptive) ----
  output$ft_secondary_plot_ui <- renderUI({
    if (!ft_secondary_applicable()) return(NULL)
    card(
      full_screen = TRUE,
      card_header(paste0(
        "Manta tow as the PREDICTED method \u2014 empirical distribution of ",
        input$ft_from, " by Manta's coral-cover category"
      )),
      plotOutput("ft_secondary_plot", height = "380px")
    )
  })

  output$ft_secondary_plot <- renderPlot({
    req(ft_secondary_applicable())
    o <- ft_current_obs()
    req(nrow(o) > 0)

    multi_level <- ft_multi_level_active()
    if (multi_level) o <- o |> mutate(level = factor(level, levels = level_choices))

    o_box <- o |>
      mutate(cat_mid = bin_to_category_mid(y) * 100, from_pct = x * 100) |>
      filter(!is.na(cat_mid))
    req(nrow(o_box) > 0)

    p <- ggplot(o_box, aes(x = cat_mid, y = from_pct)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4) +
      geom_boxplot(aes(group = cat_mid), width = 1, fill = "grey90", colour = "grey35",
                   alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.35, height = 0, alpha = 0.3, size = 0.8, shape = 16, colour = "grey20") +
      geom_hline(yintercept = input$ft_obs_value * 100, linetype = "dashed", colour = "grey40") +
      labs(
        x = "Manta tow coral cover category (midpoint, %)",
        y = paste(input$ft_from, "hard coral cover (%)"),
        subtitle = paste0(
          "Dashed line = current slider value (", round(input$ft_obs_value * 100, 1), "%). ",
          if (multi_level) "All spatial levels; " else paste0("Level: ", input$ft_level, "; "),
          "boxplots are empirical (raw data), not the fitted model curve."
        )
      ) +
      scale_x_continuous(labels = function(v) paste0(v, "%")) +
      scale_y_continuous(labels = function(v) paste0(v, "%")) +
      theme_csiro(12)

    if (multi_level) p <- p + facet_wrap(~ level, ncol = 3)

    p
  })

  ## ---- Table ----
  output$ft_pred_table <- renderTable({
    ft_interp_df() |>
      transmute(
        `Level`              = if_else(is.na(level), "\u2014", level),
        `Predicted value`    = if_else(is_probability, paste0(signif(pred * 100, 3), "%"), as.character(signif(pred, 3))),
        `Lower`              = if_else(is_probability, paste0(signif(lo * 100, 3), "%"), as.character(signif(lo, 3))),
        `Upper`              = if_else(is_probability, paste0(signif(hi * 100, 3), "%"), as.character(signif(hi, 3))),
        `Interval shown`     = used_interval,
        `n (model fit)`      = n,
        `Slope (beta, model scale)` = paste0(
          signif(beta_mean, 2), " [", signif(beta_lo95, 2), ", ", signif(beta_hi95, 2), "]"
        )
      )
  }, striped = TRUE, bordered = TRUE)

  ## ---- Caveats (Predicting from and to tab) ----
  output$ft_caveat <- renderText({
    msgs <- character(0)
    d_cur <- tryCatch(ft_current_curves(), error = function(e) NULL)

    if (isTRUE(input$ft_family == "eDNA")) {
      msgs <- c(msgs, msg_edna(FALSE))
    }

    if (isTRUE(input$ft_family == "Coral cover") &&
        (ft_multi_level_active() || isTRUE(input$ft_level == "Cull-site")) &&
        (grepl("Deep", input$ft_from %||% "") || grepl("Deep", input$ft_to %||% ""))) {
      msgs <- c(msgs, msg_deep_cull)
    }

    if (!is.null(d_cur) && nrow(d_cur) > 0 && isTRUE(unique(d_cur$is_probability)[1])) {
      msgs <- c(msgs, msg_prob)
    }

    if (!is.null(d_cur) && nrow(d_cur) > 0 && "pi_note" %in% names(d_cur) && any(!is.na(d_cur$pi_note))) {
      msgs <- c(msgs, unique(na.omit(d_cur$pi_note)))
    }

    interp <- tryCatch(ft_interp_df(), error = function(e) NULL)
    if (!is.null(interp) && input$ft_interval_type == "pi" && any(!interp$pi_available)) {
      msgs <- c(msgs, msg_pi_missing)
    }

    paste(msgs, collapse = " ")
  })

  ## ---- Descriptive uncertainty-comparison plot: reef level, all 4 methods, mean +/- 1 SD
  ## across sub-reef sampling units (optional reef-substring filter) ----
  uc_log1p_active <- reactive({
    isTRUE(input$uc_cots_log1p) && isTRUE(input$uc_quantity == "COTS count (per tow or survey)")
  })

  uc_grouping <- reactive({ input$uc_grouping %||% "reef" })
  uc_pooled_active <- reactive({ uc_grouping() == "all" })
  uc_band_active   <- reactive({ uc_grouping() == "band" })

  uc_filtered <- reactive({
    if (uc_pooled_active()) {
      src <- if (uc_log1p_active()) uncertainty_pooled_cots_log1p else uncertainty_pooled
      return(src |> filter(quantity == input$uc_quantity))
    }
    if (uc_band_active()) {
      src <- if (uc_log1p_active()) uncertainty_band_cots_log1p else uncertainty_band_comparison
      d <- src |> filter(quantity == input$uc_quantity)
      if (nzchar(input$uc_band_filter %||% "")) {
        d <- d |> filter(grepl(input$uc_band_filter, band_label, ignore.case = TRUE))
      }
      return(d)
    }
    src <- if (uc_log1p_active()) uncertainty_cots_log1p else uncertainty_comparison
    d <- src |> filter(quantity == input$uc_quantity)
    if (nzchar(input$uc_reef_filter %||% "")) {
      d <- d |> filter(grepl(input$uc_reef_filter, reef_label, ignore.case = TRUE))
    }
    d
  })

  output$uc_card_header <- renderText({
    if (uc_pooled_active()) {
      "Method comparison \u2013 all reefs pooled"
    } else if (uc_band_active()) {
      "Method comparison \u2013 by reef-ID band"
    } else {
      "Method comparison \u2013 by reef"
    }
  })

  output$uc_point_desc <- renderUI({
    tags$small(
      "Point = mean of sub-reef sampling units (tows, ~100 m grid cells or surveys); bars = \u00b11 SD",
      " of those units (spread within reefs, not the standard errors on the calibration tabs).",
      if (uc_pooled_active() || uc_band_active()) {
        " Units from different reefs are combined directly, so reefs with more units weigh more."
      }
    )
  })

  uc_by_tool <- reactive({
    !uc_pooled_active() && identical(input$uc_arrange, "tool")
  })

  # "Group by tool": one panel per method, units (reefs/bands) listed within it
  uc_tool_layout <- function(p) {
    p +
      facet_grid(rows = vars(method), scales = "free_y", space = "free_y") +
      guides(colour = "none") +
      theme(strip.text.y = element_text(angle = 0, hjust = 0))
  }

  output$uc_plot_ui <- renderUI({
    if (uc_pooled_active()) {
      return(plotOutput("uc_plot", height = "320px"))
    }
    if (uc_by_tool()) {
      h <- max(400, min(6000, 20 * nrow(uc_filtered()) + 140))
      return(plotOutput("uc_plot", height = paste0(h, "px")))
    }
    if (uc_band_active()) {
      n_band <- n_distinct(uc_filtered()$band_label)
      h <- max(300, min(2000, 28 * n_band))
      return(plotOutput("uc_plot", height = paste0(h, "px")))
    }
    n_reef <- n_distinct(uc_filtered()$reef_label)
    h <- max(400, min(4000, 28 * n_reef))
    plotOutput("uc_plot", height = paste0(h, "px"))
  })

  output$uc_plot <- renderPlot({
    d <- uc_filtered()
    req(nrow(d) > 0)
    log1p_mode <- uc_log1p_active()

    if (uc_pooled_active()) {
      method_order <- sort(unique(d$method), decreasing = TRUE)
      d <- d |> mutate(method = factor(method, levels = method_order))

      p <- ggplot(d, aes(y = method, x = mean, xmin = mean - sd, xmax = mean + sd, color = method)) +
        geom_pointrange(linewidth = 0.7, size = 0.55, show.legend = FALSE) +
        geom_text(aes(label = paste0("n=", n)), hjust = -0.2, vjust = -1.1, size = 3.6, colour = "grey30") +
        labs(
          y = "Method",
          x = if (log1p_mode) "COTS observed (count per tow or survey, log1p scale)" else unique(d$quantity)[1],
          subtitle = if (log1p_mode) {
            "Mean \u00b11 SD of log1p(count) per tow/survey; axis in count units."
          } else {
            "n = sampling units pooled"
          }
        ) +
        scale_colour_manual(values = method_cols) +
        theme_csiro(13)

      if (log1p_mode) {
        p <- p + scale_x_continuous(breaks = log1p(cots_log1p_breaks), labels = cots_log1p_breaks)
      }
      return(wrap_subtitle(p))
    }

    if (uc_band_active()) {
      band_order <- sort(unique(d$band_label), decreasing = TRUE)
      d <- d |> mutate(band_label = factor(band_label, levels = band_order))

      p <- ggplot(d, aes(y = band_label, x = mean, xmin = mean - sd, xmax = mean + sd, color = method)) +
        geom_pointrange(position = position_dodge(width = 0.7), linewidth = 0.5, size = 0.25) +
        labs(
          y = "Reef-ID band (e.g. \"21-XXX\")",
          x = if (log1p_mode) "COTS observed (count per tow or survey, log1p scale)" else unique(d$quantity)[1],
          color = "Method",
          subtitle = if (log1p_mode) {
            "Mean \u00b11 SD of log1p(count) per tow/survey; axis in count units."
          }
        ) +
        scale_colour_manual(values = method_cols) +
        theme_csiro(12)

      if (log1p_mode) {
        p <- p + scale_x_continuous(breaks = log1p(cots_log1p_breaks), labels = cots_log1p_breaks)
      }
      if (uc_by_tool()) p <- uc_tool_layout(p)
      return(wrap_subtitle(p))
    }

    reef_order <- sort(unique(d$reef_label), decreasing = TRUE)
    d <- d |> mutate(reef_label = factor(reef_label, levels = reef_order))

    p <- ggplot(d, aes(y = reef_label, x = mean, xmin = mean - sd, xmax = mean + sd, color = method)) +
      geom_pointrange(position = position_dodge(width = 0.7), linewidth = 0.4, size = 0.15) +
      labs(
        y = "Reef (ordered by reef ID)",
        x = if (log1p_mode) "COTS observed (count per tow or survey, log1p scale)" else unique(d$quantity)[1],
        color = "Method",
        subtitle = if (log1p_mode) {
          "Mean \u00b11 SD of log1p(count) per tow/survey; axis in count units."
        }
      ) +
      scale_colour_manual(values = method_cols) +
      theme_csiro(12)

    if (log1p_mode) {
      p <- p + scale_x_continuous(breaks = log1p(cots_log1p_breaks), labels = cots_log1p_breaks)
    }

    if (uc_by_tool()) p <- uc_tool_layout(p)
    wrap_subtitle(p)
  })

  ## ---- Table ----
  output$pred_table <- renderTable({
    interp_df() |>
      transmute(
        `Level`              = if_else(is.na(level), "\u2014", level),
        `Source method`      = from,
        `Predicted value`    = if_else(is_probability, paste0(signif(pred * 100, 3), "%"), as.character(signif(pred, 3))),
        `Lower`              = if_else(is_probability, paste0(signif(lo * 100, 3), "%"), as.character(signif(lo, 3))),
        `Upper`              = if_else(is_probability, paste0(signif(hi * 100, 3), "%"), as.character(signif(hi, 3))),
        `Interval shown`     = used_interval,
        `n (model fit)`      = n,
        `Slope (beta, model scale)` = paste0(
          signif(beta_mean, 2), " [", signif(beta_lo95, 2), ", ", signif(beta_hi95, 2), "]"
        )
      )
  }, striped = TRUE, bordered = TRUE)

  ## ---- Caveats (Calibration curves tab) ----
  output$caveat <- renderText({
    msgs <- character(0)
    d_cur <- tryCatch(current_curves(), error = function(e) NULL)

    if (isTRUE(input$family == "eDNA")) {
      msgs <- c(msgs, msg_edna(TRUE))
    }

    if (isTRUE(input$family == "Coral cover") &&
        (multi_level_active() || isTRUE(input$level == "Cull-site")) &&
        !is.null(d_cur) && any(grepl("Deep", c(d_cur$from, d_cur$to)))) {
      msgs <- c(msgs, msg_deep_cull)
    }

    if (!is.null(d_cur) && nrow(d_cur) > 0 && isTRUE(unique(d_cur$is_probability)[1])) {
      msgs <- c(msgs, msg_prob)
    }

    if (!is.null(d_cur) && nrow(d_cur) > 0 && "pi_note" %in% names(d_cur) && any(!is.na(d_cur$pi_note))) {
      msgs <- c(msgs, unique(na.omit(d_cur$pi_note)))
    }

    interp <- tryCatch(interp_df(), error = function(e) NULL)
    if (!is.null(interp) && input$interval_type == "pi" && any(!interp$pi_available)) {
      msgs <- c(msgs, msg_pi_missing)
    }

    paste(msgs, collapse = " ")
  })

  ## ---- Caveat (Uncertainty comparison tab) ----
  output$uc_caveat <- renderText({
    msgs <- c(paste(
      "Descriptive only: RHIS is not part of the pairwise calibration models",
      "(it has no tow-level grain), so it appears only on this tab. RHIS coral cover is hard coral only."
    ))

    if (uc_log1p_active()) {
      msgs <- c(msgs, paste(
        "COTS counts are right-skewed, so mean \u2212 SD can fall below zero on the log1p scale",
        "(log1p is applied per tow/survey before averaging)."
      ))
    }

    paste(msgs, collapse = " ")
  })
}

shinyApp(ui, server)
