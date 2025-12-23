# Prairie Dog Decision Tool v9.7
# - glmmTMB model with random effects (Site)
# - multi-year iteration capability with NEW TAB opening
# - economic calculator
# - summary plots
# - summary statistics
# - updated user interface
#         streamlined input uploads
#         clearer map layout
#         Pdog control and plague treatment polygon drawing capability
#           - note: neither Pdog control nor plague control carry over between years (this can be modified)
#         model runs on command rather than automatically with change of input variable
#
#         
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 1) Load Libraries -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

library(shiny)
library(shinythemes)
library(shinyBS)
library(shinyjs)
library(later)
library(sp)
library(raster)
library(landscapemetrics)
library(leaflet)
library(leaflet.extras)
library(viridis)
library(sf)
library(zip)
library(mapview)
library(terra)
library(plotrix)
library(lwgeom)
library(FNN)
library(glmmTMB)
library(lme4)
library(tidyverse)
library(digest)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 2) Check for required files and load data -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

required_files <- c("climate_complete_matrix_utm13.rds", "growth_model.rds", "plague_model.rds")
missing_files <- required_files[!file.exists(required_files)]

if(length(missing_files) > 0) {
  cat("ERROR: Missing required files:\n")
  for(file in missing_files) {
    cat("  -", file, "\n")
  }
  cat("\nPlease ensure these files are in your working directory:\n")
  cat("Current working directory:", getwd(), "\n")
  stop("Cannot proceed without required files")
} else {
  cat("All required files found!\n")
}

cat("Loading climate matrix...\n")
climate_matrix <- readRDS("climate_complete_matrix_utm13.rds")
cat("Climate matrix loaded:", nrow(climate_matrix), "rows,", ncol(climate_matrix), "columns\n")

cat("Loading glmmTMB models...\n")
growth_model <- readRDS("growth_model.rds")
plague_model <- readRDS("plague_model.rds")
cat("glmmTMB models loaded successfully!\n")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 2b) Multi-Year Iteration Parameters -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Insecticide persistence configuration
# Current setting: 1 year (no carryover between iterations)
# Future: Can be changed to 2, 3, etc. for multi-year persistence
INSECTICIDE_PERSISTENCE_YEARS <- 1

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 3) Helper Functions -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Extract climate data from matrix using nearest neighbor
extract_climate_from_matrix <- function(coords_utm, climate_matrix) {
  climate_coords <- as.matrix(climate_matrix[, c("x", "y")])
  valid_climate <- complete.cases(climate_coords)
  climate_coords <- climate_coords[valid_climate, ]
  climate_matrix_clean <- climate_matrix[valid_climate, ]
  
  nn_indices <- get.knnx(climate_coords, coords_utm, k = 1)$nn.index
  
  extracted_data <- as.data.frame(climate_matrix_clean[nn_indices, ])
  names(extracted_data) <- names(climate_matrix_clean)
  
  return(extracted_data)
}

# Extract climate data for sf polygons
extract_climate_data <- function(shapefile_sf, climate_matrix, buffer_distance = 100) {
  if(is.null(shapefile_sf) || nrow(shapefile_sf) == 0) {
    stop("Shapefile is empty or NULL")
  }
  
  if(is.null(climate_matrix) || nrow(climate_matrix) == 0) {
    stop("Climate matrix is empty or NULL")
  }
  
  target_crs <- "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs"
  if(is.na(st_crs(shapefile_sf)) || st_crs(shapefile_sf)$input != target_crs) {
    shapefile_utm <- st_transform(shapefile_sf, crs = target_crs)
  } else {
    shapefile_utm <- shapefile_sf
  }
  
  centroids <- st_centroid(shapefile_utm)
  centroid_coords <- st_coordinates(centroids)
  
  climate_coords <- as.matrix(climate_matrix[, c("x", "y")])
  valid_climate <- complete.cases(climate_coords)
  climate_coords <- climate_coords[valid_climate, ]
  climate_matrix_clean <- climate_matrix[valid_climate, ]
  
  nn_indices <- get.knnx(climate_coords, centroid_coords, k = 1)$nn.index
  extracted_data <- climate_matrix_clean[nn_indices, ]
  extracted_data$poly_id <- seq_len(nrow(centroid_coords))
  
  return(extracted_data)
}

# Create raster template from spatial extent
create_raster_template <- function(extent_sf, resolution = 100, buffer = 10000) {
  extent_utm <- st_transform(extent_sf, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
  extent_bbox <- st_bbox(extent_utm)
  
  r <- raster()
  extent(r) <- c(extent_bbox[1] - buffer, extent_bbox[3] + buffer, 
                 extent_bbox[2] - buffer, extent_bbox[4] + buffer)
  res(r) <- c(resolution, resolution)
  crs(r) <- "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs"
  
  return(r)
}

calculate_colony_size <- function(colony_raster) {
  temp_patch <- clump(colony_raster, directions = 8, gaps = FALSE)
  cells_patch <- freq(temp_patch)
  cells_patch <- as.data.frame(cells_patch[1:(nrow(cells_patch)-1), ])
  temp_psize <- subs(temp_patch, cells_patch, by = 1, which = 2)
  temp_psize[is.na(temp_psize)] <- 0
  return(temp_psize)
}

detect_plague_events <- function(colony_raster_t_minus_1, colony_raster_t) {
  raster1 <- rast(colony_raster_t_minus_1)
  raster2 <- rast(colony_raster_t)
  clumps1 <- patches(raster1, directions = 8, zeroAsNA = TRUE)
  cluster_sizes1 <- as.data.frame(freq(clumps1))
  names(cluster_sizes1) <- c("layer", "cluster_id", "size1")
  cluster_sizes1$size2 <- 0
  
  for (j in cluster_sizes1$cluster_id) {
    cluster_mask <- clumps1 == j
    remaining_pixels <- sum(raster2[cluster_mask] == 1, na.rm = TRUE)
    cluster_sizes1$size2[cluster_sizes1$cluster_id == j] <- remaining_pixels
  }
  
  cluster_sizes1$decline <- (cluster_sizes1$size1 - cluster_sizes1$size2) / cluster_sizes1$size1
  cluster_sizes1$plague <- ifelse(cluster_sizes1$decline >= 0.5, -1, 3)
  
  reclass_matrix <- as.matrix(data.frame(
    cluster_id = cluster_sizes1$cluster_id,
    plague_status = cluster_sizes1$plague
  ))
  
  plague_raster <- classify(clumps1, rcl = reclass_matrix, right = FALSE)
  plague_raster <- raster(plague_raster)
  
  return(plague_raster)
}

calculate_distance_to_plague <- function(colony_raster_current, colony_raster_previous) {
  plague_raster <- detect_plague_events(colony_raster_previous, colony_raster_current)
  chng <- colony_raster_current - colony_raster_previous
  values(chng) <- ifelse(values(chng) > 0, 2, ifelse(values(chng) < 0, -1, 0))
  maint <- colony_raster_current * colony_raster_previous
  status <- chng + maint
  
  status_new <- overlay(status, plague_raster, fun = function(x, y) {
    ifelse(x == -1 & y != -1, 3, x)
  })
  
  distnoplg <- status_new
  values(distnoplg) <- 50000
  y_r <- status_new
  values(y_r) <- ifelse(values(y_r) < 0, 1, NA)
  
  if(any(!is.na(values(y_r)))) {
    y_r2 <- distance(y_r == 1)
  } else {
    y_r2 <- distnoplg
  }
  
  return(y_r2)
}

# Calculate uncertainty using Delta Method
calculate_uncertainty_delta <- function(model, pred_data, predictions_link) {
  tryCatch({
    vcov_mat <- vcov(model)$cond
    fixed_formula <- lme4::nobars(formula(model))
    X <- model.matrix(fixed_formula[-2], data = pred_data)
    
    var_link <- diag(X %*% vcov_mat %*% t(X))
    se_link <- sqrt(var_link)
    
    pred_prob <- plogis(predictions_link)
    se_prob <- se_link * pred_prob * (1 - pred_prob)
    
    ci_lower <- plogis(predictions_link - 1.96 * se_link)
    ci_upper <- plogis(predictions_link + 1.96 * se_link)
    
    return(list(
      se = se_prob,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    ))
  }, error = function(e) {
    pred_prob <- plogis(predictions_link)
    se_prob <- pred_prob * 0.15
    list(
      se = se_prob,
      ci_lower = pmax(0, pred_prob - 1.96 * se_prob),
      ci_upper = pmin(1, pred_prob + 1.96 * se_prob)
    )
  })
}

# Economic calculator function
calculate_economic_scenario <- function(acres_poisoned, acres_pdogs_after, ppa_value, cost_per_acre, scenario_type) {
  # Fixed parameters (under the hood)
  PDOG_DENSITY_PER_HA <- 42
  EFFICACY_POISONING <- 100  # percent
  COWS_PER_PDOG <- 1/335  # prairie dogs per cow equivalent
  POUNDS_PER_COW_PER_DAY <- 0.025 * 1000  # 25 lbs
  DAYS_PER_YEAR <- 365
  
  # Calculate based on scenario type
  if (scenario_type == "consumption") {
    # Scenario A: Based on pdog consumption estimates
    
    # Forage SAVED (from poisoning)
    total_pdogs_poisoned <- PDOG_DENSITY_PER_HA * acres_poisoned
    cows_equivalent_saved <- total_pdogs_poisoned * COWS_PER_PDOG
    pounds_per_cow_per_year <- POUNDS_PER_COW_PER_DAY * DAYS_PER_YEAR
    forage_saved <- pounds_per_cow_per_year * cows_equivalent_saved
    
    # Forage LOST (from remaining pdogs)
    total_pdogs_remaining <- PDOG_DENSITY_PER_HA * acres_pdogs_after
    cows_equivalent_lost <- total_pdogs_remaining * COWS_PER_PDOG
    forage_lost <- pounds_per_cow_per_year * cows_equivalent_lost
    
  } else if (scenario_type == "no_forage") {
    # Scenario B: Pdog colonies have 0% forage value
    forage_saved <- acres_poisoned * ppa_value
    forage_lost <- acres_pdogs_after * ppa_value
    
  } else if (scenario_type == "half_forage") {
    # Scenario C: Pdog colonies have 50% forage value
    forage_saved <- acres_poisoned * ppa_value * 0.5
    forage_lost <- acres_pdogs_after * ppa_value * 0.5
    
  } else if (scenario_type == "no_increase") {
    # Scenario D: No increase in forage
    forage_saved <- 0
    forage_lost <- 0
  }
  
  # Calculate total cost
  total_cost <- cost_per_acre * acres_poisoned
  
  return(list(
    forage_saved = forage_saved,
    forage_lost = forage_lost,
    total_cost = total_cost
  ))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 4) User Interface -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

ui <- navbarPage("Prairie Dog Modeling Application", 
                 selected = "Predict", 
                 collapsible = TRUE, 
                 inverse = TRUE, 
                 theme = shinytheme("sandstone"),
                 
                 useShinyjs(),
                 
                 tags$head(
                   tags$style(HTML("
    #costBenefitBox.disabled {
      opacity: 0.5;
      pointer-events: none;
    }
    .scenario-box {
      border: 2px solid #3474A7;
      border-radius: 8px;
      padding: 10px;
      margin: 5px;
      background-color: #f8f9fa;
    }
    .scenario-title {
      color: #3474A7;
      font-weight: bold;
      text-align: center;
      margin-bottom: 10px;
      font-size: 13px;
    }
    .output-value {
      background-color: #d4edda;
      padding: 5px;
      border-radius: 4px;
      margin: 5px 0;
      font-weight: bold;
      font-size: 12px;
    }
    .economic-input {
      margin-bottom: 8px;
    }
    .economic-input label {
      font-size: 11px;
      margin-bottom: 2px;
    }
    .economic-input input {
      height: 28px;
      font-size: 11px;
    }
    /* Light yellow background for acres inputs */
    #acres_poisoned_A, #acres_pdogs_A,
    #acres_poisoned_B, #acres_pdogs_B,
    #acres_poisoned_C, #acres_pdogs_C,
    #acres_poisoned_D, #acres_pdogs_D {
      background-color: #fffacd !important;
    }
    /* Arrow connecting scenario table to economic calculator */
    .linking-arrow {
      position: absolute;
      left: -25px;
      top: 95px
      width: 0;
      height: 0;
      border-top: 40px solid transparent;
      border-bottom: 40px solid transparent;
      border-right: 20px solid #3474A7;
      z-index: 10;
    }
    
    /* Popover help icon styles */
    .help-icon {
      display: inline-block;
      width: 16px;
      height: 16px;
      line-height: 16px;
      text-align: center;
      background-color: #5A9BD4;
      color: white;
      border-radius: 50%;
      font-size: 11px;
      font-weight: bold;
      cursor: help;
      margin-left: 5px;
      position: relative;
      top: -1px;
    }
    
    .help-icon:hover {
      background-color: #4A8BC4;
    }
    
    /* Popover container */
    .popover-container {
      position: relative;
      display: inline-block;
    }
    
    /* Popover text */
    .popover-text {
      visibility: hidden;
      width: 300px;
      background-color: #333;
      color: #fff;
      text-align: left;
      border-radius: 6px;
      padding: 10px;
      position: absolute;
      z-index: 1000;
      left: 25px;
      top: -10px;
      opacity: 0;
      transition: opacity 0.3s;
      font-size: 12px;
      line-height: 1.4;
      box-shadow: 0 2px 8px rgba(0,0,0,0.3);
    }
    
    /* Popover arrow */
    .popover-text::before {
      content: '';
      position: absolute;
      top: 15px;
      left: -5px;
      border-width: 5px;
      border-style: solid;
      border-color: transparent #333 transparent transparent;
    }
    
    /* Show popover on hover */
    .popover-container:hover .popover-text {
      visibility: visible;
      opacity: 1;
    }
  "))
                 ),
                 ##############                 
                 # OVERVIEW TAB
                 ##############                 
                 tabPanel("Overview",
                          fluidPage(theme = shinytheme("sandstone"),
                                    titlePanel(p("PDOG MAPR - Prairie Dog Management and Planning Resource", style = "color:#3474A7"))),
                          mainPanel(
                            h3("Brief overview:", style = "color:#3474A7"),
                            h4(tags$b("Purpose:"), " PDOG MAPR is a decision support tool to help inform management of prairie dog ecosystems. This tool is intended to be the first step in an iterative process informed by conservation partners and managers that results in increasingly sophisticated decision support for prairie dog management."),
                            br(),
                            h4(tags$b("Platform:"), " PDOG MAPR is publicly available here on shinyapps.io, a service platform for hosting Shiny web apps. Users can upload prairie dog colony shapefile data anywhere within the range of black-tailed prairie dogs and freely use PDOG MAPR to simulate colony dynamics under various weather and management scenarios by visiting ", 
                               tags$a(href = "https://moped.shinyapps.io/PDOGMAPR/", "https://moped.shinyapps.io/PDOGMAPR/", target = "_blank"), "."),
                            br(),
                            h4(tags$b("User tutorial:"), " A video tutorial describing use of PDOG MAPR in greater detail can be found at ", 
                               tags$a(href = "https://youtu.be/jnRZ6dRA6mg", "https://youtu.be/jnRZ6dRA6mg", target = "_blank"), 
                               ". We highly recommend watching this brief video prior to using the application. See below for more detailed written instructions."),
                            br(),
                            h3("Questions that PDOG MAPR can versus cannot help users address:", style = "color:#3474A7"),
                            h4(tags$b("Questions PDOG MAPR can help answer at the user's grassland/site of interest:")),
                            tags$ul(
                              tags$li("What is the overall chance of plague outbreak on prairie dog colonies in the following year?"),
                              tags$li("How much colony expansion is expected in the following year?"),
                              tags$li("What is the predicted size and distribution of the colony complex in the following year?"),
                              tags$li("How does climate (temperature and precipitation) influence colony expansion and the chance of plague in the following year?"),
                              tags$li("How does plague mitigation and/or lethal control influence colony expansion and the chance of plague in the following year?")
                            ),
                            br(),
                            h4(tags$b("Application Example 1; forecasting colony dynamics under no management:"), " Users will typically upload a shapefile of prairie dog colonies that were mapped in the current year (e.g., 2025) and use the tool to obtain predictions regarding colony dynamics in the following year (e.g., 2026). Thus, users often will not know what the weather (precipitation and temperature) will be like in the following year, but using PDOG MAPR will be able to simulate what the colony complex is likely to do (colony sizes, distribution) in the following year under different weather scenarios. For instance, users can vary upcoming spring weather conditions in the application to produce scenarios that range from (1) lowest colony growth and highest chance of plague to (2) highest colony growth and lowest chance of plague. This exercise could provide managers with an idea of the most-likely range of possibilities for colony distribution across the range of potential weather scenarios."),
                            br(),
                            h4(tags$b("Application Example 2; assess impacts of management scenarios on colony complex:"), " Users can also upload shapefiles of management actions to simulate how different management tools (and the extent and distribution of those tools) can be applied across the colony complex to achieve management goals of interest, such as how to best mitigate and lessen the likelihood of a plague outbreak."),
                            br(),
                            h4(tags$b("Application Example 3; explore past colony dynamics to inform future management decisions:"), " Users can learn about their system(s) by inputting past colony and climate data. They can then vary management tools to determine how such applications might have influenced colony growth and/or the chance of plague outbreak."),
                            br(),
                            h4(tags$b("Questions PDOG MAPR cannot help answer:")),
                            h4("PDOGMAPR was designed to provide a broad picture of when large-scale plague outbreaks (i.e., epizootic plague) are most likely to occur. PDOGMAPR was not designed to provide fine-scale detail on colony dynamics. Thus, the following questions will not be addressed accurately by the application:"),
                            tags$ul(
                              tags$li("Where exactly will colonies recover on the landscape after an epizootic plague event?"),
                              tags$li("How does smaller-scale enzootic plague restrict colony growth?"),
                              tags$li("How exactly will management tactics like plague mitigation and lethal control influence colony dynamics?")
                            ),
                            br(),
                            h4("Overall, PDOG MAPR is intended to be used as a guideline, along with expert knowledge of the species and system, to predict large-scale patterns in colony growth and epizootics and explore potential effects from the addition of management practices (i.e., it is not intended for precise, small-scale predictions). PDOG MAPR was designed such that if the colonies were to be changed via management in a certain way (however the user/manager believes that could be accomplished on their grassland), then the resulting colony distribution would be expected to have a certain chance of plague outbreak, as predicted by the model."),
                            br(),
                            h3("What PDOG MAPR does:", style = "color:#3474A7"),
                            h4("PDOG MAPR provides predictions of colony distribution and chance of plague outbreak for the year following the year of the user-provided colony shapefiles. Users can vary weather variables in the application to determine the expected effect of environmental conditions on colony dynamics. Users also can simulate the effect of plague mitigation (i.e., insecticide/vaccine application) and lethal control (i.e., poisoning colonies) on plague dynamics and prairie dog colony size and distribution."),
                            br(),
                            h3("How PDOG MAPR works:", style = "color:#3474A7"),
                            h4("Our research team developed a model to predict plague outbreaks and subsequent colony recovery in the prairie dog–plague disease system (Barrile et al., 2023). To make the model user-friendly, we integrated it into this interactive web-based, interactive decision support tool in shinyapps: PDOG MAPR."),
                            br(),
                            h4("When starting PDOG MAPR, the user is first prompted to upload a shapefile (via clicking the browse button to locate all of the components of the shapefile on their local machine) of prairie dog colony boundaries. Users are then prompted to upload two additional shapefiles related to on-the-ground management: one shapefile of areas where prairie dogs were lethally controlled in year t (i.e., the year of the uploaded colony shapefiles), and another shapefile of areas treated with plague mitigation in year t. Year t refers to the year of prairie dog colonies for which the user uploaded (e.g., if the user uploaded colonies from 2025, then year t is 2025). If users do not have these mitigation/management shapefiles, they can simply re-upload the colony shapefile from step 1 to proceed with the application. A fourth map is also plotted such that users can visualize how management influenced colony distribution (e.g., lethally controlled areas are subtracted from the original colony data). From this new spatial configuration of prairie dog colonies, the server side of PDOG MAPR recalculates colony connectedness variables (e.g., mean nearest neighbor, dispersion) that were important covariates in our predictive model."),
                            br(),
                            h4("Next, users are prompted to upload colony data from the same site but for the year prior (year t - 1). The server side of PDOG MAPR uses colony data from year t - 1 to calculate model covariates related to changes in colony distribution between years. For instance, colonies that suffered greater than 50% loss in area between years are assumed to have experienced a plague outbreak. If users do not have colony data from year t - 1, they can simply reupload their colony data from year t, in which case the model assumes there were no plague outbreaks between years. The server side of PDOG MAPR then runs the machine learning model and produces a predictive map of colony distribution in year t + 1. The final step in PDOG MAPR prompts the user to again upload colony data from year t - 1 (or reupload year t). The server side of the app then runs the model for a second time and calculates the average probability of plague outbreak (based on model predictions) for each colony on the landscape. Users can alter climate covariates in this step, as described above. Finally, PDOG MAPR presents a summary table displaying the total net change in colony area (in hectares) between the actual colonies in year t and the predicted colonies during t + 1, as well as the mean and maximum percent chance of plague outbreak across prairie dog colonies as predicted by the model."),
                            br(),
                            h4("At a minimum, PDOG MAPR can be used with simply a shapefile of colony data for a given year, describing colony dynamics under no management; however, PDOG MAPR can also incorporate management in the forms of lethal control or plague mitigation by uploading separate shapefiles for the desired management strategy. Preset management strategies were not incorporated into the application as management differs greatly between grasslands and for the most realistic outcomes, managers should upload scenarios relevant to their particular regions. PDOG MAPR also can be used to explore the influence of climate on prairie dog dynamics. Users have the option to vary temperature and precipitation data to test the effects of these covariates on model predictions of colony distribution in year t + 1. Specifically, users can choose to simulate average, warm, or cool temperatures and average, wet, or dry precipitation. Here climate variables represent the long-term average, min, and max during 2000–2020, which mirror the model covariates used in Barrile et al. (2023). Users can also vary the probability thresholds for colony growth and plague die-offs to examine how these user-defined thresholds affect colony and plague dynamics."),
                            br(),
                            h3("Note on management:", style = "color:#3474A7"),
                            h4("PDOG MAPR considers management practices to be 100% effective (i.e., poison removes all prairie dogs and plague mitigation ensures persistence); however, these management practices generally are not 100% effective and likely depend on a variety of factors such as climate and population density. While we acknowledge the intricacies missing from the application, PDOG MAPR was designed such that if the colonies were to be changed via management in a certain way (however the user/manager believes that could be accomplished on their grassland), then the resulting colony distribution would be expected to have a certain chance of plague outbreak, as predicted by the model."),
                            br(),
                            h3("Note on thresholds:", style = "color:#3474A7"),
                            h4("We recommend using the thresholds identified (0.90 for colony growth and 0.75 for plague die-off) as these most closely matched realistic growth and collapse in our modeling exercise (Barrile et al. 2023). Users are able to adjust these parameters as needed based on risk tolerance and their own familiarity with their system (for example, lowering the colony growth threshold if colonies expand at a high rate within the system). PDOG MAPR is intended to be used as a guideline, along with expert knowledge of the species and system, to predict large-scale patterns in colony growth and epizootics and explore potential effects from the addition of management practices (i.e., it is not intended for precise, small-scale predictions). Indeed, from personal communications with grassland managers, PDOG MAPR already is being used for its intended purposes. For instance, during winter of 2024, the Thunder Basin Working Group varied upcoming spring weather conditions in the application to produce scenarios that ranged from (1) lowest colony growth and highest chance of plague to (2) highest colony growth and lowest chance of plague. This exercise provided an idea of the range of possibilities for colony distribution across the range of potential weather scenarios, which helped in management planning."),
                            br(),
                            h3("Future directions:", style = "color:#3474A7"),
                            h4("Our current tool is intended to be the first step in an iterative process informed by conservation partners and managers that results in increasingly sophisticated decision support for prairie dog management. Future extensions of the decision support tool will include adding in probability of effectiveness for management techniques, as well as incorporating future climate scenarios, weather extremes, and interactions between management and climate. Another valuable extension of the tool will be to include the cost of proposed management practices. Prairie dog management would benefit from knowing how limited funding should be used to optimize biodiversity conservation alongside co-existence with livestock production. Overall, interactive, web-based tools such as PDOG MAPR deliver science to the fingertips of managers and decision makers. We hope that PDOG MAPR inspires others to provide web applications as part of their research process.")
                          )
                 ),
                 
                 tabPanel("Predict",
                          fluidPage(theme = shinytheme("sandstone"),
                                    titlePanel(p("Predict the distribution of prairie dog colonies and the chance of plague outbreak", style = "color:#3474A7")),
                                    
                                    sidebarLayout(
                                      sidebarPanel(
                                        # Data source indicator for multi-year iterations
                                        uiOutput("data_source_indicator"),
                                        
                                        # Step 0: Only show if this is Year 1 (no URL parameters)
                                        uiOutput("step0_year_selector"),
                                        
                                        # STEP 1 with popover
                                        div(
                                          fileInput("colonies", 
                                                    label = div(
                                                      uiOutput("colonies_label", inline = TRUE),
                                                      div(class = "popover-container", style = "display: inline-block;",
                                                          span(class = "help-icon", "?"),
                                                          span(class = "popover-text", 
                                                               "Click the browse button below and locate the shapefile (.shp) that represents the extent of your prairie dog colonies. For a successful upload you will need to select all four parts that make up a shape file: .dbf, .prj, .shp, and .shx). Click \"Open\" at the bottom right of the file search window, and the app will add your shapefile to the \"Your Inputs\" map.")
                                                      )
                                                    ),
                                                    multiple = TRUE,
                                                    accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj"))
                                        ),
                                        
                                        # STEP 2 with popover
                                        div(
                                          fileInput("poison",
                                                    label = div(
                                                      h5("Step 2. Upload shapefile of Prarie Dog control areas (Optional)", style = "display: inline;"),
                                                      div(class = "popover-container", style = "display: inline-block;",
                                                          span(class = "help-icon", "?"),
                                                          span(class = "popover-text", 
                                                               "Follow the same instructions as Step 1, except upload a shapefile that represents areas treated for prairie dog control, if any.")
                                                      )
                                                    ),
                                                    multiple = TRUE,
                                                    accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj"))
                                        ),
                                        
                                        # STEP 3 with popover
                                        div(
                                          fileInput("insect",
                                                    label = div(
                                                      h5("Step 3. Upload shapefile of plague treatment areas (Optional)", style = "display: inline;"),
                                                      div(class = "popover-container", style = "display: inline-block;",
                                                          span(class = "help-icon", "?"),
                                                          span(class = "popover-text", 
                                                               "Follow the same instructions as in Step 1, except upload a shapefile that represents areas treated for plague control, if any.")
                                                      )
                                                    ),
                                                    multiple = TRUE,
                                                    accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj"))
                                        ),
                                        
                                        # STEP 3b with popover
                                        div(
                                          h5(tags$b("Step 3b. Draw management areas (Optional)"),
                                             div(class = "popover-container", style = "display: inline-block;",
                                                 span(class = "help-icon", "?"),
                                                 span(class = "popover-text", 
                                                      "This step allows you to add any prairie dog control areas or plague control areas that you do not already have a shapefile for. Select the button corresponding to the type of control you wish to draw and then use the drawing tools on the \"Your Inputs\" map. You can also reset any drawings with the \"Clear all drawn polygons\" button here.")
                                             )
                                          ),
                                          radioButtons("drawMode", "Select Drawing Mode:",
                                                       choices = c("PDog Control Areas" = "poison", 
                                                                   "Plague Control Areas" = "insecticide"),
                                                       selected = "poison",
                                                       inline = TRUE),
                                          p(tags$i("Use polygon tool on map to draw areas"), style = "font-size: 12px; color: #666;"),
                                          actionButton("clearDrawn", "Clear All Drawn Polygons", 
                                                       class = "btn-warning btn-sm"),
                                          br(), br()
                                        ),
                                        
                                        # STEP 4 with popover
                                        div(
                                          fileInput("yrprior",
                                                    label = div(
                                                      uiOutput("yrprior_label", inline = TRUE),
                                                      div(class = "popover-container", style = "display: inline-block;",
                                                          span(class = "help-icon", "?"),
                                                          span(class = "popover-text", 
                                                               "Follow the same instructions as Step 1, except upload a shapefile that represents the extent of your prairie dog colonies for the year prior. If you do not have this data it is acceptable to upload the same colonies shape file as used in step 1 (current year).")
                                                      )
                                                    ),
                                                    multiple = TRUE,
                                                    accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj"))
                                        ),
                                        
                                        # STEP 5 with popover
                                        div(
                                          selectInput("region", 
                                                      label = div(
                                                        h5("Step 5. Select Study Region", style = "display: inline;"),
                                                        div(class = "popover-container", style = "display: inline-block;",
                                                            span(class = "help-icon", "?"),
                                                            span(class = "popover-text", 
                                                                 "The model is calibrated to take into account the region you are working in. Selecting the region that most closely relates to your study area will give the best results. Selecting \"Average\" will run the model using average paramaters that represent the model across the entire central US and may be appropriate where your study area is well away from any of the identified regions, or spans multiple regions. Selecting a region that is not close to where your data comes from may result in unreliable results.")
                                                        )
                                                      ),
                                                      choices = c("Population Average" = "population_avg",
                                                                  "Cimarron" = "Cimarron", 
                                                                  "CMRfull" = "CMRfull",
                                                                  "ComancheNW" = "ComancheNW",
                                                                  "ComancheSE" = "ComancheSE", 
                                                                  "Kiowa" = "Kiowa", 
                                                                  "PawneeE" = "PawneeE", 
                                                                  "PawneeW" = "PawneeW", 
                                                                  "RitaBlanca" = "RitaBlanca", 
                                                                  "ThunderBasin" = "ThunderBasin"),
                                                      selected = "population_avg"),
                                          
                                          helpText("⚠️ WARNING: Select the site that matches your study area!",
                                                   "Using a different site may create biologically implausible predictions.",
                                                   style = "color: #FF6347; font-weight: bold;")
                                        ),
                                        
                                        h4("Climate Parameters & Thresholds", style = "color: #3474A7;"),
                                        
                                        h5(tags$b("Climate Variables"), "(applied to both growth and plague models):"),
                                        
                                        selectInput("ZZTEMP", 
                                                    label = textOutput("zztemp_label"),
                                                    choices = c("Average","Cold","Warm")),
                                        
                                        selectInput("WholeYear", 
                                                    label = textOutput("wholeyear_label"),
                                                    choices = c("Average","Wet","Dry")),
                                        
                                        selectInput("WinterSpring", 
                                                    label = textOutput("winterspring_label"),
                                                    choices = c("Average","Wet","Dry")),
                                        
                                        selectInput("WSSF", 
                                                    label = textOutput("wssf_label"),
                                                    choices = c("Average","Dry Fall Wet Spring","Wet Fall Dry Spring")),
                                        
                                        h5(tags$b("Prediction Thresholds:")),
                                        
                                        sliderInput("threshold", "Threshold for Colony Growth:",
                                                    min = 0.5, max = 1, value = 0.90, step = 0.05),
                                        
                                        sliderInput("thresholdp", "Threshold for Plague Die-off:",
                                                    min = 0.5, max = 1, value = 0.75, step = 0.05),
                                        
                                        helpText("Note: Default plague threshold (0.75) is calibrated for glmmTMB predictions.",
                                                 "Lower values = more permissive (more colonies persist)."),
                                        
                                        br(),
                                        actionButton("runModels", "Run Models", 
                                                     class = "btn-success btn-lg", 
                                                     style = "width: 100%; margin-top: 20px;"),
                                        
                                        br(), br(),
                                        tags$button(
                                          id = "ShapeExport",
                                          class = "btn btn-secondary btn-lg",
                                          style = "width: 100%;",
                                          disabled = "disabled",
                                          tags$i(class = "fa fa-download"), " Download Predicted Colonies"
                                        ),
                                        
                                        width = 3
                                      ),
                                      
                                      mainPanel(
                                        # Year indicator and iteration button
                                        div(style = "text-align: center; margin-bottom: 15px; padding: 10px; background-color: #e8f4f8; border-radius: 8px;",
                                            h4(uiOutput("year_indicator"), style = "color: #3474A7; margin: 5px 0;"),
                                            uiOutput("carry_forward_button")
                                        ),
                                        
                                        div(style = "margin-bottom: 20px;",
                                            uiOutput("inputs_title"),
                                            leafletOutput("main_overlay_map", height = "400px")
                                        ),
                                        
                                        div(style = "margin-bottom: 20px;",
                                            uiOutput("prediction_title"),
                                            leafletOutput("prediction_results", height = "350px")
                                        ),
                                        
                                        div(style = "margin-bottom: 20px;",
                                            uiOutput("plague_title"),
                                            leafletOutput("plgpred", height = "350px")
                                        ),
                                        
                                        div(style = "display: flex; align-items: flex-start; gap: 20px;",
                                            div(style = "flex: 0.6;",
                                                h4("Summary Statistics", style = "color: #3474A7; text-align: center;"),
                                                tableOutput("irish"),
                                                
                                                br(),
                                                
                                                # ENHANCED ECONOMIC CALCULATOR
                                                div(id = "costBenefitBox",
                                                    style = "border: 2px solid #3474A7; border-radius: 8px; padding: 15px; background-color: #f8f9fa;",
                                                    h5("Economic Analysis - 4 Scenarios", style = "color: #3474A7; text-align: center; margin-top: 0; margin-bottom: 10px;"),
                                                    
                                                    actionButton("calcCosts", "Calculate Economic Analysis", 
                                                                 class = "btn-primary btn-sm", 
                                                                 style = "width: 100%; margin-bottom: 15px;"),
                                                    
                                                    fluidRow(
                                                      column(6,
                                                             div(class = "scenario-box",
                                                                 div(class = "scenario-title", "A: Average Pdog consumption values"),
                                                                 div(class = "economic-input",
                                                                     numericInput("acres_poisoned_A", 
                                                                                  label = "Acres Pdog Control:",
                                                                                  value = 0, min = 0, step = 0.1)),
                                                                 div(class = "economic-input",
                                                                     numericInput("acres_pdogs_A", 
                                                                                  label = "Acres Pdogs after colony prediction:",
                                                                                  value = 0, min = 0, step = 0.1)),
                                                                 div(class = "economic-input",
                                                                     numericInput("cost_per_acre_A", 
                                                                                  label = "Control cost/Acre ($):",
                                                                                  value = 25, min = 0, step = 1)),
                                                                 hr(style = "margin: 8px 0;"),
                                                                 div(class = "output-value",
                                                                     strong("Change in forage (lbs):"),
                                                                     br(),
                                                                     textOutput("forage_lost_A")),
                                                                 div(class = "output-value",
                                                                     strong("Expenditure on Pdog control ($):"),
                                                                     br(),
                                                                     textOutput("total_cost_A"))
                                                             )
                                                      ),
                                                      column(6,
                                                             div(class = "scenario-box",
                                                                 div(class = "scenario-title", "B: Pdogs have no influence on forage"),
                                                                 div(class = "economic-input",
                                                                     numericInput("acres_poisoned_B", 
                                                                                  label = "Acres Pdog Control:",
                                                                                  value = 0, min = 0, step = 0.1)),
                                                                 div(class = "economic-input",
                                                                     numericInput("acres_pdogs_B", 
                                                                                  label = "Acres Pdogs after colony prediction:",
                                                                                  value = 0, min = 0, step = 0.1)),
                                                                 div(class = "economic-input",
                                                                     numericInput("cost_per_acre_B", 
                                                                                  label = "Pdog control cost/Acre ($):",
                                                                                  value = 25, min = 0, step = 1)),
                                                                 hr(style = "margin: 8px 0;"),
                                                                 div(class = "output-value",
                                                                     strong("Change in forage (lbs):"),
                                                                     br(),
                                                                     textOutput("forage_lost_B")),
                                                                 div(class = "output-value",
                                                                     strong("Expenditure on Pdog control ($):"),
                                                                     br(),
                                                                     textOutput("total_cost_B"))
                                                             )
                                                      )
                                                    ),
                                                    
                                                    fluidRow(
                                                      column(6,
                                                             div(class = "scenario-box",
                                                                 div(class = "scenario-title", "C: Pdog colonies yield 50% Forage"),
                                                                 div(class = "economic-input",
                                                                     numericInput("acres_poisoned_C", 
                                                                                  label = "Acres Pdog control:",
                                                                                  value = 0, min = 0, step = 0.1)),
                                                                 div(class = "economic-input",
                                                                     numericInput("acres_pdogs_C", 
                                                                                  label = "Acres Pdogs after colony prediction:",
                                                                                  value = 0, min = 0, step = 0.1)),
                                                                 div(class = "economic-input",
                                                                     numericInput("ppa_C", 
                                                                                  label = "PPA w/o Pdogs:",
                                                                                  value = 2000, min = 0, step = 10)),
                                                                 div(class = "economic-input",
                                                                     numericInput("cost_per_acre_C", 
                                                                                  label = "Control cost/Acre ($):",
                                                                                  value = 25, min = 0, step = 1)),
                                                                 hr(style = "margin: 8px 0;"),
                                                                 div(class = "output-value",
                                                                     strong("Change in forage (lbs):"),
                                                                     br(),
                                                                     textOutput("forage_lost_C")),
                                                                 div(class = "output-value",
                                                                     strong("Expenditure on Pdog control ($):"),
                                                                     br(),
                                                                     textOutput("total_cost_C"))
                                                             )
                                                      ),
                                                      column(6,
                                                             div(class = "scenario-box",
                                                                 div(class = "scenario-title", "D: Pdog colonies yield zero forage"),
                                                                 div(class = "economic-input",
                                                                     numericInput("acres_poisoned_D", 
                                                                                  label = "Acres Pdog control:",
                                                                                  value = 0, min = 0, step = 0.1)),
                                                                 div(class = "economic-input",
                                                                     numericInput("acres_pdogs_D", 
                                                                                  label = "Acres Pdogs after colony prediction:",
                                                                                  value = 0, min = 0, step = 0.1)),
                                                                 div(class = "economic-input",
                                                                     numericInput("ppa_D", 
                                                                                  label = "PPA w/o Pdogs:",
                                                                                  value = 2000, min = 0, step = 10)),
                                                                 div(class = "economic-input",
                                                                     numericInput("cost_per_acre_D", 
                                                                                  label = "Control cost/Acre ($):",
                                                                                  value = 25, min = 0, step = 1)),
                                                                 hr(style = "margin: 8px 0;"),
                                                                 div(class = "output-value",
                                                                     strong("Change in forage (lbs):"),
                                                                     br(),
                                                                     textOutput("forage_lost_D")),
                                                                 div(class = "output-value",
                                                                     strong("Expenditure on Pdog control ($):"),
                                                                     br(),
                                                                     textOutput("total_cost_D"))
                                                             )
                                                      )
                                                    )
                                                )
                                            ),
                                            div(style = "flex: 1.4;",
                                                h4("Plague Analysis", style = "color: #3474A7; text-align: center;"),
                                                div(style = "display: flex; gap: 10px; margin-bottom: 20px;",
                                                    div(style = "flex: 1;",
                                                        h5("Colony Count by Probability", style = "text-align: center; margin-bottom: 5px;"),
                                                        plotOutput("plague_distribution", height = "450px")
                                                    ),
                                                    div(style = "flex: 1;",
                                                        h5("Area by Probability", style = "text-align: center; margin-bottom: 5px;"),
                                                        plotOutput("plague_area_distribution", height = "450px")
                                                    )
                                                ),
                                                
                                                # SCENARIO DESCRIPTIONS TABLE
                                                div(style = "border: 2px solid #3474A7; border-radius: 8px; padding: 15px; background-color: #f8f9fa; margin-bottom: 20px; position: relative;",
                                                    # Arrow pointing left to economic calculator
                                                    div(class = "linking-arrow"),             h5("Economic Scenario Descriptions", style = "color: #3474A7; text-align: center; margin-top: 0; margin-bottom: 15px;"),
                                                    
                                                    tags$table(style = "width: 100%; border-collapse: collapse; font-size: 12px;",
                                                               tags$thead(
                                                                 tags$tr(
                                                                   tags$th(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; text-align: left;", ""),
                                                                   tags$th(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; text-align: center; font-weight: bold;", "SCENARIO A: Average Prairie Dog Consumption Values"),
                                                                   tags$th(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; text-align: center; font-weight: bold;", "SCENARIO B: Prairie Dog Colonies Yield 100% Forage"),
                                                                   tags$th(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; text-align: center; font-weight: bold;", "SCENARIO C: Prairie Dog Colonies Yield 50% Forage"),
                                                                   tags$th(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; text-align: center; font-weight: bold;", "SCENARIO D: Prairie Dog Colonies Yield 0% Forage")
                                                                 )
                                                               ),
                                                               tags$tbody(
                                                                 # Row 1: Scenario (light blue header)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Scenario"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8;", "Average values for prairie dog/stock interactions are assumed."),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8;", "It is assumed that areas occupied by prairie dogs yield the same forage as land unoccupied by prairie dogs (any change in colony area will result in no gain or loss of forage)."),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8;", "It is assumed that areas occupied by prairie dogs provide 50% forage value to stock, thus area modeled as being released from prairie dogs will result in a 50% gain in forage value."),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8;", "It is assumed that area occupied by prairie dogs provide 0% forage value to stock, thus area modeled as being released from prairie dogs will result in an 100% gain in forage value.")
                                                                 ),
                                                                 # Row 2: Pdog control efficiency (light blue header, light grey data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Pdog control efficiency"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "100%"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "100%"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "100%"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "100%")
                                                                 ),
                                                                 # Row 3: Pdog density (light blue header, light grey data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Pdog density"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "42/acre"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "null"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "null"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "null")
                                                                 ),
                                                                 # Row 4: Pdogs per head of cattle (light blue header, light grey data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Pdogs per head of cattle"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "335:1"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "null"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "null"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "null")
                                                                 ),
                                                                 # Row 5: lb forage per head of cattle per day (light blue header, light grey data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "lb forage per head of cattle per day"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "25 (9,125/year)"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "null"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "null"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "null")
                                                                 ),
                                                                 # Row 6: Forage yield multiplier (light blue header, light grey data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Forage yield multiplier"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "100%"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "100%"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "50%"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8e8e8; text-align: center;", "0%")
                                                                 ),
                                                                 # Row 7: Acres of pdog control (light blue header, light yellow data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Acres of pdog control"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #fffacd; text-align: center;", "calculated in app"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #fffacd; text-align: center;", "calculated in app"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #fffacd; text-align: center;", "calculated in app"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #fffacd; text-align: center;", "calculated in app")
                                                                 ),
                                                                 # Row 8: Acres of pdog colony after model predictions (light blue header, light yellow data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Acres of pdog colony after model predictions"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #fffacd; text-align: center;", "calculated in app"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #fffacd; text-align: center;", "calculated in app"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #fffacd; text-align: center;", "calculated in app"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #fffacd; text-align: center;", "calculated in app")
                                                                 ),
                                                                 # Row 9: Production without pdogs (light blue header, white data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Production without pdogs (lbs/acre)"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #ffffff; text-align: center;", "Null"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #ffffff; text-align: center;", "null"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #ffffff; text-align: center;", "User input"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #ffffff; text-align: center;", "User input")
                                                                 ),
                                                                 # Row 10: Cost of pdog control per acre (light blue header, white data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Cost of pdog control per acre"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #ffffff; text-align: center;", "user input"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #ffffff; text-align: center;", "user input"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #ffffff; text-align: center;", "user input"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #ffffff; text-align: center;", "user input")
                                                                 ),
                                                                 # Row 11: Change in forage (light blue header, light green data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Change in forage (lbs/acre)"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #d4edda; text-align: center;", "output"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #d4edda; text-align: center;", "output"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #d4edda; text-align: center;", "output"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #d4edda; text-align: center;", "output")
                                                                 ),
                                                                 # Row 12: Expenditure on pdog control measures (light blue header, light green data)
                                                                 tags$tr(
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #e8f4f8; font-weight: bold;", "Expenditure on pdog control measures"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #d4edda; text-align: center;", "output"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #d4edda; text-align: center;", "output"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #d4edda; text-align: center;", "output"),
                                                                   tags$td(style = "border: 1px solid #3474A7; padding: 8px; background-color: #d4edda; text-align: center;", "output")
                                                                 )
                                                               )
                                                    )
                                                ),
                                                
                                                # BACKGROUND AND RESOURCES
                                                div(style = "border: 2px solid #3474A7; border-radius: 8px; padding: 15px; background-color: #f8f9fa;",
                                                    h5("Additional Information", style = "color: #3474A7; text-align: center; margin-top: 0; margin-bottom: 15px;"),
                                                    
                                                    div(style = "margin-bottom: 15px;",
                                                        h6(strong("Scenario Background:"), style = "color: #3474A7; margin-bottom: 8px;"),
                                                        p("The four scenarios presented here are simple examples, designed to provide a general sense of how forage yield may vary under different circumstances given the prairie dog demographics predicted by the model. Scenario A uses typical values (as determined by the Prairie Dog Mapper development team) to present an 'average case'. Scenario B is effectively a 'null hypothesis' which assumes that prairie dog presence makes no change to forage values. Scenarios C and D present examples for where areas occupied by prairie dogs allow forage cover of either 50%, or 0% (barren). These scenarios are not expected to exactly match any real life situations, but offer a range and bounds within which the user can calibrate their expectations.",
                                                          style = "font-size: 13px; line-height: 1.5; margin: 0;")
                                                    ),
                                                    
                                                    div(style = "margin-bottom: 15px;",
                                                        h6(strong("Assumed Values:"), style = "color: #3474A7; margin-bottom: 8px;"),
                                                        p("Assumed values used in calculations (grey fields in the Economic Scenario Descriptions table) are typical values as determined by the Prairie Dog Mapper development team. Actual values for these variables (and other potentially relevant ones) are highly complex, inter-correlated, and likely to vary widely over space and time. The 100% efficacy value for prairie dog control is likely unrealistic, however determining realistic case-specific values would require sophisticated modeling outside the scope of this tool. We encourage users with more specific knowledge of the dynamics particular to their area(s) of interest to do their own calculations using the changes in colony area provided by the model. Links are provided below for resources that may be helpful in understanding some of these variables (also the averages chosen in Scenario A) and which may be useful in applying to your own calculations.", 
                                                          style = "font-size: 13px; line-height: 1.5; margin: 0;")
                                                    ),
                                                    
                                                    div(style = "margin-bottom: 15px;",
                                                        h6(strong("Data Resources:"), style = "color: #3474A7; margin-bottom: 8px;"),
                                                        p(strong("USDA WebSoil Survey"), " (Soils, Range Production, Yield, Land Classifications, etc.)", 
                                                          style = "font-size: 13px; line-height: 1.5; margin: 5px 0;"),
                                                        p(tags$a(href = "https://websoilsurvey.nrcs.usda.gov/app/WebSoilSurvey.aspx", 
                                                                 target = "_blank",
                                                                 "https://websoilsurvey.nrcs.usda.gov/app/WebSoilSurvey.aspx"),
                                                          style = "font-size: 13px; margin: 0 0 10px 20px;"),
                                                        p(strong("Rangeland Analysis Platform and Production Explorer:"), 
                                                          style = "font-size: 13px; line-height: 1.5; margin: 5px 0;"),
                                                        p(tags$a(href = "https://rangelands.app/", 
                                                                 target = "_blank",
                                                                 "https://rangelands.app/"),
                                                          style = "font-size: 13px; margin: 0 0 0 20px;")
                                                    ),
                                                    
                                                    div(style = "margin-bottom: 0;",
                                                        h6(strong("Further Reading:"), style = "color: #3474A7; margin-bottom: 8px;"),
                                                        p("Augustine & Derner, 2021: \"Long-Term Effects of Black-Tailed Prairie Dogs on Livestock Grazing Distribution and Mass Gain\"; The Journal of Wildlife Management, 85(7):1332-1343", 
                                                          style = "font-size: 13px; line-height: 1.5; margin: 5px 0;"),
                                                        p("Augustine et al., 2024: \"Does Drought Intensify the Effects of Black-Tailed Prairie Dogs on Livestock Production and Net Revenue in Semiarid Rangelands?\"; Rangeland Ecology and Management, 10.1016/j.rama.2024.04.011", 
                                                          style = "font-size: 13px; line-height: 1.5; margin: 5px 0;"),
                                                        p("Buehler et al., 2025; Exploring the Efficacy of Prairie Dog Boundary Management and its Application Toward Density Control; Rangeland Ecology and Management, 99(2):66-76.", 
                                                          style = "font-size: 13px; line-height: 1.5; margin: 5px 0;"),
                                                        p("Crow et al., 2022; Evaluating Prairie Dog-Cattle Competition from the Perspective of a Ranching Enterprise in the Western Great Plains: Economic Analysis of Potential Effects on Long-Term Profitability; Rangeland Ecology and Management, 85: 56-65", 
                                                          style = "font-size: 13px; line-height: 1.5; margin: 5px 0 0 0;")
                                                    )
                                                )
                                            )
                                        ),
                                        
                                        width = 9
                                      )
                                    )
                          )
                 )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 5) Server Function -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

server <- function(input, output, session) {
  # Create session-specific temporary directory for carrying forward data
  session_id <- paste0("session_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", 
                       substr(digest(runif(1)), 1, 8))
  temp_dir <- file.path(tempdir(), "prairie_dog_carryover", session_id)
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Clean up temp files when session ends
  session$onSessionEnded(function() {
    unlink(temp_dir, recursive = TRUE)
  })
  
  map_trigger <- reactiveVal(0)
  map_update_trigger <- reactiveVal(0)
  
  # Parse URL parameters to check if this is a carried-over session
  url_params <- reactiveValues(
    has_carryover = FALSE,
    carryover_id = NULL,
    year_num = 1,
    base_year = NULL
  )
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query$carryover_id)) {
      url_params$has_carryover <- TRUE
      url_params$carryover_id <- query$carryover_id
      url_params$year_num <- as.numeric(query$year_num)
      url_params$base_year <- as.numeric(query$base_year)
      cat("Loaded carried-over data: Year", url_params$year_num, "\n")
    }
  })
  
  values <- reactiveValues(
    current_colonies = NULL,
    poison_areas = NULL,
    insecticide_areas = NULL,
    processed_colonies = NULL,
    predicted_colonies = NULL,
    predicted_plague = NULL,
    models_completed = FALSE,
    last_results = NULL,
    button_clicked = FALSE,
    models_running = FALSE,
    drawn_poison_polygons = list(),
    drawn_insecticide_polygons = list(),
    poison_counter = 0,
    insecticide_counter = 0,
    poison_upload_year = NULL,
    insecticide_upload_year = NULL,
    economic_results = list(
      A = list(forage_lost = 0, total_cost = 0),
      B = list(forage_lost = 0, total_cost = 0),
      C = list(forage_lost = 0, total_cost = 0),
      D = list(forage_lost = 0, total_cost = 0)
    ),
    # Multi-year iteration state
    current_year = 1,
    base_year = NULL,  # The year selected in Step 0 (or NULL for t notation)
    carried_colonies_current = NULL,
    carried_colonies_prior = NULL
  )
  
  # Initialize year from URL if carried over
  observe({
    if (url_params$has_carryover) {
      values$current_year <- url_params$year_num
      values$base_year <- url_params$base_year
    }
  })
  
  observe({
    runjs("$('#costBenefitBox').addClass('disabled');")
  })
  
  # ============================================================================
  # LOAD CARRIED-OVER DATA FROM TEMP FILES
  # ============================================================================
  
  observe({
    if (url_params$has_carryover && !is.null(url_params$carryover_id)) {
      carryover_path <- file.path(tempdir(), "prairie_dog_carryover", url_params$carryover_id)
      
      current_path <- file.path(carryover_path, "current_colonies.rds")
      prior_path <- file.path(carryover_path, "prior_colonies.rds")
      
      if (file.exists(current_path) && file.exists(prior_path)) {
        values$carried_colonies_current <- readRDS(current_path)
        values$carried_colonies_prior <- readRDS(prior_path)
        
        showNotification(
          paste("Loaded carried-over data for Year", values$current_year),
          type = "message",
          duration = 5
        )
        
        cat("Successfully loaded carried colonies for Year", values$current_year, "\n")
      } else {
        showNotification(
          "Error: Could not load carried-over data. Files may have expired.",
          type = "error",
          duration = 10
        )
      }
    }
  })
  
  # ============================================================================
  # UI OUTPUTS
  # ============================================================================
  
  # Only show Step 0 if this is Year 1 (no carryover)
  output$step0_year_selector <- renderUI({
    if (values$current_year == 1 && !url_params$has_carryover) {
      selectInput("currentYear",
                  label = h5("Optional Step 0. What year is your data for? (just makes prompts clearer)"),
                  choices = c("Use t/t+1/t-1 notation" = "", 2000:2100),
                  selected = "")
    } else {
      NULL
    }
  })
  
  # Data source indicator (sidebar)
  output$data_source_indicator <- renderUI({
    if (url_params$has_carryover) {
      div(
        style = "background-color: #d1ecf1; border: 1px solid #bee5eb; border-radius: 5px; padding: 10px; margin-bottom: 15px;",
        h5(icon("info-circle"), " Multi-Year Iteration Mode", style = "color: #0c5460; margin-top: 0;"),
        p(paste("Using predicted colonies from Year", values$current_year - 1, "as starting point"),
          style = "color: #0c5460; margin-bottom: 5px; font-size: 13px;"),
        p("Colony uploads disabled - data carried forward automatically",
          style = "color: #0c5460; margin-bottom: 0; font-size: 12px; font-style: italic;")
      )
    } else {
      NULL
    }
  })
  
  # Helper function to get year label
  get_year_label <- function(year_offset = 0) {
    if (is.null(values$base_year)) {
      # Use t notation
      if (year_offset == 0) return("t")
      else if (year_offset == 1) return("t + 1")
      else if (year_offset == -1) return("t - 1")
      else return(paste0("t ", ifelse(year_offset > 0, "+ ", "- "), abs(year_offset)))
    } else {
      # Use actual years
      target_year <- values$base_year + (values$current_year - 1) + year_offset
      return(as.character(target_year))
    }
  }
  
  # Set base_year from Step 0 selection
  observe({
    if (!is.null(input$currentYear) && input$currentYear != "" && values$current_year == 1) {
      values$base_year <- as.numeric(input$currentYear)
    }
  })
  
  # Dynamic labels
  output$colonies_label <- renderUI({
    current_year <- get_year_label(0)
    h5(paste("Step 1. Upload colony shapefile for", current_year), style = "display: inline;")
  })
  
  output$yrprior_label <- renderUI({
    prior_year <- get_year_label(-1)
    h5(paste("Step 4. Upload colony shapefile for", prior_year), style = "display: inline;")
  })
  
  output$zztemp_label <- renderText({
    prior_year <- get_year_label(-1)
    paste("Choose temperature in summer of", prior_year)
  })
  
  output$wholeyear_label <- renderText({
    current_year <- get_year_label(0)
    paste("Choose precipitation in year", current_year)
  })
  
  output$winterspring_label <- renderText({
    next_year <- get_year_label(1)
    paste("Choose precip in spring of", next_year)
  })
  
  output$wssf_label <- renderText({
    current_year <- get_year_label(0)
    next_year <- get_year_label(1)
    paste("Choose change in precip btwn", current_year, "and", next_year)
  })
  
  output$prediction_title <- renderUI({
    next_year <- get_year_label(1)
    h3(paste("Predicted Colonies (", next_year, ")", sep=""), style = "text-align: center; color: #3474A7;")
  })
  
  # ADD THESE TWO NEW OUTPUTS:
  output$inputs_title <- renderUI({
    current_year <- get_year_label(0)
    h3(paste("Your Inputs (", current_year, ")", sep=""), style = "text-align: center; color: #3474A7;")
  })
  
  output$plague_title <- renderUI({
    next_year <- get_year_label(1)
    h3(paste("Predicted Plague Probability (", next_year, ")", sep=""), style = "text-align: center; color: #3474A7;")
  })
  
  # Year indicator output (main panel)
  output$year_indicator <- renderUI({
    current_display <- get_year_label(0)
    next_display <- get_year_label(1)
    
    if (values$current_year == 1) {
      span(paste("Year", current_display, "→", next_display, "Prediction"), style = "font-weight: bold;")
    } else {
      span(paste("Year", current_display, "→", next_display, "Prediction (Iteration", values$current_year - 1, ")"),
           style = "font-weight: bold;")
    }
  })
  
  # Carry forward button - opens NEW TAB
  output$carry_forward_button <- renderUI({
    req(values$models_completed)
    req(values$predicted_colonies)
    
    actionButton(
      inputId = "carry_forward",
      label = sprintf("Carry Forward to Year %s (Opens New Tab) ➔", get_year_label(1)),
      class = "btn-primary btn-lg",
      style = "margin-top: 10px; font-size: 16px;"
    )
  })
  
  # ============================================================================
  # CARRY FORWARD TO NEW TAB
  # ============================================================================
  
  observeEvent(input$carry_forward, {
    req(values$models_completed)
    req(values$predicted_colonies)
    req(values$current_colonies)
    
    # Create unique carryover ID
    carryover_id <- paste0("carry_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", 
                           substr(digest(runif(1)), 1, 8))
    carryover_path <- file.path(tempdir(), "prairie_dog_carryover", carryover_id)
    dir.create(carryover_path, showWarnings = FALSE, recursive = TRUE)
    
    # Save the colonies to temp files
    # Predicted colonies become "current" for next year
    # Current colonies become "prior" for next year
    saveRDS(values$predicted_colonies, file.path(carryover_path, "current_colonies.rds"))
    saveRDS(values$current_colonies, file.path(carryover_path, "prior_colonies.rds"))
    
    # Build URL for new tab
    next_year_num <- values$current_year + 1
    base_url <- session$clientData$url_pathname
    
    # Build query string
    query_params <- paste0("?carryover_id=", carryover_id,
                           "&year_num=", next_year_num,
                           "&base_year=", ifelse(is.null(values$base_year), "", values$base_year))
    
    new_url <- paste0(base_url, query_params)
    
    # Open new tab with JavaScript
    runjs(sprintf("window.open('%s', '_blank');", new_url))
    
    showNotification(
      paste("Opening new tab for Year", get_year_label(1), "- Data will load automatically"),
      type = "message",
      duration = 5
    )
    
    cat("Carried forward to new tab: Year", next_year_num, "\n")
  })
  
  # ============================================================================
  # DRAWING AND MANAGEMENT AREAS
  # ============================================================================
  
  convert_drawn_to_sf <- function(feature, crs = 4326) {
    if (is.null(feature) || length(feature) == 0) return(NULL)
    
    tryCatch({
      coords <- feature$geometry$coordinates[[1]]
      if (length(coords) == 0) return(NULL)
      
      coord_matrix <- do.call(rbind, lapply(coords, function(x) c(x[[1]], x[[2]])))
      
      if (!all(coord_matrix[1,] == coord_matrix[nrow(coord_matrix),])) {
        coord_matrix <- rbind(coord_matrix, coord_matrix[1,])
      }
      
      poly <- st_polygon(list(coord_matrix))
      sf_obj <- st_sfc(poly, crs = crs)
      sf_obj <- st_sf(geometry = sf_obj)
      
      return(sf_obj)
    }, error = function(e) {
      cat("Error converting drawn feature:", e$message, "\n")
      return(NULL)
    })
  }
  
  observeEvent(input$main_overlay_map_draw_new_feature, {
    feature <- input$main_overlay_map_draw_new_feature
    current_mode <- input$drawMode
    
    sf_poly <- convert_drawn_to_sf(feature)
    
    if (!is.null(sf_poly)) {
      if (current_mode == "poison") {
        values$poison_counter <- values$poison_counter + 1
        poly_id <- paste0("poison_", values$poison_counter)
        values$drawn_poison_polygons[[poly_id]] <- sf_poly
        cat("Added poison polygon:", poly_id, "\n")
      } else if (current_mode == "insecticide") {
        values$insecticide_counter <- values$insecticide_counter + 1
        poly_id <- paste0("insecticide_", values$insecticide_counter)
        values$drawn_insecticide_polygons[[poly_id]] <- sf_poly
        cat("Added insecticide polygon:", poly_id, "\n")
      }
    }
  })
  
  observeEvent(input$clearDrawn, {
    values$drawn_poison_polygons <- list()
    values$drawn_insecticide_polygons <- list()
    values$poison_counter <- 0
    values$insecticide_counter <- 0
    
    showNotification("All drawn polygons cleared", type = "message", duration = 2)
  })
  
  observeEvent(input$poison, {
    if (!is.null(input$poison)) {
      values$poison_upload_year <- values$current_year
      cat("Poison file uploaded in Year", values$current_year, "\n")
    }
  })
  
  observeEvent(input$insect, {
    if (!is.null(input$insect)) {
      values$insecticide_upload_year <- values$current_year
      cat("Insecticide file uploaded in Year", values$current_year, "\n")
    }
  })
  
  # Trigger map update when files are uploaded
  observeEvent(input$poison, {
    if (!is.null(input$poison)) {
      map_update_trigger(map_update_trigger() + 1)
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$insect, {
    if (!is.null(input$insect)) {
      map_update_trigger(map_update_trigger() + 1)
    }
  }, ignoreInit = TRUE)
  
  # ============================================================================
  # COLONY DATA LOADING
  # ============================================================================
  
  current_colonies <- reactive({
    # If carried over, use that data
    if (!is.null(values$carried_colonies_current)) {
      col <- values$carried_colonies_current
      values$current_colonies <- col
      return(col)
    }
    
    # Otherwise require file upload
    req(input$colonies)
    
    shpdf <- input$colonies
    tempdirname <- dirname(shpdf$datapath[1])
    
    for (i in 1:nrow(shpdf)) {
      new_path <- paste0(tempdirname, "/", shpdf$name[i])
      if (!file.exists(new_path)) {
        file.copy(shpdf$datapath[i], new_path, overwrite = TRUE)
      }
    }
    
    col <- st_read(paste(tempdirname, shpdf$name[grep(pattern = "*.shp$", shpdf$name)], sep = "/"))
    col <- st_transform(col, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
    col <- st_buffer(col, 0)
    col <- st_transform(col, crs = "+proj=longlat +datum=WGS84")
    
    values$current_colonies <- col
    return(col)
  })
  
  poison_areas <- reactive({
    poison_sf <- NULL
    
    if (!is.null(input$poison) && 
        !is.null(values$poison_upload_year) && 
        values$poison_upload_year == values$current_year) {
      
      shpdf <- input$poison
      tempdirname <- dirname(shpdf$datapath[1])
      
      for (i in 1:nrow(shpdf)) {
        new_path <- paste0(tempdirname, "/", shpdf$name[i])
        if (!file.exists(new_path)) {
          file.copy(shpdf$datapath[i], new_path, overwrite = TRUE)
        }
      }
      
      pois <- st_read(paste(tempdirname, shpdf$name[grep(pattern = "*.shp$", shpdf$name)], sep = "/"), quiet = TRUE)
      pois <- st_transform(pois, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      pois <- st_buffer(pois, 0)
      pois <- st_transform(pois, crs = "+proj=longlat +datum=WGS84")
      
      colcheck <- current_colonies()
      if (identical(colcheck, pois)) {
        poisPT <- st_sample(pois, size = 1, type = "random", by_polygon = FALSE)
        poisPT <- st_transform(poisPT, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
        poisPT <- st_buffer(poisPT, 1)
        pois <- st_transform(poisPT, crs = "+proj=longlat +datum=WGS84")
      }
      
      poison_sf <- st_geometry(pois)
      poison_sf <- st_sf(geometry = poison_sf)
    }
    
    if (length(values$drawn_poison_polygons) > 0) {
      drawn_sf <- do.call(rbind, values$drawn_poison_polygons)
      drawn_sf <- st_transform(drawn_sf, crs = "+proj=longlat +datum=WGS84 +no_defs")
      
      if (!is.null(poison_sf)) {
        combined_geom <- st_union(st_geometry(poison_sf), st_geometry(drawn_sf))
        poison_sf <- st_sf(geometry = combined_geom)
      } else {
        poison_sf <- drawn_sf
      }
    }
    
    if (!is.null(poison_sf)) {
      values$poison_areas <- poison_sf
    }
    return(poison_sf)
  })
  
  insecticide_areas <- reactive({
    insect_sf <- NULL
    
    if (!is.null(input$insect) && 
        !is.null(values$insecticide_upload_year) && 
        values$insecticide_upload_year == values$current_year) {
      
      shpdf <- input$insect
      tempdirname <- dirname(shpdf$datapath[1])
      
      for (i in 1:nrow(shpdf)) {
        new_path <- paste0(tempdirname, "/", shpdf$name[i])
        if (!file.exists(new_path)) {
          file.copy(shpdf$datapath[i], new_path, overwrite = TRUE)
        }
      }
      
      mitig <- st_read(paste(tempdirname, shpdf$name[grep(pattern = "*.shp$", shpdf$name)], sep = "/"), quiet = TRUE)
      mitig <- st_transform(mitig, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      mitig <- st_buffer(mitig, 0)
      mitig <- st_transform(mitig, crs = "+proj=longlat +datum=WGS84")
      
      colcheck <- current_colonies()
      if (identical(colcheck, mitig)) {
        insectPT <- st_sample(mitig, size = 1, type = "random", by_polygon = FALSE)
        insectPT <- st_transform(insectPT, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
        insectPT <- st_buffer(insectPT, 1)
        mitig <- st_transform(insectPT, crs = "+proj=longlat +datum=WGS84")
      }
      
      insect_sf <- st_geometry(mitig)
      insect_sf <- st_sf(geometry = insect_sf)
    }
    
    if (length(values$drawn_insecticide_polygons) > 0) {
      drawn_sf <- do.call(rbind, values$drawn_insecticide_polygons)
      drawn_sf <- st_transform(drawn_sf, crs = "+proj=longlat +datum=WGS84 +no_defs")
      
      if (!is.null(insect_sf)) {
        combined_geom <- st_union(st_geometry(insect_sf), st_geometry(drawn_sf))
        insect_sf <- st_sf(geometry = combined_geom)
      } else {
        insect_sf <- drawn_sf
      }
    }
    
    if (!is.null(insect_sf)) {
      values$insecticide_areas <- insect_sf
    }
    return(insect_sf)
  })
  
  # ============================================================================
  # MAIN OVERLAY MAP
  # ============================================================================
  
  output$main_overlay_map <- renderLeaflet({
    map_update_trigger()
    req(current_colonies())
    
    col_data <- current_colonies()
    
    map <- leaflet(col_data) %>%
      addTiles() %>%
      addPolygons(
        data = col_data,
        group = "Current Colonies",
        color = "black",
        weight = 2,
        fillColor = "brown",
        fillOpacity = 0.6,
        popup = "Current Year Colonies"
      )
    
    # Show uploaded poison if from current year
    if (!is.null(input$poison) &&
        !is.null(values$poison_upload_year) &&
        values$poison_upload_year == values$current_year) {
      pois <- poison_areas()
      if (!is.null(pois) && nrow(pois) > 0) {
        shpdf <- input$poison
        tempdirname <- dirname(shpdf$datapath[1])
        for (i in 1:nrow(shpdf)) {
          new_path <- paste0(tempdirname, "/", shpdf$name[i])
          if (!file.exists(new_path)) {
            file.copy(shpdf$datapath[i], new_path, overwrite = TRUE)
          }
        }
        uploaded_pois <- st_read(paste(tempdirname, shpdf$name[grep(pattern = "*.shp$", shpdf$name)], sep = "/"), quiet = TRUE)
        uploaded_pois <- st_transform(uploaded_pois, crs = "+proj=longlat +datum=WGS84")
        
        map <- map %>%
          addPolygons(
            data = uploaded_pois,
            group = "Pdog Control Areas",
            color = "darkred",
            weight = 2,
            fillColor = "red",
            fillOpacity = 0.6,
            popup = "Uploaded Pdog Control Areas"
          )
      }
    }
    
    # Show uploaded insecticide if from current year
    if (!is.null(input$insect) &&
        !is.null(values$insecticide_upload_year) &&
        values$insecticide_upload_year == values$current_year) {
      insect <- insecticide_areas()
      if (!is.null(insect) && nrow(insect) > 0) {
        shpdf <- input$insect
        tempdirname <- dirname(shpdf$datapath[1])
        for (i in 1:nrow(shpdf)) {
          new_path <- paste0(tempdirname, "/", shpdf$name[i])
          if (!file.exists(new_path)) {
            file.copy(shpdf$datapath[i], new_path, overwrite = TRUE)
          }
        }
        uploaded_insect <- st_read(paste(tempdirname, shpdf$name[grep(pattern = "*.shp$", shpdf$name)], sep = "/"), quiet = TRUE)
        uploaded_insect <- st_transform(uploaded_insect, crs = "+proj=longlat +datum=WGS84")
        
        map <- map %>%
          addPolygons(
            data = uploaded_insect,
            group = "Plague Treatment Areas",
            color = "darkblue",
            weight = 2,
            fillColor = "blue",
            fillOpacity = 0.6,
            popup = "Uploaded Plague Treatment Areas"
          )
      }
    }
    
    # Show drawn poison
    if (length(values$drawn_poison_polygons) > 0) {
      drawn_poison <- do.call(rbind, values$drawn_poison_polygons)
      drawn_poison <- st_transform(drawn_poison, crs = 4326)
      
      map <- map %>%
        addPolygons(
          data = drawn_poison,
          color = "#FF6666",
          weight = 2,
          fillColor = "#FFB3B3",
          fillOpacity = 0.5,
          popup = "Drawn Poisoned Area",
          group = "Drawn Poison"
        )
    }
    
    # Show drawn insecticide
    if (length(values$drawn_insecticide_polygons) > 0) {
      drawn_insect <- do.call(rbind, values$drawn_insecticide_polygons)
      drawn_insect <- st_transform(drawn_insect, crs = 4326)
      
      map <- map %>%
        addPolygons(
          data = drawn_insect,
          color = "#6666FF",
          weight = 2,
          fillColor = "#B3B3FF",
          fillOpacity = 0.5,
          popup = "Drawn Insecticide Area",
          group = "Drawn Insecticide"
        )
    }
    
    current_mode <- input$drawMode
    
    if (current_mode == "poison") {
      map <- map %>%
        addDrawToolbar(
          targetGroup = "Drawn Poison",
          polylineOptions = FALSE,
          circleOptions = FALSE,
          rectangleOptions = FALSE,
          markerOptions = FALSE,
          circleMarkerOptions = FALSE,
          polygonOptions = drawPolygonOptions(
            shapeOptions = drawShapeOptions(
              fillColor = "#FFB3B3",
              color = "#FF6666",
              weight = 2,
              fillOpacity = 0.5
            )
          ),
          editOptions = editToolbarOptions(
            selectedPathOptions = selectedPathOptions()
          )
        )
    } else {
      map <- map %>%
        addDrawToolbar(
          targetGroup = "Drawn Insecticide",
          polylineOptions = FALSE,
          circleOptions = FALSE,
          rectangleOptions = FALSE,
          markerOptions = FALSE,
          circleMarkerOptions = FALSE,
          polygonOptions = drawPolygonOptions(
            shapeOptions = drawShapeOptions(
              fillColor = "#B3B3FF",
              color = "#6666FF",
              weight = 2,
              fillOpacity = 0.5
            )
          ),
          editOptions = editToolbarOptions(
            selectedPathOptions = selectedPathOptions()
          )
        )
    }
    
    overlay_groups <- c("Current Colonies")
    if (!is.null(input$poison)) {
      overlay_groups <- c(overlay_groups, "Pdog Control Areas")
    }
    if (!is.null(input$insect)) {
      overlay_groups <- c(overlay_groups, "Plague Treatment Areas")
    }
    
    map <- map %>%
      addLayersControl(
        overlayGroups = overlay_groups,
        options = layersControlOptions(collapsed = FALSE)
      )
    
    return(map)
  })
  
  # ============================================================================
  # PROCESS COLONIES AND MODEL PREDICTIONS
  # ============================================================================
  
  process_colonies <- function() {
    req(current_colonies())
    
    mapCol <- current_colonies()
    mapCol <- st_transform(mapCol, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
    mapCol <- st_buffer(mapCol, 0)
    
    mapPois <- poison_areas()
    if (!is.null(mapPois) && nrow(mapPois) > 0) {
      mapPois <- st_transform(mapPois, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      mapPois <- st_buffer(mapPois, 0)
      mapPois <- st_union(mapPois)
    }
    
    mapInsect <- insecticide_areas()
    if (!is.null(mapInsect) && nrow(mapInsect) > 0) {
      mapInsect <- st_transform(mapInsect, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      mapInsect <- st_buffer(mapInsect, 0)
      mapInsect <- st_union(mapInsect)
    }
    
    mapCol <- st_union(mapCol)
    
    colsub <- mapCol
    if (!is.null(mapPois)) {
      colsub <- st_difference(colsub, mapPois)
      colsub <- st_union(colsub)
    }
    
    if (!is.null(mapInsect)) {
      colsub <- st_difference(colsub, mapInsect)
    }
    
    colsub <- st_transform(colsub, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
    colsub <- st_buffer(colsub, 0)
    colsub <- st_transform(colsub, crs = "+proj=longlat +datum=WGS84")
    
    values$processed_colonies <- colsub
    
    return(colsub)
  }
  
  # RASTER-BASED Colony predictions - PRODUCTION VERSION
  run_colony_predictions <- function() {
    showNotification("Running colony prediction model...", type = "message", duration = 3)
    
    tryCatch({
      # Use carried-over prior colonies if available
      if (!is.null(values$carried_colonies_prior)) {
        colpred_t1 <- values$carried_colonies_prior
        colpred_t1 <- st_transform(colpred_t1, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
        colpred_t1 <- st_buffer(colpred_t1, 0)
        cat("Using carried colonies from previous year\n")
      } else {
        # Load from uploaded file
        req(input$yrprior)
        
        shpdf2 <- input$yrprior
        tempdirname2 <- dirname(shpdf2$datapath[1])
        
        for (i in 1:nrow(shpdf2)) {
          new_path <- paste0(tempdirname2, "/", shpdf2$name[i])
          if (!file.exists(new_path)) {
            file.copy(shpdf2$datapath[i], new_path, overwrite = TRUE)
          }
        }
        
        shp_file <- paste(tempdirname2, shpdf2$name[grep(pattern = "*.shp$", shpdf2$name)], sep = "/")
        if (!file.exists(shp_file)) {
          shp_file <- shpdf2$datapath[grep(pattern = "*.shp$", shpdf2$name)]
        }
        
        colpred_t1 <- st_read(shp_file)
        colpred_t1 <- st_transform(colpred_t1, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
        colpred_t1 <- st_buffer(colpred_t1, 0)
      }
      
      mapSUB <- values$processed_colonies
      mapSUB <- st_transform(mapSUB, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      mapSUB <- st_buffer(mapSUB, 0)
      mapSUB <- st_as_sf(mapSUB)
      
      r_template <- create_raster_template(mapSUB, resolution = 100)
      
      colpred_raster <- rasterize(colpred_t1, r_template)
      colpred_raster[is.na(colpred_raster)] <- 0
      values(colpred_raster) <- ifelse(values(colpred_raster) > 0, 1, 0)
      
      tempOG <- current_colonies()
      tempOG <- st_transform(tempOG, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      tempOG <- st_buffer(tempOG, 0)
      tempOG <- st_as_sf(tempOG)
      
      tempPois <- poison_areas()
      if (!is.null(tempPois) && nrow(tempPois) > 0) {
        tempPois <- st_transform(tempPois, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
        tempPois <- st_buffer(tempPois, 0)
        tempPois <- st_as_sf(tempPois)
        
        tempOG <- st_union(tempOG)
        tempPois <- st_union(tempPois)
        
        mapSUBpois <- st_difference(tempOG, tempPois)
      } else {
        mapSUBpois <- st_union(tempOG)
      }
      
      mapSUBpois <- st_transform(mapSUBpois, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      mapSUBpois <- st_buffer(mapSUBpois, 0)
      mapSUBpois <- st_as_sf(mapSUBpois)
      
      mapSUBpois_raster <- rasterize(mapSUBpois, r_template)
      mapSUBpois_raster[is.na(mapSUBpois_raster)] <- 0
      values(mapSUBpois_raster) <- ifelse(values(mapSUBpois_raster) > 0, 1, 0)
      
      mapSUB_raster <- rasterize(mapSUB, r_template)
      mapSUB_raster[is.na(mapSUB_raster)] <- 0
      values(mapSUB_raster) <- ifelse(values(mapSUB_raster) > 0, 1, 0)
      
      tryCatch({
        dz <- data.frame(
          year = c(1, 2),
          mnn = c(as.data.frame(lsm_c_enn_mn(colpred_raster))$value[2],
                  as.data.frame(lsm_c_enn_mn(mapSUB_raster))$value[2]),
          cohes = c(as.data.frame(lsm_c_cohesion(colpred_raster))$value[2],
                    as.data.frame(lsm_c_cohesion(mapSUB_raster))$value[2]),
          contag = c(as.data.frame(lsm_l_contag(colpred_raster))$value[1],
                     as.data.frame(lsm_l_contag(mapSUBpois_raster))$value[1])
        )
      }, error = function(e) {
        dz <- data.frame(
          year = c(1, 2),
          mnn = c(1000, 1000),
          cohes = c(50, 50),
          contag = c(80, 80)
        )
      })
      
      dz$mnn[nrow(dz)] <- 1178.6
      
      dfc <- mapSUBpois_raster
      values(dfc) <- ifelse(values(dfc) > 0.5, 1, NA)
      dfc0 <- distance(dfc)
      
      mapOG <- current_colonies()
      mapOG <- st_transform(mapOG, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      mapOG <- st_buffer(mapOG, 0)
      mapOG <- st_as_sf(mapOG)
      
      mapOG_raster <- rasterize(mapOG, r_template)
      mapOG_raster[is.na(mapOG_raster)] <- 0
      values(mapOG_raster) <- ifelse(values(mapOG_raster) > 0, 1, 0)
      
      dp_raster <- calculate_distance_to_plague(mapOG_raster, colpred_raster)
      
      wyType <- input$WholeYear
      wy_var <- switch(wyType, "Wet" = "wywet", "Dry" = "wydry", "wyavg")
      
      wsType <- input$WinterSpring
      ws_var <- switch(wsType, "Wet" = "wswet", "Dry" = "wsdry", "wsavg")
      
      zztmaxType <- input$ZZTEMP
      zztmax_var <- switch(zztmaxType, "Warm" = "zztmaxmax", "Cold" = "zztmaxmin", "zztmaxavg")
      
      wssfType <- input$WSSF
      zwssf_var <- switch(wssfType, "Dry Fall Wet Spring" = "zwssfdrywet", 
                          "Wet Fall Dry Spring" = "zwssfwetdry", "zwssfavg")
      
      pixel_coords <- coordinates(r_template)
      climate_data <- extract_climate_from_matrix(pixel_coords, climate_matrix)
      
      growth_pred_data <- data.frame(
        x = climate_data$x,
        y = climate_data$y,
        dfc = as.numeric(values(dfc0)),
        hotr = climate_data$hotr,
        contagion = dz$contag[nrow(dz)],
        stringsAsFactors = FALSE
      )
      
      growth_pred_data$ws <- as.numeric(climate_data[[ws_var]])
      growth_pred_data$wy <- as.numeric(climate_data[[wy_var]])
      
      growth_predictions_link <- predict(growth_model, newdata = growth_pred_data, 
                                         type = "link", re.form = NA)
      
      if(input$region != "population_avg") {
        growth_site_re <- ranef(growth_model)$cond$site[input$region, "(Intercept)"]
        growth_predictions_link <- growth_predictions_link + growth_site_re
      }
      
      growth_predictions <- plogis(growth_predictions_link)
      
      Coloniz <- r_template
      values(Coloniz) <- growth_predictions
      
      Coloniz_masked <- Coloniz * (1 - mapSUBpois_raster)
      Coloniz_masked[is.na(Coloniz_masked)] <- 0
      
      plague_pred_data <- data.frame(
        x = climate_data$x,
        y = climate_data$y,
        mnn = dz$mnn[nrow(dz)],
        mnn_sq = dz$mnn[nrow(dz)]^2,
        cohesion = dz$cohes[nrow(dz)],
        dp = as.numeric(values(dp_raster)),
        dp_sq = as.numeric(values(dp_raster))^2,
        psize = 1,
        stringsAsFactors = FALSE
      )
      
      plague_pred_data$zwssf <- as.numeric(climate_data[[zwssf_var]])
      plague_pred_data$zztmax <- as.numeric(climate_data[[zztmax_var]])
      plague_pred_data$zztmax_sq <- as.numeric(climate_data[[zztmax_var]])^2
      
      plague_predictions_link <- predict(plague_model, newdata = plague_pred_data, 
                                         type = "link", re.form = NA)
      
      if(input$region != "population_avg") {
        plague_site_re <- ranef(plague_model)$cond$site[input$region, "(Intercept)"]
        plague_predictions_link <- plague_predictions_link + plague_site_re
      }
      
      plague_predictions <- plogis(plague_predictions_link)
      
      Persist <- r_template
      values(Persist) <- 1 - plague_predictions
      
      Persist_masked <- Persist * mapSUBpois_raster
      Persist_masked[is.na(Persist_masked)] <- 0
      
      tempInsect2 <- insecticide_areas()
      
      if (!is.null(tempInsect2) && nrow(tempInsect2) > 0) {
        tempOG2 <- current_colonies()
        tempOG2 <- st_transform(tempOG2, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
        tempOG2 <- st_buffer(tempOG2, 0)
        tempOG2 <- st_as_sf(tempOG2)
        
        tempInsect2 <- st_transform(tempInsect2, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
        tempInsect2 <- st_buffer(tempInsect2, 0)
        tempInsect2 <- st_as_sf(tempInsect2)
        
        tempOG2 <- st_union(tempOG2)
        tempInsect2 <- st_union(tempInsect2)
        
        insect_in_cols <- tryCatch({
          st_intersection(tempInsect2, tempOG2)
        }, error = function(e) {
          NULL
        })
        
        if(!is.null(insect_in_cols) && length(insect_in_cols) > 0) {
          insect_in_cols <- st_as_sf(insect_in_cols)
          
          insect_raster <- rasterize(insect_in_cols, r_template)
          insect_raster[is.na(insect_raster)] <- 0
          values(insect_raster) <- ifelse(values(insect_raster) > 0, 1, 0)
          
          Persist_masked <- Persist_masked + insect_raster
        }
      }
      
      persist_threshold <- 1 - input$thresholdp
      
      values(Persist_masked) <- ifelse(values(Persist_masked) > persist_threshold, 1, 0)
      values(Coloniz_masked) <- ifelse(values(Coloniz_masked) > input$threshold, 1, 0)
      
      Persist_masked[is.na(Persist_masked)] <- 0
      values(Coloniz_masked)[is.na(values(Coloniz_masked))] <- 0
      
      nxt <- Coloniz_masked + Persist_masked
      values(nxt) <- ifelse(values(nxt) > 0.5, 1, NA)
      
      crs(nxt) <- "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs"
      
      nxt_terra <- as(nxt, "SpatRaster")
      nxt_clumped <- patches(nxt_terra, directions = 8, zeroAsNA = TRUE)
      
      colpred <- as.polygons(nxt_clumped, dissolve = TRUE)
      colpred <- st_as_sf(colpred)
      
      if(nrow(colpred) == 0) {
        return(NULL)
      }
      
      colpred_utm <- st_transform(colpred, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      climate_data_colonies <- extract_climate_data(colpred_utm, climate_matrix)
      
      centroids <- st_centroid(colpred_utm)
      centroid_coords <- st_coordinates(centroids)
      dfc_values <- raster::extract(dfc0, centroid_coords)
      
      pred_data_uncertainty <- data.frame(
        x = climate_data_colonies$x,
        y = climate_data_colonies$y,
        dfc = ifelse(is.na(dfc_values), mean(values(dfc0), na.rm = TRUE), dfc_values),
        hotr = climate_data_colonies$hotr,
        contagion = dz$contag[nrow(dz)],
        stringsAsFactors = FALSE
      )
      
      pred_data_uncertainty$ws <- as.numeric(climate_data_colonies[[ws_var]])
      pred_data_uncertainty$wy <- as.numeric(climate_data_colonies[[wy_var]])
      
      predictions_link_uncertainty <- predict(growth_model, newdata = pred_data_uncertainty, 
                                              type = "link", re.form = NA)
      
      if(input$region != "population_avg") {
        growth_site_re_unc <- ranef(growth_model)$cond$site[input$region, "(Intercept)"]
        predictions_link_uncertainty <- predictions_link_uncertainty + growth_site_re_unc
      }
      
      colony_uncertainty <- calculate_uncertainty_delta(
        growth_model, 
        pred_data_uncertainty, 
        predictions_link_uncertainty
      )
      
      colony_predictions <- plogis(predictions_link_uncertainty)
      
      colpred$pred_prob <- round(colony_predictions * 100, 1)
      colpred$ci_lower <- round(colony_uncertainty$ci_lower * 100, 1)
      colpred$ci_upper <- round(colony_uncertainty$ci_upper * 100, 1)
      
      colpred <- st_buffer(colpred, 0)
      colpred <- st_transform(colpred, crs = "+proj=longlat +datum=WGS84")
      
      return(colpred)
      
    }, error = function(e) {
      cat("ERROR in colony predictions:", e$message, "\n")
      showNotification(paste("Error in colony predictions:", e$message), type = "error")
      return(NULL)
    })
  }
  
  # Plague predictions (visualization)
  run_plague_predictions <- function() {
    showNotification("Running plague prediction model for visualization...", type = "message", duration = 3)
    
    tryCatch({
      if (!is.null(values$carried_colonies_prior)) {
        plgpred_t1 <- values$carried_colonies_prior
        plgpred_t1 <- st_transform(plgpred_t1, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
        plgpred_t1 <- st_buffer(plgpred_t1, 0)
      } else {
        req(input$yrprior)
        
        shpdf2 <- input$yrprior
        tempdirname2 <- dirname(shpdf2$datapath[1])
        
        for (i in 1:nrow(shpdf2)) {
          new_path <- paste0(tempdirname2, "/", shpdf2$name[i])
          if (!file.exists(new_path)) {
            file.copy(shpdf2$datapath[i], new_path, overwrite = TRUE)
          }
        }
        
        shp_file <- paste(tempdirname2, shpdf2$name[grep(pattern = "*.shp$", shpdf2$name)], sep = "/")
        if (!file.exists(shp_file)) {
          shp_file <- shpdf2$datapath[grep(pattern = "*.shp$", shpdf2$name)]
        }
        
        plgpred_t1 <- st_read(shp_file)
        plgpred_t1 <- st_transform(plgpred_t1, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
        plgpred_t1 <- st_buffer(plgpred_t1, 0)
      }
      
      mapSUBZ <- values$processed_colonies
      mapSUBZ <- st_transform(mapSUBZ, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      mapSUBZ <- st_buffer(mapSUBZ, 0)
      mapSUBZ <- st_as_sf(mapSUBZ)
      
      r_template <- create_raster_template(mapSUBZ, resolution = 100)
      
      plgpred_raster <- rasterize(plgpred_t1, r_template)
      plgpred_raster[is.na(plgpred_raster)] <- 0
      values(plgpred_raster) <- ifelse(values(plgpred_raster) > 0, 1, 0)
      
      mapSUBZ_raster <- rasterize(mapSUBZ, r_template)
      mapSUBZ_raster[is.na(mapSUBZ_raster)] <- 0
      values(mapSUBZ_raster) <- ifelse(values(mapSUBZ_raster) > 0, 1, 0)
      
      tryCatch({
        dz <- data.frame(
          year = c(1, 2),
          mnn = c(as.data.frame(lsm_c_enn_mn(plgpred_raster))$value[2],
                  as.data.frame(lsm_c_enn_mn(mapSUBZ_raster))$value[2]),
          cohes = c(as.data.frame(lsm_c_cohesion(plgpred_raster))$value[2],
                    as.data.frame(lsm_c_cohesion(mapSUBZ_raster))$value[2])
        )
      }, error = function(e) {
        dz <- data.frame(
          year = c(1, 2),
          mnn = c(1000, 1000),
          cohes = c(50, 50)
        )
      })
      
      dz$mnn[nrow(dz)] <- 1178.6
      
      mapOGZ <- current_colonies()
      mapOGZ <- st_transform(mapOGZ, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      mapOGZ <- st_buffer(mapOGZ, 0)
      mapOGZ <- st_as_sf(mapOGZ)
      
      mapOGZ_raster <- rasterize(mapOGZ, r_template)
      mapOGZ_raster[is.na(mapOGZ_raster)] <- 0
      values(mapOGZ_raster) <- ifelse(values(mapOGZ_raster) > 0, 1, 0)
      
      dp_raster <- calculate_distance_to_plague(mapOGZ_raster, plgpred_raster)
      
      zztmaxType <- input$ZZTEMP
      zztmax_var <- switch(zztmaxType, "Warm" = "zztmaxmax", "Cold" = "zztmaxmin", "zztmaxavg")
      
      wssfType <- input$WSSF
      zwssf_var <- switch(wssfType, "Dry Fall Wet Spring" = "zwssfdrywet", 
                          "Wet Fall Dry Spring" = "zwssfwetdry", "zwssfavg")
      
      shape3 <- values$processed_colonies
      shape3 <- st_transform(shape3, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      shape3 <- st_buffer(shape3, 0)
      shape3 <- st_as_sf(shape3)
      shape3 <- st_cast(shape3, "POLYGON")
      
      pixel_coords <- coordinates(r_template)
      climate_data <- extract_climate_from_matrix(pixel_coords, climate_matrix)
      
      plague_pred_data <- data.frame(
        x = climate_data$x,
        y = climate_data$y,
        mnn = dz$mnn[nrow(dz)],
        mnn_sq = dz$mnn[nrow(dz)]^2,
        cohesion = dz$cohes[nrow(dz)],
        dp = as.numeric(values(dp_raster)),
        dp_sq = as.numeric(values(dp_raster))^2,
        psize = 1,
        stringsAsFactors = FALSE
      )
      
      plague_pred_data$zwssf <- as.numeric(climate_data[[zwssf_var]])
      plague_pred_data$zztmax <- as.numeric(climate_data[[zztmax_var]])
      plague_pred_data$zztmax_sq <- as.numeric(climate_data[[zztmax_var]])^2
      
      plague_predictions_link <- predict(plague_model, newdata = plague_pred_data, 
                                         type = "link", re.form = NA)
      
      if(input$region != "population_avg") {
        plague_site_re <- ranef(plague_model)$cond$site[input$region, "(Intercept)"]
        plague_predictions_link <- plague_predictions_link + plague_site_re
      }
      
      plague_predictions <- plogis(plague_predictions_link)
      
      ext_raster <- r_template
      values(ext_raster) <- plague_predictions
      
      extmask <- mask(ext_raster, shape3, updatevalue = 0)
      
      r.vals <- raster::extract(extmask, shape3)
      r.mean <- lapply(r.vals, FUN=mean, na.rm=TRUE)
      
      shape3$val <- unlist(r.mean)
      shape3$val[is.na(shape3$val)] <- 0
      shape3$val <- round(shape3$val * 100, 2)
      
      shape3$col <- "#4B0000"
      shape3$col[shape3$val < 90] <- "#8B0000"
      shape3$col[shape3$val < 80] <- "#B22222"
      shape3$col[shape3$val < 70] <- "#DC143C"
      shape3$col[shape3$val < 60] <- "#FF4500"
      shape3$col[shape3$val < 50] <- "#FF6347"
      shape3$col[shape3$val < 40] <- "#FFA500"
      shape3$col[shape3$val < 30] <- "#FFD700"
      shape3$col[shape3$val < 20] <- "#32CD32"
      shape3$col[shape3$val < 10] <- "#2E8B57"
      
      shape3$new <- paste0(shape3$val, "% plague risk")
      shape3$new[shape3$val == 0] <- "No plague risk"
      
      plgpred <- st_transform(shape3, crs = "+proj=longlat +datum=WGS84")
      
      return(plgpred)
      
    }, error = function(e) {
      cat("ERROR in plague predictions:", e$message, "\n")
      showNotification(paste("Error in plague predictions:", e$message), type = "error")
      return(NULL)
    })
  }
  
  # ============================================================================
  # MAP OUTPUTS
  # ============================================================================
  
  output$prediction_results <- renderLeaflet({
    map_trigger()
    
    isolate({
      if (!is.null(values$predicted_colonies) && 
          inherits(values$predicted_colonies, "sf") &&
          nrow(values$predicted_colonies) > 0 &&
          values$models_completed) {
        
        tryCatch({
          popup_text <- character(nrow(values$predicted_colonies))
          
          for(i in 1:nrow(values$predicted_colonies)) {
            area_m2 <- as.numeric(st_area(values$predicted_colonies[i,]))
            area_acres <- area_m2 / 4046.86
            
            if("pred_prob" %in% names(values$predicted_colonies) && 
               !is.na(values$predicted_colonies$pred_prob[i])) {
              popup_text[i] <- paste0(
                "<b>Predicted Colony</b><br>",
                "Area: ", round(area_acres, 2), " acres<br>",
                "Growth Probability: ", values$predicted_colonies$pred_prob[i], "%<br>",
                "95% CI: [", values$predicted_colonies$ci_lower[i], "%, ",
                values$predicted_colonies$ci_upper[i], "%]"
              )
            } else {
              popup_text[i] <- paste0(
                "<b>Predicted Colony</b><br>",
                "Area: ", round(area_acres, 2), " acres"
              )
            }
          }
          
          leaflet(values$predicted_colonies) %>%
            addTiles() %>%
            addPolygons(
              data = values$predicted_colonies,
              color = "darkgreen",
              weight = 2,
              fillColor = "lightgreen",
              fillOpacity = 0.7,
              popup = popup_text,
              group = "Predicted Colonies"
            ) %>%
            addLayersControl(
              overlayGroups = "Predicted Colonies",
              options = layersControlOptions(collapsed = FALSE)
            )
        }, error = function(e) {
          cat("Error in prediction_results:", e$message, "\n")
          leaflet() %>%
            addTiles() %>%
            setView(lng = -104.5, lat = 42.0, zoom = 8)
        })
      } else {
        leaflet() %>%
          addTiles() %>%
          setView(lng = -104.5, lat = 42.0, zoom = 8) %>%
          addControl(HTML("<div style='background: white; padding: 10px; border-radius: 5px;'>Click 'Run Models' to generate colony predictions</div>"), 
                     position = "topright")
      }
    })
  })
  
  output$plgpred <- renderLeaflet({
    map_trigger()
    
    isolate({
      if (!is.null(values$predicted_plague) && 
          inherits(values$predicted_plague, "sf") &&
          nrow(values$predicted_plague) > 0 &&
          "col" %in% names(values$predicted_plague) &&
          "new" %in% names(values$predicted_plague) &&
          values$models_completed) {
        
        tryCatch({
          leaflet(values$predicted_plague) %>%
            addTiles() %>%
            addPolygons(
              data = values$predicted_plague,
              color = "black",
              weight = 2,
              fillColor = values$predicted_plague$col,
              fillOpacity = 0.7,
              popup = values$predicted_plague$new,
              group = "Plague Probability"
            ) %>%
            addLayersControl(
              overlayGroups = "Plague Probability",
              options = layersControlOptions(collapsed = FALSE)
            )
        }, error = function(e) {
          cat("Error in plgpred:", e$message, "\n")
          leaflet() %>%
            addTiles() %>%
            setView(lng = -104.5, lat = 42.0, zoom = 8)
        })
      } else {
        leaflet() %>%
          addTiles() %>%
          setView(lng = -104.5, lat = 42.0, zoom = 8) %>%
          addControl(HTML("<div style='background: white; padding: 10px; border-radius: 5px;'>Click 'Run Models' to generate plague probability</div>"), 
                     position = "topright")
      }
    })
  })
  
  # ============================================================================
  # RUN MODELS BUTTON
  # ============================================================================
  
  observeEvent(input$runModels, {
    
    if (values$models_running) {
      return()
    }
    
    # Validation
    if (is.null(values$carried_colonies_current)) {
      if (is.null(input$colonies)) {
        showNotification("Please upload current year colony shapefile", type = "error", duration = 5)
        return()
      }
    }
    
    if (is.null(values$carried_colonies_prior)) {
      if (is.null(input$yrprior)) {
        showNotification("Please upload prior year colony shapefile", type = "error", duration = 5)
        return()
      }
    }
    
    values$models_running <- TRUE
    values$models_completed <- FALSE
    
    runjs("
      var btn = document.getElementById('runModels');
      btn.innerHTML = 'Models Running...';
      btn.className = 'btn btn-danger btn-lg action-button';
      btn.disabled = true;
    ")
    
    showNotification("Starting glmmTMB models...", type = "message", duration = 3)
    
    process_colonies()
    values$predicted_colonies <- run_colony_predictions()
    values$predicted_plague <- run_plague_predictions()
    
    if (!is.null(values$predicted_colonies) && nrow(values$predicted_colonies) > 0) {
      areamap <- current_colonies()
      areamap <- st_transform(areamap, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      areamap <- as.numeric(st_area(st_union(areamap))) / 4046.86
      
      areamap2 <- values$predicted_colonies
      areamap2 <- st_transform(areamap2, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      areamap2 <- as.numeric(st_area(st_union(areamap2))) / 4046.86
      
      area_change <- round(areamap2 - areamap, 2)
      area_change_formatted <- ifelse(area_change > 0, paste0("+", area_change), as.character(area_change))
      
      num_current <- tryCatch({
        nrow(suppressWarnings(st_cast(current_colonies(), "POLYGON")))
      }, error = function(e) {
        length(st_geometry(current_colonies()))
      })
      
      num_predicted <- tryCatch({
        nrow(suppressWarnings(st_cast(values$predicted_colonies, "POLYGON")))
      }, error = function(e) {
        length(st_geometry(values$predicted_colonies))
      })
      
      if (!is.null(values$predicted_plague) && "val" %in% names(values$predicted_plague)) {
        valid_vals <- values$predicted_plague$val[!is.na(values$predicted_plague$val)]
        if (length(valid_vals) > 0) {
          max_plague <- round(max(valid_vals), 2)
          mean_plague <- round(mean(valid_vals), 2)
          
          next_year_label <- get_year_label(1)
          
          values$last_results <- data.frame(
            Metric = c("Current Colonies", paste("Predicted Colonies (", next_year_label, ")", sep=""), 
                       "Area Change (acres)", "Max Plague Probability (%)", "Mean Plague Probability (%)"),
            Value = c(num_current, num_predicted, area_change_formatted, max_plague, mean_plague),
            stringsAsFactors = FALSE
          )
        } else {
          next_year_label <- get_year_label(1)
          
          values$last_results <- data.frame(
            Metric = c("Current Colonies", paste("Predicted Colonies (", next_year_label, ")", sep=""), 
                       "Area Change (acres)", "Plague Status"),
            Value = c(num_current, num_predicted, area_change_formatted, "No plague data"),
            stringsAsFactors = FALSE
          )
        }
      } else {
        next_year_label <- get_year_label(1)
        
        values$last_results <- data.frame(
          Metric = c("Current Colonies", paste("Predicted Colonies (", next_year_label, ")", sep=""),"Area Change (acres)", "Plague Status"),
          Value = c(num_current, num_predicted, area_change_formatted, "Plague data unavailable"),
          stringsAsFactors = FALSE
        )
      }
    } else {
      values$last_results <- data.frame(
        Metric = c("Status"),
        Value = c("No colonies predicted - try lowering thresholds"),
        stringsAsFactors = FALSE
      )
    }
    
    values$models_completed <- TRUE
    values$models_running <- FALSE
    
    runjs("
      var btn = document.getElementById('runModels');
      btn.innerHTML = 'Run Models';
      btn.className = 'btn btn-success btn-lg action-button';
      btn.disabled = false;
    ")
    
    showNotification("Models completed!", type = "message", duration = 3)
    
    runjs("$('#costBenefitBox').removeClass('disabled');")
    
    # Auto-populate economic inputs
    poison_area_acres <- 0
    if (!is.null(poison_areas()) && nrow(poison_areas()) > 0) {
      poison_utm <- st_transform(poison_areas(), crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      poison_area_acres <- round(as.numeric(st_area(st_union(poison_utm))) / 4046.86, 2)
    }
    
    pdogs_area_acres <- 0
    if (!is.null(values$predicted_colonies) && nrow(values$predicted_colonies) > 0) {
      pdogs_utm <- st_transform(values$predicted_colonies, crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
      pdogs_area_acres <- round(as.numeric(st_area(st_union(pdogs_utm))) / 4046.86, 2)
    }
    
    updateNumericInput(session, "acres_poisoned_A", value = poison_area_acres)
    updateNumericInput(session, "acres_pdogs_A", value = pdogs_area_acres)
    
    updateNumericInput(session, "acres_poisoned_B", value = poison_area_acres)
    updateNumericInput(session, "acres_pdogs_B", value = pdogs_area_acres)
    
    updateNumericInput(session, "acres_poisoned_C", value = poison_area_acres)
    updateNumericInput(session, "acres_pdogs_C", value = pdogs_area_acres)
    
    updateNumericInput(session, "acres_poisoned_D", value = poison_area_acres)
    updateNumericInput(session, "acres_pdogs_D", value = pdogs_area_acres)
    
    map_trigger(map_trigger() + 1)
  })
  
  # ============================================================================
  # PLOTS AND TABLES
  # ============================================================================
  
  output$plague_distribution <- renderPlot({
    if (values$models_completed && !is.null(values$predicted_plague) && 
        nrow(values$predicted_plague) > 0 && "val" %in% names(values$predicted_plague)) {
      
      plague_probs <- values$predicted_plague$val[!is.na(values$predicted_plague$val)]
      
      if (length(plague_probs) > 0) {
        bins <- seq(0, 100, by = 10)
        bin_labels <- paste0(bins[-length(bins)], "-", bins[-1], "%")
        
        plague_binned <- cut(plague_probs, breaks = bins, labels = bin_labels, include.lowest = TRUE)
        bin_counts <- table(plague_binned)
        
        plot_data <- data.frame(
          Probability_Range = names(bin_counts),
          Count = as.numeric(bin_counts)
        )
        
        max_count <- max(plot_data$Count)
        y_limit <- max_count * 1.15
        
        par(mar = c(6, 4, 2, 2))
        
        barplot(plot_data$Count, 
                names.arg = plot_data$Probability_Range,
                main = "",
                xlab = "",
                ylab = "Number of Colonies",
                ylim = c(0, y_limit),
                col = c("#2E8B57", "#32CD32", "#FFD700", "#FFA500", "#FF6347", 
                        "#FF4500", "#DC143C", "#B22222", "#8B0000", "#4B0000"),
                las = 2,
                cex.names = 1.0,
                cex.axis = 1.1,
                cex.lab = 1.2,
                border = "white")
        
        mtext("Plague Probability", side = 1, line = 4.5, cex = 1.2)
        grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
        
        text(x = seq_along(plot_data$Count), 
             y = plot_data$Count + max_count * 0.03, 
             labels = plot_data$Count, 
             cex = 1.0, 
             pos = 3)
        
      } else {
        plot.new()
        text(0.5, 0.5, "No valid plague probability data", 
             cex = 1, col = "gray50", font = 2)
      }
      
    } else {
      plot.new()
      text(0.5, 0.5, "Run models to see plague\nprobability distribution", 
           cex = 1, col = "gray50", font = 2)
    }
  })
  
  output$plague_area_distribution <- renderPlot({
    if (values$models_completed && !is.null(values$predicted_plague) && 
        nrow(values$predicted_plague) > 0 && "val" %in% names(values$predicted_plague)) {
      
      valid_indices <- !is.na(values$predicted_plague$val)
      plague_probs <- values$predicted_plague$val[valid_indices]
      
      if (length(plague_probs) > 0) {
        plague_utm <- st_transform(values$predicted_plague[valid_indices, ], 
                                   crs = "+proj=utm +zone=13 +datum=NAD83 +units=m +no_defs")
        
        areas_acres <- as.numeric(st_area(plague_utm)) / 4046.86
        
        bins <- seq(0, 100, by = 10)
        bin_labels <- paste0(bins[-length(bins)], "-", bins[-1], "%")
        
        plague_binned <- cut(plague_probs, breaks = bins, labels = bin_labels, include.lowest = TRUE)
        
        bin_areas <- tapply(areas_acres, plague_binned, sum, na.rm = TRUE)
        bin_areas[is.na(bin_areas)] <- 0
        
        plot_data <- data.frame(
          Probability_Range = names(bin_areas),
          Area = as.numeric(bin_areas)
        )
        
        max_area <- max(plot_data$Area)
        y_limit <- max_area * 1.15
        
        par(mar = c(6, 4, 2, 2))
        
        barplot(plot_data$Area, 
                names.arg = plot_data$Probability_Range,
                main = "",
                xlab = "", 
                ylab = "Total Area (acres)",
                ylim = c(0, y_limit),
                col = c("#2E8B57", "#32CD32", "#FFD700", "#FFA500", "#FF6347", 
                        "#FF4500", "#DC143C", "#B22222", "#8B0000", "#4B0000"),
                las = 2,
                cex.names = 1.0,
                cex.axis = 1.1,
                cex.lab = 1.2,
                border = "white")
        
        mtext("Plague Probability", side = 1, line = 4.5, cex = 1.2)
        grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
        
        text(x = seq_along(plot_data$Area), 
             y = plot_data$Area + max_area * 0.03, 
             labels = round(plot_data$Area, 1), 
             cex = 1.0, 
             pos = 3)
        
      } else {
        plot.new()
        text(0.5, 0.5, "No valid plague probability data", 
             cex = 1, col = "gray50", font = 2)
      }
      
    } else {
      plot.new()
      text(0.5, 0.5, "Run models to see plague\narea distribution", 
           cex = 1, col = "gray50", font = 2)
    }
  })
  
  output$irish <- renderTable({
    if (values$models_completed && !is.null(values$last_results)) {
      return(values$last_results)
    } else {
      data.frame(
        Metric = c("Status"),
        Value = c("Click 'Run Models' to generate statistics"),
        stringsAsFactors = FALSE
      )
    }
  })
  
  # ============================================================================
  # DOWNLOAD HANDLER
  # ============================================================================
  
  observe({
    if (values$models_completed && !is.null(values$predicted_colonies)) {
      runjs("
        var btn = document.getElementById('ShapeExport');
        btn.disabled = false;
        btn.className = 'btn btn-primary btn-lg';
        btn.style.width = '100%';
        btn.onclick = function() { 
          Shiny.setInputValue('download_trigger', Math.random()); 
        };
      ")
    } else {
      runjs("
        var btn = document.getElementById('ShapeExport');
        btn.disabled = true;
        btn.className = 'btn btn-secondary btn-lg';
        btn.style.width = '100%';
        btn.onclick = null;
      ")
    }
  })
  
  observeEvent(input$download_trigger, {
    if (values$models_completed && !is.null(values$predicted_colonies)) {
      showNotification("Download functionality ready for implementation", type = "message", duration = 3)
    }
  })
  
  output$ShapeExport <- downloadHandler(
    filename = function() { "predicted_colonies.txt" },
    content = function(file) {
      writeLines("Download functionality ready", file)
    }
  )
  
  # ============================================================================
  # ECONOMIC CALCULATOR
  # ============================================================================
  
  observeEvent(input$calcCosts, {
    req(values$models_completed)
    
    # Scenario A
    result_A <- calculate_economic_scenario(
      acres_poisoned = input$acres_poisoned_A,
      acres_pdogs_after = input$acres_pdogs_A,
      ppa_value = NULL,
      cost_per_acre = input$cost_per_acre_A,
      scenario_type = "consumption"
    )
    values$economic_results$A <- result_A
    
    # Scenario B
    result_B <- calculate_economic_scenario(
      acres_poisoned = input$acres_poisoned_B,
      acres_pdogs_after = input$acres_pdogs_B,
      ppa_value = NULL,
      cost_per_acre = input$cost_per_acre_B,
      scenario_type = "no_increase"
    )
    values$economic_results$B <- result_B
    
    # Scenario C
    result_C <- calculate_economic_scenario(
      acres_poisoned = input$acres_poisoned_C,
      acres_pdogs_after = input$acres_pdogs_C,
      ppa_value = input$ppa_C,
      cost_per_acre = input$cost_per_acre_C,
      scenario_type = "half_forage"
    )
    values$economic_results$C <- result_C
    
    # Scenario D
    result_D <- calculate_economic_scenario(
      acres_poisoned = input$acres_poisoned_D,
      acres_pdogs_after = input$acres_pdogs_D,
      ppa_value = input$ppa_D,
      cost_per_acre = input$cost_per_acre_D,
      scenario_type = "no_forage"
    )
    values$economic_results$D <- result_D
    
    showNotification("Economic analysis calculated for all 4 scenarios!", type = "message", duration = 3)
  })
  
  output$forage_lost_A <- renderText({
    val <- round(values$economic_results$A$forage_lost, 0)
    sign <- ifelse(val >= 0, "+", "")
    paste0(sign, format(val, big.mark = ","))
  })
  output$total_cost_A <- renderText({
    paste0("$", format(round(values$economic_results$A$total_cost, 2), big.mark = ","))
  })
  
  output$forage_lost_B <- renderText({
    val <- round(values$economic_results$B$forage_lost, 0)
    sign <- ifelse(val >= 0, "+", "")
    paste0(sign, format(val, big.mark = ","))
  })
  output$total_cost_B <- renderText({
    paste0("$", format(round(values$economic_results$B$total_cost, 2), big.mark = ","))
  })
  
  output$forage_lost_C <- renderText({
    val <- round(values$economic_results$C$forage_lost, 0)
    sign <- ifelse(val >= 0, "+", "")
    paste0(sign, format(val, big.mark = ","))
  })
  output$total_cost_C <- renderText({
    paste0("$", format(round(values$economic_results$C$total_cost, 2), big.mark = ","))
  })
  
  output$forage_lost_D <- renderText({
    val <- round(values$economic_results$D$forage_lost, 0)
    sign <- ifelse(val >= 0, "+", "")
    paste0(sign, format(val, big.mark = ","))
  })
  output$total_cost_D <- renderText({
    paste0("$", format(round(values$economic_results$D$total_cost, 2), big.mark = ","))
  })
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 6) Run App -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

app <- shinyApp(ui = ui, server = server)
runApp(app)