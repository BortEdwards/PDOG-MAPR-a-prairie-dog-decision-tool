# PDOG MAPR - Prairie Dog Management and Planning Resource

[![Version](https://img.shields.io/badge/version-9.7-blue.svg)](https://github.com/yourusername/pdog-mapr)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)

A decision support tool for managing prairie dog ecosystems, forecasting colony dynamics, and evaluating management scenarios.

## Key Features

- **Colony Growth Predictions** using glmmTMB models with site-specific random effects
- **Plague Outbreak Forecasting** based on landscape metrics and climate data
- **Multi-Year Iterations** with automatic data carryover between years
- **Economic Calculator** with 4 scenarios for cost-benefit analysis
- **Interactive Maps** with polygon drawing for management planning
- **Climate Scenario Testing** to explore impacts of different weather conditions

## Quick Start

### Prerequisites

```r
# Required R packages
install.packages(c("shiny", "shinythemes", "shinyjs", "sf", "raster", 
                   "terra", "leaflet", "leaflet.extras", "glmmTMB", 
                   "lme4", "landscapemetrics", "tidyverse", "FNN"))
```

### Required Data Files

The application requires three data files (contact authors for access):
- `climate_complete_matrix_utm13.rds` - Climate data matrix (2000-2020)
- `growth_model.rds` - Fitted colony growth model
- `plague_model.rds` - Fitted plague prediction model

### Running the Application

```r
# Clone the repository
git clone https://github.com/yourusername/pdog-mapr.git
cd pdog-mapr

# Ensure data files are in the working directory
# Then run the app
source("prairie_dog_app_glm_stablev9_7.R")
```

## Documentation

- **[Web Documentation](https://bortedwards.github.io/PDOG-MAPR-a-prairie-dog-decision-tool/)** - Web-based documentation
- **[Full Documentation](docs/prairie_dog_app_documentation_v9_7.Rmd)** - Comprehensive R Markdown documentation
- **[Video Tutorial](https://youtu.be/jnRZ6dRA6mg)** - Step-by-step walkthrough

## Use Cases

1. **Forecast colony dynamics** under different climate scenarios
2. **Assess plague outbreak risk** across your prairie dog complex
3. **Evaluate management strategies** (colony control, plague mitigation)
4. **Perform cost-benefit analysis** of different management approaches
5. **Plan multi-year interventions** with iterative modeling

## Model Details

**Growth Model (glmmTMB)**
- Predicts colony expansion based on climate, distance from edge, and landscape metrics
- Default threshold: 0.90 (adjustable)

**Plague Model (glmmTMB)**
- Predicts plague outbreak probability based on colony connectivity and climate
- Default threshold: 0.75 (adjustable)

**Regional Calibration Sites:**
- Thunder Basin, Pawnee (E/W), Comanche (NW/SE), Cimarron, Kiowa, Rita Blanca, CMR

## Citation

If you use this tool in your research, please cite:

```
Barrile et al. (2023). A big data–model integration approach for predicting epizootics and population recovery in a keystone species. Ecol. Appl. 33, e2827
```

## Contact

- **Research Team**: [Ana.Davidson@colostate.edu]
- **Issues**: [GitHub Issues](https://github.com/yourusername/pdog-mapr/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/pdog-mapr/discussions)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This tool was developed through collaborative efforts of conservation partners, grassland managers, research scientists, statistical modelers, and software developers.

## Additional Resources

### Key Publications
- Augustine & Derner (2021) - *Journal of Wildlife Management*, 85(7):1332-1343
- Augustine et al. (2024) - *Rangeland Ecology and Management*
- Buehler et al. (2025) - *Rangeland Ecology and Management*, 99(2):66-76
- Crow et al. (2022) - *Rangeland Ecology and Management*, 85:56-65

### Data Resources
- [USDA WebSoil Survey](https://websoilsurvey.nrcs.usda.gov/)
- [Rangeland Analysis Platform](https://rangelands.app/)

---

**Version 9.7** | Last Updated: December 2025
