# Boundary Line Analysis of Soil Organic Carbon

This repository demonstrates separate boundary line models for clay
and electrical conductivity (EC) within the Summer 350-650 mm
rainfall region.

## Files

- BLA_example.R: complete R script.
- example_data.csv: input dataset.
- outputs/: generated models, predictions, summaries and figures.

## Requirements

Install R and the following packages:

```r
install.packages(c("MASS", "ggplot2", "dplyr"))
```

## How to run

1. Download this repository using Code > Download ZIP.
2. Extract the ZIP file.
3. Open RStudio.
4. Select File > New Project > Existing Directory.
5. Browse to the extracted folder containing BLA_example.R.
6. Click Create Project.
7. Open BLA_example.R and click Source.

Results are saved in the outputs folder.

## Input data

The CSV must contain:

- bom_class: rainfall-region label.
- OC: soil organic carbon (%).
- Clay: clay content (%).
- EC: electrical conductivity.

Use the same EC measurement unit as the supplied example data.
The EC unit should be documented here before publication.

## Method

Each predictor is analysed separately using:

- Classical Mahalanobis screening at the empirical 98th percentile.
- Two-dimensional kernel density estimation.
- Fifteen overlapping predictor windows.
- Exclusion of observations below each window's fifth density
  percentile from boundary-point selection.
- Selection of the maximum eligible SOC in each window.
- LOESS boundary fitting with span = 2.
