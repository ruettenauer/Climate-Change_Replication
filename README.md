
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Replication-Material: More Talk, No Action? The Link between Exposure to Extreme Weather Events, Climate Change Belief and Pro-Environmental Behaviour.

This repository provides all materials for replication of results in

> Rüttenauer T. (2023) “More Talk, No Action? The Link between Exposure
> to Extreme Weather Events, Climate Change Belief and Pro-Environmental
> Behaviour”, *European Societies*, Forthcoming.

Published manuscript available here:
<https://doi.org/10.1080/14616696.2023.2277281>

<!-- Preprint available here: https://doi.org/10.31235/osf.io/574uf -->

Date: 2023-06-06

## Set up

The code for replication of the results requires the following folders:
“01_Script”, “02_Data”, “03_Output”. All R Scripts are required in
folder “01_Script”, all data will be save in “02_Data”. Original input
data (such as Understanding Society, Flood records and MET office
Weather data are required in “00_Input_data”.)

To reproduce the results of the paper, the scripts need to be executed
in order. Note that data preparation and analyses are done separately
for floods and heatwaves. Results for floods are stored and common
figures are then produced in the heatwave analyse script.

The following packages are necessary for reproduction of main results:

``` r
library(CBPS)
library(GISTools)
library(MatchIt)
library(Matching)
library(OpenStreetMap)
library(cleangeo)
library(doParallel)
library(dplyr)
library(extrafont)
library(foreign)
library(future.apply)
library(fuzzyjoin)
library(ggplot2)
library(gridExtra)
library(lfe)
library(lme4)
library(lmtest)
library(pglm)
library(plm)
library(plyr)
library(purrr)
library(rgdal)
library(rgeos)
library(sandwich)
library(sf)
library(spdep)
library(stargazer)
library(texreg)
library(tibble)
library(tmap)
library(tmaptools)
```

### Scripts:

- *01_Long-BHPS-UKHLS.do*: Prepares the longitudinal BHPS and UKHLS data
  from the origin input data.

- *02_Floods_Floods-to-LSOA*: Prepares the flood data and merges it to
  the LSOAs.

- *03_Floods_Add-UKHLS*: Merges floods per LSOA and UKHLS data.

- *04_Floods_Data-Preparation*: Recodes necessary variables in the UKHLS
  flood data.

- *05_Floods_Descriptives*: Performs matching and produces the
  descriptive figures around the climate event floods.

- *06_Floods_Analyse*: Performs the analsysis of the main text with the
  flood incidences.

- *07_Heat_Prepare-weather*: Prepares the weather data and merges it to
  the LSOAs.

- *08_Heat_Combine-UKHLS*: Merges heatwaves per LSOA and UKHLS data.

- *09_Heat_Data-Preparation*: Recodes necessary variables in the UKHLS
  heatwaves data.

- *10_Heat_Descriptives*: Performs matching and produces the descriptive
  figures around the climate event heatwaves.

- *11_Heat_Analyse*: Performs the analsysis of the main text with the
  heat incidences. Produces the final coefficient plots.

## Input Data:

Data availability: The Understanding Society and Harmonised BHPS Special
Licence Access data (12th Edition) is not publicly available because of
privacy restrictions, but access can be acquired via the UK Data Service
(SN: 6931, ) after application. The Recorded Flood Outlines data is
publicly available via GOV.UK open data, . Data was downloaded on
2020-06-23. The HadUK-Grid weather data is publicly available from the
UK MET office after registration, . Data was downloaded 2020-12-04. R
4.2.2 was used for data analysis, see Methods section for relevant
packages.

Shapefiles on LSOAs and their centroids can be found online at ONS Open
Geography Portal: <https://geoportal.statistics.gov.uk/>

The following data is required in the “00_Input_data” folder:

- *UKDA-6931-stata*: Folder containing original Stata data of the BHPS
  UKHLS data.

- *UKDA-6670-stata*: Folder containing the geo identifiers for the the
  BHPS UKHLS data.

- *EA_RecordedFloodOutlines_SHP_Full-2020-06-23*: Folder containing the
  shape files of the flood outlines (used as of 2020-06-23).

- \*Lower_Layer_Super_Output_Areas\_\_December_2001\_\_Boundaries_EW_BGC-shp\*:
  Shapefiles of the 2001 LSOAs

- *MET-Climate*: Folder containing pre-processed data on weather across
  the UK. The R Script to query original data via API and processing
  this information into a RData file is available inside this folder.

- \*Lower_Layer_Super_Output_Areas\_\_December_2001\_\_Boundaries_EW_BGC-shp\*:
  Shapefiles of the 2001 LSOAs.

- *soa*: Shapefiles of the 2001 Super Output Areas in Norther Ireland

- *SG_DataZoneCent_2001*: Shapefiles of the 2001 Scottish Datazone
  Centroids

- \*Countries\_\_December_2016\_\_Boundaries_GB_BUC-shp\*: Shapefile of
  GB boundaries (only for maps)

- \*OSNI_Open_Data\_-*Largescale_Boundaries*-\_NI_Outline-shp\*:
  Shapefile of Norther Ireland boundaries (only for maps)

## System and version information

Platform: x86_64-w64-mingw32/x64 (64-bit)

Version: R version 4.2.2
