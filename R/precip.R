#' Hourly station precipitation data
#'
#' Station details
#' * Name: Philadelphia International Airport
#' * Network ID: COOP:366889
#' * Latitude/Longitude:	39.87327°, -75.22678°
#' * Elevation: 3m
#'
#' more details can be found \href{https://www.ncdc.noaa.gov/cdo-web/datasets/PRECIP_HLY/stations/COOP:366889/detail}{here}.
#'
#' @format A data.table with 79633 rows and 2 variables:
#' \describe{
#'   \item{date}{POSIXct format date/time}
#'   \item{value}{precipitation totals}
#' }
#' @source The National Oceanic and Atmospheric Administration (NOAA)
"precip"
