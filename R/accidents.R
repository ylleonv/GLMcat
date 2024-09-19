#' Accidents Dataset
#'
#' This dataset contains information about various accidents, including details such as accident severity, road and weather conditions, light conditions, and the number of casualties.
#'
#' @format A data frame with 109,577 rows and 12 variables:
#' \describe{
#'   \item{accident_severity}{Factor with levels \code{Slight}, \code{Serious}, \code{Fatal}}
#'   \item{road_type}{Factor with levels \code{Dual carriageway}, \code{One way street}, \code{Roundabout}, \code{Single carriageway}, \code{Slip road}}
#'   \item{weather_conditions}{Factor with levels \code{Fine + high winds}, \code{Fine no high winds}, \code{Fog or mist}, \code{Raining + high winds}, \code{Raining no high winds}, \code{Snowing}}
#'   \item{light_conditions}{Factor with levels \code{Darkness}, \code{Daylight}}
#'   \item{day_of_week}{Factor with levels \code{Monday}, \code{Tuesday}, \code{Wednesday}, \code{Thursday}, \code{Friday}, \code{Saturday}, \code{Sunday}}
#'   \item{number_of_casualties}{Numeric, number of casualties in the accident}
#'   \item{urban_or_rural_area}{Factor with levels \code{Urban}, \code{Rural}}
#'   \item{speed_limit}{Numeric, speed limit at the accident location}
#'   \item{junction_detail}{Factor with levels \code{Not at junction or within 20 metres}, \code{T or staggered junction}, \code{Crossroads}, \code{Roundabout}, \code{Other junction}, \code{Private drive or entrance}}
#'   \item{carriageway_hazards}{Factor with levels \code{Any animal in carriageway (except ridden horse)}, \code{Data missing or out of range}, \code{None}, \code{Other object on road}, \code{Pedestrian in carriageway - not injured}, \code{Previous accident}, \code{Vehicle load on road}}
#'   \item{weather}{Factor with levels \code{Fine + high winds}, \code{Fine no high winds}, \code{Fog or mist}, \code{Raining + high winds}, \code{Raining no high winds}, \code{Snowing}}
#'   \item{road}{Factor with levels \code{Dual carriageway}, \code{One way street}, \code{Roundabout}, \code{Single carriageway}, \code{Slip road}}
#' }
#'
#' @keywords datasets
#'
#' @source Data from 2019, openly available at \url{https://data.gov.uk}, accessed in September 2023.
#' @examples
#' data(accidents)
"accidents"
