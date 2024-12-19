#' Add Date SIFT Function
#'
#' @description This function adds a date column (dttm) to SIFT data from the scan start time (dttm) and time (numeric) columns, and can set the timezone of date and start time columns.
#'
#' @param data A data frame (from other read_sift functions) containing the SIFT data to be processed.
#' @param time_units A character string specifying the units for time in the unstitched files. Can be "s" for seconds or "ms" for milliseconds. This argument is required.
#' @param time_input A character string specifying the column name for the time input (e.g., "time_s" or "time_ms").
#' @param start_time_input A character string specifying the column name for the start time input. Defaults to "start_time".
#' @param force_tz A character string containing the time zone to convert to, or `FALSE` to not adjust the time zone. Defaults to `FALSE`. For use in WACL, use "Europe/London".
#' @param output A character string specifying the name of the output column for the calculated date. Defaults to "date".
#' @param round_date A logical value indicating whether to round the output to the nearest second. Defaults to TRUE.
#'
#' @return A data frame with an added output column (default: `date`), and time zone set for the specified columns.
#'
#' @import dplyr
#' @import lubridate
#' @import rlang
#'
#' @export
#'

add_date_sift <- function(data, time_units, time_input, start_time_input = "start_time", force_tz = FALSE, output = "date", round_date = TRUE) {

  if (!time_units %in% c("s", "ms")) {
    stop("Invalid value for 'time_units'. Please use 's' or 'ms'. Check unstitched files to find which to use.")
  }

  if (!(time_input %in% colnames(data))) {
    stop("The specified time_input column does not exist in the data.")
  }

  if (!(start_time_input %in% colnames(data))) {
    stop("The specified start_time_input column does not exist in the data.")
  }

  if (!is.numeric(data[[time_input]])) {
    stop("The time_input column must be numeric.")
  }


  time_col <- if (time_units == "ms") data[[time_input]] / 1000 else data[[time_input]]

  data <- data %>%
    dplyr::mutate(
      !!rlang::sym(output) := .data[[start_time_input]] + time_col
    )

  if (round_date) {
    data <- data %>%
      dplyr::mutate(
        !!rlang::sym(output) := lubridate::round_date(.data[[output]], unit = "second")
      )
  }

  if (force_tz != FALSE) {
    data <- data %>%
      dplyr::mutate(
        !!rlang::sym(output) := lubridate::force_tz(.data[[output]], tzone = force_tz),
        !!rlang::sym(start_time_input) := lubridate::force_tz(.data[[start_time_input]], tzone = force_tz)
      )
  } else {
    warning("Timezone not specified. Set a timezone using the force_tz parameter. Use ?add_date_sift to see function documentation.")
  }

  data <- data %>% dplyr::relocate(.data[[output]])

  return(data)
}
