#' Add Date SIFT Function
#'
#' @description This function adds a date column to SIFT data, sets the timezone of date and start_time columns, and optionally filters by phase.
#'
#' @param data A data frame containing the SIFT data to be processed.
#' @param units A character string specifying the units for time in the unstitched files. Can be "s" for seconds or "ms" for milliseconds. This argument is required.
#' @param filter_phase A character string specifying the phase to filter by. Options are "SAMPLE", "BACKGROUND", "CALIBRATION", or `FALSE` to not filter. Defaults to 'FALSE'. Check the 'Phase' column in SIFT files to find correct option; for most users this will be "SAMPLE" for raw SIM or fullscan data and 'FALSE' for processed data.
#' @param force_tz A character string containing the time zone to convert to, or `FALSE` to not adjust the time zone. R must recognize the name contained in the string as a time zone on your system. Defaults to "Europe/London".
#'
#' @return A data frame with an added `date` column, time zone set for date and start time columns, and optionally filtered by phase.
#'
#' @import dplyr
#' @import lubridate
#'
#' @export
#'

add_date_sift <- function(data, units, filter_phase = FALSE, force_tz = "Europe/London") {
  # Validate the units argument
  if (!units %in% c("s", "ms")) {
    stop("Invalid value for 'units'. Please use 's' or 'ms'. Check unstitched files to find which to use.")
  }

  # Validate the filter_phase argument
  if (!filter_phase %in% c("SAMPLE", "BACKGROUND", "CALIBRATION", FALSE)) {
    stop("Invalid value for 'filter_phase'. Please use 'SAMPLE', 'BACKGROUND', 'CALIBRATION', or FALSE.")
  }

  # Select the appropriate time column based on units
  if (units == "ms") {
    time_col <- data$time_ms / 1000
  } else {
    time_col <- data$time_s
  }

  # Add date (end_time)
  data <- data %>%
    mutate(
      date = (start_time + time_col)) %>%
    mutate(
      date = round_date(
        date,
        unit = "second"))

   # Adjust time zone if force_tz is not FALSE
  if (force_tz != FALSE) {
    data <- data %>%
      mutate(
        date = force_tz(date, tzone = force_tz),
        start_time = force_tz(start_time, tzone = force_tz)
      )
  }

  # Relocate date column to the front
  data <- data %>% relocate(date)

  # Apply filter if filter_phase is not FALSE
  if (filter_phase != FALSE) {
    data <- data %>% filter(phase == filter_phase)
  }

  return(data)
}
