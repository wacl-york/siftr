#' Read SIFT Full Scan Data
#'
#' This function reads extracts the "Time vs Mass" table from a SIFT-MS file. Useful for fullscan SIFT-MS data which is missing several tables usually found in a SIM SIFT-MS file and causes an error in read_sift.
#'
#' @param file The file path for the SIFT .csv file.
#' @param drop_prep Is the PREPARATION sub-data missing? TRUE/FALSE.
#' @param chatty \code{TRUE}/\code{FALSE}. Should the function communicate what it is doing? Useful for debugging.
#' @param warn logical. if a section of the file cannot be read, should a warning or error be thrown? Default FALSE leads to an error being produced. For use with \code{read_many_sift()} to skip bad files
#'
#' @return A tibble containing processed time vs mass data.
#'
#' @export
#'

read_fullscan_sift <- function(file, drop_prep = F, chatty = T, warn = FALSE) {
  # Time vs mass
  flags <- c("Mass Vs Time", "Cycle vs Product", "Time vs Mass",
             "Detailed Compound Concentrations", "Analyte vs Time",
             "Summary,")
  lines <- tibble::tibble(line = readLines(file)) %>% dplyr::mutate(start = dplyr::row_number())

  start_ends <- lines %>% dplyr::filter(stringr::str_detect(line,
                                                            paste(flags, collapse = "|"))) %>% rbind(tibble::tibble(line = "start",
                                                                                                                    start = 0)) %>% dplyr::arrange(start) %>% dplyr::mutate(end = dplyr::lead(start)) %>%
    tidyr::replace_na(list(end = max(lines$start))) %>%
    dplyr::mutate(start = dplyr::if_else(line == "Summary,:",
                                         start + 1, start))
  read_data <- function(start, end) {
    if (start == 0) {
      df <- suppressMessages(suppressWarnings(readr::read_csv(file,
                                                              col_types = readr::cols(), na = c(":", " ",
                                                                                                ""), skip = start, n_max = end - start, col_names = F)))
    }
    else {
      df <- suppressMessages(suppressWarnings(readr::read_csv(file,
                                                              col_types = readr::cols(), na = c(":", " ",
                                                                                                ""), skip = start, n_max = end - (start +
                                                                                                                                    2))))
    }
    return(df)
  }
  raw <- purrr::map2(.x = start_ends$start, .y = start_ends$end,
                     .f = ~read_data(start = .x, end = .y))
  n <- 0



  # Meta
  meta <- raw[[1]] %>% janitor::remove_empty(which = c("rows",
                                                       "cols")) %>% tidyr::drop_na(X2) %>% dplyr::select(-X3) %>%
    tidyr::pivot_wider(names_from = "X1", values_from = "X2") %>%
    janitor::clean_names()
  start_time <- meta$job_start_date %>% lubridate::ymd_hms()

  # Time vs mass
  if (drop_prep) {
    n <- n - 1
  }
  time_vs_mass <- raw[[n + 6]] %>%
    janitor::remove_empty(which = c("rows", "cols")) %>%
    tidyr::pivot_longer(-(1:10), names_to = "ion") %>%
    janitor::clean_names()

  cbind(start_time, time_vs_mass) |> as_tibble()
}
