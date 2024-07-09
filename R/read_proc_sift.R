#' Read Processed SIFT Data
#'
#' @description This is a function to read data that has already been processed
#'   by the SYFT into three tables - analyte concentrations, concentrations per
#'   reagent and concentrations per product.
#'
#' @param file The file path for the pre-processed SIFT .csv file.
#'
#' @return A list.
#' @export
#'

read_proc_sift <- function(file) {
  lines <- tibble::tibble(line = readLines(file)) %>%
    dplyr::mutate(start = dplyr::row_number())

  # Calculate start and end lines for each section
  start_ends <- lines %>%
    dplyr::filter(stringr::str_detect(line, ".xml")) %>%
    dplyr::mutate(end = dplyr::lead(start)-2) %>%
    tidyr::replace_na(list(end = max(lines$start)))

  # Subfunction: Read in each bit of the data individually by mapping over our starts and ends
  read_data <- function(start, end) {
    title <- suppressWarnings(readr::read_csv(
      file = file, col_types = readr::cols(),
      skip = start - 1, n_max = 1, col_names = F
    )) %>%
      dplyr::pull(1) %>%
      gsub("\\..*", "", .) %>%
      janitor::make_clean_names()

    df <- suppressWarnings(readr::read_csv(file,
      col_types = readr::cols(),
      skip = start, n_max = end - start,
      locale = readr::locale(encoding = "ISO-8859-1")
    )) %>%
      dplyr::mutate(
        table = title,
        start_time = stringr::word(table, start = -2, end = -1, sep = "_") %>% lubridate::ymd_hms(),
        table = stringr::word(table, start = 1, end = -3, sep = "_")
      ) %>%
      dplyr::relocate(table, start_time) %>%
      tidyr::pivot_longer(-(1:3), values_to = "conc") %>%
      dplyr::mutate(name = stringr::str_replace_all(name, "\\(R\\)", "R")) %>%
      dplyr::mutate(name = stringr::str_replace_all(name, "\\(S\\)", "S")) %>%
      janitor::clean_names() %>%
      tidyr::drop_na()

    return(df)
  }

  # Map using the above function - returns a list.
  raw <- purrr::map2(.x = start_ends$start, .y = start_ends$end, .f = ~ read_data(start = .x, end = .y))

  analyte_conc <- raw[[1]] %>%
    tidyr::separate(name, into = c("name", "numbers", "unit"), sep = " \\(") %>%
    dplyr::mutate(dplyr::across(where(is.character), stringr::str_remove, "\\)"),
      unit = stringr::str_replace_all(unit, "\U00B3", "3")
    )

  conc_per_reagent <- raw[[2]] %>%
    tidyr::separate(name, into = c("reagent_ion", "compound", "numbers", "unit"), sep = " \\(| / ") %>%
    dplyr::mutate(dplyr::across(where(is.character), stringr::str_remove, "\\)"))

  conc_per_product <- raw[[3]] %>%
    tidyr::separate(name,
      into = c("product_ion", "product_ion_mass", "reagent_ion", "compound", "numbers", "unit"),
      sep = " \\[|  / | / | \\("
    ) %>%
    dplyr::mutate(dplyr::across(where(is.character), stringr::str_remove, "\\)|\\]"))

  sift <- list(
    analyte_conc = analyte_conc,
    conc_per_reagent = conc_per_reagent,
    conc_per_product = conc_per_product
  )

  return(sift)
}
