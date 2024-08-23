#' filter_proc_sift: Filter and Process Processed SIFT-MS Data Tables
#'
#' @description
#' This function is designed to filter processed SIFT-MS data outputted from the
#' read_proc_sift function. It allows you to filter the data frame for one of six
#' tables, can filter to remove incomplete scans and sets the units for
#' concentration.
#'
#' @param df A data frame outputted from the read_proc_sift function containing
#'           Processed SIFT-MS data, which includes columns `table`, `time_s`,
#'           `compound_identity`, and `conc`.
#'
#' @param table_type A character string indicating the type of concentration data
#'        to filter. Acceptable values include:
#'        - "analyte_conc": Analyte concentrations
#'        - "conc_per_reagent": Analyte concentration per reagent ion
#'        - "conc_per_product": Analyte concentration per product ion
#'        - "raw_analyte_conc": Raw analyte concentrations
#'        - "raw_conc_per_reagent": Raw analyte concentration per reagent ion
#'        - "raw_conc_per_product": Raw analyte concentration per product ion
#' @param remove_incomplete_scans A logical value. If TRUE (default), scans
#'        which are incomplete will be removed. Incomplete scans can result when
#'        a scan is ended.
#' @param conc_unit A character string specifying the desired concentration unit for
#'        conversion. Acceptable values are "ppm", "ppb", or "ppt".
#'
#' @return A data frame
#' @export

filter_proc_sift <- function(df, table_type, remove_incomplete_scans = TRUE, conc_unit) {

  valid_units <- c("ppm", "ppb", "ppt")
  if (!conc_unit %in% valid_units) {
    stop("Invalid conc_unit provided. Choose from 'ppm', 'ppb', or 'ppt'.")
  }

  df <- df |>
    separate(col = table,
             into = c("table", "data_type", "table_type"),
             sep = ":") |>
    mutate(data_type = str_trim(data_type),
           table_type = str_trim(table_type),
           start_time = str_extract(table, "\\d{8}-\\d{6}") %>%
             ymd_hms(tz = "UTC", truncated = 3))

  if (remove_incomplete_scans) {
    df <- df %>%
      group_by(time_s, table) %>%
      filter(!any(is.na(conc))) %>%
      ungroup()
  }

  lookup <- list(
    analyte_conc = "Analyte concentrations",
    conc_per_reagent = "Analyte concentration per reagent ion",
    conc_per_product = "Analyte concentration per product ion",
    raw_analyte_conc = "Raw analyte concentrations",
    raw_conc_per_reagent = "Raw analyte concentration per reagent ion",
    raw_conc_per_product = "Raw analyte concentration per product ion"
  )

  if(!table_type %in% names(lookup)) {
    stop("Invalid table_type provided.")
  }

  filtered_df <- df[df$table_type == lookup[[table_type]], ]

  if (table_type %in% c("analyte_conc", "raw_analyte_conc")) {
    filtered_df <- filtered_df %>%
      mutate(
        unit = str_extract(compound_identity, "(?<=\\()ppb|ppm|ppt(?=\\))"),
        cas_number = str_extract(compound_identity, "(?<=\\()\\d{1,7}-\\d{2}-\\d{1,2}(?=\\))"),
        compound = str_remove(compound_identity, "\\s*\\([^()]*\\)\\s*$") %>%
          str_remove_all("\\s*\\([^()]*\\)\\s*$") %>%
          str_trim()
      ) %>%
      select(start_time, table, table_type, time_s, compound_identity, compound, cas_number, conc, unit)
  }

  if (table_type %in% c("conc_per_reagent", "raw_conc_per_reagent")) {
    filtered_df <- filtered_df %>%
      separate(compound_identity, into = c("reagent_ion", "compound_identity_2"), sep = " / ") %>%
      mutate(
        reagent_ion = str_trim(reagent_ion, side = "right"),
        compound_identity_2 = str_trim(compound_identity_2, side = "left")
      ) %>%
      mutate(compound_identity = filtered_df$compound_identity) %>%
      mutate(
        unit = str_extract(compound_identity_2, "(?<=\\()ppb|ppm|ppt(?=\\))"),
        cas_number = str_extract(compound_identity_2, "(?<=\\()\\d{1,7}-\\d{2}-\\d{1,2}(?=\\))"),
        compound = str_remove(compound_identity_2, "\\s*\\([^()]*\\)\\s*$") %>%
          str_remove_all("\\s*\\([^()]*\\)\\s*$") %>%
          str_trim()
      ) %>%
      select(start_time, table, table_type, time_s, compound_identity, reagent_ion, compound, cas_number, conc, unit)
  }

  if (table_type %in% c("conc_per_product", "raw_conc_per_product")) {
    filtered_df <- filtered_df %>%
      separate(compound_identity, into = c("product_ion_2", "reagent_ion", "compound_identity_2"), sep = " / ") %>%
      mutate(compound_identity = filtered_df$compound_identity) %>%
      separate(product_ion_2, into = c("product_ion", "product_mass"), sep = "\\s+\\[", remove = FALSE) %>%
      mutate(
        product_mass = gsub("]", "", product_mass),
        product_mass = trimws(product_mass)
      ) %>%
      mutate(
        unit = str_extract(compound_identity_2, "(?<=\\()ppb|ppm|ppt(?=\\))"),
        cas_number = str_extract(compound_identity_2, "(?<=\\()\\d{1,7}-\\d{2}-\\d{1,2}(?=\\))"),
        compound = str_remove(compound_identity_2, "\\s*\\([^()]*\\)\\s*$") %>%
          str_remove_all("\\s*\\([^()]*\\)\\s*$") %>%
          str_trim()
      ) %>%
      select(start_time, table, table_type, time_s, compound_identity, reagent_ion, product_ion, product_mass, compound, cas_number, conc, unit)
  }

  if (conc_unit == "ppb") {
    filtered_df <- filtered_df %>%
      mutate(
        conc = if_else(unit == "ppm", conc * 1000, conc),
        unit = if_else(unit == "ppm", "ppb", unit),
        compound_identity = str_replace_all(compound_identity, "\\(ppm\\)", "(ppb)"),
        conc = if_else(unit == "ppt", conc / 1000, conc),
        unit = if_else(unit == "ppt", "ppb", unit),
        compound_identity = str_replace_all(compound_identity, "\\(ppt\\)", "(ppb)")
      )
  }

  if (conc_unit == "ppm") {
    filtered_df <- filtered_df %>%
      mutate(
        conc = if_else(unit == "ppb", conc / 1000, conc),
        unit = if_else(unit == "ppb", "ppm", unit),
        compound_identity = str_replace_all(compound_identity, "\\(ppb\\)", "(ppm)"),
        conc = if_else(unit == "ppt", conc / 1000000, conc),
        unit = if_else(unit == "ppt", "ppm", unit),
        compound_identity = str_replace_all(compound_identity, "\\(ppt\\)", "(ppm)")
      )
  }

  if (conc_unit == "ppt") {
    filtered_df <- filtered_df %>%
      mutate(
        conc = if_else(unit == "ppb", conc * 1000, conc),
        unit = if_else(unit == "ppb", "ppt", unit),
        compound_identity = str_replace_all(compound_identity, "\\(ppb\\)", "(ppt)"),
        conc = if_else(unit == "ppm", conc * 1000000, conc),
        unit = if_else(unit == "ppm", "ppt", unit),
        compound_identity = str_replace_all(compound_identity, "\\(ppm\\)", "(ppt)")
      )
  }

  unique_units <- unique(filtered_df$unit)
  if (length(unique_units) > 1) {
    warning("Not all units are the same after conversion. Multiple units found: ", paste(unique_units, collapse = ", "))
  }

  return(filtered_df)
}
