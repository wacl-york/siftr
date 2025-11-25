#' Filter Processed SIFT-MS Data
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
#' @param table_type A character string indicating the type of table to filter.
#'        Acceptable values include: `"analyte_conc"`, `"conc_per_reagent"`,
#'        `"conc_per_product"`, `"raw_analyte_conc"`, `"raw_conc_per_reagent"`, and
#'        `"raw_conc_per_product"`. Raw tables do not have the ICF correction applied.
#' @param remove_incomplete_scans A logical value. If TRUE (default), scans
#'        which are incomplete will be removed. Incomplete scans can result when
#'        a scan is ended.
#' @param conc_unit A character string specifying the desired concentration unit for
#'        conversion. Acceptable volumetric units are `"ppm"`, `"ppb"`, and `"ppt"`; and
#'        mass based units `"mg/m³"`, `"µg/m³"`, and `"ng/m³"`. The function does not convert
#'        between mass and volumetric - ensure the correct type when exporting from LabSyft.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{data}{A data frame with processed and converted concentration values in the requested unit.}
#'   \item{removed_conflicts}{A data frame of rows removed due to volumetric/mass unit conflicts.}
#' }
#'
#' @export
#'

filter_proc_sift <- function(df, table_type, remove_incomplete_scans = TRUE, conc_unit) {

  volumetric_units <- c("ppm", "ppb", "ppt")
  mass_units       <- c("mg/m³", "µg/m³", "ng/m³")
  valid_units      <- c(volumetric_units, mass_units)

  if (!conc_unit %in% valid_units) {
    stop("Invalid conc_unit. Choose from: ppm, ppb, ppt, mg/m³, µg/m³, ng/m³")
  }

  if (remove_incomplete_scans) {
    df <- df |>
      dplyr::group_by(.data$time_s, .data$table) |>
      dplyr::filter(!any(is.na(.data$conc))) |>
      dplyr::ungroup()
  }

  df <- df |>
    tidyr::separate(
      col = "table",
      into = c("table", "data_type", "table_type"),
      sep = ":"
    ) |>
    dplyr::mutate(
      data_type = stringr::str_trim(.data$data_type),
      table_type = stringr::str_trim(.data$table_type),
      start_time = stringr::str_extract(.data$table, "\\d{8}-\\d{6}") |>
        lubridate::ymd_hms(tz = "UTC", truncated = 3)
    )

  lookup <- list(
    analyte_conc = "Analyte concentrations",
    conc_per_reagent = "Analyte concentration per reagent ion",
    conc_per_product = "Analyte concentration per product ion",
    raw_analyte_conc = "Raw analyte concentrations",
    raw_conc_per_reagent = "Raw analyte concentration per reagent ion",
    raw_conc_per_product = "Raw analyte concentration per product ion"
  )

  if(!table_type %in% names(lookup)) stop("Invalid table_type provided.")

  filtered_df <- df[df$table_type == lookup[[table_type]], ]

  if (table_type %in% c("analyte_conc", "raw_analyte_conc")) {

    filtered_df <- filtered_df |>
      dplyr::mutate(
        unit = stringr::str_extract(
          .data$compound_identity,
          "(?<=\\()(ppb|ppm|ppt|mg/m³|µg/m³|ng/m³)(?=\\))"
        ),
        cas_number = stringr::str_extract(
          .data$compound_identity,
          "(?<=\\()\\d{1,7}-\\d{2}-\\d{1,2}(?=\\))"
        ),
        compound = stringr::str_remove(.data$compound_identity, "\\s*\\([^()]*\\)\\s*$") |>
          stringr::str_remove_all("\\s*\\([^()]*\\)\\s*$") |>
          stringr::str_trim()
      ) |>
      dplyr::select(.data$start_time, .data$table, .data$table_type, .data$time_s, .data$compound_identity,
                    .data$compound, .data$cas_number, .data$conc, .data$unit)
  }

  if (table_type %in% c("conc_per_reagent", "raw_conc_per_reagent")) {

    filtered_df <- filtered_df |>
      tidyr::separate(
        .data$compound_identity,
        into = c("reagent_ion", "compound_identity_2"),
        sep = " / "
      ) |>
      dplyr::mutate(
        reagent_ion = stringr::str_trim(.data$reagent_ion, side = "right"),
        compound_identity_2 = stringr::str_trim(.data$compound_identity_2, side = "left"),
        compound_identity = .data$compound_identity
      ) |>
      dplyr::mutate(
        unit = stringr::str_extract(
          .data$compound_identity_2,
          "(?<=\\()(ppb|ppm|ppt|mg/m³|µg/m³|ng/m³)(?=\\))"
        ),
        cas_number = stringr::str_extract(
          .data$compound_identity_2,
          "(?<=\\()\\d{1,7}-\\d{2}-\\d{1,2}(?=\\))"
        ),
        compound = stringr::str_remove(.data$compound_identity_2, "\\s*\\([^()]*\\)\\s*$") |>
          stringr::str_remove_all("\\s*\\([^()]*\\)\\s*$") |>
          stringr::str_trim()
      ) |>
      dplyr::select(.data$start_time, .data$table, .data$table_type, .data$time_s, .data$compound_identity,
                    .data$reagent_ion, .data$compound, .data$cas_number, .data$conc, .data$unit)
  }

  if (table_type %in% c("conc_per_product", "raw_conc_per_product")) {

    filtered_df <- filtered_df |>
      tidyr::separate(
        .data$compound_identity,
        into = c("product_ion_2", "reagent_ion", "compound_identity_2"),
        sep = " / "
      ) |>
      dplyr::mutate(compound_identity = .data$compound_identity) |>
      tidyr::separate(.data$product_ion_2,
                      into = c("product_ion", "product_mass"),
                      sep = "\\s+\\[",
                      remove = FALSE) |>
      dplyr::mutate(
        product_mass = gsub("]", "", .data$product_mass),
        product_mass = trimws(.data$product_mass),
        unit = stringr::str_extract(
          .data$compound_identity_2,
          "(?<=\\()(ppb|ppm|ppt|mg/m³|µg/m³|ng/m³)(?=\\))"
        ),
        cas_number = stringr::str_extract(
          .data$compound_identity_2,
          "(?<=\\()\\d{1,7}-\\d{2}-\\d{1,2}(?=\\))"
        ),
        compound = stringr::str_remove(.data$compound_identity_2, "\\s*\\([^()]*\\)\\s*$") |>
          stringr::str_remove_all("\\s*\\([^()]*\\)\\s*$") |>
          stringr::str_trim()
      ) |>
      dplyr::select(.data$start_time, .data$table, .data$table_type, .data$time_s, .data$compound_identity,
                    .data$reagent_ion, .data$product_ion, .data$product_mass,
                    .data$compound, .data$cas_number, .data$conc, .data$unit)
  }

  filtered_df <- filtered_df |>
    dplyr::mutate(unit_class = dplyr::case_when(
      .data$unit %in% volumetric_units ~ "volumetric",
      .data$unit %in% mass_units       ~ "mass",
      TRUE ~ "unknown"
    ))

  target_class <- if (conc_unit %in% volumetric_units) "volumetric" else "mass"

  removed_conflicts <- filtered_df |>
    dplyr::filter(.data$unit_class != target_class)

  cleaned_df <- filtered_df |>
    dplyr::filter(.data$unit_class == target_class)

  if (target_class == "volumetric") {
    cleaned_df <- cleaned_df |>
      dplyr::mutate(
        conc = dplyr::case_when(
          .data$unit == "ppm" & conc_unit == "ppb" ~ .data$conc * 1000,
          .data$unit == "ppm" & conc_unit == "ppt" ~ .data$conc * 1e6,
          .data$unit == "ppb" & conc_unit == "ppm" ~ .data$conc / 1000,
          .data$unit == "ppb" & conc_unit == "ppt" ~ .data$conc * 1000,
          .data$unit == "ppt" & conc_unit == "ppb" ~ .data$conc / 1000,
          .data$unit == "ppt" & conc_unit == "ppm" ~ .data$conc / 1e6,
          TRUE ~ .data$conc
        ),
        unit = conc_unit
      )
  }

  if (target_class == "mass") {
    cleaned_df <- cleaned_df |>
      dplyr::mutate(
        conc = dplyr::case_when(
          .data$unit == "mg/m³" & conc_unit == "µg/m³" ~ .data$conc * 1000,
          .data$unit == "mg/m³" & conc_unit == "ng/m³" ~ .data$conc * 1e6,
          .data$unit == "µg/m³" & conc_unit == "mg/m³" ~ .data$conc / 1000,
          .data$unit == "µg/m³" & conc_unit == "ng/m³" ~ .data$conc * 1000,
          .data$unit == "ng/m³" & conc_unit == "µg/m³" ~ .data$conc / 1000,
          .data$unit == "ng/m³" & conc_unit == "mg/m³" ~ .data$conc / 1e6,
          TRUE ~ .data$conc
        ),
        unit = conc_unit
      )
  }

  if (nrow(removed_conflicts) > 0) {
    warning(
      sprintf(
        "Unit mismatch detected: %d row(s) contained volumetric ↔ mass concentration conflicts and were removed. Check conflicts table and confirm conc_unit and input data unit types.",
        nrow(removed_conflicts)
      ),
      call. = FALSE
    )
  }

  return(list(
    data = cleaned_df |> dplyr::select(-.data$unit_class),
    removed_conflicts = removed_conflicts
  ))
}
