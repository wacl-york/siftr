#' Validate SIFT data as quantitative
#'
#' This function validates that fullscan SIFT-MS data is quantitative through
#' addition of a logical column. It ensures that the sum of product ions for a
#' particular reagent ion and date (dttm) do not exceed that of the reagent ion
#' (for the H30+ ion this includes the first water cluster). This function is
#' designed for data obtained with the H3O+, NO+ and/or O2+ reagent ions.
#'
#' @param data A data frame containing fullscan data SIFT including columns
#'  `ion`, `date`, and `value`.
#'
#' @return A list of tibbles. `validated_data` is a data frame with additional
#' columns (`product_ion_mass`, `reagent_ion_mass`, `total_ion _count`,
#' `reagent_ion_count`, `sum_product_ions`, `ratio_reagent_to_product` and
#' `quantitative`) containing validation information, and there are two other
#' data frames to aid with potential warning messages.
#'
#' @export
#'

validate_as_quantitative_sift <- function(data) {

  data <- data %>%
    separate(ion, into = c("product_ion_mass", "reagent_ion_mass"), sep = " \\(", remove = FALSE) %>%
    mutate(
      product_ion_mass = as.numeric(gsub("\\+", "", product_ion_mass)),
      reagent_ion_mass = as.numeric(gsub("\\+\\)", "", reagent_ion_mass))
    )

  valid_reagents <- c(19, 30, 32)
  if (any(!data$reagent_ion_mass %in% valid_reagents)) {
    stop("Invalid reagent_ion_mass values detected. Check values are only 19, 30 or 32")
  }

  data <- data %>%
    group_by(date, reagent_ion_mass) %>%
    mutate(total_ion_count = sum(value, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(date) %>%
    mutate(
      reagent_ion_count = case_when(
        reagent_ion_mass == 19 ~ sum(value[ion %in% c("19+ (19+)", "37+ (19+)")], na.rm = TRUE),
        reagent_ion_mass == 30 ~ sum(value[ion == "30+ (30+)"], na.rm = TRUE),
        reagent_ion_mass == 32 ~ sum(value[ion == "32+ (32+)"], na.rm = TRUE),
        TRUE ~ NA_real_
      ),
      sum_product_ions = total_ion_count - reagent_ion_count,
      ratio_reagent_to_product = case_when(
        sum_product_ions != 0 ~ reagent_ion_count / sum_product_ions,
        TRUE ~ NA_real_
      ),
      quantitative = !is.na(ratio_reagent_to_product) & ratio_reagent_to_product > 1
    ) %>%
    ungroup()

  na_rows <- data %>%
    filter(is.na(sum_product_ions)) %>%
    select(date, ion) %>%
    distinct()

  na_summary <- if (nrow(na_rows) > 0) {
    na_summary <- na_rows %>%
      group_by(date) %>%
      summarise(ions = paste(unique(ion), collapse = ", "), .groups = 'drop')

    warning_msg1 <- "Some rows have NA values for sum_product_ions. Here are the details:\n"
    warning_msg1 <- paste0(warning_msg1,
                          paste0(apply(na_summary, 1, function(row) {
                            paste("Date:", row["date"], "Ions:", row["ions"])
                          }), collapse = "\n"))

    warning(warning_msg1)
    na_summary
  } else {
    NULL
  }

  non_quantitative_rows <- data %>%
    filter(!quantitative) %>%
    select(date, ion) %>%
    distinct()

  non_quantitative_summary <- if (nrow(non_quantitative_rows) > 0) {
    non_quantitative_summary <- non_quantitative_rows %>%
      group_by(date) %>%
      summarise(ions = paste(unique(ion), collapse = ", "), .groups = 'drop')

    warning_msg2 <- "Some rows have non-quantitative results. Here are the details:\n"
    warning_msg2 <- paste0(warning_msg2,
                          paste0(apply(non_quantitative_summary, 1, function(row) {
                            paste("Date:", row["date"], "Ions:", row["ions"])
                          }), collapse = "\n"))

    warning(warning_msg2)
    non_quantitative_summary
  } else {
    NULL
  }

  result <- list(
    validated_data = data,
    sum_product_ions_warning_data = na_summary,
    non_quantitative_summary_data = non_quantitative_summary
  )

  return(result)
}
