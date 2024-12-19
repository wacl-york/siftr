#' Read and Combine Multiple Fullscan XML SIFT Files
#'
#' This function reads multiple XML files using `read_fullscan_xml_sift` and combines the resulting dataframes
#' from each file into five tables: `settings`, `method`, `scan_masses`, `measurements` and `errors`
#'
#' @param xml_file_paths A character vector of XML file paths to be processed.
#' @param tzone A character string specifying the time zone for parsing the file timestamps.
#'        Default is "Europe/London".
#'
#' @return A list containing:
#' \describe{
#'   \item{settings}{Combined settings data from all successfully processed XML files.}
#'   \item{method}{Combined method data from all successfully processed XML files.}
#'   \item{scan_masses}{Combined scan mass data from all successfully processed XML files.}
#'   \item{measurements}{Combined measurements data from all successfully processed XML files.}
#'   \item{errors}{A dataframe listing files that failed to process and their error messages.}
#' }
#'
#' @export
read_many_fullscan_xml_sift <- function(xml_file_paths, tzone = "Europe/London") {
  message("Reading multiple XML files...")

  safe_read <- purrr::safely(read_fullscan_xml_sift)

  results <- purrr::imap(xml_file_paths, ~ {

    result <- safe_read(.x, tzone)

    if (!is.null(result$error)) {
      warning("Error reading file: ", .x, "\n  ", result$error$message, call. = FALSE)
    }

    result
  })

  successful_results <- purrr::map(results, "result") %>% purrr::compact()
  error_results <- purrr::map(results, "error") %>% purrr::imap_dfr(~ {
    if (!is.null(.x)) {
      tibble::tibble(file = xml_file_paths[.y], error_message = .x$message)
    }
  })

  if (nrow(error_results) == 0) {
    error_results <- tibble::tibble(file = character(), error_message = character())
  }

  if (length(successful_results) == 0) {
    stop("No valid XML files were processed. Check the `errors` table for details.")
  }

  result <- list(
    settings = dplyr::bind_rows(purrr::map(successful_results, "settings")),
    method = dplyr::bind_rows(purrr::map(successful_results, "method")),
    scan_masses = dplyr::bind_rows(purrr::map(successful_results, "scan_masses")),
    measurements = dplyr::bind_rows(purrr::map(successful_results, "measurements")),
    errors = error_results
  )

  message("All successfully processed files have been combined.")
  return(result)
}
