#' Read Full-Scan XML SIFT Data
#'
#' Parses a full-scan SIFT XML file and extracts information including settings, methods,
#' scan masses, and measurements. The data is returned as a list of data frames.
#'
#' @param xml_file_path Character. Path to the XML file to be read.
#' @param tzone Character. Time zone for parsing timestamps. Defaults to `Europe/London`.
#'
#' @return A named list containing the following elements:
#' \item{settings}{Data frame of settings extracted from the XML file.}
#' \item{method}{Data frame of method values extracted from the XML file.}
#' \item{scan_masses}{Data frame of masses scanned.}
#' \item{measurements}{Data frame containing fullscan data.}
#'
#' @export
#'

read_fullscan_xml_sift <- function(xml_file_path, tzone = 'Europe/London') {

  xml_file <- xml2::read_xml(xml_file_path)

  message("Reading XML file: ", xml_file_path)

  # Extract settings
  settings_xml <- xml2::xml_find_all(xml_file, "//settings/setting")
  settings_df <- data.frame(
    name = xml2::xml_attr(settings_xml, "name"),
    value = xml2::xml_text(settings_xml),
    units = xml2::xml_attr(settings_xml, "units")
  )

  # Derive start_time, job.name, and job.id from settings
  job_start_time <- settings_df %>%
    dplyr::filter(name == "job.start.time") %>%
    dplyr::pull(value) %>%
    as.numeric()

  start_time <- as.POSIXct(
    job_start_time / 1000,
    origin = "1970-01-01",
    tz = tzone
  ) + (job_start_time %% 1000) / 1000


  job_name <- settings_df %>%
    dplyr::filter(name == "job.name") %>%
    dplyr::pull(value)

  job_id <- settings_df %>%
    dplyr::filter(name == "job.id") %>%
    dplyr::pull(value) %>%
    as.numeric()

  # Include start_time, job_name, and job_id in settings_df
  settings_df <- settings_df %>%
    dplyr::mutate(start_time = start_time, job_name = job_name, job_id = job_id)

  # Process method
  method <- xml2::xml_find_all(xml_file, "//method/values/*")
  method_df <- data.frame(
    start_time = start_time,
    job_name = job_name,
    job_id = job_id,
    field = xml2::xml_name(method),
    value = xml2::xml_text(method)
  )

  # Process scan masses
  scan_masses_xml <- xml2::xml_find_all(xml_file, "//userSelectedMasses/scanMass")
  scan_masses_df <- data.frame(
    start_time = start_time,
    job_name = job_name,
    job_id = job_id,
    precursor = xml2::xml_find_first(scan_masses_xml, "precursor") %>% xml2::xml_text() %>% as.numeric(),
    start = xml2::xml_find_first(scan_masses_xml, "start") %>% xml2::xml_text() %>% as.numeric(),
    end = xml2::xml_find_first(scan_masses_xml, "end") %>% xml2::xml_text() %>% as.numeric(),
    step = xml2::xml_find_first(scan_masses_xml, "step") %>% xml2::xml_text() %>% as.numeric(),
    interval = xml2::xml_find_first(scan_masses_xml, "interval") %>% xml2::xml_text() %>% as.numeric(),
    count_limit = xml2::xml_find_first(scan_masses_xml, "count_limit") %>% xml2::xml_text() %>% as.numeric()
  )

  # Process measurements
  measurements_xml <- xml2::xml_find_all(xml_file, "//measurements/datum")
  measurements_df <- data.frame(
    start_time = start_time,
    job_name = job_name,
    job_id = job_id,
    ion_time_ms = xml2::xml_attr(measurements_xml, "time") %>% as.numeric(),
    reagent_ion = xml2::xml_attr(measurements_xml, "reagent"),
    product_mass = xml2::xml_attr(measurements_xml, "product") %>% as.numeric(),
    period = xml2::xml_attr(measurements_xml, "period") %>% as.numeric(),
    ups_intensity_nA = xml2::xml_attr(measurements_xml, "ups") %>% as.numeric() * 1e9,
    dws_intensity_pA = xml2::xml_attr(measurements_xml, "dws") %>% as.numeric() * 1e12,
    flowtube_temperature_C = xml2::xml_attr(measurements_xml, "flowtubeTemperature") %>% as.numeric(),
    flowtube_pressure_torr = xml2::xml_attr(measurements_xml, "flowtubePressure") %>% as.numeric(),
    carrier_flow_tls = xml2::xml_attr(measurements_xml, "carrierFlow") %>% as.numeric(),
    sample_flow_tls = xml2::xml_attr(measurements_xml, "sampleFlow") %>% as.numeric(),
    reaction_time_ms = xml2::xml_attr(measurements_xml, "reactionTime") %>% as.numeric(),
    icf = xml2::xml_attr(measurements_xml, "icf") %>% as.numeric(),
    af = xml2::xml_attr(measurements_xml, "af") %>% as.numeric(),
    dilution_factor_sift = xml2::xml_attr(measurements_xml, "df") %>% as.numeric(),
    raw_count = xml2::xml_text(measurements_xml) %>% as.numeric()
  ) %>%
    dplyr::arrange(.data$ion_time_ms) %>%
    dplyr::mutate(
      reagent_mass = dplyr::case_when(
        .data$reagent_ion == "H3O+" ~ 19,
        .data$reagent_ion == "NO+" ~ 30,
        .data$reagent_ion == "O2+" ~ 32,
        TRUE ~ NA_real_
      )
    )

  first_reagent_ion <- measurements_df %>% dplyr::arrange(.data$ion_time_ms) %>% dplyr::pull(.data$reagent_ion) %>% dplyr::first()
  first_product_mass <- measurements_df %>% dplyr::arrange(.data$ion_time_ms) %>% dplyr::pull(.data$product_mass) %>% dplyr::first()

  time_ms <- measurements_df %>%
    dplyr::filter(.data$reagent_ion == first_reagent_ion & .data$product_mass == first_product_mass) %>%
    dplyr::pull(.data$ion_time_ms)

  scan_intervals <- c(0, time_ms)

  measurements_df <- measurements_df %>%
    dplyr::mutate(count_cps = .data$raw_count / (period / 1e6)) %>%
    dplyr::mutate(
      count_ICF_corrected_cps = (.data$count_cps * .data$icf),
      time_ms = sapply(.data$ion_time_ms, function(x) {
        idx <- max(which(x >= scan_intervals))
        scan_intervals[idx]
      })
    ) %>%
    dplyr::select(c(start_time, job_name, job_id, .data$time_ms, .data$ion_time_ms, .data$reagent_ion, .data$reagent_mass, dplyr::everything()))

  method_df <- method_df %>%
    tidyr::pivot_wider(names_from = .data$field, values_from = .data$value)

  settings_df <- settings_df %>%
    tidyr::pivot_wider(names_from = .data$name, values_from = c(.data$value, units), names_glue = "{name}{ifelse(.value == 'value', '', paste0('_', .value))}")

  result <- list(
    settings = settings_df,
    method = method_df,
    scan_masses = scan_masses_df,
    measurements = measurements_df
  )

  message("  File read successfully")

  return(result)
}
