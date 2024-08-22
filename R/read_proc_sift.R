#' Read Processed SIFT Data
#'
#' @description This is a function to read data that has already been processed
#'   by LabSyft (into up to six different tables - Analyte concentrations,
#'   Analyte concentrations per reagent ion, Analyte concentrations per product
#'   ion, and their raw equivalents) into a single data frame. Each .csv file
#'   can contain multiple times. Use lapply() and bind_rows() to read multiple
#'   sheets at once.
#'
#' @param file_path A character string specifying the path to the Processed SIFT .csv file.
#'
#' @return A data frame with columns: table, time_s, compound_identity, conc.
#' @export
#'

read_proc_sift <- function(file_path) {
  # Read the entire CSV file
  raw_data <- readLines(file_path)

  # Initialize an empty data frame to store the final result
  result_df <- data.frame(
    table = character(),
    time_s = numeric(),
    compound_identity = character(),
    conc = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop through the lines and find headers
  i <- 1
  while (i <= length(raw_data)) {
    line <- raw_data[i]

    # Check if the line contains ".xml" which indicates the start of a new table
    if (grepl(".xml", line)) {
      table <- line  # Capture the entire .xml line as the table name

      # Ensure there is a line for column names
      if ((i + 1) > length(raw_data)) {
        warning(paste("No column names found after .xml header at line", i))
        break
      }

      # Use read.csv to correctly parse column names, which may contain commas
      colnames_line <- read.csv(text = raw_data[i + 1], header = FALSE, stringsAsFactors = FALSE)
      colnames_line <- unlist(colnames_line)

      # Find the index of "Time (s)"
      time_col_idx <- which(trimws(colnames_line) == "Time (s)")

      if (length(time_col_idx) == 0) {
        warning(paste('"Time (s)" column not found after .xml header at line', i))
        break
      }

      # Extract the table data (assuming data follows the column names)
      data_start <- i + 2
      j <- data_start

      while (j <= length(raw_data) && !grepl(".xml", raw_data[j])) {
        j <- j + 1
      }

      table_data <- read.csv(
        text = paste(raw_data[data_start:(j-1)], collapse = "\n"),
        header = FALSE, stringsAsFactors = FALSE
      )

      # Remove trailing empty rows if present
      table_data <- table_data[complete.cases(table_data[, time_col_idx]), ]

      # Loop through each column (except the Time (s) column)
      for (col_idx in (time_col_idx + 1):ncol(table_data)) {
        # Ensure that the concentration column has the same length as the time column
        conc_values <- as.numeric(table_data[[col_idx]])
        conc_values <- ifelse(is.na(conc_values), NA, conc_values)

        # Construct the data frame for the current compound
        temp_df <- data.frame(
          table = rep(table, length(conc_values)),
          time_s = as.numeric(table_data[[time_col_idx]]),
          compound_identity = rep(colnames_line[col_idx], length(conc_values)),
          conc = conc_values,
          stringsAsFactors = FALSE
        )

        # Explicitly set row names to NULL to avoid warnings
        rownames(temp_df) <- NULL

        # Bind this to the result dataframe
        result_df <- rbind(result_df, temp_df)
      }

      # Move to the next table
      i <- j
    } else {
      # Just move to the next line
      i <- i + 1
    }
  }

  return(result_df)
}
