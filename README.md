# SIFT-R {siftr} 📦 <img src="man/figures/logo.png" align="right" height="138" />

Tools for streamlining the use of SIFT-MS data in WACL.

To install, just write the following into your R console.

```
devtools::install_github("wacl-york/siftr")
```

# Reading Data

Currently, the only feature of the package is simplifying the reading-in of SIFT data. 

## Raw SIFT data

The key function is `read_sift`, which simply reads in raw SIFT data and returns a list of data frames. By "raw" I refer to data laid out in nine sub-tables, beginning with a meta-data section, followed by the PREPARATION phase, then SAMPLE phase, etc., concluding with the summary section.

```
siftr::read_sift("~/my_sift_data.csv")
```

Depending on the mode the SIFT runs in, the structure of its .csv can change. For example, the preparation phase may not be included. If this is the case, use the `drop_prep` argument.

```
siftr::read_sift("~/my_sift_data_batched.csv", drop_prep = T)
```

`siftr` also defaults to being very chatty, letting you know what it is doing to keep you updated and help debug any issues. If you want it to shut up, set the `chatty` argument to `FALSE`.

If you have multiple SIFT files and are only interested in one of the tables, you can use the `read_many_sift` function with a vector of file paths and a table (one of "time", "meta", "prep_phase", "sample_phase", "phase_mean_values", "intensity_corrected", "time_vs_mass", "concentrations", "analytes" or "summary", defaulting to "concentrations").

```
siftr::read_many_sift(my_files, table = "summary")
```

`siftr::read_many_sift()` can also be passed the `drop_prep` argument.

## Pre-Processed SIFT data

Depending on the set-up of the SYFT, data output can be pre-processed into a more readable format with just three sub-tables; analyte concentrations, concentrations per reagents and concentrations per products. In this case, the function `read_proc_sift` can be used to quickly read these into R.

```
siftr::read_proc_sift("~/my_sift_lab_data.csv")
```

There is an equivalent function for reading a single table from multiple files. Once again, provide a vector of file paths and a table (one of "analyte_conc", "conc_per_reagent" or "conc_per_product", defaulting to "analyte_conc".)

```
siftr::read_many_proc_sift(my_files, table = "conc_per_product")
```
