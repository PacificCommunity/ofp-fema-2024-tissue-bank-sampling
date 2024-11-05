## Make a folder if it doesn't already exist
make_folder <- function(path, ...) if(!dir.exists(path)) dir.create(path, ...)

## Load in an .RData and rename the object
rename_RData <- function(...) {
    obj_name <- load(...)
    get(obj_name)
}

## Function to get folders at a location
get_folders <- function(fld_path) {
    list.dirs(fld_path, recursive = FALSE)
}

## Convenience function to set NAs to zero
na_to_zero <- function(x) {
    x[is.na(x)] <- 0
    x
}

## Safe implementation of sample - see ?sample
resample <- function(x, ...) x[sample.int(length(x), ...)]

## Convert from longitude -180 - 180 to 0 - 360
lon_reformat <- function(lon) {
    ndx <- lon < 0 & !is.na(lon)
    lon[ndx] <- lon[ndx] + 360
    lon
}

## Colour blind friendly palette by Okabe & Ito
ito_cols <- c("#56B4E9", "#E69F00", "#009E73", "#F0E442", "#CC79A7", "#0072B2", "#D55E00", "#000000")
names(ito_cols) <- c("light blue", "orange", "green", "yellow", "pink", "dark blue", "red", "black")
