# heat

## Installation

### For Positron Users

If you're using Positron IDE, you need to install the GitHub version of the `progress` package for progress bars to work correctly:

``` r
# Install pak if not available
if (!requireNamespace('pak', quietly = TRUE)) {
  install.packages('pak')
}

# Install progress from GitHub (required for Positron compatibility)
pak::pak('r-lib/progress')

# Then install heat
devtools::install()
```

### For RStudio Users

Standard installation works fine:

``` r
devtools::install()
```

The package will automatically install the required GitHub version of `progress` when you use `devtools::install()` or `remotes::install_github()` thanks to the `Remotes` field in DESCRIPTION.

## Progress Bars

The package displays progress bars during batch processing showing: - Visual progress indicator - Percentage complete - Current batch out of total batches - Estimated time remaining (ETA)

Progress bars are automatically displayed for operations with multiple batches and work in both RStudio and Positron.