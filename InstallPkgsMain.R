#### Introduction #####
# This script calls functions to install R-packages and to obtain information
# about them. It was created by Jesse Alderliesten, see
# https://github.com/JesseAlderliesten


#### Required user input ####
# Library path to directory where R-packages are stored or should be installed
lib <- "C:/Program Files/R/R-4.1.3/library"

# Countries to select mirrors from, in order of decreasing priority
mirror_countries <- c("Belgium", "Germany")

#### Preparations to use this script ####

# Load required functions and package lists
source("./InstallPkgsFuncs.R")
source("./InstallPkgsLists.R")

# Perform various checks to ensure the rest of the script can run. It is
# sufficient to run this only once after (re)installing R
prepare_install()

# Select mirrors
choose_mirrors(countries = mirror_countries,
               databases = c("BioConductor", "CRAN"))

# Package lists
pkgs_lists <- list(high_prio_pkgs = high_prio_pkgs, used_pkgs_UU = used_pkgs_UU,
                   used_pkgs_UvA = used_pkgs_UvA)

pkgs_lists <- list(high_prio_pkgs = high_prio_pkgs, 
                   coding_data_handling = coding_data_handling,
                   coding_documentation = coding_documentation,
                   coding_package_management = coding_package_management,
                   coding_performance = coding_performance,
                   graphics_palettes = graphics_palettes,
                   graphics_combine_plots = graphics_combine_plots,
                   graphics_pkgs = graphics_pkgs,
                   math_distr_random = math_distr_random,
                   math_multidimensional_networks = math_multidimensional_networks,
                   math_num_optim = math_num_optim,
                   math_stats_sets = math_stats_sets,
                   math_stats_MCMC_bayes = math_stats_MCMC_bayes,
                   math_stats_regression = math_stats_regression,
                   math_stats = math_stats,
                   microbio_eco = microbio_eco,
                   ODE_epi_growth = ODE_epi_growth,
                   ODE_handling = ODE_handling,
                   ODE_PBPK = ODE_PBPK,
                   taxonomy_phylogeny = taxonomy_phylogeny
)


#### Check if package lists contain duplicates #### 
check_duplicates(pkgs_lists = pkgs_lists, neglect_repos = TRUE)


#### Check if all packages from the lists are installed and functional #### 
list_nonfunctional_pkgs(pkgs = pkgs_lists, save_file = TRUE, sort = TRUE,
                        quietly = FALSE, verbose = FALSE)


#### Check if all installed packages are up-to-date ####
check_status(checkBuilt = TRUE, type = "binary", save_file = TRUE,
             print_output = "both")


#### Install new packages #### 
# The 'update' argument is set to FALSE to prevent inadvertently changing the
# version of already installed packages when installing new packages. However,
# updating out-of-date packages might be preferable to prevent compatibility
# issues between already installed packages and newly installed packages. If
# packages are not functional after updating, re-install them using the argument
# force = TRUE.
# Note:
#   Run R as administrator to install or update packages!
#   'pkgs' should be a charactor vector, not a list of character vectors.
BiocManager::install(pkgs = new_pkgs, lib = lib, verbose = FALSE, type = "binary",
                     update = FALSE, ask = FALSE, checkBuilt = TRUE, force = FALSE)


#### List dependencies of specific packages ####
list_dependencies(pkgs, deps_type = "strong", recursive = TRUE, name_per_pkg = FALSE,
                  number_per_pkg = TRUE, name_total = TRUE,
                  add_pkgs_to_total = FALSE, exclude_high_prio = FALSE,
                  exclude_pkgs = NULL, sort_ndeps_by = "ndeps")


#### Save details of installed packages to a .csv file ####
save_details(PC_name = "desktop")
