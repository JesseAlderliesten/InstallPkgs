#### Introduction #####
# This script can be used to check if R-packages are installed, get information
# about installed R-packages, and to install new R-packages. It was created by
# Jesse Alderliesten, see https://github.com/JesseAlderliesten

# It is assumed R has been installed, see the file 'InstallR.txt' for
# instructions. It is also assumed the R-scripts "InstallPkgsFuncs.R" and
# "InstallPkgsLists.R" are in the same folder as the current R-script.
# See also the file 'InstructionsPkgs.txt' for details.


#### To do ####
# - Add remark about questions 'update packages' and 'Do you want to install
#   from sources the package which needs compilation'? that will probably arise
#   during installation.


#### Wishlist ####
# - Various functions (e.g., list_dependencies()), fail if no internet
#   connection is present. How do other packages handle this?


#### Preparations to use this script ####
# Load required functions and package lists
source(file.path(".", "InstallPkgsFuncs.R"))
source(file.path(".", "InstallPkgsLists.R"))

# Get the path where R-packages are (or will be) installed from the first
# non-empty element of .libPaths() (if present one that contains the current R
# version number), see annotation of the function 'get_paths()' for details.
lib_path <- get_paths(path = character(0), quietly = FALSE)$first_path

# Perform various checks to ensure the rest of the script can run.
prepare_install()

# Package lists
pkgs_lists <- list(high_prio_pkgs = high_prio_pkgs,
                   used_pkgs_UU = used_pkgs_UU,
                   used_pkgs_UvA = used_pkgs_UvA)

pkgs_lists_UvA <- list(used_pkgs_UvA, UvA_data, UvA_data_val, UvA_EnzKin, 
                       UvA_ExpDesign, UvA_Fitting, UvA_misc, UvA_ODEs, UvA_plot, 
                       UvA_simul, UvA_stat, UvA_task_views)

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
                   taxonomy_phylogeny = taxonomy_phylogeny)


#### Check if package lists contain duplicates #### 
check_duplicates(pkgs_lists = pkgs_lists, neglect_repos = TRUE, quietly = FALSE)


#### Check if all packages from the lists are installed and functional #### 
# Note:
# - Input for argument 'pkgs' can be a character vector or a list of character
#   vectors.
find_nonfunctional_pkgs(pkgs = pkgs_lists, lib = lib_path, save_file = TRUE,
                        sort = TRUE, quietly = FALSE, verbose = FALSE)


#### Check if all installed packages are up-to-date ####
check_status(lib = lib_path, checkBuilt = TRUE, type = "both", save_file = TRUE,
             print_output = "both")


#### Save details of installed packages to a CSV file ####
save_details(PC_name = "desktop")


#### List dependencies of specific packages ####
# Note:
# - Input for argument 'pkgs' can be a character vector or a list of character
#   vectors.
list_dependencies(pkgs = pkgs_lists, deps_type = "strong", recursive = TRUE,
                  name_per_pkg = FALSE, number_per_pkg = TRUE,
                  name_total = TRUE, add_pkgs_to_total = FALSE,
                  exclude_high_prio = FALSE, exclude_pkgs = NULL,
                  sort_ndeps_by = "ndeps")


#### Install new packages #### 
# Notes:
# - Run RStudio as administrator to install or update packages! If the warning
#   'lib = ... is not writeable' was issued, you most likely forgot to run
#   RStudio as administrator: close RStudio and open it with administrator
#   rights by right-clicking on the RStudio icon and selecting 'Run as
#   administrator'. Within RStudio, open the R-project file 'InstallPkgs.Rproj',
#   and try installing again.
# - Argument 'update' is set to FALSE to prevent inadvertently changing the
#   version of already-installed packages when installing new packages. However,
#   updating out-of-date packages might be preferable to prevent compatibility
#   issues between already-installed packages and newly-installed packages.
#   To update packages in base-R (i.e., not using BioConductor), use 
#   update.packages(instlib = lib, ask = FALSE, checkBuilt = TRUE).
# - If packages are not functional after updating, re-install them using the
#   argument force = TRUE.
# - For an overview which Bioconductor release corresponds to which R version,
#   see http://bioconductor.org/about/release-announcements/#release-versions
BiocManager::install(pkgs = unlist(new_pkgs, use.names = FALSE),
                     lib.loc = lib, lib = lib, verbose = FALSE,
                     type = "both", update = FALSE, ask = FALSE,
                     checkBuilt = TRUE, force = FALSE, version = 3.16)
