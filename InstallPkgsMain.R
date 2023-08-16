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
# - 'list' is used both in the strict sense as 'an object of class (mode?) list'
#   and in the loose sense as 'a character vector'. Distinguish between these
#   two uses.
# - 'pkgs' should be a character vector, not a list of character vectors, so it
#   is needed to use 'new_pkgs <- unique(unlist(pkgs_lists))'. Let that be
#   handled inside the functions if necessary.


#### Required user input ####
# The directory where R-packages are or will be installed is taken from the
# first element of .libPaths() that contains the string 'Program Files'. See
# ?libPaths and ?Startup for details, and the function 'select_libpath()' for
# alternative options.
lib <- grep("Program Files", x = .libPaths(), value = TRUE, fixed = TRUE)[1]
message(lib, " will be used as library path")


#### Preparations to use this script ####
# Load required functions and package lists
source(file.path(".", "InstallPkgsFuncs.R"))
source(file.path(".", "InstallPkgsLists.R"))

# Perform various checks to ensure the rest of the script can run. It is
# sufficient to run this only once after (re)installing R.
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
list_nonfunctional_pkgs(pkgs = pkgs_lists, save_file = TRUE, sort = TRUE,
                        quietly = FALSE, verbose = FALSE)


#### Check if all installed packages are up-to-date ####
check_status(checkBuilt = TRUE, type = "both", save_file = TRUE,
             print_output = "both")


#### Save details of installed packages to a CSV file ####
save_details(PC_name = "desktop")


#### List dependencies of specific packages ####
list_dependencies(pkgs, deps_type = "strong", recursive = TRUE,
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
# - 'pkgs' should be a character vector, not a list of character vectors.
#   So you probably first want to run 'new_pkgs <- unique(unlist(pkgs_lists))'.
# - Argument 'update' is set to FALSE to prevent inadvertently changing the
#   version of already-installed packages when installing new packages. However,
#   updating out-of-date packages might be preferable to prevent compatibility
#   issues between already-installed packages and newly-installed packages.
# - If packages are not functional after updating, re-install them using the
#   argument force = TRUE.
# - For an overview which Bioconductor release corresponds to which R version,
#   see http://bioconductor.org/about/release-announcements/#release-versions
BiocManager::install(pkgs = new_pkgs, lib.loc = lib, lib = lib, verbose = FALSE,
                     type = "both", update = FALSE, ask = FALSE,
                     checkBuilt = TRUE, force = FALSE, version = 3.16)
