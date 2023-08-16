#### Introduction #####
# Functions to install R-packages and obtain information about installed packages.


#### To do ####
# - Change the documentation on installing RTools. See the point on Rtools in
#   the wishlist for additional info.
#   'RTools ('the C++ Toolchain) has to be downloaded from
#   https://cran.r-project.org/bin/windows/Rtools/ and has to be installed.
#   In addition, RTools has to be configured by putting the location of Rtools
#   utilities (bash, make, etc) on the search path if it is not there yet.
#   devtools::find_rtools(debug = TRUE) checks if this is done, and gives
#   instructions how to perform it if it is not yet done.
#   For R versions 4.0.0 to 4.1.3, putting Rtools on the search path can be
#   achieved by running the following line:
#   write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron",
#   append = TRUE)'


#### Wishlist ####
# - Write a utility function to check if an 'output' folder exists and to create
#   one if it does not yet exist, write results to file, and print message on
#   how to read data back in R. That is now coded repeatedly.
# - Should I provide a function to install and / or check the status of RTools?
#   But the benefit of providing a link is that I don't need to update the info!
#   See devtools::find_rtools(debug = TRUE), cmdstanr::check_cmdstan_toolchain(),
#   and documentation at https://cran.r-project.org/bin/windows/Rtools/,
#   https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows,
#   https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#installing-rstan)


#### Utility functions ####
# Utility functions to write more succinct code in stopifnot()

is_logical <- function(x) {
  is.logical(x) && length(x) == 1L
}

# Includes test of nchar(x) > 0 to return FALSE in case of ""
all_characters <- function(x) {
  is.character(x) && length(x) > 0 && all(nchar(x) > 0)
}


#### select_libpath ####
# To do:
#   Check select_libpath() and implement it into prepare_install().
# Wishlist:
#   Let selection of library path depend on currently used R-version
select_libpath <- function() {
  warning("This function has not yet been checked.")
  if(grepl("windows", tolower(Sys.info()["sysname"]), fixed = TRUE)) {
    # The location to install R packages is chosen based on existing library
    # paths. A path containing 'Program Files' (case-insensitive) is used, or,
    # if such path does not exist, the path where R is installed. If that also
    # does not exist, the current working directory is used. Only the first
    # match is selected if multiple matches are present, to make downstream code
    # work.
    Rstring <- paste0("R-", paste(R.Version()[c("major", "minor")],
                                  collapse = "."))
    
    indices_path_current_R <- grep(tolower(Rstring), tolower(.libPaths()),
                                   fixed = TRUE)
    indices_path_program_files <- grep("program files", tolower(.libPaths()),
                                       fixed = TRUE)
    same_indices <- indices_path_current_R %in% indices_path_program_files
    if(any(same_indices)) {
      # Use the first path that contains both the current R version and "program
      # files"
      use_path <- .libPaths()[indices_path_current_R[same_indices[1]]]
    } else {
      # Use first path that contains 'program files'
      if(length(indices_path_program_files) > 0) {
        use_path <- .libPaths()[indices_path_program_files]
      } else {
        warning("Packages will be installed at ", lib, ".\nThis library path",
                " does not point to 'Program Files' because no such library",
                " path was found.\nIf the library path points to a network",
                " drive, installation of R packages and cmdstan might fail.")
      }
    }
    
    message("R packages will be installed at ", lib)
    
    index_first_OK_path <- index_path_program_files[1]
    if(length(index_first_OK_path) > 0) {
      lib <- .libPaths()[index_first_OK_path]
      
    } else {
      index_first_OK_path <- grep("/R/", .libPaths(), fixed = TRUE,
                                  ignore.case = FALSE)[1]
      if(length(index_first_OK_path) > 0) {
        
      } else {
        lib <- getwd()
        warning("Packages will be installed in the current working directory:\n",
                lib, ".\nThis library path does not point to 'Program Files'",
                " or a parent folder of an R installation\nbecause no such",
                " library paths were found.\nTherefor installation of R packages",
                " might fail.")
      }
    }
  } else {
    warning("This script was written using Windows, but you appear to be using",
            " another operating system.\nTherefore this script may fail. Check",
            " the installation instructions at https://cran.r-project.org/.")
  }
}


#### prepare_install ####
# Perform various checks to ensure the rest of the script can run: (1) if Rtools
#   utilities have been put on the search path, (2) if the library path 'lib' is
#   specified and matches the currently used version of R, and (3) if the
#   BiocManager package is installed and functional.
# Input:
#   None
# Return:
#   invisible NULL
# Side-effects:
#   The location of Rtools utilities is put on the search path if it is not
#     there yet.
#   The BiocManager package is installed if it is not installed and functional
# Notes:
#   It is sufficient to run this only once after (re)installing R.
prepare_install <- function() {
  if(grepl("windows", tolower(Sys.info()["sysname"]), fixed = TRUE) == FALSE) {
    warning("This script might fail because it is meant to be used on Windows,",
            " whereas you appear to be using ", Sys.info()["sysname"], 
            ".\nPlatform-specific installation instructions are available at",
            " https://cran.r-project.org/index.html.")
  } else {
    if(nchar(Sys.which("make")) == 0) {
      if(as.numeric(substr(R.version$minor, start = 1, stop = 1)) < 2) {
        # Put the location of Rtools utilities (bash, make, etc) on the search
        # path if it is not there yet.
        write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron",
              append = TRUE)
      }
      stop("Download Rtools from https://cran.r-project.org/bin/windows/Rtools/",
           " if you have not yet installed Rtools.\nYou are using ", R.version.string,
           ".\nThen restart R and run this function 'prepare_install()' again.",
           "\nSee https://cran.r-project.org/bin/windows/Rtools/",
           " for help on Rtools if this error keeps occuring.")
    }
  }
  
  # Check if the user-supplied library path 'lib' is specified, points to an R
  # library path, and matches the currently used version of R.
  rversion <- paste0("R-", as.character(getRversion()))
  msg_lib <- paste0("Specify global variable 'lib' containing a character",
                    " string giving the\nlibrary path where packages are or",
                    " should be installed by running the following line:\n",
                    "lib <- file.path(\"C:\", \"Program Files\", \"R\",",
                    " paste0(\"R-\", paste(R.Version()[c(\"major\", \"minor\")],",
                    " collapse = \".\")), \"library\")")
  if(!exists("lib")) {
    stop(msg_lib)
  } else {
    rversionpath <- regmatches(lib,
                               regexpr("R-[[:digit:]].[[:digit:]].[[:digit:]]",
                                       lib, ignore.case = TRUE, fixed = FALSE))
    if(length(rversionpath) == 0) {
      stop(msg_lib)
    } else {
      if(rversion != rversionpath) {
        stop("The R-version in the specified library path (", lib, ")\ndoes",
             " not correspond to the currently used R-version (", rversion,
             ")!\n", msg_lib)
      }
    }
  }
  
  # (re)install the BiocManager package from CRAN if it is not installed or not
  # functional
  if(!requireNamespace("BiocManager", lib.loc = lib, quietly = TRUE)) {
    message("Installing BiocManager package")
    install.packages("BiocManager", lib = lib, type = "binary")
    if(requireNamespace("BiocManager", quietly = TRUE)) {
      stop("Installed BiocManager package. Restart R session before proceeding.")
    } else {
      stop("Installation of the BiocManager package failed.",
           "\nIf the warning 'lib = \"", lib, "\" is not writeable' was issued,",
           "\nyou most likely forgot to run R as administrator.\nPlease close",
           " R, right-click on the R or RStudio icon, select 'Run as",
           " administrator',\nopen the 'InstallPkgs' R-project file, and try",
           " again.", call. = FALSE)
    }
  }
  
  message("Succesfully completed preparations.")
  invisible(NULL)
}


#### check_duplicates ####
# Check for duplicates within or across the various package lists.
# Input:
#   pkgs_lists: NULL, or a list containing (possibly named) character vectors of
#     package names to be checked for duplicates.
#   neglect_repos: logical indicating if the repository name should be excluded
#     when checking for duplicates, such that packages with the same name from
#     different repositories (e.g., 'pkgname' and 'repositoryname/pkgname') are
#     considered duplicates.
#   quietly: a logical. If FALSE (default), a message will be printed if no
#     duplicates were found.
# Return:
#   A list containing the names of duplicated packages if any are found 
#     (returned invisibly), and NULL otherwise.
# Side-effects:
#   Names of duplicate package are printed to the console if they are found,
#     with a warning.
check_duplicates <- function(pkgs_lists, neglect_repos = TRUE, quietly = FALSE) {
  stopifnot(is.list(pkgs_lists), all_characters(unlist(pkgs_lists)),
            is_logical(neglect_repos), is_logical(quietly))
  if(is.null(names(pkgs_lists))) {
    names(pkgs_lists) <- paste0("unnamed_list_entry_", seq_along(pkgs_lists))
  }
  
  githuburls <- grep("github", unlist(pkgs_lists, use.names = FALSE),
                     ignore.case = TRUE, value = TRUE)
  if(length(githuburls) > 0) {
    warning("Package name(s) '", paste(githuburls, collapse = "', '"),
            "'\nshould have the format username/repository instead of being",
            " full URLs pointing to GitHub")
  }
  
  checklist <- NULL
  for(index_first_element in seq_along(pkgs_lists)) {
    unlisted1 <- unlist(pkgs_lists[index_first_element], use.names = FALSE)
    if(neglect_repos == TRUE) {
      unlisted1 <- sub(".*/", "", unlisted1)
    }
    
    # Start index_second_element at index_first_element to prevent double
    # counting while including checks for duplicates within lists.
    for(index_second_element in (index_first_element:length(pkgs_lists))) {
      unlisted2 <- unlist(pkgs_lists[index_second_element], use.names = FALSE)
      if(neglect_repos == TRUE) {
        unlisted2 <- sub(".*/", "", unlisted2)
      }
      
      if(index_first_element == index_second_element) {
        # Check for duplicate names within package list
        add_to_checklist <- structure(
          list(unlisted1[anyDuplicated(unlisted1)]),
          names = paste0("duplicates_in_", names(pkgs_lists[index_first_element]))
        )
      } else {
        # Check for duplicate names across package lists
        add_to_checklist <- structure(
          list(unique(unlisted1[(unlisted1 %in% unlisted2)])),
          names = paste0(names(pkgs_lists[index_first_element]), "_in_",
                         names(pkgs_lists[index_second_element]))
        )
      }
      
      if(length(unlist(add_to_checklist, use.names = FALSE)) > 0) {
        checklist <- c(checklist, add_to_checklist)
      }
    }
  }
  
  if(length(checklist) > 0) {
    warning("Package lists contain duplicate package names.")
    print(checklist)
  } else {
    if(quietly == FALSE) {
      message("No duplicates were found in package lists.")
    }
  }
  
  invisible(checklist)
}


#### list_nonfunctional_pkgs ####
# Function to check for non-functional packages
# Input:
#   pkgs: a character vector, or list of character vectors, of package names to
#     be checked. Names of packages from GitHub can be in the format
#     username/repository or only the repository name.
#   save_file: a logical indicating if the names of non-functional packages
#     should be saved as a .txt-file, such that they can be obtained easily
#     after restarting the R-session.
#   sort: a logical indicating if the names of non-functional packages should be
#     sorted.
#   quietly: a logical indicating if printing the reason why packages are not
#     functioning should be suppressed.
#   verbose: a logical indicating if other warnings issued when packages are
#     loaded should be printed.
# Return:
#   a character vector containing the names of non-functional packages, sorted
#     if argument sort is TRUE, returned invisibly.
# Side-effects:
#   Packages are loaded, with the result that updating packages might fail if
#     they are not unloaded first. A message is printed urging to restart R
#     before continuing to prevent this.
#   The names of non-functional packages are printed to the console, and, if
#     save_file = TRUE, saved as a text-file inside the subfolder 'output'.
# Notes:
#   This function uses requireNamespace() instead of installed.packages(),
#     because installed.packages() does not check if packages are functional
#     nor if all needed dependencies are installed and functional. In addition,
#     installed.packages() can be slow such that its help page states that
#     requireNamespace() or require() should be used instead.
# Wishlist:
#   When testing if packages are functioning correctly, differentiate between
#     missing and non-functioning packages.
list_nonfunctional_pkgs <- function(pkgs, save_file = TRUE, sort = TRUE,
                                    quietly = FALSE, verbose = FALSE) {
  if(is.list(pkgs)) {
    pkgs <- unlist(pkgs, use.names = FALSE)
  }
  
  stopifnot(all_characters(pkgs), is_logical(save_file), is_logical(sort),
            is_logical(quietly), is_logical(verbose))
  
  pkgs <- unique(pkgs)
  pkgs_input <- pkgs
  pkgs <- sub(".*/", "", pkgs) # Package name is the part after the last '/'
  
  if(verbose == FALSE) {
    index_nonfunctional <- suppressWarnings(suppressPackageStartupMessages(
      which(!sapply(X = pkgs, FUN = requireNamespace, lib.loc = lib,
                    quietly = quietly, simplify = TRUE))
    ))
  } else {
    index_nonfunctional <- suppressPackageStartupMessages(
      which(!sapply(X = pkgs, FUN = requireNamespace, lib.loc = lib,
                    quietly = quietly, simplify = TRUE))
    )
  }
  
  nonfunctional_pkgs <- pkgs_input[index_nonfunctional]
  
  if(sort == TRUE) {
    nonfunctional_pkgs <- sort(nonfunctional_pkgs)
  }
  
  if(length(nonfunctional_pkgs) > 0) {
    message("Names of non-functional packages:")
    dput(nonfunctional_pkgs)
    
    if(save_file == TRUE) {
      file_name <- paste0("nonfunc_pkgs_",
                          format(Sys.time(), format = "%Y_%m_%d_%H_%M"),
                          "_", paste0("R", as.character(getRversion())), ".txt")
      dir_path <- file.path(".", "output")
      if(!dir.exists(dir_path)){
        dir.create(dir_path)
      }
      read_back_path <- file.path(dir_path, file_name)
      message_read_back <- paste0("\nTo read the package names back into R use",
                                  "\nnonfunctional_pkgs <- dget(\"",
                                  read_back_path, "\")")
      if(file.exists(read_back_path)) {
        warning("Textfile with names of non-functional packages already exists,",
                " not saved again!")
        message(message_read_back)
      } else {
        dput(nonfunctional_pkgs, file = read_back_path)
        message("Textfile with names of non-functional packages saved as ",
                file_name, "\nin ", file.path(getwd(), "output"), message_read_back)
      }
    }
  } else {
    message("All packages are functional")
  }
  message("Restart R before continuing, to prevent problems arising from",
          " updating loaded packages")
  invisible(nonfunctional_pkgs)
}


#### check_status ####
# Function to check the status of installed packages
# Input:
#   checkBuilt: a logical of length one. If TRUE, a package built under an
#     earlier major.minor version of R (e.g., 3.4) is considered to be old.
#   type: a character vector indicating the type of available package (e.g.,
#     binary, source) to check validity against.
#   save_file: a logical indicating if the details of invalid packages (i.e.,
#     packages that are outdated or too new) should be saved in a .csv-file,
#     such that they can be obtained easily after restarting the R-session.
#   print_output: a character vector indicating what information on invalid
#     packages should be printed to the console: their details, a character
#     vector with their names, both of these, or none.
# Return:
#   A list (returned invisibly) containing a character vector with the names of
#     packages that are outdated or too new, and a matrix containing their
#     details. Both entries in the list are NULL if all packages are valid.
# Side-effects:
#   Packages are loaded, with the result that updating packages might fail if
#     they are not unloaded first. A message is printed urging to restart R
#     before continuing to prevent this.
#   The details of invalid packages are saved as a .csv file inside the
#     subfolder 'output' if the argument 'save_file' is TRUE.
check_status <- function(checkBuilt = TRUE,
                         type = c("binary", "both", "source", "win.binary"),
                         save_file = TRUE,
                         print_output = c("both", "pkgs_details", "pkgs_names",
                                          "none")) {
  type <- match.arg(type, several.ok = FALSE)
  print_output <- match.arg(print_output, several.ok = TRUE)
  
  stopifnot(is_logical(checkBuilt), is_logical(save_file))
  if(length(print_output) > 1) {
    print_output <- print_output[1]
    warning("Only the first element (\'", print_output, "\') of the argument ",
            "'print_output' will be used.")
  }
  
  valid_out <- BiocManager::valid(checkBuilt = checkBuilt, type = type)
  out_of_date_names <- NULL
  too_new_names <- NULL
  invalid_details <- NULL
  # valid() returns TRUE if all packages are valid, not an empty list
  if(length(valid_out) > 1) {
    if(length(valid_out$out_of_date) > 0) {
      out_of_date_names <- dimnames(valid_out$out_of_date)[[1]]
      cols_specs <- c("Package", "Installed", "ReposVer", "Built")
      invalid_details <- valid_out$out_of_date[, cols_specs, drop = FALSE]
      if(any(print_output %in% c("both", "pkgs_details"))) {
        message("Some packages are out-of-date:")
        print(invalid_details)
      }
    }
    if(length(valid_out$too_new) > 0) {
      too_new_names <- dimnames(valid_out$too_new)[[1]]
      too_new_specs <- data.frame(Package = rownames(valid_out$too_new),
                                  Installed = valid_out$too_new[, "Version"],
                                  ReposVer = NA, Built = NA)
      invalid_details <- rbind(invalid_details, too_new_specs)
      if(any(print_output %in% c("both", "pkgs_details"))) {
        message("Some packages are too new:")
        print(valid_out$too_new[, c("Version"), drop = FALSE])
      }
    }
    
    invalid_names <- c(out_of_date_names, too_new_names)
    
    if(save_file == TRUE) {
      rownames(invalid_details) <- NULL
      file_name <- paste0("invalid_pkgs_",
                          format(Sys.time(), format = "%Y_%m_%d_%H_%M"),
                          "_", paste0("R", as.character(getRversion())), ".csv")
      dir_path <- file.path(".", "output")
      if(!dir.exists(dir_path)){
        dir.create(dir_path)
      }
      read_back_path <- file.path(dir_path, file_name)
      message_read_back <- paste0("\nTo read the information back into R use",
                                  "\ninvalid_pkgs <- read.csv(\"",
                                  read_back_path, "\")")
      if(file.exists(read_back_path)) {
        warning("File with invalid packages already exists, not saved again!")
        message(message_read_back)
      } else {
        write.csv(invalid_details, file = read_back_path, row.names = FALSE)
        message("A .csv file with details of invalid packages saved as ",
                file_name, "\nat ", file.path(getwd(), "output"), message_read_back)
      }
    }
    
    if(any(print_output %in% c("both", "pkgs_names"))) {
      message("Package names of invalid packages:")
      dput(invalid_names)
    }
    warning("Restart R before continuing, to prevent problems when updating",
            " loaded packages")
    invisible(list(names = invalid_names, details = invalid_details))
  } else {
    warning("All packages are up-to-date. Restart R before continuing,\n",
            "to prevent problems when updating loaded packages")
    invisible(list(names = NULL, details = NULL))
  }
}


#### list_dependencies ####
# Function to list package dependencies and the number of dependencies.
# Input:
#   pkgs: character vector of package names. Names of packages from GitHub can
#     be in the format username/repository or only the repository name.
#   deps_type: character vector indicating the type of dependencies. See the
#     description of the argument 'which' at help(tools::package_dependencies)
#     for details.
#   recursive: logical indicating if recursive dependencies should be included
#   name_per_pkg: logical indicating if the names of dependencies should be
#     given separately for each package in pkgs.
#   number_per_pkg: logical indicating if the number of dependencies should
#     be given separately for each package in pkgs.
#   name_total: logical indicating if the names of the total set of dependencies
#     should be given.
#   add_pkgs_to_total: logical indicating if the package names for which the
#     dependencies are requested should be added themselves to the list and
#     number of total dependencies.
#   exclude_high_prio: logical indicating if the high priority packages
#     should be omitted from the total counts and names of dependencies
#   exclude_pkgs: NULL, a character vector (or a list containing character
#     vectors) with package names to be omitted from the listed dependencies.
#     Note that dependencies of those packages will not be excluded.
#   sort_ndeps_by: character vector indicating if the number of dependencies per
#     per package should be sorted on package names ("names") or on their number 
#     of dependencies ("ndeps"). The other elements of the returned list are 
#     are sorted on package names irrespective of the value of this argument.
# Return
#   A list with for each package in 'pkgs' the names and number of dependencies,
#     the names and number of all dependencies, and the number of unique
#     packages supplied in the arguments 'exclude_pkgs' and 'exclude_high_prio'.
#     Note that the latter does NOT indicate the number of packages that were
#     actually excluded.
# WARNING:
#   Information about dependencies for non-CRAN packages is currently NOT
#     correctly handled. In those cases 'NULL' is given as dependency, whereas
#     integer(0) indicates the package does not have any dependencies. As a
#     consequence, this function returns too few dependencies for non-CRAN
#     packages.
# To do:
#   tools::package_dependencies() requires internet access? Rewrite to use
#     info from local package versions instead?
#   Separately handle packages with no dependencies (returning integer(0)) from
#     packages that were not found in the database such that no information on
#     dependencies could be obtained (returning NULL). Alternatively, use
#     utils::packageDescription() or utils::installed.packages()?
list_dependencies <- function(pkgs, deps_type = "strong", recursive = TRUE,
                              name_per_pkg = FALSE, number_per_pkg = TRUE,
                              name_total = TRUE, add_pkgs_to_total = FALSE,
                              exclude_high_prio = FALSE, exclude_pkgs = NULL,
                              sort_ndeps_by = c("ndeps", "names")) {
  sort_ndeps_by <- match.arg(sort_ndeps_by, several.ok = FALSE)
  if(is.list(pkgs)) {
    pkgs <- unlist(pkgs, use.names = FALSE)
  }
  if(is.list(exclude_pkgs)) {
    exclude_pkgs <- unlist(exclude_pkgs, use.names = FALSE)
  }
  
  stopifnot(all_characters(pkgs), all_characters(deps_type),
            is_logical(recursive), is_logical(name_per_pkg),
            is_logical(number_per_pkg), is_logical(name_total),
            is_logical(add_pkgs_to_total),
            is_logical(exclude_high_prio),
            is.null(exclude_pkgs) || all_characters(exclude_pkgs))
  
  pkgs <- sub(".*/", "", pkgs) # Package name is the part after the last '/'

  if(sort_ndeps_by == "names") {
    pkgs <- sort(pkgs)
  }
  
  if(exclude_high_prio == TRUE) {
    if(!exists("high_prio_pkgs")) {
      message("Collecting names of installed high priority packages because ",
              "the set 'high_prio_pkgs' was not loaded.\nNames of high ",
              "priority packages that are not installed will not be excluded ",
              "from the dependencies.")
      high_prio_pkgs <- unname(installed.packages(priority = "high")[, "Package"])
    }
    
    exclude_pkgs <- unique(c(high_prio_pkgs, exclude_pkgs))
  } else {
    exclude_pkgs <- unique(exclude_pkgs)
  }
  
  deps_per_pkgs <- tools::package_dependencies(pkgs, which = deps_type,
                                               recursive = recursive)
  pkgs_not_found <- NULL
  for(index_pkgs in seq_along(deps_per_pkgs)) {
    # NULL indicates the package was not found in the database such that no
    # information about dependencies is obtained, whereas integer(0) indicates
    # the package does not have any dependencies.
    pkg_deps <- unlist(deps_per_pkgs[index_pkgs], use.names = FALSE)
    if(length(pkg_deps) > 0) {
      deps_per_pkgs[index_pkgs][[1]] <- sort(pkg_deps[!(pkg_deps %in% exclude_pkgs)])
    }
  }
  
  deps_total <- unlist(deps_per_pkgs, use.names = FALSE)
  if(add_pkgs_to_total == TRUE) {
    deps_total <- c(deps_total, pkgs)
  }
  deps_total <- unique(deps_total)
  deps_total <- sort(deps_total[!(deps_total %in% exclude_pkgs)])
  
  out <- NULL
  if(name_per_pkg == TRUE) {
    out <- c(out, deps_per_pkgs)
  }
  
  if(number_per_pkg == TRUE) {
    ndeps_per_pkgs <- lengths(deps_per_pkgs)
    if(sort_ndeps_by == "ndeps") {
      ndeps_per_pkgs <- sort(ndeps_per_pkgs)
    }
    out <- c(out, list(ndeps_per_pkgs = ndeps_per_pkgs))
  }
  
  if(name_total == TRUE) {
    out <- c(out, list(deps_total = sort(deps_total)))
  }
  
  c(out, list(ndeps_total = length(deps_total),
              length_excl = length(exclude_pkgs)))
}


#### save_details ####
# Function to save details of installed packages as a .csv file, for example to
#   compare installed sets of packages across multiple machines or over time.
# Input:
#  - PC_name: NULL or character vector of length 1 giving the name of the PC to
#    be added to the name of the .csv-file.
# Return:
# - A matrix containing the details of the installed packages is returned
#   invisible.
# Side effects:
# - A .csv-file giving details of installed packages is saved inside the
#     subfolder 'output'.
save_details <- function(PC_name = "desktop") {
  stopifnot((length(PC_name) == 1L && is.character(PC_name)) || is.null(PC_name))
  PC_name <- gsub("[^[:alnum:]_]", "_", PC_name)
  
  installed_pkgs <- installed.packages()[, c("Package", "Version", "Built",
                                             "NeedsCompilation", "Priority")]
  rownames(installed_pkgs) <- NULL
  file_name <- paste0("installed_pkgs", PC_name, "_",
                      format(Sys.Date(), format = "%Y_%m_%d"),
                      "_", paste0("R", as.character(getRversion())), ".csv")
  dir_path <- file.path(".", "output")
  if(!dir.exists(dir_path)){
    dir.create(dir_path)
  }
  read_back_path <- file.path(dir_path, file_name)
  message_read_back <- paste0("\nTo read the information back into R use",
                              "\ninstalled_pkgs <- read.csv(\"", read_back_path,
                              "\")")
  if(file.exists(read_back_path)) {
    warning("File with details of installed packages already exists,",
            " not saved again!")
    message(message_read_back)
  } else {
    write.csv(installed_pkgs, file = file.path(dir_path, file_name),
              row.names = FALSE)
    message("A .csv-file with details of installed packages saved as ",
            file_name, "\nat ", file.path(getwd(), "output"), message_read_back)
  }
  invisible(installed_pkgs)
}


message("Sourced script containing the functions to check and install R-packages.")
