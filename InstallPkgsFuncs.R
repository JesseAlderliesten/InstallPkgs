#### Introduction #####
# Functions to install R-packages and obtain information about packages.


#### Utility functions ####
# Utility functions to write more succinct code in stopifnot()

is_logical <- function(x) {
  is.logical(x) && length(x) == 1L
}

# Includes test of nchar(x) > 0 to return FALSE in case of ""
all_characters <- function(x) {
  is.character(x) && length(x) > 0 && all(nchar(x) > 0)
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
# Note:
#   It is sufficient to run this only once after (re)installing R.
prepare_install <- function() {
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
  
  # Check if the user-supplied library path 'lib' is specified, points to an R
  # library path, and matches the currently used version of R.
  rversion <- paste0("R-", as.character(getRversion()))
  msg_lib <- paste0("Specify global variable 'lib' containing a character",
                    " string giving the\nlibrary path where packages are or",
                    " should be installed by running the following line:\n",
                    "lib <- file.path(\"C:\", \"Program Files\", \"R\",",
                    "paste0(\"R-\", paste(R.Version()[c(\"major\", \"minor\")],",
                    "collapse = \".\")), \"library\")")
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
  
  # (re)install the BiocManager package if it is not installed and functional
  if(!requireNamespace("BiocManager", lib.loc = lib, quietly = TRUE)) {
    message("Installing BiocManager package")
    install.packages(Bioc_missing, lib = lib, type = "source")
    stop("Installed BiocManager package. Restart R session before proceeding.")
  }
  
  message("Succesfully completed preparations.")
  invisible(NULL)
}


#### choose_mirrors ####
# Choose CRAN and BioConductor mirrors based on selected countries.
# Input:
#   countries: Countries to select mirrors from, with most relevant first
#   databases: databases to set mirrors for. "BioConductor", "CRAN", or both.
# Return:
#   invisible NULL
# Side-effects:
#   The CRAN and BioConductor mirrors are set
# Note:
#   The default mirror will be used if no mirror is available from any of the
#     provided countries, with a warning. 
choose_mirrors <- function(countries = c("Belgium", "Germany"),
                           databases = c("BioConductor", "CRAN")) {
  databases <- match.arg(databases, several.ok = TRUE)
  stopifnot(all_characters(countries))
  
  for(database in databases) {
    if(database == "BioConductor") {
      BioC <- read.csv(paste0(R.home(), "/doc/BioC_mirrors.csv"))[, "Country"]
      indices_matched <- which(BioC %in% countries)
    }
    
    if(database == "CRAN") {
      CRAN <- getCRANmirrors()[, "Country"]
      indices_matched <- which(CRAN %in% countries)
    }
    
    if(length(indices_matched) > 0) {
      mirror_index <- indices_matched[1]
    } else {
      mirror_index <- 1
      warning("No ", database, " mirror available for any of the provided",
              " countries (", paste0(countries, collapse = ", "), ").\nDefault",
              " mirror will be used instead")
    }
    
    if(database == "BioConductor") {
      chooseBioCmirror(ind = mirror_index, local.only = TRUE)
    }
    
    if(database == "CRAN") {
      chooseCRANmirror(ind = mirror_index, local.only = TRUE)
    }
  }
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
#     (returned invisible), and NULL otherwise.
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
  for(index1 in seq_along(pkgs_lists)) {
    unlisted1 <- unlist(pkgs_lists[index1], use.names = FALSE)
    if(neglect_repos == TRUE) {
      unlisted1 <- sub(".*/", "", unlisted1)
    }
    
    # Start index2 at index1 to prevent double counting while including checks
    # for duplicates within lists.
    for(index2 in (index1:length(pkgs_lists))) {
      unlisted2 <- unlist(pkgs_lists[index2], use.names = FALSE)
      if(neglect_repos == TRUE) {
        unlisted2 <- sub(".*/", "", unlisted2)
      }
      
      if(index1 == index2) {
        # Check for duplicate names within package list
        add_to_checklist <- structure(
          list(unlisted1[anyDuplicated(unlisted1)]),
          names = paste0("duplicates_in_", names(pkgs_lists[index1]))
        )
      } else {
        # Check for duplicate names across package lists
        add_to_checklist <- structure(
          list(unique(unlisted1[(unlisted1 %in% unlisted2)])),
          names = paste0(names(pkgs_lists[index1]), "_in_",
                         names(pkgs_lists[index2]))
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
#   sort: a logical indicating if the names of nonfunctioning packages should be
#     sorted.
#   quietly: a logical indicating if printing the reason why packages are not
#     functioning should be suppressed, regarded to be FALSE if verbose is TRUE.
#   verbose: a logical indicating if other warnings issued when packages are
#     loaded should be printed. The argument 'quietly' is regarded to be FALSE
#     if verbose is TRUE.
# Return:
#   a character vector containing the names of non-functional packages, sorted
#     if argument sort is TRUE, returned invisible.
# Side-effects:
#   Packages are loaded, with the result that updating packages might fail if
#     they are not unloaded first. A message is printed urging to restart R
#     before continuing to prevent this.
#   The names of non-functional packages are printed to the console, and, if
#     save_file = TRUE, saved as a .txt-file.
# Note:
#   installed.packages() does not check if packages are functional, nor if all
#     needed dependencies are installed and functional. In addition, it can be
#     slow such that ?installed.packages() states that requireNamespace() or
#     require() should be used instead.
list_nonfunctional_pkgs <- function(pkgs, save_file = FALSE, sort = TRUE,
                                    quietly = FALSE, verbose = FALSE) {
  if(is.list(pkgs)) {
    pkgs <- unlist(pkgs, use.names = FALSE)
  }
  
  stopifnot(all_characters(pkgs), is_logical(save_file), is_logical(sort),
            is_logical(quietly), is_logical(verbose))
  
  pkgs <- unique(pkgs)
  pkgs_input <- pkgs
  pkgs <- sub(".*/", "", pkgs) # Package name is the part after the last '/'
  
  index_nonfunctional <- if(verbose == FALSE) {
    suppressWarnings(suppressPackageStartupMessages(
      which(!sapply(X = pkgs, FUN = requireNamespace, lib.loc = lib,
                    quietly = quietly, simplify = TRUE))
    ))
  } else {
    suppressPackageStartupMessages(
      which(!sapply(X = pkgs, FUN = requireNamespace, lib.loc = lib,
                    quietly = FALSE, simplify = TRUE))
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
      dput(nonfunctional_pkgs, file = file_name)
      message("Textfile with names of non-functional packages saved as ",
              file_name, "\nat ", getwd(), "\nUse nonfunctional_pkgs <- dget(\"",
              file_name, "\") to read the package names back into R.")
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
#   print_output: a character vector indicating what information on nonvalid
#     packages should be printed to the console: their details, a character
#     vector with their names, both of these, or none.
# Return:
#   A list (returned invisible) containing a character vector with the names of
#     packages that are outdated or too new, and a matrix containing their
#     details. Both entries in the list are NULL if all packages are valid.
# Side-effects:
#   Packages are loaded, with the result that updating packages might fail if
#     they are not unloaded first. A message is printed urging to restart R
#     before continuing to prevent this.
#   The details of invalid packages are saved in a .csv file if the argument
#     'save_file' is TRUE.
check_status <- function(checkBuilt = TRUE,
                         type = c("binary", "both", "source", "win.binary"),
                         save_file = FALSE,
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
  # NOTE: valid() returns TRUE if all packages are valid, not an empty list
  if(length(valid_out) > 1) {
    if(length(valid_out$out_of_date) > 0) {
      out_of_date_names <- dimnames(valid_out$out_of_date)[[1]]
      cols_specs <- c("Package", "Installed", "ReposVer", "Built")
      invalid_details <- valid_out$out_of_date[, cols_specs, drop = FALSE]
      if(any(print_output %in% c("both", "pkgs_details"))) {
        print("Some packages are out-of-date:")
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
        print("Some packages are too new:")
        print(valid_out$too_new[, c("Version"), drop = FALSE])
      }
    }
    
    invalid_names <- c(out_of_date_names, too_new_names)
    
    if(save_file == TRUE) {
      rownames(invalid_details) <- NULL
      file_name <- paste0("invalid_pkgs_",
                          format(Sys.time(), format = "%Y_%m_%d_%H_%M"),
                          "_", paste0("R", as.character(getRversion())), ".csv")
      write.csv(invalid_details, file = file.path(getwd(), file_name),
                row.names = FALSE)
      message("A .csv file with details of invalid packages saved as ",
              file_name, "\nat ", getwd(), "\nUse invalid_pkgs <- ",
              "read.csv(\"", file.path(getwd(), file_name), "\")\nto read the",
              " information back into R.")
    }
    
    if(any(print_output %in% c("both", "pkgs_names"))) {
      print("Package names of invalid packages:")
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
# - A .csv-file giving details of installed packages is saved.
save_details <- function(PC_name = "desktop") {
  stopifnot((length(PC_name) == 1L && is.character(PC_name)) || is.null(PC_name))
  PC_name <- gsub("[^[:alnum:]_]", "_", PC_name)
  
  installed_pkgs <- installed.packages()[, c("Package", "Version", "Built",
                                             "NeedsCompilation", "Priority")]
  rownames(installed_pkgs) <- NULL
  file_name <- paste0("installed_pkgs", PC_name, "_", format(Sys.Date(), format = "%Y_%m_%d"),
                      "_", paste0("R", as.character(getRversion())), ".csv")
  write.csv(installed_pkgs, file = file.path(getwd(), file_name),
            row.names = FALSE)
  message("A .csv-file with details of installed packages saved as ",
          file_name, "\nat ", getwd(), "\nUse installed_pkgs <- ",
          "read.csv(\"", file.path(getwd(), file_name), "\")\nto read the",
          " information back into R.")
  invisible(installed_pkgs)
}

message("The script containing the functions to install R-packages has been sourced.")
