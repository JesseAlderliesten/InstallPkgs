#### Introduction #####
# Functions to install R-packages and obtain information about installed
# packages. Also includes some functions to check the installation of R.


#### Wishlist ####
# - Write a utility function to check if an 'output' folder exists and to create
#   one if it does not yet exist, write results to file, and print message on
#   how to read data back in R. That is now coded repeatedly. See my function
#   create_directory() in utils.R from the UvA.
# - Write tests for the functions.


#### Utility functions ####
# Utility functions to write more succinct code in stopifnot()
is_logical <- function(x) {
  is.logical(x) && length(x) == 1L
}

all_characters <- function(x, allow_empty_char = FALSE, allow_zero_char = FALSE,
                           allow_NA_char = FALSE) {
  # is.character(x) is TRUE for NA_character_, but FALSE for NULL, NA and NaN.
  is.character(x) &&
    (allow_zero_char == TRUE || length(x) > 0) &&
    (allow_empty_char == TRUE || !any(x == "", na.rm = TRUE)) &&
    (allow_NA_char == TRUE || !anyNA(x))
}


#### get_paths ####
# Get paths where packages are (or should be) installed.
# Input:
# - path: character(0) (default) or a character vector indicating paths.
# - quietly: logical of length 1 (default FALSE) indicating if messages should
#   be suppressed.
# Notes:
# - It is NOT checked if paths supplied in argument 'path' are valid file paths.
# - A warning is issued if the working directory is returned as element
#   'first_path', because that implies that no path was provided in argument
#   'path' and no paths are present in '.libPath'.
# Return:
# - A list with five elements, all containing characters (possibly character(0)):
#   'first_path' with the first non-empty path, 'argument_paths' with the value
#   of argument 'path', 'Rversion_paths' with the paths from .libPath that
#   contain the current R version number, 'other_paths' with the paths from
#   .libPath that do not contain the current R version number, and 'wd_path'
#   with the working directory. Elements for which no path is found are set to
#   character(0).
# Programming notes:
# - To implement a check if paths supplied in argument 'path' are valid file
#   paths, see normalizePath(), path.expand(), checkmate::testPathForOutput(),
#   fs::path_sanitize(), https://en.wikipedia.org/wiki/Filename#In_Windows,
#   and https://learn.microsoft.com/en-us/dotnet/standard/io/file-path-formats.
get_paths <- function(path = character(0), quietly = FALSE) {
  stopifnot(all_characters(path, allow_zero_char = TRUE), is_logical(quietly))
  
  paths_libPath <- .libPaths()
  bool_paths_Rversion <- grepl(pattern = paste0("R-", as.character(getRversion())),
                               x = paths_libPath, fixed = TRUE)
  
  # Using a list because elements other than 'wd_path' can have lengths larger
  # than one. The selection of paths returns character(0) if no path is present,
  # because 'paths_libPath' is character.
  paths_possible <- list(argument_paths = path,
                         Rversion_paths = paths_libPath[bool_paths_Rversion],
                         other_paths = paths_libPath[!bool_paths_Rversion],
                         wd_path = getwd())
  
  # With these settings all_characters() only returns 'TRUE' for non-empty,
  # non-NA_character_ character strings containing more than one character. The
  # first such element in paths_possible is selected, and the first element of
  # that string is included as element 'first_path' in the returned list.
  path_first_OK <- which(vapply(X = paths_possible, FUN = all_characters,
                                FUN.VALUE = logical(1), allow_empty_char = FALSE,
                                allow_zero_char = FALSE, allow_NA_char = FALSE,
                                USE.NAMES = FALSE))[1]
  path_first_name <- names(paths_possible[path_first_OK])
  path_first <- paths_possible[path_first_OK][[1]][1]
  
  if(path_first_name == "wd_path") {
    warning("Returning the working directory ('", path_first,
            "')\nas first path because no path was supplied in argument 'path'",
            " and no path was\nfound in '.libPaths'. Installation of R",
            " packages might fail if the path points\nto a network drive.")
  } else {
    if(!quietly) {
      message("Path '", path_first,
              "' is element 'first_path' of the returned list.")
    }
  }
  c(list(first_path = path_first), paths_possible)
}


#### check_OS_is_Windows ####
# Check if the operating system is windows.
# Input:
# - on_error: character string containing 'warn', 'message' or 'quiet' that
#   determines how the message if the operating system is not Windows is
#   displayed.
# Return:
# - A logical value indicating if the operating system is Windows, returned
#   invisibly.
check_OS_is_Windows <- function(on_error = c("warn", "message", "quiet")) {
  on_error <- match.arg(on_error, several.ok = FALSE)
  # Using 'Sys.info()' which returns information about the platform R is running
  # on, not 'R.version()' which returns information about the platform R was
  # built on.
  system_name <- Sys.info()["sysname"]
  
  # Need tolower() because 'ignore_case' does not work if 'fixed' is TRUE.
  if(grepl("windows", tolower(system_name), fixed = TRUE)) {
    OS_is_Windows <- TRUE
  } else {
    OS_is_Windows <- FALSE
    text <- paste0("This script might fail because it is meant to be used on",
                   " Windows, whereas you\nare using '", system_name,
                   "' instead. See the platform-specific installation",
                   " instructions\nat https://cran.r-project.org/ and the",
                   " 'R installation and administration manual'\nat",
                   " https://cran.r-project.org/doc/manuals/r-release/R-admin.html",
                   " for help.")
    
    switch(on_error,
           message = message(text),
           quiet = NULL,
           warning(text))
  }
  invisible(OS_is_Windows)
}


#### prepare_install ####
# Perform checks to ensure the rest of the script can run. Checks (1) if Rtools
# utilities have been put on the search path (an error occurs with a message
# proposing steps to fix it if they are not yet on the path); (2) if the global
# variable 'lib_path' is specified (otherwise an error occurs) and that it
# matches the used version of R (otherwise a warning is issued); and (3) if the
# BiocManager package is installed and functional.
# Input:
# - None.
# Return:
# - invisible NULL
# Side-effects:
# - For R versions 4.0.0 - 4.1.3, the location of Rtools utilities is put on the
#   search path if it is not there yet.
# - The BiocManager package is installed if it is not installed and functional.
# Wishlist:
# - Check https://stackoverflow.com/questions/26244530
prepare_install <- function() {
  if(check_OS_is_Windows(on_error = "warn") == TRUE) {
    if(nchar(Sys.which("make")) == 0) {
      # Put the location of Rtools (the C++ Toolchain) utilities (e.g., bash,
      # make) on the search path if it is not there yet for R versions 4.0.0 to
      # 4.1.3.
      if(as.integer(R.version$major) == 4L &&
         as.integer(substr(R.version$minor, start = 1, stop = 1)) < 2L) {
        write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron",
              append = TRUE)
      }
      stop("Utilities for RTools (the toolchain bundle used for building R",
           " packages\nthat need compilation of C/C++ or Fortran code from",
           " source) were not yet on the\nsearch path. If you have not yet",
           " done so, download the version of Rtools that\ncorresponds to the",
           " R version you are using (", R.version.string, ") from",
           "\nhttps://cran.r-project.org/bin/windows/Rtools/ and install it.",
           " If you are using\nR versions 4.0.0 - 4.1.3, you have to restart R",
           " afterwards and run this\nfunction 'prepare_install()' again to",
           " put the location of Rtools utilities\n(bash, make, etc) on the",
           " search path. For instructions see the link given above,\nor",
           " install package 'pkgbuild' and run",
           " pkgbuild::check_build_tools(debug = TRUE)\nto check if Rtools is",
           " set up correctly and get instructions how to fix it if\nnot.")
    }
  }
  
  # Check if the variable 'lib_path' is specified, if any of its strings contain
  # an R version number, and if that number matches the used version of R. The
  # command to obtain all paths in the last line of the message contains
  # unique(unlist(...)) to remove potential duplicated entries.
  msg_lib <- paste0(" global variable 'lib_path' containing character strings",
                    " with the paths\nwhere packages are (or should be)",
                    " installed by running the following line:",
                    "\nlib_path <- get_paths()$first_path\nAlternatively, run",
                    " the following line to use all paths returned by",
                    " get_paths():\nlib_path <- unique(unlist(get_paths(quietly",
                    " = TRUE), use.names = FALSE))")
  if(!exists("lib_path")) {
    stop(paste0("Specify", msg_lib))
  }
  lib_path <- unlist(lib_path, use.names = FALSE)
  rversionpath <- regmatches(lib_path,
                             regexpr("R-[[:digit:]].[[:digit:]].[[:digit:]]",
                                     lib_path, ignore.case = TRUE, fixed = FALSE))
  rversion <- paste0("R-", as.character(getRversion()))
  if(length(rversionpath) == 0) {
    warning("None of the specified path(s) (", paste0(lib_path, collapse = ",\n"),
            ")\ncontains an R version number. You might want to\nspecify",
            msg_lib)
  }
  if(!any(rversionpath == rversion)) {
    warning("The version number in none of the specified paths\n(",
            paste0(lib_path, collapse = ",\n"), ")\ncorresponds to the version",
            " of R you are using (", rversion, ")!\nYou might want to\nspecify",
            msg_lib)
  }
  
  # (re)install the BiocManager package from CRAN if it is not installed or not
  # functional
  if(!requireNamespace("BiocManager", lib.loc = lib_path, quietly = TRUE)) {
    message("Trying to install package 'BiocManager' in ", lib_path)
    install.packages("BiocManager", lib.loc = lib_path, lib = lib_path,
                     type = "binary")
    if(requireNamespace("BiocManager", lib.loc = lib_path, quietly = TRUE)) {
      stop("If no error message is printed below, the BiocManager package",
           "  was\nsuccesfully installed. Restart R session before proceeding.")
    } else {
      stop("Installation of the BiocManager package failed.\nIf a warning like",
           " 'lib = \"", lib_path[[1]][1], "\" is not writeable'\nwas issued,",
           " you most likely forgot to run R as administrator, or used a wrong",
           " path (the warnings printed below might point to that).\nClose R",
           " and\nrestart R as administrator (e.g., right-click on the R or",
           " RStudio icon, select\n'Run as administrator', open the",
           " 'InstallPkgs' R-project file, and try again.")
    }
  }
  
  message("Successfully completed preparations.")
  invisible(NULL)
}


#### check_duplicates ####
# Check for duplicates within or across the various package lists.
# Input:
#   pkgs_lists: a list containing (possibly named) character vectors of package
#     names to be checked for duplicates.
#   neglect_repos: logical indicating if the repository name should be excluded
#     when checking for duplicates, such that packages with the same name from
#     different repositories (e.g., 'pkgname' and 'repositoryname/pkgname') are
#     considered duplicates.
#   quietly: a logical (default FALSE) indicating if the message that is printed
#     if no duplicates are found should be suppressed.
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
  } else {
    bool_ind_empty_names <- nchar(names(pkgs_lists), type = "width") == 0
    if(any(bool_ind_empty_names)) {
      names(pkgs_lists)[bool_ind_empty_names] <- paste0(
        "unnamed_list_entry_", which(bool_ind_empty_names))
    }
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
    warning("Package lists contain duplicate package names, see the printed",
            " output.", call. = FALSE)
    print(checklist)
  } else {
    if(quietly == FALSE) {
      message("No duplicates were found in package lists.")
    }
  }
  
  invisible(checklist)
}


#### find_nonfunctional_pkgs ####
# Function to check for non-functional packages
# Input:
#   pkgs: a character vector, or list of character vectors, of package names to
#     be checked. Names of packages from GitHub can be in the format
#     username/repository or only the repository name.
#   lib: character vector giving the paths where packages are installed.
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
#     save_file = TRUE, saved as a text-file inside the subfolder 'output' of
#     the working directory.
# Notes:
#   This function uses requireNamespace() instead of installed.packages(),
#     because installed.packages() does not check if packages are functional
#     nor if all needed dependencies are installed and functional. In addition,
#     installed.packages() can be slow such that its help page states that
#     requireNamespace() or require() should be used instead.
# Wishlist:
#   When testing if packages are functioning correctly, differentiate between
#     missing and non-functioning packages.
find_nonfunctional_pkgs <- function(pkgs, lib, save_file = TRUE, sort = TRUE,
                                    quietly = FALSE, verbose = FALSE) {
  if(is.list(pkgs)) {
    pkgs <- unlist(pkgs, use.names = FALSE)
  }
  
  stopifnot(all_characters(pkgs), all_characters(lib), is_logical(save_file),
            is_logical(sort), is_logical(quietly), is_logical(verbose))
  
  pkgs <- unique(pkgs)
  pkgs_input <- pkgs
  # Remove the last forward slash and everything before it, because the package
  # name is the part after the last forward slash in GitHub repository names.
  pkgs <- sub(pattern = ".*/", replacement = "", x = pkgs)
  
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
                          format(Sys.time(), format = "%Y_%m_%d_%H_%M_%S"),
                          "_", paste0("R", as.character(getRversion())), ".txt")
      dir_path <- file.path(".", "output")
      
      # Notes:
      # - The path is normalised such that the printed path also works to read
      #   the data back into R after the working directory has changed, for
      #   example because R is not opened through the InstallPkgs project.
      # - Using "/" instead of "\\" as winslash so the printed path can be
      #   directly used in dget().
      read_back_path <- normalizePath(file.path(dir_path, file_name),
                                      winslash = "/", mustWork = FALSE)
      
      message_read_back <- paste0("\nTo read the package names back into R",
                                  " use:\nnonfunctional_pkgs <- dget(\"",
                                  read_back_path, "\")")
      
      if(!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE)
      }
      if(file.exists(read_back_path)) {
        warning("Textfile with names of non-functional packages already exists,",
                " not saved again!")
        message(message_read_back)
      } else {
        dput(nonfunctional_pkgs, file = read_back_path)
        message("Saved textfile with names of non-functional packages.",
                message_read_back)
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
#   lib: character vector giving the paths where packages are installed.
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
check_status <- function(lib, checkBuilt = TRUE,
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
  
  valid_out <- BiocManager::valid(lib.loc = lib, checkBuilt = checkBuilt,
                                  type = type)
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
                          format(Sys.time(), format = "%Y_%m_%d_%H_%M_%S"),
                          "_", paste0("R", as.character(getRversion())), ".csv")
      dir_path <- file.path(".", "output")
      if(!dir.exists(dir_path)){
        dir.create(dir_path, recursive = TRUE)
      }
      read_back_path <- normalizePath(file.path(dir_path, file_name),
                                      winslash = "/", mustWork = FALSE)
      message_read_back <- paste0("\nTo read the information back into R use:",
                                  " invalid_pkgs <- read.csv(\n\"",
                                  read_back_path, "\")")
      if(file.exists(read_back_path)) {
        warning("File with invalid packages already exists, not saved again!")
        message(message_read_back)
      } else {
        write.csv(invalid_details, file = read_back_path, row.names = FALSE)
        message("Saved .csv-file with details of invalid packages.",
                message_read_back)
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
#   See also https://pak.r-lib.org/reference/pkg_deps_explain.html which gives
#     details about which function creates the dependency and
#     https://github.com/yihui/xfun/blob/main/R/revcheck.R
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
  stopifnot(is.null(PC_name) ||
              (length(PC_name) == 1L && is.character(PC_name)))
  PC_name <- gsub("[^[:alnum:]_]", "_", PC_name)
  
  installed_pkgs <- installed.packages()[, c("Package", "Version", "Built",
                                             "NeedsCompilation", "Priority")]
  rownames(installed_pkgs) <- NULL
  file_name <- paste0("pkgs_installed_", PC_name, "_",
                      format(Sys.time(), format = "%Y_%m_%d_%H_%M_%S"), "_",
                      paste0("R", as.character(getRversion())), ".csv")
  dir_path <- file.path(".", "output")
  if(!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  read_back_path <- normalizePath(file.path(dir_path, file_name),
                                  winslash = "/", mustWork = FALSE)
  message_read_back <- paste0("\nTo read the information back into R use:",
                              " installed_pkgs <- read.csv(\n\"",
                              read_back_path, "\")")
  if(file.exists(read_back_path)) {
    warning("File with details of installed packages already exists,",
            " not saved again!")
    message(message_read_back)
  } else {
    write.csv(installed_pkgs, file = read_back_path, row.names = FALSE)
    message("Saved .csv-file with details of installed packages.",
            message_read_back)
  }
  invisible(installed_pkgs)
}


message("Sourced script containing the functions to check and install R-packages.")
