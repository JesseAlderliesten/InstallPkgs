# InstallPkgs

Repository containing utility functions to install R packages and obtain
information about installed packages.
Also contains information on installing R, RStudio, Rtools and GitHub.
The information is intended for Windows-users and the code was run with R
version 4.1.3 on Windows 10.

Overview of the files:

- InstallR.txt: text file with installation instructions for R, RStudio, and
  Rtools. The main script assumes R is installed. If Rtools is not installed or
  not set up, it will try to do so.
- InstallPkgsMain.R: main R script to be run to install R packages and obtain
  information about installed packages, such as the status of packages and
  package dependencies.
- InstallPkgsFuncs.R: R script containing the functions used by the main script,
  with their documentation.
- InstallPkgsLists.R: R script containing thematic collections of packages I
  use, find useful, or want to store for potential later use. These collections
  are used by the main script.
- InstructionsPkgs.txt: text file with information on installing and getting
  information about R packages. Includes sections about installing alternative
  package versions and task views, lists of CRAN mirrors and package
  repositories, and a selection of books and websites that discuss R packages.
- InstructionsGitGithub.text: text file that provides information on setting up
  Git and GitHub.
