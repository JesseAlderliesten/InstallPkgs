# InstallPkgs

Repository containing information on installing and configuring R, RStudio, and
Rtools; setting up and using Git and GitHub; obtaining information about
to-be-installed and already-installed packages. The repository also contains
utility functions to install R packages. The information is intended for
Windows-users and has been tested with R versions 4.1.2 and 4.2.3 on Windows 10.

The steps in InstallR.txt should be followed before using InstallPkgsMain.R,
because that R script assumes R is installed.

Overview of the files:

- InstallR.txt: text file with installation instructions on the installation and
  configuration of R, RStudio, and Rtools on Windows. It also provides an
  overview of functions to get information about the used R installation and
  links to documentation and help on R and RStudio.
- InstallPkgsMain.R: This script can be used to check if R-packages are
  installed, get information about installed R-packages, and install new
  R-packages.
- InstallPkgsFuncs.R: R script with the functions used by InstallPkgsMain.R to
  install R-packages and obtain information about installed packages. It also
  includes some functions to check the installation of R.
- InstallPkgsLists.R: R script containing thematic collections of packages I
  use or want to store for potential later use. These package collections are
  used by InstallPkgsMain.R.
- InstructionsPkgs.txt: text file containing information on installing and
  getting information about R packages that are not yet installed or are already
  installed. It also includes a list of package repositories and details methods
  to obtain the source code of R functions.
- InstructionsGitGithub.txt: text file that provides information on setting up
  and using Git and GitHub on Windows.
