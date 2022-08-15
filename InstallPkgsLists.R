#### Introduction #####
# This script contains thematic collections of packages I use, find useful, or
#   want to store for potential later use.
# The set 'high_prio_pkgs' contains the high-priority packages, consisting of
#   the 'base' R-packages that are installed if R is installed and the
#   'recommended' packages which are usually installed with R as well
#   (https://cran.r-project.org/doc/FAQ/R-FAQ.html#R-Add_002dOn-Packages;
#   https://cran.r-project.org/doc/FAQ/R-FAQ.html#Add_002don-packages-from-CRAN).
#   To see which high-priority packages are currently installed, run
#   dput(unname(installed.packages(priority = "high")[, "Package"]))
# The sets 'used_pkgs_UU' and 'used_pkgs_UvA' contain packages I explicitly load
#   in my scripts.
# The remaining sets are thematic collections of packages


# Note:
#   the 'translations' package is not a recommended package, but will be
#   installed with R if that option is set during the installation.
high_prio_pkgs <- c("base", "boot", "class", "cluster", "codetools", "compiler",
                    "datasets", "foreign", "graphics", "grDevices", "grid",
                    "KernSmooth", "lattice", "MASS", "Matrix", "methods",
                    "mgcv", "nlme", "nnet", "parallel", "rpart", "spatial",
                    "splines", "stats", "stats4", "survival", "tcltk", "tools",
                    "utils")

used_pkgs_UU <- c(
  ## Package management
  "BiocManager", # Access the Bioconductor Project Package Repository
  "devtools",    # Tools to Make Developing R Packages Easier
  "remotes",     # R Package Installation from Remote Repositories, Including
                 #   'GitHub'

  ## Data handling
  "dplyr",       # A Grammar of Data Manipulation
  "purrr",       # Functional Programming Tools
  "readr",       # Read Rectangular Text Data
  "reshape",     # (RETIRED! Use tidyr) Flexibly Reshape Data
  "reshape2",    # (RETIRED! Use tidyr) Flexibly Reshape Data: A Reboot of the
                 #   Reshape Package
  "tidyr",       # Tidy Messy Data
  "tidyverse",   # Easily Install and Load the 'Tidyverse'
  
  ## Graphics
  "cowplot",     # Streamlined Plot Theme and Plot Annotations for 'ggplot2'
  "ggalt",       # Extra Coordinate Systems, 'Geoms', Statistical Transformations,
                 #   Scales and Fonts for 'ggplot2'
  "ggdist",      # Visualizations of distributions and uncertainty
  "ggplot2",     # Create Elegant Data Visualisations Using the Grammar of Graphics
  "ggpubr",      # 'ggplot2' Based Publication Ready Plots
  "ggrepel",     # Automatically Position Non-Overlapping Text Labels with 'ggplot2'
  "gridExtra",   # Miscellaneous Functions for "Grid" Graphics
  "plot3D",      # Plotting Multi-Dimensional Data
  "RColorBrewer",# ColorBrewer Palettes (I prefer 'viridis' for continuous scales)
  "viridis",     # Colorblind-Friendly Color Maps for R
  
  ## Math
  "bayesplot",   # Plotting for Bayesian Models
  "bayestestR",  # Understand and Describe Bayesian Models and Posterior
                 #   Distributions
  "dunn.test",   # Dunn's Test of Multiple Comparisons Using Rank Sums
  "fitdistrplus",# Help to Fit of a Parametric Distribution to Non-Censored or
                 #   Censored Data
  "insight",     # Easy Access to Model Information for Various Model Objects
  "lme4",        # Linear Mixed-Effects Models using 'Eigen' and S4
  "lmtest",      # Testing Linear Regression Models
  "psych",       # Procedures for Psychological, Psychometric, and Personality
                 #   Research
  "ranger",      # A Fast Implementation of Random Forests
  "TruncatedNormal", # Truncated Multivariate Normal and Student Distributions
  
  ## Microbiology, ecology
  "ape",         # Analyses of Phylogenetics and Evolution
  "Biostrings",  # Efficient manipulation of biological strings
  "dada2",       # Accurate, high-resolution sample inference from amplicon
                 #   sequencing data
  "microbiome",  # Microbiome Analytics
  "phyloseq",    # Handling and analysis of high-throughput microbiome census data
  "vegan",       # Community Ecology Package
  
  ## ODEs
  "deSolve",     # Solvers for Initial Value Problems of Differential Equations
                 #   ('ODE', 'DAE', 'DDE')
  "FME",         # A Flexible Modelling Environment for Inverse Modelling,
                 #   Sensitivity, Identifiability and Monte Carlo Analysis
  "rootSolve",   # Nonlinear Root Finding, Equilibrium and Steady-State Analysis
                 #   of Ordinary Differential Equations
  
  ## Others
  "knitr",       # A General-Purpose Package for Dynamic Report Generation in R
  "microbenchmark", # Accurate Timing Functions
  "rmarkdown"    # Dynamic Documents for R
)

used_pkgs_UvA <- c(
  "BiocManager", # Access the Bioconductor Project Package Repository
  "deSolve",     # Solvers for Initial Value Problems of Differential Equations
                 #   ('ODE', 'DAE', 'DDE')
  "devtools",    # Tools to Make Developing R Packages Easier
  "FME",         # A Flexible Modelling Environment for Inverse Modelling,
                 #   Sensitivity, Identifiability and Monte Carlo Analysis
  "knitr",       # A General-Purpose Package for Dynamic Report Generation in R
  "matrixStats", # Functions that Apply to Rows and Columns of Matrices (and to
                 #   Vectors)
  "minpack.lm",  # R Interface to the Levenberg-Marquardt Nonlinear
                 #   Least-Squares Algorithm Found in MINPACK, Plus Support for
                 #   Bounds
  "remotes",     # R Package Installation from Remote Repositories, Including
                 #   'GitHub'
  "rmarkdown",   # Dynamic Documents for R
  "rootSolve",   # Nonlinear Root Finding, Equilibrium and Steady-State Analysis
                 #   of Ordinary Differential Equations
  "tidyr",       # Tidy Messy Data
  "TruncatedNormal" # Truncated Multivariate Normal and Student Distributions
)

coding_data_handling <- c(
  "abind",       # Combine Multidimensional Arrays
  "Biostrings",  # Efficient manipulation of biological strings
  "data.table",  # Extension of 'data.frame'
  "dplyr",       # A Grammar of Data Manipulation
  "gdata",       # Various R Programming Tools for Data Manipulation
  "purrr",       # Functional Programming Tools
  "readr",       # Read Rectangular Text Data
  "rlang",       # Functions for Base Types and Core R and 'Tidyverse' Features
  "runner",      # Running Operations for Vectors
  "sloop",       # Helpers for 'OOP' in R
  "tidyr",       # Tidy Messy Data
  "tidytable",   # Tidy Interface to 'data.table'
  "tidyverse"    # Easily Install and Load the 'Tidyverse'
)

coding_documentation <- c(
  "bookdown",    # Authoring Books and Technical Documents with R Markdown
  "knitr",       # A General-Purpose Package for Dynamic Report Generation in R
  "rmarkdown",   # Dynamic Documents for R
  "roxygen2",    # In-Line Documentation for R
  "shiny"        # Web Application Framework for R
)

coding_package_management <- c(
  "BiocManager", # Access the Bioconductor Project Package Repository
  "credentials", # Tools for Managing SSH and Git Credentials
  "ctv",         # CRAN Task Views
  "devtools",    # Tools to Make Developing R Packages Easier
  "needs",       # Attaches and Installs Packages
  "remotes",     # R Package Installation from Remote Repositories, Including
                 #   'GitHub'
  "sessioninfo"  # R Session Information
)

coding_performance <- c(
  ## Data inspection
  "lobstr",      # Visualize R Data Structures with Trees
  
  ## Defensive programming
  "assertthat",  # Easy Pre and Post Assertions
  "testthat",    # Unit Testing for R
  
  ## Code timing
  "bench",       # High Precision Timing of R Expressions
  "microbenchmark", # Accurate Timing Functions
  "proftools",   # Profile Output Processing Tools for R
  "profvis",     # Interactive Visualizations for Profiling R Code
  
  ## Speeding up
  "furrr",       # Apply Mapping Functions in Parallel using Futures
  "future",      # Unified Parallel and Distributed Processing in R for Everyone
  "future.apply",# Apply Function to Elements in Parallel using Futures
  "memoise",     # 'Memoisation' of Functions (caching)
  "parallel",    # Support for parallel computation, including by forking and by
                 #   sockets, and random-number generation (note: part of base R)
  "Rcpp",        # Seamless R and C++ Integration
  "RcppParallel" # Parallel Programming Tools for 'Rcpp'
)

graphics_palettes <- c(
  "khroma",      # Colour Schemes for Scientific Data Visualization
  "paletteer",   # Comprehensive Collection of Color Palettes
  "RColorBrewer",# ColorBrewer Palettes
  "viridis"      # Colorblind-Friendly Color Maps for R
)

graphics_combine_plots <- c(
  "cowplot",     # Streamlined Plot Theme and Plot Annotations for 'ggplot2'
  "gclus",       # Clustering Graphics
  "patchwork"    # The Composer of Plots
)

graphics_pkgs <- c(   
  "aRtsy",       # Generative Art with 'ggplot2'
  "beanplot",    # Visualization via Beanplots (like Boxplot/Stripchart/Violin Plot)
  "corrplot",    # Visualization of a Correlation Matrix
  "GGally",      # Extension to 'ggplot2'
  "ggalt",       # Extra Coordinate Systems, 'Geoms', Statistical Transformations,
                 #   Scales and Fonts for 'ggplot2'
  "gghighlight", # Highlight Lines and Points in 'ggplot2'
  "ggmuller",    # Creates Muller plots for visualizing evolutionary dynamics
  "ggplot2",     # Create Elegant Data Visualisations Using the Grammar of Graphics
  "ggpubr",      # 'ggplot2' Based Publication Ready Plots
  "ggrepel",     # Automatically Position Non-Overlapping Text Labels with 'ggplot2'
  "ggthemes",    # Extra Themes, Scales and Geoms for 'ggplot2'
  "labeling",    # Axis Labeling
  "pheatmap",    # Pretty Heatmaps
  "plotrix",     # Various Plotting Functions
  "png",         # Read and write PNG images
  "rgl",         # 3D Visualization Using OpenGL
  "scales",      # Scale Functions for Visualization
  "vcd",         # Visualizing Categorical Data
  "vcdExtra",    # 'vcd' Extensions and Additions
  "xtable"       # Export Tables to LaTeX or HTML
)

math_distr_random <- c(
  "distributional", # Vectorised Probability Distributions
  "fitdistrplus",# Help to Fit of a Parametric Distribution to Non-Censored or
                 #   Censored Data
  "ggdist",      # Visualizations of distributions and uncertainty
  "mvnTest",     # Goodness of Fit Tests for Multivariate Normality
  "randtoolbox", # Toolbox for Pseudo and Quasi Random Number Generation and
                 #   Random Generator Tests
  "RDieHarder",  # R Interface to the 'DieHarder' RNG Test Suite
  "rlecuyer",    # R Interface to RNG with Multiple Streams
  "TruncatedNormal",# Truncated Multivariate Normal and Student Distributions
  "truncdist"    # Truncated Random Variables
)

math_multidimensional_networks <- c(
  "ade4",        # Analysis of Ecological Data: Exploratory and Euclidean
                 #   Methods in Environmental Sciences
  "adegraphics", # An S4 Lattice-Based Package for the Representation of
                 #   Multivariate Data
  "BioNet",      # Routines for the functional analysis of biological networks
  "corrr",       # Correlations in R
  "ggmulti",     # High Dimensional Data Visualization
  "ggnetwork",   # Geometries to Plot Networks with 'ggplot2'
  "ggtern",      # An Extension to 'ggplot2', for the Creation of Ternary Diagrams
  "igraph",      # Network Analysis and Visualization
  "labdsv",      # Ordination and Multivariate Analysis for Ecology
  "mvabund",     # Statistical Methods for Analysing Multivariate Abundance Data
  "MVar",        # Multivariate Analysis
  "network",     # Classes for Relational Data
  "plot3D",      # Plotting Multi-Dimensional Data
  "plot3Drgl",   # Plotting Multi-Dimensional Data - Using 'rgl'
  "PMA",         # Penalized Multivariate Analysis
  "rrcov"        # Scalable Robust Estimators with High Breakdown Point
                 #   (principal component analysis, among others)
)

math_num_optim <- c(
  "caracas",     # Computer Algebra (symbolically solve equations, find integrals)
  "DEoptim",     # Global Optimization by Differential Evolution
  "DEoptimR",    # Differential Evolution Optimization in Pure R
  "dfoptim",     # Derivative-Free Optimization
  "ecr",         # Evolutionary Computation in R
  "expm",        # Matrix Exponential, Log, 'etc'
  "GA",          # Genetic Algorithms
  "GenSA",       # Generalized Simulated Annealing
  "irace",       # Iterated Racing for Automatic Algorithm Configuration
  "MaOEA",       # Many Objective Evolutionary Algorithm
  "mcga",        # Machine Coded Genetic Algorithms for Real-Valued Optimization
                 #   Problems
  "mco",         # Multiple Criteria Optimization Algorithms and Related Functions
  "minpack.lm",  # R Interface to the Levenberg-Marquardt Nonlinear Least-Squares
                 #   Algorithm Found in MINPACK, Plus Support for Bounds
  "nleqslv",     # Solve Systems of Nonlinear Equations
  "NMOF",        # Numerical Methods and Optimization in Finance
  "numDeriv",    # Accurate Numerical Derivatives
  "pracma",      # Practical Numerical Math Functions
  "pso",         # Particle Swarm Optimization
  "wnl"          # Minimization Tool for Pharmacokinetic-Pharmacodynamic Data
                 #   Analysis
)

math_stats_sets <- c(
  "bbmle",       # Tools for General Maximum Likelihood Estimation
  "combinat",    # combinatorics utilities
  "DescTools",   # Tools for Descriptive Statistics
  "e1071",       # Misc Functions of the Department of Statistics, Probability
                 #   Theory Group
  "EnvStats",    # Package for Environmental Statistics, Including US EPA Guidance  
  "epiR",        # Tools for the Analysis of Epidemiological Data
  "HH",          # Statistical Analysis and Data Display
  "Hmisc",       # Functions useful for data analysis, high-level graphics,
                 #   utility operations, computing sample size and power,
                 #   simulation, importing and annotating datasets..."
  "matrixStats", # Functions that Apply to Rows and Columns of Matrices (and to
                 #   Vectors)
  "misty",       # Miscellaneous Functions 'T. Yanagida' (Statistics in psychology)
  "parameters",  # Processing of Model Parameters
  "psych",       # Procedures for Psychological, Psychometric, and Personality
                 #   Research
  "robustbase",  # Basic Robust Statistics
  "see",         # Model Visualisation Toolbox for 'easystats' and 'ggplot2'
  "SHT",         # Statistical Hypothesis Testing Toolbox
  "statpsych"    # Statistical Methods for Psychologists
)

math_stats_MCMC_bayes <- c(
  "abc",         # Tools for Approximate Bayesian Computation (ABC)
  "BayesFactor", # Computation of Bayes Factors for Common Designs
  "bayesplot",   # Plotting for Bayesian Models
  "bayestestR",  # Understand and Describe Bayesian Models and Posterior
                 #   Distributions
  "coda",        # Output Analysis and Diagnostics for MCMC
  "ggmcmc",      # Tools for Analyzing MCMC Simulations from Bayesian Inference
  "loo",         # Efficient Leave-One-Out Cross-Validation and WAIC for
                 #   Bayesian Models
  "MCMCvis",     # Tools to Visualize, Manipulate, and Summarize MCMC Output
  "posterior",   # Tools for Working with Posterior Distributions
  "SimDesign",   # Structure for Organizing Monte Carlo Simulation Designs
  "simhelpers",  # Helper Functions for Simulation Studies (a.o. parallel running
                 #   over each row in tibble)
  "simsalapar",  # Tools for Simulation Studies in Parallel
  "simTool",     # Conduct Simulation Studies with a Minimal Amount of Source Code
  "tidybayes",   # Tidy Data and 'Geoms' for Bayesian Models
  "vizdraws"     # Visualize Draws from the Prior and Posterior Distribution
)

math_stats_regression <- c(
  "brglm2",      # Bias Reduction in Generalized Linear Models
  "broom",       # Convert Statistical Objects into Tidy Tibbles
  "car",         # Companion to Applied Regression
  "caret",       # Classification and Regression Training
  "dunn.test",   # Dunn's Test of Multiple Comparisons Using Rank Sums
  "gglm",        # Grammar of Graphics for Linear Model Diagnostic Plots
  "lme4",        # Linear Mixed-Effects Models using 'Eigen' and S4
  "lmerTest",    # Tests in Linear Mixed Effects Models
  "lmtest",      # Testing Linear Regression Models
  "merTools",    # Tools for Analyzing Mixed Effect Regression Models
  "nlmixr2",     # Nonlinear Mixed Effects Models in Population PK/PD
  "nlsMicrobio", # Nonlinear Regression in Predictive Microbiology
  "nlstools",    # Tools for Nonlinear Regression Analysis
  "olsrr",       # Tools for Building OLS Regression Models (ordinary least squares) 
  "ordinal",     # Regression Models for Ordinal Data
  "performance", # Assessment of Regression Models Performance
  "ranger",      # A Fast Implementation of Random Forests
  "standardize", # Tools for Standardizing Variables for Regression in R
  "visreg"       # Visualization of Regression Models
)

math_stats <- c(
  "binom",       # Binomial Confidence Intervals For Several Parameterizations
  "cAIC4",       # Conditional Akaike Information Criterion for 'lme4' and 'nlme'
  "coin",        # Conditional Inference Procedures in a Permutation Test Framework
  "cocorresp",   # Co-Correspondence Analysis Methods
  "confintr",    # Confidence Intervals
  "emmeans",     # Estimated Marginal Means, aka Least-Squares Means
  "fdrtool",     # Estimation of (Local) False Discovery Rates and Higher Criticism
  "flashlight",  # Shed Light on Black Box Machine Learning Models
  "fpc",         # Flexible Procedures for Clustering
  "glmnet",      # Lasso and Elastic-Net Regularized Generalized Linear Models
  "kernlab",     # Kernel-Based Machine Learning Lab
  "limma",       # Linear Models for Microarray Data
  "marginaleffects", # Marginal Effects, Marginal Means, Predictions, and Contrasts
  "MKinfer",     # Inferential Statistics
  "nlts",        # Nonlinear Time Series Analysis
  "pbkrtest",    # Parametric Bootstrap, Kenward-Roger and Satterthwaite Based
                 #   Methods for Test in Mixed Models
  "PMCMRplus",   # Calculate Pairwise Multiple Comparisons of Mean Rank Sums
                 #   Extended
  "pvclust",     # Hierarchical Clustering with P-Values via Multiscale
                 #   Bootstrap Resampling
  "pwr",         # Basic Functions for Power Analysis
  "randomForest",# Breiman and Cutler's Random Forests for Classification and
                 #   Regression
  "simglm",      # Simulate Models Based on the Generalized Linear Model
  "simpr",       # Flexible 'Tidyverse'-Friendly Simulations
  "simstudy",    # Simulation of Study Data
  "sva",         # Surrogate Variable Analysis
  "tidymodels",  # Easily Install and Load the 'Tidymodels' Packages
  "vsn",         # Variance stabilization and calibration for microarray data
  "workboots"    # Generate Bootstrap Prediction Intervals from a 'tidymodels'
                 # Workflow
)

microbio_eco <- c(
  # See also the sets 'graphics_multidimensional_networks' and
  #   'taxonomy_phylogeny'
  "KasperSkytte/ampvis2",# Tools for visualising amplicon data
  "ANCOMBC",     # Analysis of compositions of microbiomes with bias correction
  "biovizBase",  # Basic graphic utilities for visualization of genomic data
  "compositions",# Compositional Data Analysis
  "dada2",       # Accurate, high-resolution sample inference from amplicon
                 #   sequencing data
  "DECIPHER",    # Tools for curating, analyzing, and manipulating biological
                 #   sequences
  "DESeq2",      # Differential gene expression analysis based on the negative
                 #   binomial distribution [MSMB]
  "dimensio",    # Multivariate Data Analysis (Simple Principal Components Analysis
                 #   (PCA) and Correspondence Analysis (CA))
  "DirichletMultinomial",  # Dirichlet-Multinomial Mixture Model Machine Learning
                 #   for Microbiome Data
  "ecolMod",     # A practical guide to ecological modelling - using R as a
                 #   simulation platform (K. Soetaert)
  "genefilter",  # Methods for filtering genes from high-throughput experiments
  "HMP",         # Hypothesis Testing and Power Calculations for Comparing
                 #   Metagenomic Samples from HMP
  "microbiome",  # Microbiome Analytics
  "nlsMicrobio", # Nonlinear Regression in Predictive Microbiology
  "phyloseq",    # Handling and analysis of high-throughput microbiome census data
  "phyloseqGraphTest",  # Graph-Based Permutation Tests for Microbiome Data
  "primer",      # Functions and Data for the Book, a Primer of Ecology with R
  "CenterForStatistics-UGent/RCM",  # A model-based visualization method for
                 #   microbiome data
  "seqTools",    # Analysis of nucleotide, sequence and quality content on fastq
                 #   files
  "zdk123/SpiecEasi",  # Sparse InversE Covariance estimation for Ecological
                 #   Association and Statistical Inference
  "vegan",       # Community Ecology Package
  "vegan3d"      # Static and Dynamic 3D Plots for the 'vegan' Package
)

ODE_epi_growth <- c(
  "biogrowth",   # Modelling of Population Growth
  "EpiDynamics", # Dynamic Models in Epidemiology (replicating scripts for
                 #   Keeling and Rohani 'Modeling Infectious Diseases in Humans
                 #   and Animals')
  "EpiEstim",    # Estimate Time Varying Reproduction Numbers from Epidemic Curves
  "epimdr",      # Functions and Data for "Epidemics: Models and Data in R"
  "EpiModel",    # Mathematical Modeling of Infectious Disease Dynamics
  "gauseR",      # Lotka-Volterra Models for Gause's 'Struggle for Existence'
  "genSEIR",     # Predict Epidemic Curves with Generalized SEIR Modeling (not
                 #   based on deSolve, but on own implementation of RK4)
  "growthrates", # Estimate Growth Rates from Experimental Data
  "idmodelr",    # Infectious Disease Model Library and Utilities
  "microPop",    # Process-Based Modelling of Microbial Populations
  "outbreaker2", # Bayesian Reconstruction of Disease Outbreaks by Combining
                 #   Epidemiologic and Genomic Data
  "streambugs"   # Parametric Ordinary Differential Equations Model of Growth,
                 #   Death, and Respiration of Macroinvertebrate and Algae Taxa
)

ODE_handling <- c(
  # See also the set 'math_num_optim'
  "cOde",        # Automated C Code Generation for 'deSolve', 'bvpSolve'
  "Deriv",       # Symbolic Differentiation
  "deSolve",     # Solvers for Initial Value Problems of Differential Equations
                 #   ('ODE', 'DAE', 'DDE')
  "deTestSet",   # Test Set for Differential Equations
  "diffeqr",     # Solving Differential Equations (ODEs, SDEs, DDEs, DAEs)
                 #   (based on Julia instead of deSolve)
  "ecolMod",     # A practical guide to ecological modelling - using R as a
                 #   simulation platform (K. Soetaert)
  "FME",         # A Flexible Modelling Environment for Inverse Modelling,
                 #   Sensitivity, Identifiability and Monte Carlo Analysis
  "GillespieSSA2", # Gillespie's Stochastic Simulation Algorithm for Impatient
                 #   People
  "mkin",        # Kinetic Evaluation of Chemical Degradation Data
  "mrgsim.parallel", # Simulate with 'mrgsolve' in Parallel
  "mrgsim.sa",   # Sensitivity Analysis with 'mrgsolve'
  "mrgsolve",    # Simulate from ODE-Based Models (not based on deSolve, C++ code)
  "odeintr",     # C++ ODE Solvers Compiled on-Demand
  "ODEsensitivity", # Sensitivity Analysis of Ordinary Differential Equations
  "odin",        # ODE Generation and Integration (uses deSolve, but automatically
                 #   translates to C)
  "phaseR",      # Phase Plane Analysis of One- And Two-Dimensional Autonomous
                 #   ODE Systems
  "PSPManalysis",# Analysis of Physiologically Structured Population Models
  "r2sundials",  # Wrapper for 'SUNDIALS' Solving ODE and Sensitivity Problem
                 #   (not based on deSolve)
  "rodeo",       # A Code Generator for ODE-Based Models
  "sde",         # Simulation and Inference for Stochastic Differential Equations
  "tpetzoldt/rodeoExt", # Extensions to the 'rodeo' Package (stoichiometry notation)
  "rootSolve",   # Nonlinear Root Finding, Equilibrium and Steady-State Analysis
                 #   of Ordinary Differential Equations
  "RxODE",       # Facilities for Simulating from ODE-Based Models (developed for
                 #   PKPD analyses; automatically translates to C)
  "SimInf",      # A Framework for Data-Driven Stochastic Disease Spread Simulations
  "simlandr",    # Simulation-Based Landscape Construction for Dynamical Systems
  "sundialr"     # An Interface to 'SUNDIALS' Ordinary Differential Equation
                 #   (ODE) Solvers (not using deSolve)
)

ODE_PBPK <- c(
  "ecotox",      # Analysis of Ecotoxicology
  "nlmixr2",     # Nonlinear Mixed Effects Models in Population PK/PD 
  "PKconverter", # The Parameter Converter of the Pharmacokinetic Models
  "ssdtools",    # Species Sensitivity Distributions
  "stanette",    # R Interface to Stan"
  "wnl"          # Minimization Tool for Pharmacokinetic-Pharmacodynamic Data
                 #   Analysis  
)

taxonomy_phylogeny <- c(
  # See also the sets 'graphics_multidimensional_networks' and 'microbio_eco'
  "ape",         # Analyses of Phylogenetics and Evolution
  "BAMMtools",   # Analysis and Visualization of Macroevolutionary Dynamics on
                 #   Phylogenetic Trees
  "beastier",    # Call 'BEAST2'. NOTE: BEAST2 needs to be installed, see
                 #   https://www.beast2.org/.
  "ggdendro",    # Create Dendrograms and Tree Diagrams Using 'ggplot2'
  "ggtree",      # An R package for visualization of tree and annotation data
  "ggtreeExtra", # An R Package To Add Geometric Layers On Circular Or Other
                 # Layout Tree Of "ggtree"
  "metacoder",   # Tools for Parsing, Manipulating, and Graphing Taxonomic
                 #   Abundance Data
  "phangorn",    # Phylogenetic Reconstruction and Analysis
  "donkeyshot/phybreak", # Outbreak reconstruction with sequence data
  "phylobase",   # Base Package for Phylogenetic Structures and Comparative Data
  "picante",     # Integrating Phylogenies and Ecology
  "rRDP",        # Interface to the RDP Classifier
  "taxize",      # Taxonomic Information from Around the Web
  "TDbook",      # Companion Package for the Book "Data Integration,
                 #   Manipulation and Visualization of Phylogenetic Trees"
  "tidytree",    # A Tidy Tool for Phylogenetic Tree Data Manipulation
  "TipDatingBeast", # Using Tip Dates with Phylogenetic Trees in BEAST
  "TiPS",        # Trajectories and Phylogenies Simulator (also ODE)
                 #   https://gitlab.in2p3.fr/ete/tips
  "TreeDist"     # Calculate and Map Distances Between Phylogenetic Trees
)

message("The script containing the lists of R-packages has been sourced.")
