

# ===========================================================
# 📌 FIRST-TIME REQUIREMENT:
# ===========================================================

#  Install missing packages (only needed once)

	# This section should ONLY be run the first time you execute the script on a new R installation.
	# Once the packages are installed, you don't need to run this part again.Required Libraries and Installation Order

	# List of required packages
	required_packages <- c("readr", "ggplot2", "dplyr", "tidyr", "tibble", "matrixStats", "optimx")

	# Check for missing packages and install them if necessary
	new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
	if(length(new_packages)) install.packages(new_packages)

# LOAD LIBRARIES: This should always be executed before running the main script
	# Even though installation is only needed once, libraries must be loaded in every session.

	lapply(required_packages, library, character.only = TRUE)

# Now the main script can run without issues