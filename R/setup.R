.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "neonMicrobe relies on a pre-generated directory structure. ",
    "If this is your first time using neonMicrobe, or you want ",
    "to create a new directory structure, first set your working ",
    "directory to the location where you would like to create ",
    "this structure, and then run makeDataDirectories(). Learn ",
    "more in the vignette 'Download NEON Data'.",
    "\n\n",
    "neonMicrobe generates default output directories based ",
    "on the current working directory. To hold this constant, ",
    "set a 'base directory' using setBaseDirectory().",
    "\n"
  )
}

# Create custom package environment for storing parameters
neonmicrobe_env <- list2env(getDadaOpt())

# Create custom package environment for storing batch-specific parameters
batch_env <- new.env(parent=neonmicrobe_env)
