.onLoad <- function(lib, pkg) {
  library.dynam("flexrsurv", pkg, lib)
}#end of .onLoad


.onUnload <- function(lib)
{
    library.dynam.unload("flexrsurv", lib)
}


.onAttach <- function(lib, pkg) {
  vv <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), c("Version", "Date"))
  packageStartupMessage('Successfully loaded package flexrsurv version ',
                        as.character(packageVersion("flexrsurv")),'. For help type ?flexrsurv')
}

