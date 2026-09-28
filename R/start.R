.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "DEA version ", utils::packageVersion("DEA"), "\n",
    "Type citation('DEA') for citing this package in publications."
  )
}