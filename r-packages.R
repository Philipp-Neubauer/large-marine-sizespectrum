old <- getOption("defaultPackages"); 
r <- getOption("repos")
r["CRAN"] <- "http://cran.stat.auckland.ac.nz"
options(defaultPackages = c(old, "MASS"), repos = r)

packages <- c(
  "mizer","Hmisc","RODBC"
)

for (p in packages) {
  if (!require(p, character.only=TRUE)) {
    install.packages(p)
  }
}