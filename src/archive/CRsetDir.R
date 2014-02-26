############################### PREAMBLE ################################
# PROGRAM: CRsetDir
# CODER: Cameron Roach
# DATE: 16/7/2013
# DESCRIPTION: Make sure the correct directory is chosen

switch(Sys.info()["sysname"],
       Darwin =  {setwd("/Volumes/Data/Users/chroach/Dropbox/Uni/Research/")},
       Windows = {dirBootcamp <- "C:/Users/chroach/Dropbox/Uni/Research"
                  dirWindows <- "C:/Users/Cameron/Dropbox/Uni/Research"
                  if (file.exists(dirBootcamp)) {
                    setwd(dirBootcamp)
                  } else if (file.exists(dirWindows)) {
                    setwd(dirWindows)
                  } else {
                    cat("Directories are all messed up.")
                  }}
)