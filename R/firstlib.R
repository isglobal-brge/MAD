
############ First.lib ###############

.onLoad <- function(lib, pkg){
   library.dynam("mad", pkg, lib)
}


.onUnload <- function(libpath) {
    library.dynam.unload("mad", libpath)
}

############ End of .First.lib ###############


