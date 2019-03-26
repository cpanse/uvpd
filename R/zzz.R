#R

.onAttach <- function(lib, pkg){
	if(interactive()){
		version <- packageVersion('uvpd')
		packageStartupMessage("Package 'uvpd' version ", version)
	  invisible()
	}
}
