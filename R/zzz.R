#R

.onAttach <- function(lib, pkg){
	if(interactive()){
		version <- packageVersion('UVPD')
		packageStartupMessage("Package 'UVPD' version ", version)
	  invisible()
	}
}
