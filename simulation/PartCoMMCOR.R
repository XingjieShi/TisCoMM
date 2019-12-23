
source("PartType1ErrCOR.R")


assig = function(n_args){
    # example:  
	# n_args = c(2, 3, 4)
	cargs <- vector("list", length(n_args))
	for(i in 1:length(n_args)) cargs[[i]] = 1:n_args[i]
	t(expand.grid(cargs))
}
n_args = c(3, 5, 3, 2, 2)
jobs = assig(n_args)

h_z <- 0.01

main <- function(number0){
	number <- as.numeric(number0)	
	a0 <- 1
	rhoX <- c(0.2, 0.5, 0.8)[jobs[1, number]]
	h_c <- c(0.025, 0.05, 0.1, 0.2, 0.4)[jobs[2, number]]
	s = c(0.1, 0.5, 0.9)[1]
	rhoW = c(0.2, 0.5, 0.8)[jobs[3, number]]
	Ti <- 3
	nz_ti = c(1, 2, 3)[jobs[4, number]]
 	batch <- (1:10)[jobs[5, number]]

 	rhoE = 0.5
	Rep = 500
	set.seed(20102014*batch)
	case = c(rhoX = rhoX, rhoW = rhoW, h_c = h_c, s = s, nz_ti = nz_ti, batch = batch)
	print(case)

	PartType1Err(h_z, h_c, s, rhoX, rhoE, rhoW, nz_ti, Ti, a0, batch, Rep)
 }

args <- commandArgs(TRUE)
main(args[1])
