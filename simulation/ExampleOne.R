
source("evalType1Err.R")


assig = function(n_args){
    # example:  
	# n_args = c(2, 3, 4)
	cargs <- vector("list", length(n_args))
	for(i in 1:length(n_args)) cargs[[i]] = 1:n_args[i]
	t(expand.grid(cargs))
}
n_args = c(2, 5, 3, 3, 5)
jobs = assig(n_args)


main <- function(number0){
	number <- as.numeric(number0)	
	h_z <- c(0, 0.01)[jobs[1, number]]
	h_c <- c(0.025, 0.05, 0.1, 0.2, 0.4)[jobs[2, number]]
	s = c(0.1, 0.5, 0.9)[jobs[3, number]]
	rhoX = c(0.2, 0.5, 0.8)[jobs[4, number]]
 	batch <- (1:10)[jobs[5, number]]

 	rhoE = 0.5
	Ti = 10
	Rep = 10
	set.seed(20132014*batch)
	case = c(h_z = h_z, h_c = h_c, s = s, rhoX = rhoX, batch = batch)
	print(case)

	evalType1Err(h_z, h_c, s, rhoX, rhoE, Ti, batch, Rep)
 }

args <- commandArgs(TRUE)
main(args[1])
