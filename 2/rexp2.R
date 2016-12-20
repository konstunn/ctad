#!/usr/bin/Rscript

# generate sample from two parameters exponential distribution
rexp2 = function(n,rate=1,loc=0) {
	rexp(n,rate=rate)+loc
}

