#!/bin/csh -f

# ======================================================================
# first generate some synthetic data using brownian dynamics
#


### EQUILIBRIUM RUN WITH NO UMBRELLA
#======================================================================
# This uses brumbrella to run a 1D Brownian dynamics simulation on 
# a bistable potential - just to generate some diffusive data to fit
# 
# Normally, one would fit more interesting data...
# 
# uniform 1D diffusion coefficient 
set D = 0.1	# (D)
set dt = 0.001	# (L**2/D) - standard
set nprint = 200
# barrier height of bistable potential:
set eps = 3	# (kBT)

set dtsave = `python -c "print float($nprint)*$dt"`

./brumbrella -o br_dt${dt}_p${nprint}.dat -D $D -d $dt \
-e $eps -x 0.0 -k 0.0 -n 100000000 -p ${nprint} -s 27041994 

## ======================================================================
## BIN TRAJECTORY INTO TRANSITION MATRIX
##
# number of bins for discretizing 1D coordinate 
set nbin = 30

# This is making separate transition matrices with different lag times t

foreach t (  2  4  6  8  10 )
./thist -c 2 -L -1.65 -H 1.65 -n $nbin -o br_tmat_dt${dt}_lag${t}.dat \
	-t $t -s $dtsave -k 0 -q 0 \
	br_dt${dt}_p${nprint}.dat

end
#
## ======================================================================
# do bayesian diffusion analysis
#
# This is the fun part... and also most time consuming calculation. 

# First optimizing without any smoothening prior (diffusion coefficients
# may end up being a little noisy)
#
set stiff = `python -c "print ${D}*0.1"`
#
#foreach t (  2  4  6  8  10 )
foreach t ( 2 )

set deltat = `python -c "print ${dtsave}*${t}"`

# "equilibration" to find near optimal parameters
# Check log file to make sure the likelihood function has reached
# stable value
./diffit -o br_diff_lag${t}_eq.log -D 0.25 -d 0.1 -f 0.1 -n 50000 \
	-p 500 br_tmat_dt${dt}_lag${t}.dat -R br_diff_lag${t}_eq_rst.txt
#
# "production" where we sample from the Bayesian posterior
# likelihood in log file should just be fluctuating about mean value
./diffit -o br_diff_lag${t}.log -D 0.25 -d 0.1 -f 0.1 -n 50000 \
	-p 500 br_tmat_dt${dt}_lag${t}.dat -r br_diff_lag${t}_eq_rst.txt

# analyze log file to get F(x) and D(x):
./diffit_cred.py br_diff_lag${t}.log br_lag${t}_fq.dat br_lag${t}_dq.dat 

end

# now, optimization with smoothing
foreach t ( 2 )

set deltat = `python -c "print ${dtsave}*${t}"`

# "equilibration" to find near optimal parameters
# Check log file to make sure the likelihood function has reached
# stable value
./diffit -o br_diff_lag${t}_S0.01_eq.log -Q 0.1 -D 0.25 -d 0.1 -f 0.1 -n 50000 \
	-p 500 br_tmat_dt${dt}_lag${t}.dat -R br_diff_lag${t}_S0.01_eq_rst.txt

# "production" where we sample from the Bayesian posterior
# likelihood in log file should just be fluctuating about mean value
./diffit -o br_diff_lag${t}_S0.01.log -Q 0.01 -D 0.25 -d 0.1 -f 0.1 -n 50000 \
	-p 500 br_tmat_dt${dt}_lag${t}.dat -r br_diff_lag${t}_S0.01_eq_rst.txt

# analyze log file to get F(x) and D(x):
./diffit_cred.py br_diff_lag${t}_S0.01.log br_lag${t}_S0.01_fq.dat br_lag${t}_S0.01_dq.dat 

end

