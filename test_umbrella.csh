#!/bin/csh -f

# ======================================================================
# first generate some synthetic data using brownian dynamics
#

# SOME PARAMETERS...
#======================================================================
#set x1 = -1.0 	# (L)
set D = 0.1	# (D)
#set dt = 0.1	# (L**2/D) 
#set nprint = 2
#
#set dt = 0.01	# (L**2/D) 
#set nprint = 20
#
set dt = 0.001	# (L**2/D) - standard
set nprint = 200
#
#set dt = 0.0001	# (L**2/D)
#set nprint = 2000
set eps = 3	# (kBT)
set k1 = 30	# (kBT/L)
#set nprint = 20
set nbin = 30
set dtsave = `python -c "print float($nprint)*$dt"`

### EQUILIBRIUM RUN WITH NO UMBRELLA
#======================================================================
# 
## for dt = 0.0001
#./brumbrella -o X/br_noum_dt${dt}_p${nprint}.dat -D $D -d $dt \
#	-e $eps -x 0.0 -k 0.0 -n 1000000000 -p ${nprint} -s 27041994 
## for dt = 0.001
#./brumbrella -o X/br_noum_dt${dt}_p${nprint}.dat -D $D -d $dt \
#-e $eps -x 0.0 -k 0.0 -n 100000000 -p ${nprint} -s 27041994 
# for dt = 0.01
#./brumbrella -o X/br_noum_dt${dt}_p${nprint}.dat -D $D -d $dt \
#-e $eps -x 0.0 -k 0.0 -n 10000000 -p ${nprint} -s 27041994 
# for dt = 0.1
#./brumbrella -o X/br_noum_dt${dt}_p${nprint}.dat -D $D -d $dt \
#-e $eps -x 0.0 -k 0.0 -n 1000000 -p ${nprint} -s 27041994 

#
#
##
# RUNS WITH UMBRELLA SAMPLING (DT=0.001)
#======================================================================
#foreach x1 ( -1.5 -1.25 -1.0 -0.75 -0.5 -0.25 0.0 0.25 0.5 0.75 1.0 1.25 1.5 )
##foreach x1 ( 0.0  )
#
#	echo "umbrella sampling with x1 = ${x1}"
#	set seed = `python -c "print int($x1 * 10000 + 100000)"`
#	./brumbrella -o X/brum_x${x1}.dat -D $D -d $dt \
#		-e $eps -x $x1 -k $k1 -n 100000000 -p ${nprint} -s $seed
#
#end

#
#
## ======================================================================
## BIN TRAJECTORY INTO TRANSITION MATRIX
##
# EQUILIBRIUM (NO UMBRELLA)
#foreach t (  2  4  6  8  10 )
## new data format
#./thist -c 2 -L -1.65 -H 1.65 -n $nbin -o tmatX/br_noum_tmat_dt${dt}_lag${t}.dat \
#	-t $t -s $dtsave -k 0 -q 0 \
#	X/br_noum_dt${dt}_p${nprint}.dat
#
### format for GH program
##./thist -c 2 -L -1.65 -H 1.65 -n $nbin -o tmatX/br_noum_tmatg_dt${dt}_lag${t}.dat \
##	-t $t -s $dtsave -k 0 -q 0 -g \
##	X/br_noum_dt${dt}_p${nprint}.dat
#end

#cat tmatX/br_noum_tmatg_dt${dt}_lag{2,4,6,8,10}.dat > \
#	tmatX/br_noum_tmatg_dt${dt}_all.dat
##
### UMBRELLA SAMPLING
#foreach x1 ( -1.5 -1.25 -1.0 -0.75 -0.5 -0.25 0.0 0.25 0.5 0.75 1.0 1.25 1.5 )
#foreach t (  2  4  6  8  10 )
##foreach x1 ( 0.0  )
##foreach t ( 6 )
### new data format
#./thist -c 2 -L -1.65 -H 1.65 -n $nbin -o tmatX/brum_x${x1}_tmat_dt${dt}_lag${t}.dat \
#	-t $t -s $dtsave -k ${k1} -q ${x1} \
#	X/brum_x${x1}.dat
#
#### format for GH program
###./thist -c 2 -L -1.65 -H 1.65 -n $nbin -o tmatX/br_noum_tmatg_dt${dt}_lag${t}.dat \
###	-t $t -s $dtsave -k 0 -q 0 -g \
###	X/br_noum_dt${dt}_p${nprint}.dat
#end
#end
##
##cat tmatX/br_noum_tmatg_dt${dt}_lag{2,4,6,8,10}.dat > \
##	tmatX/br_noum_tmatg_dt${dt}_all.dat
##
## ======================================================================
## do bayesian diffusion analysis
##
#
set FC = "gfortran -O3"
set stiff = `python -c "print ${D}*0.1"`
#
#foreach t (  2  4  6  8  10 )
#
#set deltat = `python -c "print ${dtsave}*${t}"`
#
## format for GH program
#
#${FC} -o mc_eq maximum_likelihood_MC_diffusion09.F \
#	-D_N_=${nbin} -D_TIME_=${deltat} -D_NTIM_=1 -D_NMC_=50000 \
#	-D_DR_=0.1 -D_DP_=0.1 -D_TEMP0_=1.0 -D_TEMP1_=1.0 \
#	-D_DELTA_=0.11 -DLINEAR -DPLOG -D_DIFF_=${D} -D_DIFFSTIFF=$stiff
#${FC} -o mc_pr maximum_likelihood_MC_diffusion09.F \
#	-D_N_=${nbin} -D_TIME_=${deltat} -D_NTIM_=1 -D_NMC_=30000 \
#	-D_DR_=0.1 -D_DP_=0.1 -D_TEMP0_=1.0 -D_TEMP1_=1.0 \
#	-D_DELTA_=0.11 -DLINEAR -DPLOG -D_DIFF_=${D} -D_DIFFSTIFF=$stiff
#
#./mc_eq tmatX/br_noum_tmatg_dt${dt}_lag${t}.dat
#
#cp maximum_likelihood_MC_diffusion.log  diffg/br_noum_diff_dt${dt}_lag${t}_eq.log
#
#cp maximum_likelihood_MC_diffusion.{save,in}
#./mc_pr tmatX/br_noum_tmatg_dt${dt}_lag${t}.dat r > junk.dat
#
#./credibility.py -1.65 1.65 maximum_likelihood_MC_diffusion.log \
#	> diffg/br_noum_diff_dt${dt}_lag${t}_cred.dat
#
#cp maximum_likelihood_MC_diffusion.log  diffg/br_noum_diff_dt${dt}_lag${t}.log
#
#end

## GH program
#set tmin = `python -c "print ${dtsave}*2"`
#
#${FC} -o mc_eq maximum_likelihood_MC_diffusion09.F \
#	-D_N_=${nbin} -D_TIME_=${tmin} -D_NTIM_=5 -D_NMC_=50000 \
#	-D_DR_=0.1 -D_DP_=0.1 -D_TEMP0_=1.0 -D_TEMP1_=1.0 \
#	-D_DELTA_=0.11 -DLINEAR -DPLOG -D_DIFF_=${D} -D_DIFFSTIFF=$stiff
#${FC} -o mc_pr maximum_likelihood_MC_diffusion09.F \
#	-D_N_=${nbin} -D_TIME_=${tmin} -D_NTIM_=5 -D_NMC_=30000 \
#	-D_DR_=0.1 -D_DP_=0.1 -D_TEMP0_=1.0 -D_TEMP1_=1.0 \
#	-D_DELTA_=0.11 -DLINEAR -DPLOG -D_DIFF_=${D} -D_DIFFSTIFF=$stiff
#
#./mc_eq tmatX/br_noum_tmatg_dt${dt}_all.dat
#
#cp maximum_likelihood_MC_diffusion.log  diffg/br_noum_diff_dt${dt}_all_eq.log
#
#cp maximum_likelihood_MC_diffusion.{save,in}
#./mc_pr tmatX/br_noum_tmatg_dt${dt}_all.dat r > junk.dat
#
#./credibility.py -1.65 1.65 maximum_likelihood_MC_diffusion.log \
#	> diffg/br_noum_diff_dt${dt}_all_cred.dat
#
#cp maximum_likelihood_MC_diffusion.log  diffg/br_noum_diff_dt${dt}_all.log

# new program

# no umbrella ... seems to work
#./diffit -o diffn/br_noum_diff_eq.log -D 0.25 -d 0.1 -f 0.1 -n 50000 \
#	-p 500 tmatX/br_noum_tmat_dt${dt}_lag6.dat -R rst/br_noum_diff_eq_rst.txt
#
#./diffit -o diffn/br_noum_diff.log -D 0.25 -d 0.1 -f 0.1 -n 50000 \
#	-p 500 tmatX/br_noum_tmat_dt${dt}_lag6.dat -r rst/br_noum_diff_eq_rst.txt

#foreach x1 ( 0.0  )
#foreach t ( 6 )
### new data format
#./diffit -o diffn/brum_x${x1}_diff_eq.log -D 0.25 -d 0.1 -f 0.1 -n 50000 \
#	-p 500 tmatX/brum_x${x1}_tmat_dt${dt}_lag${t}.dat -R rst/brum_x${x1}_diff_eq_rst.txt

#
#./diffit -o diffn/brum_x${x1}_diff.log -D 0.25 -d 0.1 -f 0.1 -n 50000 \
#	-p 500 tmatX/brum_x${x1}_tmat_dt${dt}_lag${t}.dat -r rst/brum_x${x1}_diff_eq_rst.txt
#
#end


foreach t ( 6 )
#./diffit -o diffn/brum_all_lag${t}_diff_anneal.log -D 0.25 -d 0.1 -f 0.1 -n 30000 \
#	-R rst/brum_all_lag${t}_diff_anneal_rst.txt -S 0.01 -A 10 -T 1 \
#	-p 100 tmatX/brum_x{-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0,1.25,1.5}_tmat_dt${dt}_lag${t}.dat 
#
#./diffit -o diffn/brum_all_lag${t}_diff_eq.log -D 0.25 -d 0.1 -f 0.1 -n 50000 \
#	-r rst/brum_all_lag${t}_diff_anneal_rst.txt \
#	-R rst/brum_all_lag${t}_diff_eq_rst.txt -S 0.01 \
#	-p 100 tmatX/brum_x{-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0,1.25,1.5}_tmat_dt${dt}_lag${t}.dat 
#
#./diffit -o diffn/brum_all_lag${t}_diff.log -D 0.25 -d 0.1 -f 0.1 -n 50000 \
#	-r rst/brum_all_lag${t}_diff_eq_rst.txt -S 0.01 \
#	-p 100 tmatX/brum_x{-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0,1.25,1.5}_tmat_dt${dt}_lag${t}.dat 
./diffit_cred.py diffn/brum_all_lag${t}_diff.log cred/brum_all_lag${t}_fq.dat cred/brum_all_lag${t}_dq.dat 
end

