#!/usr/bin/env python3

import sys,math,functools

def block_err(data,nblock):
	ndat = len(data) 
	blen = ndat/nblock

	SUM=0.0
	SUM2=0.0

	for b in range(nblock):
		mean = reduce(lambda x,y:x+y,data[b*blen:(b+1)*blen])/float(blen)
		SUM+=mean
		SUM2+=mean**2
	
	mean = SUM/float(nblock)
	#print abs(SUM2/float(nblock)-mean**2)
	if abs(SUM2/float(nblock)-mean**2) < 1.e-10:
		sdev,sdev_mean=0.,0.
	else:
		sdev = math.sqrt((SUM2/float(nblock)-mean**2)/float(nblock))
		sdev_mean = sdev/math.sqrt(float(nblock))
	return mean,sdev,sdev_mean

diffit_dat = open(sys.argv[1]).readlines()

f_qaxis = map(float,diffit_dat[1].split()[1:])
d_qaxis = map(float,diffit_dat[3].split()[1:])
nf = len(f_qaxis)
F = {}
D = {}
rate = []
for i in range(nf):
	F[i] = []
nd = len(d_qaxis)
for i in range(nf):
	D[i] = []

for d in diffit_dat[4:]:
	ds= d.split()
	rate.append(float(ds[2]))
	for i in range(nf):
		F[i].append(float(ds[3+i]))
	for i in range(nd):
		D[i].append(float(ds[3+nf+i]))

fq_out = open(sys.argv[2],"w")
dq_out = open(sys.argv[3],"w")
mu_rate,sd_rate,sdm_rate = block_err(rate,10)
dq_out.write("# rate = %12.6e (%12.6e)\n"%(mu_rate,sdm_rate))

for i in range(nf):
	q = f_qaxis[i]
	fq,sd,sdm = block_err(F[i],10)
	fq_out.write("%12.6f %12.6f %12.6f\n"%(q,fq,sdm))

for i in range(nd):
	q = d_qaxis[i]
	dq,sd,sdm = block_err(D[i],10)
	dq_out.write("%14.6e %14.6e %14.6e\n"%(q,dq,sdm))

