#!/usr/bin/env python

import sys,math,numarray

def read_dat(file):
	# first three fields: time, log_like, rate
	inp = map(lambda x: map(float,x.split()[3:]), open(file).readlines())
	data = numarray.array(inp, numarray.Float64)
	data.transpose()
	return data

def cumu(d):
	fndat = float(len(d))
	mean = (reduce(lambda x,y:x+y, d))/fndat
	var = (reduce(lambda x,y:x+y*y, d)-d[0]+d[0]**2)
	var = (var/fndat-mean**2)#/max((fndat-1.0,1.0))
	return mean, var

lo = float(sys.argv[1])
hi = float(sys.argv[2])
dat = read_dat(sys.argv[3])
nbin = len(dat)/2

peq = dat[:nbin]
diff = dat[nbin:]
dx = (hi-lo)/float(nbin)

for i in range(nbin):
	fm,fv = cumu(map(lambda x: -math.log(x/dx), peq[i]))
	dm,dv = cumu(diff[i])
	f_x = lo+(float(i)+0.5)*dx
	d_x = f_x+0.5*dx
	if fv >= 0:
		sqfv = math.sqrt(fv)
	else:
		sqfv = 0
	if dv >= 0:
		sqdv = math.sqrt(dv)
	else:
		sqdv = 0
	print f_x, d_x, fm, sqfv, dm, sqdv



