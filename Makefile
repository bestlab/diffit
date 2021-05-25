
TARGETS=brumbrella thist diffit propagators 1d_diff diffit_cred.py
# gogo
#CFLAGS= 
#LDFLAGS=-lgsl -lgslcblas
# bigmac
CFLAGS=-I/opt/local/include
LDFLAGS=-lgsl -L/opt/local/lib -lgsl -lgslcblas
OPT=-O3

all: $(TARGETS)

diff_model.o: diff_model.cc diff_model.h
	$(CXX) $(OPT) -c diff_model.cc $(CFLAGS)

brumbrella: brumbrella.cc
	$(CXX) $(OPT) -o brumbrella brumbrella.cc $(CFLAGS) $(LDFLAGS)

diffit: diffit.cc diff_model.o
	$(CXX) $(OPT) -o diffit diffit.cc diff_model.o $(CFLAGS) $(LDFLAGS)

propagators: propagators.cc diff_model.o
	$(CXX) $(OPT) -o propagators propagators.cc diff_model.o $(CFLAGS) $(LDFLAGS)

thist: thist.cc
	$(CXX) $(OPT) -o thist thist.cc

1d_diff: 1d_diff.cc ap.o spline3.o
	$(CXX) $(OPT) -o 1d_diff 1d_diff.cc ap.o spline3.o $(CFLAGS) $(LDFLAGS)

install: all
	cp $(TARGETS) $(HOME)/bin

clean:
	-rm $(TARGETS)

distclean: clean
	-rm *.o
