
TARGETS=brumbrella thist diffit propagators
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

install: 
	cp $(TARGETS) $(HOME)/bin

clean:
	-rm $(TARGETS)

distclean: clean
	-rm *.o
