PROGNAME=prog

IDIR =./src/include
ODIR=.
SRCDIR=.

CC=g++ -std=c++17
C=gcc

CFLAGST=-I$(IDIR) `pkg-config --cflags gsl` -pthread -ggdb -fexceptions -Wall -pg # for testing
CFLAGS=-I$(IDIR) `pkg-config --cflags gsl` -O3 -pthread # for stuff with RNAfold 2.7.0

LIBS=-lm `pkg-config --libs gsl` -fno-lto -Wl,-fno-lto -lRNA -fopenmp -lgsl -lgslcblas -lpthread -lstdc++ -fopenmp # for RNAlib 2.7.0

_DEPS =  
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = list.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	@mkdir -p ${ODIR}
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	@mkdir -p ${ODIR}
	$(C) -c -o $@ $< $(CFLAGS)

$(PROGNAME): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: gdb
gdb: debug
gdb: CFLAGS=$(CFLAGST)
debug: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGST) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  

.PHONY: run

run:
	./$(PROGNAME)  

