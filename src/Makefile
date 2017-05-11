# list of executable files to produce
MSD = msd
ENR = enr

# .o files necessary to build the executables
OBJ_MSD = AATree.o EvolutionModel.o EvolutionModelProt.o EvolutionModelProtStored.o Utils.o Tree.o NLOpt.o Alignment.o ConvergenceExpectation.o ColumnLikelihood.o DensityUtils.o AlignmentLikelihood.o MolSigDetect.o Character.o
OBJ_ENR = StatAnnotation.o Annotation.o  Utils.o Enrichment.o

########### MODIFY ONLY TO CHANGE OPTIONS ############

# compiler and its options
CC = gcc
CFLAGS = -I/usr/local/include -I/usr/include -Wall -Wno-char-subscripts -D_POSIX_SOURCE -std=c99 -Wall -pedantic -march=native -O3

# linker and its options
LD = $(CC)
LDFLAGS =  -rdynamic -lm /usr/local/lib/libnlopt.a -lgsl -lgslcblas 
LDFLAGS_TH =  -rdynamic -lm /usr/local/lib/libnlopt.a -lgsl -lgslcblas -D_REENTRANT -L/usr/local/lib -L/usr/lib/ -lpthread  

############ LIST OF EXECUTABLE TARGETS (MODIFY ONLY TO ADD AN EXECUTABLE) ##############

all: Makefile.d $(MSD) $(ENR)


# build the executable
$(MSD): $(OBJ_MSD)
	$(LD) $^ -o $@ $(LDFLAGS_TH)


# build the executable
$(ENR): $(OBJ_ENR)
	$(LD) $^ -o $@ $(LDFLAGS)


############ DO NOT MODIFY ANYTHING BELOW THIS LINE ##############

# create .o from .c
.c.o:
	$(CC) $(CFLAGS) -c $<

# remove non essential files
clean:
	$(RM) *.o *~ *.log Makefile.d

# clean everything but sources
distclean: clean
	$(RM) $(EXE)

# dependencies
Makefile.d:
	$(CC) -MM $(CFLAGS) *.c > Makefile.d

# only real files can be non phony targets
.PHONY: all clean distclean debug release

-include Makefile.d
