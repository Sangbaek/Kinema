# location of the Geant and other libraries
CERN_PRO = /usr/local/cern/pro
CERN_LIB = $(CERN_PRO)/lib
GLIBS = -lmathlib -lkernlib

#Kinematics Directory
KINDIR = $(HOME)/kinema
KIN    = $(KINDIR)/kin
XS     = $(KINDIR)/xs
PHOTO  = $(KINDIR)/photo

G77    = g77
CXX    = g++ -I$(HOME)/lib/include -I$(KIN) -I$(XS)


OBJ    = $(PHOTO)/PhProd.o $(PHOTO)/FF.o $(PHOTO)/calcXS.o \
	$(KIN)/relativistic.o $(KIN)/calc.o $(KIN)/Lab2CM.o $(KIN)/Math.o \
	$(XS)/Phys.o

#Includes
INCLUDE= $(HOME)/lib/include 
FLAG   = -I$(INCLUDE)

BINARY = strong PhProd Fermi
BINDIR = $(HOME)/bin

all: strong PhProd Fermi

strong: $(PHOTO)/strong.f
	$(G77) $^ -L$(CERN_LIB) $(GLIBS) -o $@

PhProd:$(OBJ)
	$(CXX) $(FLAG) $^ -o $@

Fermi:$(PHOTO)/Fermi.cc
	$(CXX) $(FLAG) $^ -o $@

install: $(BINARY)
	install -c $^ $(BINDIR)

clean:
	rm -f $(BINARY) $(BINDIR)/$(BINARY) *.o core *~


