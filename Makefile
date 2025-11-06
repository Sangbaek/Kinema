
CXX    = g++ -I../lib/include -I../XSisolator


SRC    = kinema.cc dkinema.cc calc.cc Lab2CM.cc relativistic.cc xs.cc Math.cc gkinema.cc gkinema2pi.cc\

OBJ    = calc.o Lab2CM.o relativistic.o xs.o Math.o \

INCLUDE= $(HOME)/include
FLAG   = -I$(INCLUDE) -Wall

BINARY = kinema dkinema gkinema gkinema2pi
BINDIR = $(HOME)/bin

all: kinema dkinema smite gkinema gkinema2pi

clean:
	rm -f *.o core 

kinema: $(OBJ) kinema.o
	$(CXX) $^ -o $@

gkinema: $(OBJ) gkinema.o
	$(CXX) $^ -o $@

gkinema2pi: $(OBJ) gkinema2pi.o
	$(CXX) $^ -o $@

dkinema: $(OBJ) dkinema.o
	$(CXX) $(FLAG) $^ -o $@

smite: $(OBJ) smite.o
	$(CXX) $(FLAG) $^ -o $@

install: $(BINARY)
	install -c $^ $(BINDIR)

