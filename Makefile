.PHONY: sasa install clean
PREFIX ?= /usr
BINDIR = ${DESTDIR}${PREFIX}/bin

sasa: 
	@mkdir -p bin include
	@gfortran src/main.f90 -lgmxfort -I/usr/include -J./include -ljsonfortran -fopenmp -Wall -o bin/sasa

install: sasa
	@install -Dm755 bin/* -t ${BINDIR}

clean:
	@rm -r bin include
