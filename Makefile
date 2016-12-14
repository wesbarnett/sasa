.PHONY: sasa install clean
PROGRAM = sasa
PREFIX ?= /usr/local
BINDIR = ${DESTDIR}${PREFIX}/bin
CFLAGS += -I/usr/include -Jinclude `pkg-config --cflags libgmxfort` -fopenmp -Wall
LDFLAGS += -ljsonfortran `pkg-config --libs libgmxfort`

${PROGRAM}: 
	@mkdir -p bin include
	@gfortran src/main.f90 -o bin/$@ ${CFLAGS} ${LDFLAGS}

install: sasa
	@install -Dm755 bin/* -t ${BINDIR}

clean:
	@rm -r bin include
