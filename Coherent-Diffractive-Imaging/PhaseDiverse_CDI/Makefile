#include make.inc
DO_IDL=TRUE

all: 
	make -C src
	@echo "done making the NADIA libraries"
	make -C tools
	@echo "done making the tools"
	make -C examples
	@echo "done making the examples"
ifeq ($(DO_IDL), TRUE)
	make -C interfaces/idl
	@echo "done making the IDL interface"
endif
	@echo ""
	@echo "\"make\" completed successfully."
	@echo ""
	@echo "Look at the code in the examples directory to see how to get started"
	@echo ""



#install:
#	mkdir -p $(PREFIX)/include
#	cp src/*.h $(PREFIX)/include
#	mkdir -p $(PREFIX)/lib
#	cp src/lib* $(PREFIX)/lib/
#	mkdir -p bin $(PREFIX)/bin/
#	cp tools/*.exe $(PREFIX)/bin/

clean: 
	make clean -C src
	make clean -C examples
	make clean -C tools
	make clean -C doc
	rm -f interfaces/*~
ifeq ($(DO_IDL), TRUE)
	make clean -C interfaces/idl
endif

clobber:
	make clean
	rm -f lib/libNADIA.so lib/libNADIA.a
ifeq ($(DO_IDL), TRUE)
	make clobber -C interfaces/idl
endif

#use with caution
bare:
	make clean
	make clobber
	rm -f examples/*.ppm tools/*.ppm
	rm -f \#*# examples/\#*# src/\#*# tools/\#*#
	rm -f src/Makefile example/Makefile tools/Makefile
	rm -f config.log config.status
ifeq ($(DO_IDL), TRUE)
	rm -f interfaces/idl/Makefile 
endif
