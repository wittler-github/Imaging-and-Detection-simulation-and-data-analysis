export LD_RUN_PATH=:/home/17526019/Work/NCDI/PDCODE/NADIA_A/lib

IDL_SRC=IDL_interface.c
IDL_O=IDL_interface.o
LIB=libIDLNADIA.so

all: $(LIB)

$(LIB): $(IDL_O)
	gcc -shared -fPIC -O3 -ffast-math  -fPIC -o $(LIB) $(IDL_O) \
	-L/home/17526019/Work/NCDI/PDCODE/NADIA_A/lib/ -lNADIA \
	 -ltiff -lmfhdf -ldf -llapack -lz -ljpeg -lm -lfftw3f -lstdc++ 

$(IDL_O): $(IDL_SRC)
	g++ -fPIC -O3 -ffast-math  -fPIC -c \
	-I/home/17526019/Work/NCDI/PDCODE/NADIA_A/include \
	$(IDL_SRC)

clean:                                                     
	rm -f *~ $(IDL_O)

clobber: clean
	rm -f $(LIB)
