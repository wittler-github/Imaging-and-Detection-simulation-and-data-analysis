export LD_RUN_PATH=:/home/hwittler/TW/CODE/cxs_software/lib

IDL_SRC=IDL_interface.c
IDL_O=IDL_interface.o
LIB=libIDLCOECXS.so

all: $(LIB)

$(LIB): $(IDL_O)
	gcc -shared -O3  -fPIC -o $(LIB) $(IDL_O) \
	-L/home/hwittler/TW/CODE/cxs_software/lib/ -lCOECXS \
	 -ltiff -lmfhdf -ldf -lz -ljpeg -lm -lfftw3 -lstdc++ 

$(IDL_O): $(IDL_SRC)
	g++ -O3  -fPIC -c \
	-I/home/hwittler/TW/CODE/cxs_software/include \
	$(IDL_SRC)

clean:                                                     
	rm -f *~ $(IDL_O)

clobber: clean
	rm -f $(LIB)
