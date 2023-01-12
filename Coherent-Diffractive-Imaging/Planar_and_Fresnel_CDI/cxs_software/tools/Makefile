export LD_RUN_PATH=:/home/hwittler/TW/CODE/cxs_software/lib

SOURCE_TOOLS=hdf2ppm.c tiff2ppm.c CDI_reconstruction.c dbin2ppm.c cplx2ppm.c cplx2dbin.c dbin2tiff.c
EXEC_TOOLS=$(SOURCE_TOOLS:.c=.exe)

LIB_SO=../lib/libCOECXS.so
LIB_A=../lib/libCOECXS.a

all: $(LIB_A) $(LIB_SO) $(EXEC_TOOLS)

$(LIB_A): ../src/*.c++ ../src/*.h ../src/*.c
	make -C ../src

$(LIB_SO): ../src/*.c++ ../src/*.h ../src/*.c
	make -C ../src

%.exe: %.c ../src/*.c++ ../src/*.h ../src/*.c
	g++ -O3  -fPIC -I/home/hwittler/TW/CODE/cxs_software/include $< -o $@ \
	-L/home/hwittler/TW/CODE/cxs_software/lib -lCOECXS \
	 -ltiff -lmfhdf -ldf -lz -ljpeg -lm -lfftw3 -lstdc++ 
	cp $@ /home/hwittler/TW/CODE/cxs_software/bin/$@

clean:                                                     
	rm -f *.o *~ 

clobber: clean
	rm -f *.exe 
