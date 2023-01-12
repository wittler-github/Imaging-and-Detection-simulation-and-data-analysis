export LD_RUN_PATH=:/home/17526019/Work/NCDI/PDCODE/NADIA_A/lib

SOURCE_TOOLS=hdf2ppm.c tiff2ppm.c CDI_reconstruction.c dbin2ppm.c \
             cplx2ppm.c cplx2dbin.c dbin2tiff.c PhaseDiverseFresnelRec.c
EXEC_TOOLS=$(SOURCE_TOOLS:.c=.exe)

LIB_SO=../lib/libNADIA.so
LIB_A=../lib/libNADIA.a

all: $(LIB_A) $(LIB_SO) $(EXEC_TOOLS)

$(LIB_A): ../src/*.c++ ../src/*.c
	make -C ../src

$(LIB_SO): ../src/*.c++ ../src/*.c
	make -C ../src

%.exe: %.c ../src/*.c++ ../src/*.c
	g++ -fPIC -O3 -ffast-math  -fPIC -I/home/17526019/Work/NCDI/PDCODE/NADIA_A/include $< -o $@ \
	-L/home/17526019/Work/NCDI/PDCODE/NADIA_A/lib -lNADIA \
	 -ltiff -lmfhdf -ldf -llapack -lz -ljpeg -lm -lfftw3f -lstdc++ 
	cp $@ /home/17526019/Work/NCDI/PDCODE/NADIA_A/bin/$@

clean:                                                     
	rm -f *.o *.exe  *~ 

