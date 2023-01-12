
SOURCE_FILES_CXX=Complex_2D.c++ BaseCDI.c++ PlanarCDI.c++ \
		 Config.c++ FresnelCDI_WF.c++ FresnelCDI.c++ \
		 TransmissionConstraint.c++ PhaseDiverseCDI.c++ \
		 PartialCharCDI.c++ PartialCDI.c++ PolyCDI.c++

SOURCE_FILES_C=io_hdf.c io_ppm.c io_tiff.c io_dbin.c \
	       io_cplx.c utils.c io_spec.c

OBJECT_FILES=$(SOURCE_FILES_CXX:.c++=.o) $(SOURCE_FILES_C:.c=.o)
HEADER_FILES=$(SOURCE_FILES_CXX:.c++=.h) io.h utils.h Double_2D.h

LIB_A= /home/17526019/Work/NCDI/PDCODE/NADIA_A/lib/libNADIA.a
LIB_SO= /home/17526019/Work/NCDI/PDCODE/NADIA_A/lib/libNADIA.so

all: $(LIB_A) $(LIB_SO) 

$(LIB_A): $(OBJECT_FILES) 
	ar cr $(LIB_A)$  $(OBJECT_FILES)

$(LIB_SO): $(OBJECT_FILES)
	g++ -I/home/17526019/Work/NCDI/PDCODE/NADIA_A/include -shared  -o $(LIB_SO) $(OBJECT_FILES) -ltiff -lmfhdf -ldf -llapack -lz -ljpeg -lm -lfftw3f -lstdc++ 

#install_headers: $(HEADER_FILES)
#	cp $(HEADER_FILES) /home/17526019/Work/NCDI/PDCODE/NADIA_A/include

%.o: %.c++
	g++ -fPIC -O3 -ffast-math  -fPIC -I/home/17526019/Work/NCDI/PDCODE/NADIA_A/include -c $< -o $@


%.o: %.c
	g++ -fPIC -O3 -ffast-math  -fPIC -I/home/17526019/Work/NCDI/PDCODE/NADIA_A/include -c $< -o $@

clean:                                                     
	rm -f *.o *.exe *~ 

clobber: clean
	rm -f *.exe 
