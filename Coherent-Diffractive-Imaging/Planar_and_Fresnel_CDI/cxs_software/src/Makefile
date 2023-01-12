
SOURCE_FILES_CXX=Complex_2D.c++ Double_2D.c++ FourierT.c++ \
	PlanarCDI.c++ Config.c++ FresnelCDI_WF.c++ FresnelCDI.c++
SOURCE_FILES_C=io_hdf.c io_ppm.c io_tiff.c io_dbin.c \
	io_cplx.c
OBJECT_FILES=$(SOURCE_FILES_CXX:.c++=.o) $(SOURCE_FILES_C:.c=.o)
HEADER_FILES=$(SOURCE_FILES_CXX:.c++=.h) io.h 

LIB_A= /home/hwittler/TW/CODE/cxs_software/lib/libCOECXS.a
LIB_SO= /home/hwittler/TW/CODE/cxs_software/lib/libCOECXS.so

all: $(LIB_A) $(LIB_SO) install_headers 

$(LIB_A): $(OBJECT_FILES) *.h
	ar cr $(LIB_A)$  $(OBJECT_FILES)

$(LIB_SO): $(OBJECT_FILES) *.h
	g++ -shared  -o $(LIB_SO) $(OBJECT_FILES) -ltiff -lmfhdf -ldf -lz -ljpeg -lm -lfftw3 -lstdc++ 

install_headers: $(HEADER_FILES)
	cp $(HEADER_FILES) /home/hwittler/TW/CODE/cxs_software/include

%.o: %.c++
	g++ -O3  -fPIC -c $< -o $@


%.o: %.c
	g++ -O3  -fPIC -DPACKAGE_NAME=\"CXS_Software\" -DPACKAGE_TARNAME=\"cxs_software\" -DPACKAGE_VERSION=\"0.0\" -DPACKAGE_STRING=\"CXS_Software\ 0.0\" -DPACKAGE_BUGREPORT=\"Nadia\ Davidson\ \<nadiamd@unimelb.edu.au\>\" -DPACKAGE_URL=\"\" -DHAVE_LIBSTDC__=1 -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIBFFTW3=1 -DHAVE_LIBM=1 -DHAVE_HDF_MFHDF_H=1 -DHAVE_LIBJPEG=1 -DHAVE_LIBZ=1 -DHAVE_LIBDF=1 -DHAVE_LIBMFHDF=1 -DHAVE_LIBTIFF=1 -DHAVE_POW=1 -DHAVE_SQRT=1 -c $< -o $@

clean:                                                     
	rm -f *.o *~ 

clobber: clean
	rm -f *.exe 
