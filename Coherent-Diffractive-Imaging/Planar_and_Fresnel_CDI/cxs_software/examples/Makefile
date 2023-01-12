export LD_RUN_PATH=:/home/hwittler/TW/CODE/cxs_software/lib

EXAMPLE_SRC=PlanarCDI_example.c PlanarCDI_simulation_example.c	\
FresnelCDI_example.c FresnelCDI_WF_example.c \
FresnelCDI_simulation_example.c Fresnel_NAS.c WF_Fresnel_NAS.c 
EXAMPLE_EXEC=$(EXAMPLE_SRC:.c=.exe)

all: $(EXAMPLE_EXEC)

%.exe: %.c
	g++ -O3  -fPIC -I/home/hwittler/TW/CODE/cxs_software/include $< -o $@ \
	-L/home/hwittler/TW/CODE/cxs_software/lib -lCOECXS \
	`pkg-config --cflags --libs gsl` -ltiff -lmfhdf -ldf -lz -ljpeg -lm -lfftw3 -lstdc++ 


clean:                                                     
	rm -f *~

clobber: clean
	rm -f *.exe
