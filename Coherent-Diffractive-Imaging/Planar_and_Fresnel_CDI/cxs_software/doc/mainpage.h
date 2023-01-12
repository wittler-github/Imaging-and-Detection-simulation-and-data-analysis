/** 

 @mainpage The CXS Software Project
 
 Instruction on how to install and begin using the software package
 are provided on the <a
 href=http://www.ph.unimelb.edu.au/~ndavidson/cxs/index.html>main
 website</a>.  The doxygen documentation provided here is aimed at
 more advanced users who are comfortable with c++ and this style of
 reference documentation.

 <hr>
 @section examples Examples

 All examples and there associated data files are provided in the
 examples directory. A brief description of each is listed below.

 <ul>
 <li> @copydoc PlanarCDI_example.c
 <li> @copydoc PlanarCDI_simulation_example.c
 <li> @copydoc FresnelCDI_WF_example.c
 <li> @copydoc FresnelCDI_example.c
 <li> \a PlanarCDI_example.sh - Some of the same reconstuctions as above will
 be performed using the CDI_reconstruction.exe command line program. A
 simple bash script and configuration file (planar_example.config)
 shows how this tool can be used. Other parts of the bash script can
 be uncommented to run Fresnel reconstruction or to run multiple times
 with a different starting seed.
 </ul>

 <hr>

 @section clt Command Line Tools
 <hr>
 @copydoc CDI_reconstruction.c
 <hr>
 @copydoc dbin2ppm.c
 <hr>
 @copydoc hdf2ppm.c
 <hr>
 @copydoc tiff2ppm.c
 <hr> 
 @copydetails cplx2dbin.c
 <hr> 
 @copydoc cplx2ppm.c
 <hr> 


*/
