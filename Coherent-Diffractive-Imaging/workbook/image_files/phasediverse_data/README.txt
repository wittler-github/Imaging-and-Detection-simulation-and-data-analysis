This directory contains Au/Cr X-ray data ready to reconstruct. Once you have compiled the code to do phase diverse CDI, go into the reconstruction directory and type (make sure the binary is in your path somewhere or you know how to execute it):

cd recon
phase_diverse_cdi ../dlist_20100617.txt 2500 1 1 5.0

This will run a reconstruction using 7 datasets which include both longitudinal and lateral phase diversity. Each has its own diffraction data, its own reconstructed whitefield, and its own params.txt. All the diffraction data and wf data is in the ./data/ directory. Diffraction data is in double binary format (1024x1024). And white-field data is in double cplx format in accordance with fftw_complex. There is no headers in these files.

The command above runs the reconstruction for 2500 major iterations, with 1 FCDI iteration at each step. It uses gamma=5.0 which gives a nice combination of the data. Generally I always use beta=1.0.
