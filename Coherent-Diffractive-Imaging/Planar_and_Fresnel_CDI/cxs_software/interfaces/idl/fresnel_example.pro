; This file demonstrates the steps in performing Fresnel CDI
; reconstruction using the library in IDL. 

; Load up the module containing the wrapper code
.Compile CXS_interface.pro

; Load some image files of the white-field data and 
; zone-plate support. A 2D array of doubles is returned.
; Please replace the file-name with your 
; "cxs_software/example/image_files" directory if you are not
; running this example on osiris.
s = cxs_read_tiff(1024,1024,'/data/nadia/cxs_software_rel_0/cxs_software/examples/image_files/FCDI_wf_support.tiff')
d = cxs_read_dbin(1024,1024,'/data/nadia/cxs_software_rel_0/cxs_software/examples/image_files/FCDI_wf_data.dbin')

; Set-up everything ready for reconstruction of the 
; white-field phase. You need to pass the data, support
; and experimental parameters.
cxs_init_fresnel_wf, d, s, 4.892e-10, 16.353e-3, 0.9078777, 13.5e-6

; Perform 20 iterations.
white_field = cxs_iterate(20)


; Now load the files of the sample support and data with the sample in place
s = cxs_read_tiff(1024,1024,'/data/nadia/cxs_software_rel_0/cxs_software/examples/image_files/FCDI_support.tiff')
d = cxs_read_dbin(1024,1024,'/data/nadia/cxs_software_rel_0/cxs_software/examples/image_files/FCDI_data.dbin')

; and set-up everything ready for the sample reconstruction.
cxs_init_fresnel, d, s, white_field, 4.892e-10, 0.9078777, 2.16e-3, 13.5e-6, 0.95

; Perform 20 iterations
a = cxs_iterate(20)

; Now get the transmission function based on the result of the
; final iteration.
a = cxs_get_transmission_function()
cxs_clear_memory

; We are finished with the white-field reconstruction now, so
; let free some memory

; Now you can play with "a" however you like in IDL.

; e.g. get the phase and display it:
phase = ATAN(a, /PHASE)
window, XSIZE=512, YSIZE=512
TVSCL, rebin(phase,512,512)

; or the magnitude:
; TVSCL, rebin(abs(a),512,512)

; or save the result to a file:
; write_cplx(a , 'result_of_my_FCDI_CDI.cplx')



