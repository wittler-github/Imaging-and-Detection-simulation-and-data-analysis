; This file demonstrates the steps in performing planar CDI
; reconstruction using the library in IDL. 

; Load up the module containing the wrapper code
.Compile CXS_interface.pro

; Load some image files of the data and the support. 
; A 2D array of doubles is returned.
; Please replace the file-name with your 
; "cxs_software/example/image_files" directory if you are not
; running this example on osiris.
support = cxs_read_tiff(1024,1024,'../../examples/image_files/planar_support.tiff')
data = cxs_read_tiff(1024,1024,'../../examples/image_files/planar_data.tif')

; Set-up everything ready for planar reconstruction.
; You need to pass the image data and support.
; By default the HIO algorithm is used with a relaxation
; parameter of beta=0.9
cxs_init_planar, data, support

;Perform an iteration
a = cxs_iterate()

; Apply shrink wrap (default Gaussian width = 1.5 pixels,
; threshold = 0.1 time the maximum pixel value). And then 
; perform 50 more iteration. Do this 3 times.
FOR I=0,3 DO BEGIN cxs_apply_shrinkwrap & a = cxs_iterate(50) & ENDFOR

cxs_apply_shrinkwrap

a = cxs_iterate(49)

; Change to the error-reduction algorithm 
cxs_set_algorithm, 'ER'

;Perform an iteration
a = cxs_iterate()

FOR I=0,1 DO BEGIN cxs_apply_shrinkwrap & a = cxs_iterate(50) & ENDFOR

cxs_apply_shrinkwrap

a = cxs_iterate(49)

; Apply shrink-wrap again. Lets be less tight this time
;cxs_apply_shrinkwrap, 1, 0.05 

; Do one more iteration
;a = cxs_iterate()


; Lets get the result with the lowest error from the
; last 100 iterations:
;a = cxs_get_best_result()

; and use this to restart the reconstruction starting with this
; result:
;cxs_init_planar, data, support, a

; We have finished with the reconstruction now, so free-up the
; memory we allocated earlier (note this does not effect "a").
; You will not longer be able to call cxs_iterate or any of the
; cxs get and set methods.
cxs_clear_memory

; Now you can play with "a" however you like in IDL.

; e.g. get the phase and display it:
; phase = ATAN(a, /PHASE)
; window, XSIZE=512, YSIZE=512
; TVSCL, rebin(phase,512,512)

; or the magnitude:
; TVSCL, rebin(abs(a),512,512)

; or save the result to a file:
; write_cplx(a , 'result_of_my_planar_CDI.cplx')

