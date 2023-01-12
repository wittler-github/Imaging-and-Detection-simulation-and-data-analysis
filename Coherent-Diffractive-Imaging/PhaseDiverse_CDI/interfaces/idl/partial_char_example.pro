.Compile NADIA_interface.pro

support = nadia_read_tiff(2048,2048,'../../examples/image_files/part_support.tiff')
data = nadia_read_dbin(2048,2048,'../../examples/image_files/part_data.dbin')

nadia_init_partial_char, data, support, 0.0000135, 0.0000135, 1400.0, 1.4

nadia_set_algorithm, 'ER'

FOR I=0,4 DO BEGIN a = nadia_iterate(50) & nadia_apply_shrinkwrap & ENDFOR


nadia_set_algorithm, 'ER'

; Do another 150 iterations with shrink-wrap applied every 50 iterations
ake cl
; Do one last iteration to get the final result
a = nadia_iterate()

; Apply shrink-wrap again. Lets be less tight this time
; nadia_apply_shrinkwrap, 1, 0.05 

; Do one more iteration
;a = nadia_iterate()


; Lets get the result with the lowest error from the
; last 100 iterations:
;a = nadia_get_best_result()

; and use this to restart the reconstruction starting with this
; result:
;nadia_init_planar, data, support, a

; We have finished with the reconstruction now, so free-up the
; memory we allocated earlier (note this does not effect "a").
; You will not longer be able to call nadia_iterate or any of the
; nadia get and set methods.
nadia_clear_memory

; Now you can play with "a" however you like in IDL.

; e.g. get the phase and display it:
; phase = ATAN(a, /PHASE)
; window, XSIZE=512, YSIZE=512
; TVSCL, rebin(phase,512,512)

; or the magnitude:
; TVSCL, rebin(abs(a),512,512)

; or save the result to a file:
; write_cplx(a , 'result_of_my_planar_CDI.cplx')

