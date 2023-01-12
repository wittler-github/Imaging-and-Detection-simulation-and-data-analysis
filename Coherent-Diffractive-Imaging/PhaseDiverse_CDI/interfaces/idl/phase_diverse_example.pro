; Copyright 2011 Nadia Davidson 
; for The ARC Centre of Excellence in Coherent X-ray Science. 
;
; This program is distributed under the GNU General Public License. 
; We also ask that you cite this software in publications where you made 
; use of it for any part of the data analysis.

; This file demonstrates the steps in performing Phase Diverse
; FCDI reconstruction using the library in IDL. 
; Thanks for Corey for the data which this example replies on.

; Load up the module containing the wrapper code
.Compile NADIA_interface.pro

; Load some image files of the white-field data and 
; zone-plate support. A 2D array of doubles is returned.
; Please replace the file-name with your 
; "NADIA/example/image_files" directory if you are not
; running this example on osiris.

white_field = ["wf_A_1024.cplx", $
               "wf_B_1024.cplx", $
               "wf_C_1024.cplx", $
               "wf_D_1024.cplx", $
               "wf_E_1024.cplx", $
               "wf_F_1024.cplx", $
               "wf_G_1024.cplx"]

data = ["A.dbin", $
        "B.dbin", $
        "C.dbin", $
        "D.dbin", $
        "E.dbin", $
        "F.dbin", $
        "G.dbin" ]

wavelength = 4.892e-10 ; wavelength
fz = 16.353e-3  ; zone plate to focal distance
ps = 13.5e-6 ; pixel size

;zone plate to detector distances
fd = [0.909513, 0.909388, 0.909263, 0.909213, $
      0.909088, 0.909088, 0.909088] 

;sample to detector distances
fs = [18.513e-3, 18.388e-3, 18.263e-3, 18.213e-3, $
      18.088e-3, 18.088e-3, 18.088e-3]

;white-field normalisation
norm = [0.984729833,0.97700431,0.986270638, 0.967487825,$
        0.980945916, 0.97628279, 0.963066039 ]

;relative positions in pixels 
x_pos = [0,23,3,35,30,37,60]
y_pos = [-150,-145,-134,-123,-150,-93,210]


;make a support (use the same for all frames)
s = NADIA_GET_ROUND_SUPPORT(1024,1024,0.5)

;Make a new PhaseDiverseCDI object
; NADIA_INIT_FRESNEL
NADIA_INIT_PHASE_DIVERSE, beta=1, gamma=3


FOR I=0,6 DO BEGIN $
   d = NADIA_READ_DBIN(1024,1024,data[I]) & $
   w = NADIA_READ_CPLX(1024,1024,white_field[I]) & $

   NADIA_INIT_FRESNEL, d, s, w, wavelength, fd[I]-fz, fs[I]-fz, ps, 1 & $
   ;norm[I] & $
   
   NADIA_SET_CHARGE_FLIPPING, 1 & $
   NADIA_SET_TRANS_UNITY_CONSTRAINT, 1 & $
   
   NADIA_PHASE_DIVERSE_ADD_POSITION, x_pos[I], y_pos[I] & $

ENDFOR


a = NADIA_PHASE_DIVERSE_INIT_ESTIMATE()

NADIA_PHASE_DIVERSE_ADJUST_POSITIONS

a = NADIA_PHASE_DIVERSE_ITERATE(10) 

NADIA_PHASE_DIVERSE_ADJUST_POSITIONS, TYPE=1

b = NADIA_PHASE_DIVERSE_ITERATE(5) 


; Now get the transmission function based on the result of the
; final iteration.
; a = nadia_get_transmission_function()
; nadia_clear_memory

; We are finished with the white-field reconstruction now, so
; let free some memory

; Now you can play with "a" however you like in IDL.

; e.g. get the phase and display it:

phase = ATAN(b, /PHASE)
show, phase 

; or the magnitude:
; TVSCL, rebin(abs(a),512,512)

; or save the result to a file:
; nadia_write_cplx(a, 'result_of_my_FCDI_CDI.cplx')

NADIA_CLEAR_MEMORY

