;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Copyright 2011 Nadia Davidson, 2012 T'Mir Julius
; for The ARC Centre of Excellence in Coherent X-ray Science. 
;
; This program is distributed under the GNU General Public License. 
; We also ask that you cite this software in publications where you made 
; use of it for any part of the data analysis.
;
; This code provides wrappers to the functions in the NADIA C++
; library. It calls the methods which are defined in the
; IDL_interface.c file (and compiled into libIDLNADIA.so), thus
; allowing CDI reconstruction to be performed in IDL. Some examples
; are provided in this directory, showing how you can use these
; methods.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Couldn't work out how to make global variables
; so this is my way around it.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function lib_name
return, 'libIDLNADIA.so'
end


function pixels
return, 512
end

;dimensions are reversed.
function ny
return, call_external(lib_name(),'IDL_get_array_x_size')
end

function nx
return, call_external(lib_name(),'IDL_get_array_y_size')
end


;pro show, image
;n =  size(image)
;if ( (n[2] MOD pixels()) eq 0.0 ) and ((n[1] MOD pixels()) eq 0.0 ) then begin
;   window, XSIZE=pixels(), YSIZE=pixels()
;   TVSCL, rebin(image,pixels(),pixels())
;endif
;end

pro show, image
  n =  size(image)
  new_x = pixels()
  new_y = n[2]*new_x/n[1]
  window, XSIZE=new_x, YSIZE=new_y
  TVSCL, CONGRID( image, new_x, new_y)
end


pro check_size, array
n = size(array)
if (nx() ne n[2]) or (ny() ne n[1]) then begin
    print, 'Please provide a 2D array with dimensions ', nx, ' x ', ny
    stop
endif
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Here are the real wrappers to the C++ code
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;+
; NAME:
;       NADIA_INIT_PLANAR
;
; PURPOSE:
;       Set-up a planar CDI reconstruction. This will
;       initialise the reconstruction with the data
;       and support. Some defaults will be set and
;       memory will be allocated ready for reconstruction.
;       It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or 
;       calling NADIA_ITERATE).
;
;       Calling this procedure will initialise the reconstruction 
;       algorithm to hybrid-input-output with a relaxation parameter of 0.9.
;
;
; CALLING SEQUENCE:
;
;	NADIA_INIT_PLANAR, Data, Support [,Starting_point]
;
;
; INPUTS:
;
;	Data: 
;             The detector illumination. It should be
;             in the form of a 2D array
;
;       Support: 
;             A 2D Array giving the sample support.
;             Values or 1 or greater are considered inside
;             the support. All others are considered to be
;             outside the support.
;
;       Starting_point: 
;             As an option you may supply an initial 
;             guess of the exit-surface-wave for the sample. 
;             This maybe useful, for example, if you wish to 
;             start from the end point of a previous run. The
;             format of this parameter much be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the starting point is initialised to be zero outside
;             the support and a random number inside the support, 
;             for both the magnitude and phase.
;
; EXAMPLE:
;        An example of loading two 2D arrays from file and using
;        them to initialise the planar reconstruction:
;
;        my_support = nadia_read_tiff(1024,1024,'planar_support.tiff')
;        my_data = nadia_read_tiff(1024,1024,'planar_data.tiff')
;        NADIA_INIT_PLANAR, my_data, my_support
;
;-
pro nadia_init_planar, data, support, complex_array
n = size(data)
nx = n[2]
ny = n[1]
IF N_Params() EQ 3 THEN $
  b = call_external(lib_name(),'IDL_planar_init',nx,ny,complex_array) $
ELSE $
  b = call_external(lib_name(),'IDL_planar_init',nx,ny)
nadia_set_support, support
nadia_set_intensity, data
IF N_Params() EQ 2 THEN $
  nadia_initialise_esw
end


;+
; NAME:
;       NADIA_INIT_FRESNEL_WF
;
; PURPOSE:
;       Set-up a Fresnel white-field CDI reconstruction. This will
;       initialise the reconstruction with the white-field intensity, 
;       zone-plate support and experimental parameters. Some defaults 
;       will be set and memory will be allocated ready for 
;       reconstructing the white-field (phase and magnitude) in the
;       detector plane. It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or calling NADIA_ITERATE).
;
; CALLING SEQUENCE:
;
;	NADIA_INIT_FRESNEL_WF, data, support, beam_wavelength,
;                            zone_focal_length, focal_detector_length,
;                            pixel_size [,starting_point]
;
; INPUTS:
;
;       You may use any length units for the experimental parameters
;       below, as long as all quantities are given in the same units.
;
;	data: 
;             The white-field illumination. It should be
;             in the form of a 2D array.
;
;       support: 
;             A 2D array of integers or doubles giving the zone-plate support.
;             Values or 1 or greater are considered inside
;             the support. All others are considered to be
;             outside the support.
;
;       beam_wavelength:
;             The beam wavelength.
;
;       zone_focal_length:
;             The distance between the zone plate and the focal point.
;
;       focal_detector_length:
;             The distance between the focal point and the detector.
;
;       pixel_size:
;             The side length of one detector pixel.
;
;       starting_point: 
;             As an option you may supply an initial guess of the 
;             white-field. This may be useful, for example, if you 
;             wish to start from the end point of a previous run. The
;             format of this parameter must be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the starting point is initialised to be zero outside
;             the support, a random number inside the support 
;             for the magnitude and zero for the phase.
;
; EXAMPLE:
;
;        An example of loading two 2D arrays from file and using
;        them to initialise the white-field reconstruction for FCDI:
;
;        my_support = nadia_read_tiff(1024,1024,'support.tiff')
;        my_data = nadia_read_tiff(1024,1024,'data.tiff')
;        nadia_init_fresnel_wf, my_data, my_support, 4.892e-10, 16.353e-3, 0.9078777,13.5e-6
;-
pro nadia_init_fresnel_wf, data, $
                         support, $
                         beam_wavelength, $
                         zone_focal_length, $
                         focal_detector_length, $
                         pixel_size, $
                         complex_array
n = size(data)
nx = n[2]
ny = n[1]
IF N_Params() EQ 7 THEN $
  b = call_external(lib_name(),'IDL_fresnel_wf_init',nx,ny, $
                    double(beam_wavelength), $
                    double(zone_focal_length), $
                    double(focal_detector_length), $
                    double(pixel_size), $
                    complex_array) $
ELSE $
  b = call_external(lib_name() ,'IDL_fresnel_wf_init',nx,ny, $
                    double(beam_wavelength), $
                    double(zone_focal_length), $
                    double(focal_detector_length), $
                    double(pixel_size))

nadia_set_support, support
nadia_set_intensity, data

IF N_Params() EQ 6 THEN $
  nadia_initialise_esw

end

;+
; NAME:
;       NADIA_INIT_FRESNEL
;
; PURPOSE:
;       Set-up a Fresnel CDI reconstruction. This will
;       initialise the reconstruction using a previously reconstructed
;       white-field, detector data, sample support and experimental 
;       parameters. Some defaults will be set and memory will be
;       allocated ready for reconstructing the sample
;       exit-surface-wave. It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or calling NADIA_ITERATE).
;
;       Calling this procedure will initialise the reconstruction algorithm
;       to the error-reduction with a relaxation parameter of 0.9.
;
; CALLING SEQUENCE:
;
;	NADIA_INIT_FRESNEL, data, support, white-field, beam_wavelength,
;	                  focal_detector_length, focal_sample_length, 
;                         pixel_size [, normalisation, starting_point ]
;
; INPUTS:
;
;       You may use any length units for the experimental parameters
;       below, as long as all quantities are given in the same units.
;
;	data: 
;             The detector data with the sample in place. It should be
;             in the form of a 2D array.
;
;       support: 
;             A 2D array of integers or doubles which give the sample support.
;             Values or 1 or greater are considered inside
;             the support. All others are considered to be
;             outside the support.
;
;       white-field:
;             A COMPLEX 2D array of the reconstructed white-field in
;             the detector plane. This can be recovered using 
;             NADIA_INIT_FRESNEL_WF followed by NADIA_ITERATE.
;
;       beam_wavelength:
;             The beam wavelength.
;
;       focal_detector_length:
;             The distance between the focal point and the detector.
;
;       focal_sample_length:
;             The distance between the focal point and the sample.
;
;       pixel_size:
;             The side length of one detector pixel.
;
;       normalisation: 
;             The factor to scale the white-field before
;             performing FCDI. If this parameter is excluded, the
;             ratio of the square-root of the intensity data and the
;             white-field magnitude is used as the normalisation.
;
;       starting_point: 
;             As an option you may supply an initial 
;             guess of the exit-surface-wave for the sample. 
;             This maybe useful, for example, if you wish to 
;             start from the end point of a previous run. The
;             format of this parameter much be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the initialisation described in Harry's review paper: 
;             page 29. (in particular e.q. 137) is used.
;
; EXAMPLE:
;
;        nadia_init_fresnel, my_data, my_supports, my_white-field, $
;                          4.892e-10, 0.9078777, 2.16e-3, $
;                          13.5e-6
;-
pro nadia_init_fresnel, data, support, $
                      white_field, $
                      beam_wavelength, $
                      focal_detector_length, $
                      focal_sample_length, $
                      pixel_size, $
                      normalisation,$
                      complex_array
n = size(data)
nx = n[2]
ny = n[1]
IF N_Params() EQ 7 THEN BEGIN
  mag2_wf = abs(white_field)^2
  normalisation = total(data) / total(mag2_wf)
  print, normalisation
ENDIF
IF N_Params() EQ 9 THEN $
  b = call_external(lib_name(),'IDL_fresnel_init',long(nx),long(ny), $
                    white_field, $
                    double(beam_wavelength), $
                    double(focal_detector_length), $
                    double(focal_sample_length), $
                    double(pixel_size), $
                    double(normalisation),$
                    complex_array) $
ELSE $
  b = call_external(lib_name() ,'IDL_fresnel_init',long(nx),long(ny), $
                    white_field, $
                    double(beam_wavelength), $
                    double(focal_detector_length), $
                    double(focal_sample_length), $
                    double(pixel_size), $
                    double(normalisation))

nadia_set_support, support
nadia_set_intensity, data

IF N_Params() LT 9 THEN $
  nadia_initialise_esw

end

;+
; NAME:
;       NADIA_INIT_PARTIAL
;
; PURPOSE:
;       Set-up a Partial CDI reconstruction. This will
;       initialise the reconstruction with sample support and experimental 
;       parameters. Some defaults will be set and memory will be
;       allocated ready for reconstructing the sample
;       exit-surface-wave. It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or calling NADIA_ITERATE).
;
;       Calling this procedure will initialise the reconstruction algorithm
;       to the error-reduction with a relaxation parameter of 0.9.
;
; CALLING SEQUENCE:
;
; NADIA_INIT_PARTIAL, data, support, lcx, lcy, pxsize, pysize, beam_energy
;                   zsd [, starting_point ]
;
; INPUTS:
;
;       Please use metres for the lengths, and eV for the energy.
;
;       data: 
;             The detector data with the sample in place. It should be
;             in the form of a 2D array.
;
;       support: 
;             A 2D array of integers or doubles which give the sample support.
;             Values or 1 or greater are considered inside
;             the support. All others are considered to be
;             outside the support.
;
;       lcx:
;             The coherence length in the x direction
;             
;       lcy:
;             The coherence length in the y direction
;
;       pxsize:
;             The size of the detector pixel in the x direction
;
;       pysize:
;             The size of the detector pixel in the y direction
;
;       beam energy: 
;             The beam energy in eV
;
;	zsd:
;	      The distance between the object and the detector
;
;       starting_point: 
;             As an option you may supply an initial 
;             guess of the exit-surface-wave for the sample. 
;             This maybe useful, for example, if you wish to 
;             start from the end point of a previous run. The
;             format of this parameter much be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the initialisation described in Harry's review paper: 
;             page 29. (in particular e.q. 137) is used.
;
; EXAMPLE:
;
;        nadia_init_partial, my_data, my_supports, 13.3e-6, $
;			    40.0e-3, 13.5e-6, 13.5e-6, 1400, 1.4, 32, 7
;
;
;-
pro nadia_init_partial, data, support, $
                      lcx, lcy, $
                      pxsize, pysize, $
                      energy, $
                      zsd, $
		      nleg, $
		      nmode, $
                      complex_array
                   
n = size(data)
nx = n[2]
ny = n[1]



IF N_Params() EQ 11 THEN BEGIN

  ;  print,'nx ' & print,nx
  ;  print,'ny ' & print,ny
  ;  print,'lcx ' & print,lcx
  ;  print,'lcy ' & print,lcy
  ;  print,'pxsize ' & print,pxsize
  ;  print,'pysize ' & print,pysize
  ;  print,'energy ' & print,energy
  ;  print,'zsd ' & print,zsd
  ;  print,'nleg ' & print,nleg
  ;  print,'nmode ' & print,nmode
    b = call_external(lib_name(),'IDL_partial_init',$ 
      nx, ny, $
      double(lcx), double(lcy), $
      double(pxsize), double(pysize), $
      double(energy), $
      double(zsd), $
      nleg, $
      nmode, $
      complex_array)

ENDIF
IF N_Params() EQ 10 THEN BEGIN
  ;print,'nx ' & print,nx
  ;print,'ny ' & print,ny 
  ;print,'lcx ' & print,lcx
  ;print,'lcy ' & print,lcy
  ;print,'pxsize ' & print,pxsize
  ;print,'pysize ' & print,pysize
  ;print,'energy ' & print,energy
  ;print,'zsd ' & print,zsd
  ;print,'nleg ' & print,nleg
  ;print,'nmode ' & print,nmode
  b = call_external(lib_name(),'IDL_partial_init', $
    nx, ny, $
    double(lcx), double(lcy), $
    double(pxsize), double(pysize), $
    double(energy), $
    double(zsd), $
    double(nleg), $
    double(nmode))
ENDIF

nadia_set_support, support
nadia_set_intensity, data

IF N_Params() EQ 10 THEN $
  nadia_initialise_esw

end

;+
; NAME:
;       NADIA_INIT_PARTIAL_CHAR
;
; PURPOSE:
;       Set-up a Partial Charactrisation CDI reconstruction. This will
;       initialise the reconstruction with sample support and experimental 
;       parameters. Some defaults will be set and memory will be
;       allocated ready for reconstructing the sample
;       exit-surface-wave. It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or calling NADIA_ITERATE).
;
;       Calling this procedure will initialise the reconstruction algorithm
;       to the error-reduction with a relaxation parameter of 0.9.
;
; CALLING SEQUENCE:
;
; NADIA_INIT_PARTIAL_CHAR, data, support, zsd, beam_energy, pxsize, pysize
;                   [, starting_point ]
;
; INPUTS:
;
;       Please use metres for the lengths, and eV for the energy.
;
;       data: 
;             The detector data with the sample in place. It should be
;             in the form of a 2D array.
;
;       support: 
;             A 2D array of integers or doubles which give the sample support.
;             Values or 1 or greater are considered inside
;             the support. All others are considered to be
;             outside the support.
;
;       pxsize:
;             The size of the detector pixel in the x direction
;
;       pysize:
;             The size of the detector pixel in the y direction
;
;       beam energy: 
;             The beam energy in eV
;
;       zsd:
;             The distance between the object and the detector
;
;       starting_point: 
;             As an option you may supply an initial 
;             guess of the exit-surface-wave for the sample. 
;             This maybe useful, for example, if you wish to 
;             start from the end point of a previous run. The
;             format of this parameter much be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the initialisation described in Harry's review paper: 
;             page 29. (in particular e.q. 137) is used.
;
; EXAMPLE:
;
;        nadia_init_partial_char, my_data, my_support, 1400,$
;				     13.5e-6, 13.5e-6, 1400, 1,4
;
;
;-
pro nadia_init_partial_char, data, support, $
  pxsize, pysize, $
  energy, $
  zsd, $
  complex_array

n = size(data)
nx = n[2]
ny = n[1]

IF N_Params() EQ 7 THEN BEGIN
  b = call_external(lib_name(),'IDL_part_char_init',long(nx),long(ny), $
  double(pxsize), double(pysize), $
  double(energy), $
  double(zsd), $
  complex_array)

ENDIF

IF N_Params() EQ 6 THEN BEGIN
  b = call_external(lib_name(),'IDL_part_char_init',long(nx),long(ny), $
  double(pxsize), double(pysize), $
  double(energy), $
  double(zsd))

ENDIF

nadia_set_support, support
nadia_set_intensity, data

IF N_Params() EQ 6 THEN $
  nadia_initialise_esw
end

;+
; NAME:
;       NADIA_INIT_POLY
;
; PURPOSE:
;       Set-up a Polychromatic CDI reconstruction. This will
;       initialise the reconstruction with sample support and experimental 
;       parameters. Some defaults will be set and memory will be
;       allocated ready for reconstructing the sample
;       exit-surface-wave. It is necessary to call this procedure before
;       attempting to call any of the reconstruction methods
;       (e.g. before setting the algorithm or calling NADIA_ITERATE).
;
;       Calling this procedure will initialise the reconstruction algorithm
;       to the error-reduction with a relaxation parameter of 0.9.
;
; CALLING SEQUENCE:
;
; NADIA_INIT_PARTIAL, data, support, beta, [, starting_point ]
;
; INPUTS:
;
;       Please use metres for the lengths, and eV for the energy.
;
;       data: 
;             The detector data with the sample in place. It should be
;             in the form of a 2D array.
;
;       support: 
;             A 2D array of integers or doubles which give the sample support.
;             Values or 1 or greater are considered inside
;             the support. All others are considered to be
;             outside the support.
;
;       beta:
;             The relaxation parameter.
;	
;       starting_point: 
;             As an option you may supply an initial 
;             guess of the exit-surface-wave for the sample. 
;             This maybe useful, for example, if you wish to 
;             start from the end point of a previous run. The
;             format of this parameter much be a 2D array of
;             COMPLEX variables. If this parameter is not supplied,
;             the initialisation described in Harry's review paper: 
;             page 29. (in particular e.q. 137) is used.
;
; EXAMPLE:
;
;        nadia_init_fresnel, my_data, my_supports
;-
pro nadia_init_poly, data, support, $
  beta, $
  complex_array

n = size(data)
nx = n[2]
ny = n[1]
IF N_Params() EQ 4 THEN $
  b = call_external(lib_name(),'IDL_poly_init',long(nx),long(ny), $
  double(beta), $
  complex_array)
IF N_Params() EQ 3 THEN $
  b = call_external(lib_name(),'IDL_poly_init',long(nx),long(ny), $
  double(beta))

IF N_Params() eq 2 THEN $
  b = call_external(lib_name(),'IDL_poly_init',long(nx),long(ny))

nadia_set_support, support
nadia_set_intensity, data

IF N_Params() LT 4 THEN $
  nadia_initialise_esw



end

;+
; NAME:
;       NADIA_SET_SUPPORT
;
; PURPOSE:
;       SET THE SUPPORT SHAPE TO BE USED IN RECONSTRUCTION. THIS WITH OVERRIDE 
;       THE SUPPORT GIVEN TO ANY OF THE NADIA_INIT METHODS AND MAYBE CALLED
;       AT ANY TIME DURING THE RECONSTRUCTION.  
;
; CALLING SEQUENCE:
;
;	NADIA_SET_SUPPORT, SUPPORT
;
; INPUTS:
;
;       SUPPORT: 
;             A 2D ARRAY OF DOUBLES OR INTEGERS GIVING THE SAMPLE'S 
;             (OR ZONE-PLATE'S) SUPPORT. VALUES OF 1 OR GREATER ARE 
;             CONSIDERED INSIDE THE SUPPORT. ALL OTHERS ARE CONSIDERED 
;             TO BE OUTSIDE THE SUPPORT.
;
; EXAMPLE:
;
;       NADIA_SET_SUPPORT, MY_SUPPORT
;
;-
pro nadia_set_support, array
n = size(array)
a = call_external(lib_name(),'IDL_set_support',n[1],n[2],double(array))
end




;+
; NAME:
;       NADIA_GET_ROUND_SUPPORT
;
; PURPOSE: 
;       This function will return a simple array with a central
;       circular region where each element has the value 1.0. Outside
;       the value returned is zero. This function has been written to
;       allow each create of the support array. For example when
;       reconstructing the white-field for Fresnel CDI reconstruction.
;       Note that it does not actually call any method from the nadia
;       software library.
;
; CALLING SEQUENCE:
;
;	NADIA_GET_ROUND_SUPPORT, nx, ny, radius
;
; INPUTS:
;
;       n_x: 
;             The number of pixels in the horizontal direction of the
;             output array.
;       n_y: 
;             The number of pixels in the vertical direction of the
;             output array.
;       radius: 
;             The radius of the circle in pixels.
;
; RETURN:
;       A two dimensional array which can be used to set the support
;       for any of the CDI reconstructions.
;
; EXAMPLE:
;
;       nadia_get_round_support, 1024, 1024, 0.25*1024
;-
function nadia_get_round_support, n_x, n_y, radius
result = make_array(n_x,n_y,/DOUBLE)
b = call_external(lib_name() ,'IDL_get_round_support', $
  long(n_x),long(n_y),double(radius),result) 
show, result
return, result
end


;+
; NAME:
;       NADIA_SET_BEAM_STOP
;
; PURPOSE: 
;       For PlanarCDI, you may set the beam stop region in the
;       detector plane. This region will be left to float (i.e. left
;       unscaled) when the intensity scaling in performed. Zeros in
;       the array (which is passed) indicate the beam-stop region,
;       all other values indicate that those pixel should be scaled.
;
; CALLING SEQUENCE:
;
;	NADIA_SET_BEAM_STOP, beam_stop_region
;
; INPUTS:
;
;       beam_stop_region: 
;             A 2D array of doubles or integers. Values of 0 are
;             considered inside the beam-stop region. All others are
;             considered to be outside the beam-stop region.
;
; EXAMPLE:
;
;       nadia_set_beam_stop, my_beam_stop_region
;
;-
pro nadia_set_beam_stop, array
n = size(array)
a = call_external(lib_name(),'IDL_set_beam_stop',n[1],n[2],double(array))
end

;+
; NAME:
;       NADIA_SET_INTENSITY
;
; PURPOSE:
;       Set the detector intensity data. This will override 
;       the intensity given to any of the NADIA_INIT methods.
;       In general, users should not need to call this method.
;
; CALLING SEQUENCE:
;
;	NADIA_SET_INTENSITY, data
;
; INPUTS:
;
;       data: 
;             The detector data. It should be in the form of a 2D
;             array or doubles or integers.
;
;
; EXAMPLE:
;
;       nadia_set_intensity, my_data
;
;-
pro nadia_set_intensity, array
n = size(array)
a= call_external(lib_name(),'IDL_set_intensity',n[1],n[2],double(array))
end

;+
; NAME:
;       NADIA_INITIALISE_ESW
;
; PURPOSE:
;       Initialise the exit-surface-wave guess. The initialisation
;       will depend on the reconstruction type (see the procedures
;       which begin "NADIA_INIT_" for a description). The "INIT" procedure  
    ;       will call this procedure if no starting guess is provided. 
    ;       It is useful if you wish to run the same reconstruction
    ;       several time with a different random starting point, or if 
    ;       you wish to reset the reconstruction back to the original guess.
    ;
    ;
    ; CALLING SEQUENCE:
    ;
    ;	NADIA_INITIALISE_ESW, seed
    ;
    ; INPUTS:
    ;
    ;       seed: 
    ;             Seed for the random number generator used to initialise
    ;             the guess. It should be an integer. This is ignored in 
    ;             the case of Fresnel CDI.
    ;
    ;
    ; EXAMPLE:
    ;
    ;       nadia_initialise_esw, 6
    ;-
    pro nadia_initialise_esw, seed
    IF N_Params() EQ 0 THEN $
      seed = 0
    b = call_external(lib_name() ,'IDL_initialise_esw',long(seed)) 
    end

    ;+
    ; NAME:
    ;       NADIA_INITIALISE_MATRICES
    ;
    ; PURPOSE:
    ;       generate the S and J matrices for the decomposition 
    ;       of the partially coherent wave where JC=nSC where
    ;       H = integral(P*l(r)J(r1, r2)Pm(r2)) dr1 dr2 and 
    ;       S=integral(P*l(r)pm(r))dr where Pl is an orhtonormal
    ;       basis set, in this case, the Legendre polynomials
    ;
    ; CALLING SEQUENCE:
    ;
    ; NADIA_INITIALISE_MATRICES, nleg, nmodes
    ;
    ; INPUTS:
    ;
    ;       nleg: 
    ;             Number of orthogonal components for the light source.
    ;             
    ;       nmodes: 
    ;             Number of component modes for the light sources. nleg
    ;             must be greater than nmodes.
    ;
    ;
    ; EXAMPLE:
    ;
    ;       nadia_initialise_matrices, 5, 6
    ;-
    pro nadia_initialise_matrices, nleg, nmodes
    b = call_external(lib_name() ,'IDL_initialise_matrices', double(nleg), double(nmodes)) 
    end

    ;+
    ; NAME:
    ;       NADIA_SET_SPECTRUM
    ;
    ; PURPOSE:
    ;	Set the spectrum. This can be read in from a Spectra
    ;	file.
    ;
    ; CALLING SEQUENCE:
    ;
    ; NADIA_SET_SPECTRU, spectrum
    ;
    ; INPUTS:
    ;
    ;       spectrum:
    ;
    ; EXAMPLE:
    ;
    ;       nadia_initialise_matrices, spectrum
    ;-
    pro nadia_set_spectrum, spectrum
    b = call_external(lib_name() ,'IDL_set_spectrum', spectrum)
    end
    ;+
    ; NAME:
    ;       NADIA_ITERATE
    ;
    ; PURPOSE:
    ;       Perform the iterative reconstruction. The number of iterations
    ;       to perform should be given and the result of the final
    ;       iteration is returned. For planar and Fresnel CDI this will be
    ;       the exit surface wave of the sample. For Fresnel white-field
    ;       reconstruction it will be the white-field at the detector surface.
    ;       The magnitude of the result is also displayed on the screen as
    ;       a 512x512 pixel image. The iteration number and corresponding
    ;       error (see NADIA_GET_ERROR) will be printed on the screen.
    ;
    ;       NADIA_ITERATE may be called several time and the reconstruction
    ;       will start from where is ended. i.e. calling NADIA_ITERATE(50),
    ;       followed by a second NADIA_ITERATE(50) is equivalent to calling
    ;       NADIA_ITERATE(100).
    ;
    ;       During the reconstruction, the best (lowest error) result is
    ;       also stored. It maybe retrieved by calling NADIA_GET_BEST_RESULT.
    ; 
    ;
    ; CALLING SEQUENCE:
    ;
    ;	result = NADIA_ITERATE([iterations])
    ;
    ; INPUTS:
    ;
    ;       iterations: 
    ;             The number of iterations to perform (an integer). 
    ;             If left empty, one iteration is performed.
    ;
    ; OUTPUTS:
    ;
    ;       result:
    ;             A COMPLEX 2D array. For planar and Fresnel CDI this will be
    ;             the exit surface wave of the sample. For Fresnel white-field
    ;             reconstruction it will be the white-field at the detector surface.
    ;
    ;
    ; EXAMPLE:
    ;
    ;       my_result = nadia_iterate(100)
    ;-
    function nadia_iterate, iterations
    nx = call_external(lib_name(),'IDL_get_array_x_size')
    ny = call_external(lib_name(),'IDL_get_array_y_size')
    result = make_array(nx,ny,/COMPLEX)
    IF N_Params() EQ 0 THEN $
      iterations = 1
    b = call_external(lib_name() ,'IDL_iterate',long(iterations),result) 
    show, ABS(result)
    return, result
    end

    ;+
    ; NAME:
    ;       NADIA_SET_RELAXATION_PARAMETER
    ;
    ; PURPOSE:
    ;       Set the relaxation parameter. The default relaxation
    ;       parameter used in Planar and Fresnel reconstruction is 0.9.
    ;       In Fresnel white-field reconstruction, this parameter is 
    ;       not used.
    ;
    ; CALLING SEQUENCE:
    ;
    ;	NADIA_SET_RELAXATION_PARAMETER, beta
    ;
    ; INPUTS:
    ;
    ;       beta: 
    ;             The relaxation parameter.
    ;
    ; EXAMPLE:
    ;
    ;       nadia_set_relaxation_parameter, 0.9
    ;-
    pro nadia_set_relaxation_parameter, beta
    b = call_external(lib_name() ,'IDL_set_relaxation_parameter',double(beta)) 
    end

    ;+
    ; NAME:
    ;       NADIA_APPLY_SHRINKWRAP
    ;
    ; PURPOSE:
    ;       Apply the shrinkwrap algorithm. The current exit-surface-wave
    ;       magnitude is used to update the support; it is convoluted with
    ;       a Gaussian and then a threshold is applied. You can use the
    ;       nadia_get_support function to see how the support have been
    ;       modified after calling this procedure.
    ;
    ; CALLING SEQUENCE:
    ;
    ;	NADIA_APPLY_SHRINKWRAP [,gauss_width, threshold ]
    ;
    ; INPUTS:
    ;
    ;       gauss_width: 
    ;             The width (1-standard deviation. in pixels) of the Gaussian
    ;             used for smearing. If this parameter is not passed, 
    ;             a width of 1.5 pixels is used.
    ;
    ;       threshold:
    ;             All pixels which are below the threshold are set to zero 
    ;             (outside the support). The threshold should be given as
    ;             a fraction of the maximum pixel value (in double format). 
    ;             If this parameter is not passed a threshold of 0.1 (10%)
    ;             is used.
    ;
    ; EXAMPLE:
    ;       Perform 1000 iterations in total, applying shrink-wrap at the 400th iteration:
    ;
    ;       ....
    ;       a = NADIA_ITERATE(400)
    ;       NADIA_APPLY_SHRINKWRAP
    ;       a = NADIA_ITERATE(600)
    ;-
    pro nadia_apply_shrinkwrap, gauss_width, threshold
    if N_Params() lt 2 then threshold = 0.1 
    if N_Params() lt 1 then gauss_width = 1.5
    b = call_external(lib_name() ,'IDL_apply_shrinkwrap',double(gauss_width),double(threshold))
    end

    ;+
    ; NAME:
    ;       NADIA_SET_ALGORITHM
    ;
    ; PURPOSE:
    ;       Select the reconstruction algorithm to use. Options are:
    ;       'ER' - error reduction 
    ;       'BIO' - basic input-output 
    ;       'BOO' - basic output-output 
    ;       'HIO' - hybrid input-output 
    ;       'DM' - difference map 
    ;       'SF' - solvent-flipping 
    ;       'ASR' - averaged successive reflections 
    ;       'HPR' - hybrid projection reflection 
    ;       'RAAR' - relaxed averaged alternating reflectors
    ;
    ; CALLING SEQUENCE:
    ;
    ;       NADIA_SET_ALGORITHM, algorithm
    ;
    ; INPUTS:
    ;
    ;       algorithm:
    ;             Set the reconstruction algorithm. It should be
    ;             one of the options listed above. Note that it is
    ;             passed as a string.
    ;
    ; EXAMPLE:
    ;       Perform 1000 iterations in total, changing to error-reduction at the 400th iteration:
    ;       ....
    ;       a = NADIA_ITERATE(400)
    ;       NADIA_SET_ALGORITHM, 'ER'
    ;       a = NADIA_ITERATE(600)
    ;-
    pro nadia_set_algorithm, algorithm
    b = call_external(lib_name() ,'IDL_set_algorithm',algorithm)
    end

    ;+
    ; NAME:
    ;       NADIA_SET_CUSTOM_ALGORITHM
    ;
    ; PURPOSE:
    ;       Set a custom reconstruction algorithm.
    ;       
    ;       Iterative reconstruction algorithms can be expressed as a
    ;       combination of several operators. The parameters to this
    ;       procedure set the coefficients for combinations of these
    ;       operators. For a description, please see page 22 of the
    ;       H.M. Quiney review: TUTORIAL REVIEW, Coherent diffractive
    ;       imaging using short wavelength light sources, Journal of
    ;       Modern Optics, 2010, DOI: 10.1080/09500340.2010.495459.  Some
    ;       description is also given in the C++ doxygen documentation
    ;       for PlanarCDI::set_custom_algorithm.
    ;
    ;
    ; CALLING SEQUENCE:
    ;
    ;       NADIA_SET_CUSTOM_ALGORITHM, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10
    ;
    ; INPUTS:
    ;
    ;       m1-m10:
    ;             Coefficients to the operator combinations.
    ;
    ; EXAMPLE:
    ;       Perform 1000 iterations in total, changing to a custom
    ;       algorithm at the 400th iteration and print the algorithm to screen:

    ;       ....
    ;       a = NADIA_ITERATE(400)
    ;       NADIA_SET_ALGORITHM, 0.5, 0.5, 0.5, 0.5, 0.5, $ 
      ;                          0.5, 0.5, 0.5, 0.5, 0.5
    ;       NADIA_PRINT_ALGORITHM                    
    ;       a = NADIA_ITERATE(600)
    ;-
    pro nadia_set_custom_algorithm, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10
    b = call_external(lib_name() ,'IDL_set_custom_algorithm', $
      double(m1), $
      double(m2), $
      double(m3), $ 
      double(m4), $ 
      double(m5), $
      double(m6), $
      double(m7), $
      double(m8), $
      double(m9), $
      double(m10) )
    end


    ;+
    ; NAME:
    ;       NADIA_GET_BEST_RESULT
    ;
    ; PURPOSE:
    ;       Get the best (lowest error) result found during the reconstruction.
    ;       The result is returned and the magnitude of the result is
    ;       displayed as a 512x512 pixel image on the screen.
    ;
    ; CALLING SEQUENCE:
    ;
    ;       result = NADIA_GET_BEST_RESULT()
    ;
    ; OUTPUTS:
    ;
    ;       result:
    ;             A 2D array of COMPLEX variables. For planar and Fresnel
    ;             reconstruction this will be the  exit-surface-wave for
    ;             the sample. For Fresnel white-field reconstruction, this
    ;             will be the complex white-field in the detector plane.
    ;
    ;
    ; EXAMPLE:
    ;       Perform 1000 iterations and get the lowest error result.
    ;       ....
    ;       a = NADIA_ITERATE(1000)
    ;       a = NADIA_GET_BEST_RESULT()
    ;-
    function nadia_get_best_result 
    result = make_array(nx(),ny(),/COMPLEX)
    b = call_external(lib_name() ,'IDL_get_best_result',result)
    show, ABS(result)
    return, result
    end

    ;+
    ; NAME:
    ;       NADIA_GET_SUPPORT
    ;
    ; PURPOSE:
    ;       Get the support. This maybe useful to see how shrinkwrap has
    ;       effected the support. The result is displayed as a 512x512 
    ;       pixel image on the screen.
    ;
    ; CALLING SEQUENCE:
    ;
    ;       support = NADIA_GET_SUPPORT()
    ;
    ; OUTPUTS:
    ;
    ;       support:
    ;             The support. A 2D array of doubles. A pixel with value below 1
    ;             is outside the support and 1 and above is inside the support.
    ;
    ; EXAMPLE:
    ;       View the support after applying shrink-wrap.
    ; 
    ;       NADIA_APPLY_SHRINKWRAP
    ;       a = NADIA_GET_SUPPORT()
    ;-
    function nadia_get_support
    result = make_array(nx(),ny(),/DOUBLE)
    b = call_external(lib_name() ,'IDL_get_support',result) 
    show, result
    return, result
    end

    ;+
    ; NAME:
    ;       NADIA_GET_ERROR
    ;
    ; PURPOSE:
    ;       Get the error metric. This is defined as the difference
    ;       between the estimated diffraction and the actual diffraction
    ;       pattern and is calculated as:
    ;           sum over pixels of ( M - sqrt(I) ) ^2 / sum(I) 
    ;       Where I is the detector data and M is the magnitude of the
    ;       estimate in the detector plane. Note that due to the way this
    ;       quantity is calculated, it actually corresponds to the
    ;       previous estimate rather than the current iteration.
    ;
    ;
    ; CALLING SEQUENCE:
    ;
    ;       error = NADIA_GET_ERROR()
    ;
    ; OUTPUTS:
    ;
    ;       error:
    ;             The error. 
    ;
    ; EXAMPLE:
    ;
    ;       a = NADIA_GET_ERROR()
    ;- 
    function nadia_get_error
    result = double(0.0)
    b = call_external(lib_name(),'IDL_get_error',result)
    return, result
    end

    ;+
    ; NAME:
    ;       NADIA_SET_TRANSMISSION
    ; 
    ; PURPOSE:
    ;       SET THE TRANSMISSION FUNCTION TO BE BEGIN THE RECONSTRUCTION. 
      ;       THIS WILL OVERRIDE THE TRANSMISSION FUNCTION AND MAY BE CALLED
      ;       AT ANY TIME DURING THE RECONSTRUCTION. NOTE THAT THIS FUNCTION 
      ;       IS ONLY AVAILIBLE FOR THE FRESNEL CDI RECONSTRUCTION. 
      ;
      ; CALLING SEQUENCE:
      ;
      ; NADIA_SET_TRANSMISSION, TRANSMISSION
      ;
      ; INPUTS:
      ;
      ;       TRANSMISSION_FUNCTION: 
      ;             The transmission function. It should be in 
      ;             the form of a 2D array.
      ;
      ; EXAMPLE:
      ;
      ;      nadia_transmission, transmission
      ;
      ;-
      pro nadia_set_transmission, array
      n = size(array)
      b = call_external(lib_name() ,'IDL_set_transmission', n[1],n[2],double(array)) 
      end

      ;+
      ; NAME:
      ;       NADIA_GET_TRANSMISSION_FUNCTION
      ;
      ; PURPOSE:
      ;       Get the transmission function from the current estimate of the
      ;       exit-surface-wave. The magnitude of the result will be
      ;       displayed on the screen. Note that this functions
      ;       is only available for Fresnel CDI reconstruction.
      ;
      ; CALLING SEQUENCE:
      ;
      ;       result = NADIA_GET_TRANSMISSION_FUNCTION()
      ;
      ; OUTPUTS:
      ;
      ;       result:
      ;             The transmission function for the sample. A 2D array of
      ;             COMPLEX variables is returned. 
      ;
      ; EXAMPLE:
      ;       a = NADIA_GET_TRANSMISSION_FUNCTION()
      ;-
      function nadia_get_transmission_function
      result = make_array(nx(),ny(),/COMPLEX)
      b = call_external(lib_name() ,'IDL_get_transmission_function',result) 
      show, abs(result)
      return, result
      end

      ;+
      ; NAME:
      ;       NADIA_GET_TRANSMISSION
      ;
      ; PURPOSE:
      ;       Get the transmission function from the current estimate of the
      ;       exit-surface-wave. The magnitude of the result will be
      ;       displayed on the screen. Note that this functions
      ;       is only available for Fresnel CDI reconstruction.
      ;
      ; CALLING SEQUENCE:
      ;
      ;       result = NADIA_GET_TRANSMISSION()
      ;
      ; OUTPUTS:
      ;
      ;       result:
      ;             The transmission function for the sample. A 2D array of
      ;             COMPLEX variables is returned. 
      ;
      ; EXAMPLE:
      ;       a = NADIA_GET_TRANSMISSION()
      ;-
      function nadia_get_transmission
      result = make_array(nx(),ny(),/COMPLEX)
      b = call_external(lib_name() ,'IDL_get_transmission',result) 
      show, abs(result)
      return, result
      end

      ;+
      ; NAME:
      ;       NADIA_GET_MODE
      ;
      ; PURPOSE:
      ;       Get the nth mode. This maybe useful to see the ways the 
      ;       modes have been divided for the partialCDI
      ;
      ; CALLING SEQUENCE:
      ;
      ;       mode = NADIA_GET_MODE(n)
      ;
      ;INPUTS
      ;       n: 
      ;             The number of the mode retrieved. 0 is the most
      ;             dominant mode. If the number is greater than the 
      ;             number of availible modes it will return the 
      ;             highest indexed modes
      ;             
      ; OUTPUTS:
      ;
      ;       mode:
      ;             The nth mode. A 2D array of doubles. 
      ;             
      ; EXAMPLE:
      ;       View the nth mode.
      ; 
      ;       a = nadia_get_mode(1) 
      ;-
      function nadia_get_mode, n
      result = make_array(nx(),ny(),/COMPLEX)
      b = call_external(lib_name() ,'IDL_get_mode',result, double(n)) 
      show, abs(result)
      return, result
      end

      ;+;+
      ; NAME:
      ;       NADIA_SET_THRESHOLD
      ; 
      ; PURPOSE:
      ;       SET THE THRESHOLD OF THE PROMINENCE OF THE LOWEST MODE 
      ;       TO BE INCLUDED IN THE RECONSTRUCTION. THIS IS TO ELIMINATE
      ;       MODES IN ORDER TO SPEED UP PROCESSING
      ;       
      ; CALLING SEQUENCE:
      ;
      ; NADIA_SET_THRESHOLD, THRESHOLD_VALUE
      ;
      ; INPUTS:
      ;
      ;       THRESHOLD_VALUE: 
      ;             The threshold value. It should be a double
      ;
      ; EXAMPLE:
      ;
      ;      nadia_set_threshold, 0.0001
      ;
      ;
      pro nadia_set_threshold, threshold_value
      b = call_external(lib_name() ,'IDL_set_transmission', threshold_value) 
      end

      ;+;+
      ; NAME:
      ;       NADIA_SET_INITIAL_COHERENCE_GUESS 
      ; 
      ; PURPOSE:
      ;       SET AN INITIAL GUESS FOR THE BEAM COHERENCE LENGTH IN PIXELS
      ;	    NOTE: WHILE BAD VALUES FOR THIS GUESS WON'T STOP CONVERGENCE,
      ;	    THEY WILL SLOW IT DOWN CONSIDERABLY
      ;
      ;       TO BE INCLUDED IN THE RECONSTRUCTION. THIS IS TO ELIMINATE
      ;       MODES IN ORDER TO SPEED UP PROCESSING
      ;       
      ; CALLING SEQUENCE:
      ;
      ; NADIA_SET_INITIAL_COHERENCE_GUESS, GUESS_X, GUESS_Y
      ;
      ; INPUTS:
      ;
      ;       GUESS_X: 
      ;             The initial guess for the coherence length in pixels in the x direction
      ;
      ;	    GUESS_Y: 
      ;             The initial guess for the coherence length in pixels in the y direction
      ;
      ; EXAMPLE:
      ;
      ;      nadia_set_initial_coherence_guess, 0.7, 0.7
      ;
      ;
      pro nadia_set_initial_coherence_guess, guess_x, guess_y
      b = call_external(lib_name() ,'IDL_set_initial_coherence_guess', guess_x, guess_y)
      end

      ;+;+
      ; NAME:
      ;       NADIA_SET_INITIAL_COHERENCE_GUESS_IN_M 
      ; 
      ; PURPOSE:
      ;       SET AN INITIAL GUESS FOR THE BEAM COHERENCE LENGTH IN METRES
      ;       NOTE: WHILE BAD VALUES FOR THIS GUESS WON'T STOP CONVERGENCE,
      ;       THEY WILL SLOW IT DOWN CONSIDERABLY
      ;
      ;       TO BE INCLUDED IN THE RECONSTRUCTION. THIS IS TO ELIMINATE
      ;       MODES IN ORDER TO SPEED UP PROCESSING
      ;       
      ; CALLING SEQUENCE:
      ;
      ; NADIA_SET_INITIAL_COHERENCE_GUESS_IN_M, GUESS_X, GUESS_Y
      ;
      ; INPUTS:
      ;
      ;       GUESS_X: 
      ;             The initial guess for the coherence length in pixels in the x direction
      ;
      ;       GUESS_Y: 
      ;             The initial guess for the coherence length in pixels in the y direction
      ;
      ; EXAMPLE:
      ;
      ;      nadia_set_initial_coherence_guess_in_m, 9e-6, 9e-6
      ;
      ;
      pro nadia_set_initial_coherence_guess_in_m, guess_x, guess_y
      b = call_external(lib_name() ,'IDL_set_initial_coherence_guess_in_m', guess_x, guess_y)
      end

      ;+;+
      ; NAME:
      ;       NADIA_SET_MINIMA_SEARCH_BOUNDS_COEFFICIENT
      ; 
      ; PURPOSE:
      ;       SET THE FACTOR THAT WILL BE USED TO DETERMINE THE SEARCH
      ;	    REGION. THE SEARCH WILL GO FROM 1/COEF TO COEF TIMES THE 
      ;	    CURRENT ESTIMATE.
      ;       
      ; CALLING SEQUENCE:
      ;
      ; NADIA_SET_MINIMA_SEARCH_BOUNDS_COEFFICIENT, COEF
      ;
      ; INPUTS:
      ;
      ;       COEF:
      ;             The coefficient that will determine the search region
      ;
      ; EXAMPLE:
      ;
      ;      nadia_set_minima_search_bounds_coefficient, 3.0
      ;
      ;
      pro nadia_set_minima_search_bounds_coefficient, coef
      b = call_external(lib_name() ,'IDL_set_minima_search_bounds_coefficient', coef)
      end

      ;+;+
      ; NAME:
      ;       NADIA_SET_MINIMA_SEARCH_TOLERANCE
      ;                                       
      ; PURPOSE:                                  
      ;       SET THE TOLERANCE OF THE LX/LY MINIMA IN PIXELS.
      ;       THE SEARCH WILL TERMNATE WHEN THE UNCERTAINTY IS 
      ;       LESS THAN THIS VALUE
      ;       
      ; CALLING SEQUENCE:
      ;
      ; NADIA_SET_MINIMA_SEARCH_TOLERANCE, TOL
      ;
      ; INPUTS:
      ;
      ;       TOL:
      ;             The tolerance of the lx/ly minima in pixels. It should be a double
      ;
      ; EXAMPLE:
      ;
      ;      nadia_set_minima_search_bounds_tolerance, 0.1
      ;
      pro nadia_set_minima_search_bounds_tolerance, tol
      b = call_external(lib_name() ,'IDL_set_minima_search_bounds_tolerance', tol)
      end

      ;+;+
      ; NAME:
      ;       NADIA_SET_MINIMA_SEARCH_TOLERANCE_IN_M
      ;                                       
      ; PURPOSE:                                  
      ;       SET THE TOLERANCE OF THE LX/LY MINIMA IN METRES.
      ;       THE SEARCH WILL TERMNATE WHEN THE UNCERTAINTY IS 
      ;       LESS THAN THIS VALUE
      ;       
      ; CALLING SEQUENCE:
      ;
      ; NADIA_SET_MINIMA_SEARCH_TOLERANCE_IN_M, TOL
      ;
      ; INPUTS:
      ;
      ;       TOL:
      ;             The tolerance of the lx/ly minima in metres. It should be a double
      ;
      ; EXAMPLE:
      ;
      ;      nadia_set_minima_search_bounds_tolerance_in_m, 1e-6
      ;
      pro nadia_set_minima_search_bounds_tolerance_in_m, tol
      b = call_external(lib_name() ,'IDL_set_minima_search_bounds_tolerance_in_m', tol)
      end

      ;+;+
      ; NAME:
      ;       NADIA_SET_MINIMA_MOVING_AVERAGE_WEIGHT
      ;                                       
      ; PURPOSE:                                  
      ;	    LX AND LY VALUES WILL BE TIME AVERAGED WITH THE PREVIOUS ESTIMATE IN
      ;	    ORDER TO SMOOTH THE SLIGHTLY ERRATIC CHANGES FROM ITERATION TO ITERATION
      ;	    W IS THE WEIGHT OF THE NEWLY CALCULATED VALUE. THE CURRENT ESTIMATE 
      ;	    WILL BE WEIGHTED AT (1-W) TO FORM A NEW CURRENT ESTIMATE. THIS IS AN 
      ;	    EXPONENTIAL MOVING AVERAGE.
      ;       
      ; CALLING SEQUENCE:
      ;	    NADIA_SET_MINIMA_MOVING_AVERAGE_WEIGHT, W
      ;
      ; INPUTS:
      ;
      ;       W:
      ;           W is the weight of the newly calculated value. the current estimate 
      ;		will be weighted at (1-w) to form a new current estimate. This is an
      ;		exponential moving average.
      ;
      ; EXAMPLE:
      ;
      ;      nadia_set_minima_moving_average_weight, 0.01
      ;
      pro nadia_set_minima_moving_average_weight, w
      b = call_external(lib_name() ,'IDL_set_minima_moving_average_weight', w)
      end

      ;+;+
      ; NAME:
      ;       NADIA_SET_MINIMA_RECALCULATION_INTERVAL
      ;                                       
      ; PURPOSE:                                  
      ;	    THE VALUES OF LX AND LY WILL ONLY BE UPADATED EVERY IVAL ITERATIONS
      ;	    (AS THE PROCESS IS SOMEWHAT TIME CONSUMING).
      ;	    THERE ISN'T MUCH ADVANTAGE TO INCREASING THE UPDATE FREQUENCY UNLESS 
      ;	    THERE ARE VERY FEW ITERATIONS BEING CALCULATED OVERALL.
      ;       
      ; CALLING SEQUENCE:
      ;       NADIA_SET_MINIMA_RECALCULATION_INTERVAL, IVAL
      ;
      ; INPUTS:
      ;
      ;       IVAL:
      ;           The number of iterations between updating lx and ly.
      ;
      ; EXAMPLE:
      ;
      ;      nadia_set_minima_recalculation_interval, 5
      ;
      pro nadia_set_minima_recalculation_interval, ival
      b = call_external(lib_name() ,'IDL_set_minima_recalculation_interval', ival)
      end

      ;+;+
      ; NAME:                                                 
      ;       NADIA_GET_X_COHERENCE_LENGTH
      ;                                                               
      ; PURPOSE:                                                          
      ;       RETRIEVE THE CALCULATED COHERENCE LENGTH OF THE BEAM (DEFINED AS THE 
      ;	    STD. DEVIATION OF THE FITTED GAUSSIAN IN THE X DIRECTION IN METRES
      ;	    AT THE OBJECT PLANE.
      ;	    THIS SHOULD BE CALLED ONLY AFTER ITERATING UNTIL A SATISFACTORY 
      ;	    IMAGE HAS BEEN PRODUCED, AND THEN IT WILL ONLY BE ACCURATE TO
      ;	    WITHIN THE TOLERANCES SET USING NADIA_SET_MINIMA_SEARCH_TOLERANCE_IN_M
      ;	    (Z*WL*NADIA_SET_MINIMA_SEARCH_TOLERANCE/(2*AVERAGE(P_SIZE_X, P_SIZE_Y)
      ;
      ; CALLING SEQUENCE:
      ;
      ;       result = NADIA_GET_X_COHERENCE_LENGTH
      ;
      ; OUTPUTS:
      ;
      ;       result:
      ;		the calculated coherence length of the beam (the std. deviation 
      ;		of the fitted gaussian) in the x direction, in meters, in the 
      ;		object plane.
      ;
      ;
      ; EXAMPLE:
      ;       Perform enough iterations to produce a satisfactory image.
      ;	    ..... 
      ;       a = NADIA_ITERATE(100)
      ;       coh_l_x = NADIA_GET_X_COHERENCE_LENGTH
      ;-
      function nadia_get_x_coherence_length
      return,  call_external(lib_name() ,'IDL_get_x_coherence_length')
      end

      ;+;+
      ; NAME:                                                 
      ;       NADIA_GET_Y_COHERENCE_LENGTH
      ;                                                               
      ; PURPOSE:                                                          
      ;       RETRIEVE THE CALCULATED COHERENCE LENGTH OF THE BEAM (DEFINED AS THE 
      ;       STD. DEVIATION OF THE FITTED GAUSSIAN IN THE Y DIRECTION IN METRES
      ;       AT THE OBJECT PLANE.
      ;       THIS SHOULD BE CALLED ONLY AFTER ITERATING UNTIL A SATISFACTORY 
      ;       IMAGE HAS BEEN PRODUCED, AND THEN IT WILL ONLY BE ACCURATE TO
      ;       WITHIN THE TOLERANCES SET USING NADIA_SET_MINIMA_SEARCH_TOLERANCE_IN_M
      ;       (Z*WL*NADIA_SET_MINIMA_SEARCH_TOLERANCE/(2*AVERAGE(P_SIZE_X, P_SIZE_Y)
      ;
      ; CALLING SEQUENCE:
      ;
      ;       result = NADIA_GET_Y_COHERENCE_LENGTH
      ;
      ; OUTPUTS:
      ;
      ;       result:
      ;           the calculated coherence length of the beam (the std. deviation 
      ;           of the fitted gaussian) in the y direction, in meters, in the 
      ;           object plane.
      ;
      ;
      ; EXAMPLE:
      ;       Perform enough iterations to produce a satisfactory image.
      ;       ..... 
      ;       a = NADIA_ITERATE(100)
      ;       coh_l_y = NADIA_GET_Y_COHERENCE_LENGTH
      ;-
      function nadia_get_y_coherence_length
      return,  call_external(lib_name() ,'IDL_get_y_coherence_length')
      end


      ;+;+
      ; NAME:                                                 
      ;       NADIA_GET_X_COHERENCE_LENGTH_IN_PIXELS
      ;                                                               
      ; PURPOSE:                                                          
      ;       RETRIEVE THE CALCULATED COHERENCE LENGTH OF THE BEAM (DEFINED AS THE 
      ;       STD. DEVIATION OF THE FITTED GAUSSIAN IN THE X DIRECTION IN PIXELS
      ;       AT THE OBJECT PLANE.
      ;       THIS SHOULD BE CALLED ONLY AFTER ITERATING UNTIL A SATISFACTORY 
      ;       IMAGE HAS BEEN PRODUCED, AND THEN IT WILL ONLY BE ACCURATE TO
      ;       WITHIN THE TOLERANCES SET USING NADIA_SET_MINIMA_SEARCH_TOLERANCE
      ;
      ; CALLING SEQUENCE:
      ;
      ;       result = NADIA_GET_X_COHERENCE_LENGTH_IN_PIXELS
      ;
      ; OUTPUTS:
      ;
      ;       result:
      ;           the calculated coherence length of the beam (the std. deviation 
      ;           of the fitted gaussian) in the x direction, in pixels, in the 
      ;           object plane.
      ;
      ;
      ; EXAMPLE:
      ;       Perform enough iterations to produce a satisfactory image.
      ;       ..... 
      ;       a = NADIA_ITERATE(100)
      ;       coh_l_x = NADIA_GET_X_COHERENCE_LENGTH_IN_PIXELS
      ;-
      function nadia_get_x_coherence_length_in_pixels
      return,  call_external(lib_name() ,'IDL_get_x_coherence_length')
      end

      ;+;+
      ; NAME:                                                 
      ;       NADIA_GET_Y_COHERENCE_LENGTH_IN_PIXELS
      ;                                                               
      ; PURPOSE:                                                          
      ;       RETRIEVE THE CALCULATED COHERENCE LENGTH OF THE BEAM (DEFINED AS THE 
      ;       STD. DEVIATION OF THE FITTED GAUSSIAN IN THE Y DIRECTION IN PIXELS
      ;       AT THE OBJECT PLANE.
      ;       THIS SHOULD BE CALLED ONLY AFTER ITERATING UNTIL A SATISFACTORY 
      ;       IMAGE HAS BEEN PRODUCED, AND THEN IT WILL ONLY BE ACCURATE TO
      ;       WITHIN THE TOLERANCES SET USING NADIA_SET_MINIMA_SEARCH_TOLERANCE
      ;                                   
      ; CALLING SEQUENCE:
      ;                                       
      ;       result = NADIA_GET_Y_COHERENCE_LENGTH_IN_PIXELS
      ;                                           
      ; OUTPUTS:
      ;                                               
      ;       result:                                     
      ;           the calculated coherence length of the beam (the std. deviation 
      ;           of the fitted gaussian) in the y direction, in pixels, in the 
      ;           object plane.
      ;                                                               
      ;                                                                   
      ; EXAMPLE:                                                              
      ;       Perform enough iterations to produce a satisfactory image.          
      ;       ..... 
      ;       a = NADIA_ITERATE(100)                                                  
      ;       coh_l_y = NADIA_GET_Y_COHERENCE_LENGTH_IN_PIXELS
      ;-
      function nadia_get_y_coherence_length_in_pixels
      return,  call_external(lib_name() ,'IDL_get_y_coherence_length')
      end



      ;
      ;       
      ; CALLING SEQUENCE:
      ;       NADIA_GET_MINIMA_RECALCULATION_INTERVAL, IVAL
      ;
      ; INPUTS:
      ;
      ;       IVAL:
      ;           The number of iterations between updating lx and ly.
      ;
      ; EXAMPLE:
      ;
      ;      nadia_set_minima_recalculation_interval, 5
      ;
      pro nadia_set_minima_recalculation_interval, ival
      b = call_external(lib_name() ,'IDL_set_minima_recalculation_interval', ival)
      end

      ; NAME:
      ;       NADIA_CLEAR_MEMORY
      ;
      ; PURPOSE:
      ;       Clean-up after a reconstruction has been performed. This 
      ;       procedure should be called at the very end of a program.
      ;       It will free up the memory that was allocated when one of 
      ;       the "NADIA_INIT_.." methods was called.
      ;
      ; CALLING SEQUENCE:
      ;
      ;       NADIA_CLEAR_MEMORY
      ;
      ; EXAMPLE:
      ;       NADIA_CLEAR_MEMORY
      ;-
      pro nadia_clear_memory
      b = call_external(lib_name() ,'IDL_deallocate_memory')
      end


      ;+
      ; NAME:
      ;       NADIA_GET_INTENSITY_AUTOCORRELATION
      ;
      ; PURPOSE:
      ;       Get the autocorrelation function from the intensity data.
      ;       This method is only useful for planar CDI reconstruction.
      ;
      ; CALLING SEQUENCE:
      ;
      ;       a = NADIA_GET_INTENSITY_AUTOCORRELATION()
      ;
      ; OUTPUTS:
      ;       a:
      ;             The autocorrelation function (a 2D array of doubles).
      ;
      ; EXAMPLE:
      ;       a = NADIA_GET_INTENSITY_AUTOCORRELATION()
      ;-
      function nadia_get_intensity_autocorrelation
      result = make_array(nx(),ny(),/DOUBLE)
      b = call_external(lib_name() ,'IDL_get_intensity_autocorrelation',result) 
      show, result
      return, result
      end



      ;+
      ; NAME:
      ;       NADIA_PRINT_ALGORITHM
      ;
      ; PURPOSE:
      ;       Output the form of the current algorithm to the screen. 
      ;       It will be written in terms of the support and intensity
      ;       projection operators.
      ;
      ;
      ; CALLING SEQUENCE:
      ;       NADIA_PRINT_ALGORITHM
      ;
      ;
      ; EXAMPLE:
      ;       NADIA_SET_ALGORITHM, 'RAAR'
      ;       NADIA_PRINT_ALGORITHM
      ;-
      pro nadia_print_algorithm
      b = call_external(lib_name() ,'IDL_print_algorithm')
      end

      ;+
      ; NAME:
      ;       NADIA_PROPAGATE_FROM_DETECTOR
      ;
      ; PURPOSE:
      ;       Propagate the given wave field from the detector plane to the 
      ;       sample plane (planar and Fresnel) or zone-plate plane (Fresnel
      ;       white-field reconstruction). The result will be returned and
      ;       displayed on the screen by default.
      ;
      ;       This function is called when using (NADIA_ITERATE) so generally
      ;       won't need to be call explicitly. An exception to this is if
      ;       the user wishes to extend the reconstruction with an addition
      ;       constraint (see the example below).
      ;
      ;
      ; CALLING SEQUENCE:
      ;
      ;	result = NADIA_PROPAGATE_FROM_DETECTOR( complex_array [,/SUPPRESS_DISPLAY] )
      ;
      ; INPUTS:
      ;
      ;       complex_array:
      ;             A 2D array of COMPLEX values which represents a wave in
      ;             the detector plane.
      ;
      ; KEYWORD PARAMETERS:
      ;
      ;       /SUPPRESS_DISPLAY:
      ;             Do not display the result on the screen. This maybe useful
      ;             if this function is used within a for loop. 
      ;
      ; OUTPUTS:
      ;
      ;       result:
      ;             A COMPLEX 2D array which is either the field in the
      ;             sample plane (for planar and Fresnel reconstruction) or
      ;             the zone-plate plane for Fresnel white-field reconstruction.
      ;
      ; EXAMPLE:
      ;       Performing the reconstruction with a new constraint applied in the
      ;       sample plane (e.g. called "NEW_SUPPORT")
      ;
      ;              FOR K = 0, 100 DO BEGIN 
	;                  a = NADIA_PROPAGATE_TO_DETECTOR(a,/SUPPRESS_DISPLAY)
	;                  a = NADIA_SCALE_INTENSITY(a,/SUPPRESS_DISPLAY)
	;                  a = NADIA_PROPAGATE_FROM_DETECTOR(a,/SUPPRESS_DISPLAY)
	;                  a = NADIA_APPLY_SUPPORT(a,/SUPPRESS_DISPLAY)
	;                  a = NEW_SUPPORT(a)
	;              ENDFOR 
	;
	;-
	function nadia_propagate_from_detector, complex_array, $
	  SUPPRESS_DISPLAY=suppress_display
	check_size, complex_array
	result = make_array(nx(),ny(),/COMPLEX)
	b = call_external(lib_name() ,'IDL_propagate_from_detector',complex_array,result) 
	if not KEYWORD_SET(suppress_display) THEN $
	  show, abs(result)
	return, result
	end

	;+
	; NAME:
	;       NADIA_PROPAGATE_TO_DETECTOR
	;
	; PURPOSE:
	;       Propagate the given wave field to the detector plane from the 
	;       sample plane (planar and Fresnel) or zone-plate plane (Fresnel
	;       white-field reconstruction). The result will be returned and
	;       displayed on the screen by default.
	;
	;       This function is called when using (NADIA_ITERATE) so generally
	;       won't need to be call explicitly. An exception to this is if
	;       the user wishes to perform simulation or to extend their
	;       reconstruction with an addition constraint (see the previous
	;       example).
	;
	;
	; CALLING SEQUENCE:
	;
	;	result = NADIA_PROPAGATE_TO_DETECTOR( complex_array [,/SUPPRESS_DISPLAY] )
	;
	; INPUTS:
	;
	;       complex_array:
	;             A 2D array of COMPLEX values which represents a wave in
	;             either the sample plane (for planar and Fresnel
	;             reconstruction) or the zone-plate plane for Fresnel 
	;             white-field reconstruction.
	;
	; KEYWORD PARAMETERS:
	;
	;       /SUPPRESS_DISPLAY:
	;             Do not display the result on the screen. This maybe useful
	;             if this function is used within a for loop. 
	;
	; OUTPUTS:
	;
	;       result:
	;             A COMPLEX 2D array which is the field in the detector plane.
	;
	; EXAMPLE:
	;       Performing Fresnel simulation (assuming you already have a
	;       complex white-field called "white_field" and the transmission
	;       function of the object, "trans").
	;
	;       white_field_at_sample = NADIA_PROPAGATE_FROM_DETECTOR(white_field)
	;       esw_at_sample = trans*white_field_at_sample
	;       esw_at_detector = NADIA_PROPAGATE_TO_DETECTOR(esw_at_sample)
	;       diffraction_pattern = abs(esw_at_detector)^2
	;-
	function nadia_propagate_to_detector, complex_array, $
	  SUPPRESS_DISPLY=suppress_display
	check_size, complex_array
	result = make_array(nx(),ny(),/COMPLEX)
	b = call_external(lib_name() ,'IDL_propagate_to_detector',complex_array,result) 
	if not KEYWORD_SET(suppress_display) THEN $
	  show, abs(result)
	return, result
	end


	;+
	; NAME:
	;       NADIA_APPLY_SUPPORT
	;
	; PURPOSE:
	;       Apply the support constraint to the given complex array. All
	;       elements outside the support with be reset to zero. Elements
	;       within the support will be left as they are. The support must
	;       have been previously set using either one of the NADIA_INIT 
	;       functions or NADIA_SET_SUPPORT.
	;
	;       This function is called when using (NADIA_ITERATE) so generally
	;       won't need to be call explicitly. An exception to this is if
	;       the user wishes to extend their reconstruction with an
	;       addition constraint.
	;
	;
	; CALLING SEQUENCE:
	;
	;	result = NADIA_APPLY_SUPPORT( complex_array [,/SUPPRESS_DISPLAY] )
	;
	; INPUTS:
	;
	;       complex_array:
	;             A 2D array of COMPLEX values which represents a wave in
	;             either the sample plane (for planar and Fresnel
	;             reconstruction) or the zone-plate plane for Fresnel 
	;             white-field reconstruction.
	;
	; KEYWORD PARAMETERS:
	;
	;       /SUPPRESS_DISPLAY:
	;             Do not display the result on the screen. This maybe useful
	;             if this function is used within a for loop. 
	;
	; OUTPUTS:
	;
	;       result:
	;             A COMPLEX 2D array after the support is applied.
	;
	; EXAMPLE:
	;       See the example for NADIA_PROPAGATE_FROM_DETECTOR
	;-
	function nadia_apply_support, complex_array, $
	  SUPPRESS_DISPLY=suppress_display
	check_size, complex_array
	result = make_array(nx(),ny(),/COMPLEX)
	b = call_external(lib_name() ,'IDL_apply_support',complex_array,result) 
	if not KEYWORD_SET(suppress_display) THEN $
	  show, abs(result)
	return, result
	end

	;+
	; NAME:
	;       NADIA_SCALE_INTENSITY
	;
	; PURPOSE:
	;       Scale the intensity (magnitude squared) of the given complex
	;       array to the data. The intensity data must have been
	;       previously set using either one of the NADIA_INIT functions or
	;       NADIA_SET_INTENSITY.
	;
	;       If Fresnel reconstruction is being done, the white-field will
	;       automatically be added to the complex array prior to scaling,
	;       and will be subtracted afterward.
	;       
	;       This function is called when using (NADIA_ITERATE) so generally
	;       won't need to be call explicitly. An exception to this is if
	;       the user wishes to extend their reconstruction with an
	;       addition constraint.
	;
	;
	; CALLING SEQUENCE:
	;
	;	result = NADIA_SCALE_INTENSITY( complex_array [,/SUPPRESS_DISPLAY] )
	;
	; INPUTS:
	;
	;       complex_array:
	;             A 2D array of COMPLEX values which represents a wave in
	;             the detector plane.
	;
	; KEYWORD PARAMETERS:
	;
	;       /SUPPRESS_DISPLAY:
	;             Do not display the result on the screen. This maybe useful
	;             if this function is used within a for loop. 
	;
	; OUTPUTS:
	;
	;       result:
	;             A COMPLEX 2D array after the intensity has been scaled
	;             to data.
	;
	; EXAMPLE:
	;       See the example for NADIA_PROPAGATE_FROM_DETECTOR
	;-
	function nadia_scale_intensity, complex_array, $
	  SUPPRESS_DISPLY=suppress_display
	check_size, complex_array
	result = make_array(nx(),ny(),/COMPLEX)
	b = call_external(lib_name() ,'IDL_scale_intensity',complex_array,result) 
	if not KEYWORD_SET(suppress_display) THEN $
	  show, abs(result)
	return, result
	end


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; io functions
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


	;+
	; NAME:
	;       NADIA_READ_PPM
	;
	; PURPOSE:
	;       Read a ppm file.
	;
	; CALLING SEQUENCE:
	;
	;       image = NADIA_READ_PPM( nx, ny, filename)
	;
	; INPUTS:
	;
	;       nx:
	;             The number of pixels horizontally (need to be checked?)
	;
	;       ny:
	;             The number of pixels vertically (need to be checked?)
	;
	;       filename:
	;             A string containing the name of the file to read.
	;
	; OUTPUTS:
	;       image:
	;             The image in the format of a 2D array of doubles.
	;
	; EXAMPLE:
	;       data = NADIA_READ_PPM(1024, 1024, 'data_file.ppm')
	;-
	function nadia_read_ppm, nx, ny, filename
	result = make_array(nx,ny, /DOUBLE)
	b = call_external(lib_name() ,'IDL_read_ppm',long(nx),long(ny),filename,result)
	show, result
	return, result
	end

	;+
	; NAME:
	;       NADIA_READ_DBIN
	;
	; PURPOSE:
	;       Read a double binary file (a binary file of doubles).
	;
	; CALLING SEQUENCE:
	;
	;       image = NADIA_READ_DBIN( nx, ny, filename)
	;
	; INPUTS:
	;
	;       nx:
	;             The number of pixels horizontally (need to be checked?)
	;
	;       ny:
	;             The number of pixels vertically (need to be checked?)
	;
	;       filename:
	;             A string containing the name of the file to read.
	;
	; OUTPUTS:
	;       image:
	;             The image in the format of a 2D array of doubles.
	;
	; EXAMPLE:
	;       data = NADIA_READ_DBIN(1024, 1024, 'data_file.dbin')
	;-
	function nadia_read_dbin, nx, ny, filename
	result = make_array(nx,ny, /DOUBLE)
	b = call_external(lib_name() ,'IDL_read_dbin',long(nx),long(ny),filename,result)
	show, result
	return, result
	end

	;+
	; NAME:
	;       NADIA_READ_TIFF
	;
	; PURPOSE:
	;       Read a tiff image file. Note that there are other IDL commands
	;       which perform this same function.
	;
	; CALLING SEQUENCE:
	;
	;       image = NADIA_READ_TIFF( nx, ny, filename)
	;
	; INPUTS:
	;
	;       nx:
	;             The number of pixels horizontally (need to be checked?)
	;
	;       ny:
	;             The number of pixels vertically (need to be checked?)
	;
	;       filename:
	;             A string containing the name of the file to read.
	;
	; OUTPUTS:
	;       image:
	;             The image in the format of a 2D array of doubles.
	;
	; EXAMPLE:
	;       data = NADIA_READ_TIFF(1024, 1024, 'data_file.tif')
	;-
	function nadia_read_tiff, nx, ny, filename
	result = make_array(nx,ny, /DOUBLE)
	b = call_external(lib_name() ,'IDL_read_tiff',long(nx),long(ny),filename,result)
	show, result
	return, result
	end

	;+
	; NAME:
	;       NADIA_READ_CPLX
	;
	; PURPOSE:
	;       Read a binary file of fftw complex numbers. This is
	;       useful for storing and restoring the result of a reconstruction.
	;
	; CALLING SEQUENCE:
	;
	;       complex_array = NADIA_READ_CPLX( nx, ny, filename)
	;
	; INPUTS:
	;
	;       nx:
	;             The number of pixels horizontally (need to be checked?)
	;
	;       ny:
	;             The number of pixels vertically (need to be checked?)
	;
	;       filename:
	;             A string containing the name of the file to read.
	;
	; OUTPUTS:
	;       complex_array:
	;             A 2D array of COMPLEX numbers.
	;
	; EXAMPLE:
	;       white_field = NADIA_READ_CPLX(1024, 1024, 'white_field_file.cplx')
	;-
	function nadia_read_cplx, nx, ny, filename
	result = make_array(nx,ny, /COMPLEX)
	b = call_external(lib_name() ,'IDL_read_cplx',long(nx),long(ny),filename,result)
	show, abs(result)
	return, result
	end

	;+
	; NAME:
	;       NADIA_WRITE_CPLX
	;
	; PURPOSE:
	;       Write a binary file of fftw complex numbers. This is
	;       useful for storing and restoring the result of a reconstruction.
	;
	; CALLING SEQUENCE:
	;
	;       NADIA_WRITE_CPLX( complex_array, filename )
	;
	; INPUTS:
	;
	;       complex_array:
	;             A 2D array of COMPLEX numbers
	;
	;       filename:
	;             A string containing the name of the file to write.
	;
	; EXAMPLE:
	;       NADIA_WRITE_CPLX, white_field, 'white_field_file.cplx'
	;-
	pro nadia_write_cplx, complex_array, filename
	n = size(complex_array)
	b = call_external(lib_name() ,'IDL_write_cplx',n[1],n[2],complex_array, filename)
	end

	;+
	; NAME:
	;       NADIA_WRITE_DBIN
	;
	; PURPOSE:
	;       Write a binary file of doubles. 
	;
	; CALLING SEQUENCE:
	;
	;       NADIA_WRITE_DBIN( image, filename )
	;
	; INPUTS:
	;
	;       image:
	;             A 2D array of numbers
	;
	;       filename:
	;             A string containing the name of the file to write.
	;
	; EXAMPLE:
	;       NADIA_WRITE_DBIN, result, 'reco_magnitude.dbin'
	;-
	pro nadia_write_dbin, array, filename
	n = size(array)
	b = call_external(lib_name() ,'IDL_write_dbin',n[1],n[2],double(array), filename)
	end

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; Complex Constraint Code
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;+
	; NAME:
	;       NADIA_SET_CHARGE_FLIPPING
	;
	; PURPOSE: For use with FresnelCDI reconstruction. Enabeling this
	;       procedure will constrain the transmission function phase to
	;       lie between -PI and 0. During each iteration (directly after
	;       applying the support constraint), the phase of the
	;       transmission function will be flipped if it is positive
	;       (i.e. if the phase, phi, lies between 0 and PI it will be
	;       reset to -phi). If this procedure is called for a PlanarCDI
	;       reconstruction, the contraint will be applied to the
	;       exit-surface-wave.
	;
	; CALLING SEQUENCE:
	;
	;       NADIA_SET_CHARGE_FLIPPING, enable
	;
	; INPUTS:
	;
	;       enable:
	;             This should be either 0 - turn off or 1 - turn on.
	;
	; EXAMPLE:
	;       NADIA_SET_CHARGE_FLIPPING, 1
	;-
	pro nadia_set_charge_flipping, enable
	b = call_external(lib_name() ,'IDL_set_charge_flipping',long(enable))
	end

	;+
	; NAME:
	;       NADIA_SET_TRANS_UNITY_CONSTRAINT
	;
	; PURPOSE: For use with FresnelCDI reconstruction. Enabeling this
	;       procedure will constrain the transmission function magnitude
	;       to lie 0 and 1. During each iteration (directly after applying
	;       the support constraint), the magnitude of the transmission
	;       function will be reset to 1 at locations where it is greater
	;       than 1.
	;
	; CALLING SEQUENCE:
	;
	;       NADIA_SET_TRANS_UNITY_CONSTRAINT, enable
	;
	; INPUTS:
	;
	;       enable:
	;             This should be either 0 - turn off or 1 - turn on.
	;
	; EXAMPLE:
	;       NADIA_SET_TRANS_UNITY_CONSTRAINT, 1
	;-
	pro nadia_set_trans_unity_constraint, enable
	b = call_external(lib_name() ,'IDL_set_trans_unity_constraint',long(enable))
	end


	;+
	; NAME:
	;       NADIA_ADD_COMPLEX_CONSTRAINT_REGION
	;
	; PURPOSE: For use with FresnelCDI reconstruction. Calling this
	;       procedure will constrain the transmission function magnitude
	;       and phase using the method from the paper "Use of a complex
	;       constraint in coherent diffractive imaging", J. N. Clark
	;       et. al., 2010. The notation used there is also used
	;       here. Users should have knowledge of this or similar papers
	;       before using the procedure.
	;
	;       A section of the reconstructed transmission function will be
	;       updated according to a restriction on the value of c (c =
	;       beta/delta). If the material is of a known element, then c can
	;       be fixed to the know value, and the phase and magnitude of the
	;       transmission function updated accordingly. Alternatively, c
	;       can be left to float, and calculated for each iteration from a
	;       mean over the defined region. The parameters alpha1 and alpha2
	;       are used to control the strength of the constraint. This
	;       procedure should be called once before each region of
	;       homogenious material.
	;
	; CALLING SEQUENCE:
	;
	;       NADIA_ADD_COMPLEX_CONSTRAINT_REGION, region, alpha1, alpha2 [, fixed_c]
	;
	; INPUTS:
	;
	;       region:
	;             A 2 array of doubles which is used to indicate which
	;             pixels the constraint should be applied to. 0 -
	;             indicated the array element does not belong to the
	;             region. Any other value indicated that it does.
	;
	;       alpha1:
	;             Constraint strength parameter for the amplitude (double
	;             type).
	;
	;       alpha2:
	;             Constraint strength parameter for the phase (double
	;             type).
	;             
	;       fixed_c:
	;             Optional parameter to fix the value of c = beta/delta
	;             (double type).
	;
	; EXAMPLE:
	;       Setting a complex constraint for two image regions.
	;       For the second, the value of c=beta/delta is fixed.
	;
	;       r1 =  nadia_read_tiff(1024,1024,'region_1.tiff')
	;       r2 =  nadia_read_tiff(1024,1024,'region_2.tiff')
	;
	;       delta = 6.45e-4;
	;       beta = 1.43e-4;
	;
	;       nadia_add_complex_constraint_region, r1, 0.5, 0.5
	;       nadia_add_complex_constraint_region, r2,   1,   0, beta/delta
	;-
	pro nadia_add_complex_constraint_region, region, alpha1, alpha2, fixed_c
	IF N_Params() EQ 3 THEN $
	  b = call_external(lib_name() ,'IDL_add_complex_constraint_region', $
	  region, double(alpha1), double(alpha2)) $
	  ELSE $
	  b = call_external(lib_name() ,'IDL_add_complex_constraint_region', $
	  region, double(alpha1), double(alpha2), double(fixed_c))  
	end


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;   Phase Diverse Code
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;

	;
	function nadia_phase_diverse_init_estimate
	;dimensions are reversed.
	nx = call_external(lib_name(),'IDL_get_phase_diverse_array_x_size')
	ny = call_external(lib_name(),'IDL_get_phase_diverse_array_y_size')
	result = make_array(nx,ny,/COMPLEX)
	b = call_external(lib_name(),'IDL_phase_diverse_init_estimate',result)
	show, ABS(result)
	return, result
	end


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;
	;
	pro nadia_init_phase_diverse, beta = beta, gamma = gamma , parallel = parallel

	; set default values
	if not keyword_set(beta) then beta = 1.0
	if not keyword_set(gamma) then gamma = 1.0
	if not keyword_set(parallel) then parallel = 0

	b = call_external(lib_name(),'IDL_phase_diverse_init', $
	  double(beta), double(gamma), long(parallel))

	end


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;
	;
	pro nadia_phase_diverse_add_position, x, y, alpha = alpha
	; set default values
	if not keyword_set(alpha) then alpha = 1.0
	b = call_external(lib_name(),'IDL_phase_diverse_add_position', $
	  double(x), double(y), double(alpha))
	end

	function nadia_phase_diverse_iterate, iterations
	; dimensions are reversed in IDL compared to the library.
	nx = call_external(lib_name(),'IDL_get_phase_diverse_array_x_size')
	ny = call_external(lib_name(),'IDL_get_phase_diverse_array_y_size')
	result = make_array(nx,ny,/COMPLEX)
	IF N_Params() EQ 0 THEN $
	  iterations = 1
	b = call_external(lib_name() ,'IDL_phase_diverse_iterate',long(iterations),result) 
	show, ABS(result)
	return, result
	end


	pro nadia_phase_diverse_iterations_per_cycle, iterations
	b = call_external(lib_name(),'IDL_phase_diverse_iterations_per_cycle', long(iterations))
	end

	pro nadia_phase_diverse_set_transmission, array
	n = size(array)
	b = call_external(lib_name(),'IDL_phase_diverse_set_transmission', long(n[2]), long(n[1]), double(array))
	end

	pro nadia_phase_diverse_adjust_positions, type=type, forward=forward, $
	  x_min=x_min, x_max=x_max, $
	  y_min=y_min, y_max=y_max, $
	  step_size=step_size

	if not keyword_set(type) then type = 0
	if not keyword_set(forward) then forward = 1
	if not keyword_set(x_min) then x_min = -50
	if not keyword_set(x_max) then x_max = 50
	if not keyword_set(y_min) then y_min = -50
	if not keyword_set(y_max) then y_max = 50
	if not keyword_set(step_size) then step_size = 4

	b = call_external(lib_name(),'IDL_phase_diverse_adjust_positions', $
	  long(type), long(forward), $
	  long(x_min),  long(x_max), $
	  long(y_min),  long(y_max), $
	  double(step_size))

	end

	function nadia_phase_diverse_get_final_x_position, nprobe
	return, call_external(lib_name(), $
	  'IDL_phase_diverse_get_final_x_position', $
	  long(nprobe))
	end

	function nadia_phase_diverse_get_final_y_position, nprobe
	return, call_external(lib_name(), $
	  'IDL_phase_diverse_get_final_y_position', $
	  long(nprobe))
	end
