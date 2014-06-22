"""
This is an attempt at re-implementing a Brand & Blitz kinematic distance tool
in Python, based on Chris Beaumont's kdist.pro and using the new fancy parameters 
from Mark Reid et al. (2014).

"""

from __future__ import division

import numpy as np

import astropy.table
from astropy.coordinates import Galactic
import astropy.units as u
import astropy.constants as c

brand_blitz_parameters = dict()

def brand_blitz_kinematic_distance(longitude, velocity, galaxy_parameters=brand_blitz_parameters):
	"""
	Computes a kdist for a given (l,b) pair.

	- notation:
	 v : rotation velocity
	 vo: Solar rotation velocity
	 vr: Radial velocity wrt LSR
	 r : Galactocentric radius
	 ro: Solar galactocentric radius

	rotation curve parameters from Brand and Blitz 1993
	 v/vo = a1 * (R/Ro)^a2 + a3
	 vr = sin(l) * vo * [ a1 * (r/ro)^(a2-1) + a3*(r/ro)^-1 - 1]

	"""
	l = np.radians(longitude) % 2*np.pi

	a1=1
	a2=0
	a3=0

	#- the following are from Reid et al. (2014)
	ro=8.34
	vo=240 # aka theta_0

	r = (np.arange(2000) + 1) / 2000 * 2 * ro
	root = np.sin(l) * vo * (a1 * (r/ro)**(a2-1) + a3 * (r/ro)**(-1) - 1) - velocity

	#find the zero crossing
	root *= np.roll(root,1)
	root[0] = 1
	root[-1] = 1

	hit = np.where(root < 0)

	if len(hit[0]) < 1:
		print 'Cannot determine galactocentric distance'
		print 'returned 0 (NaN) solutions'
		return [np.nan, np.nan]
	r=r[hit[0]][0]

	warnmsg = 'Warning: The galactocentric distance extrapolates the measured galactic roation curve'
	if (r < 2 or r > 17):
		print warnmsg

	# Having determined "r", we now need to solve the triangle and find
	# "d". If "l" is acute, there are (usually) two possible solutions; if
	# "l" is obtuse, then there should be one unique solution.

	# Case I: obtuse l - 1 solution

	if (l >= np.pi/2) and (l < 3*np.pi/2): 

	   # "Gamma" and "Beta" refer to the angles opposite "ro" and "d", 
	   # respectively. (I made a triangle to sketch this all out.)
	   gamma = np.arcsin( ro / r * np.sin(l))
	   beta = np.pi - ( l + gamma)
	   d = r * np.sin(beta) / np.sin( l )

	   print 'returned one solution'
	   return [d, d]

	# Case II: acute l - 2 or 0 solutions

	else:

	   rmin = ro*np.cos(l)
	   dr = np.sqrt(r**2-(ro*np.sin(l))**2)
	   if ~np.isfinite(dr):
	      print 'motion cannot be reproduced via galactic rotation'
	      print 'returned 0 (NaN) solutions'
	      return [np.nan, np.nan]
	   else:
	      print 'returned 2 solutions'
	      return [rmin-dr,rmin+dr]




''';+
; NAME:
;  KDIST2
;
; PURPOSE:
;  This function calculates the near and far kinematic distances for a
;  given galactic longitude and radial velocity. The calculation uses
;  the Galactic rotaion model of Brand and Blitz 1993, A&A, 275 : 67.
;
; CATEGORY:
;  coordinate systems
;
; CALLING SEQUENCE:
;  result=KDIST2( longitude, velocity, [/RADIANS])
;
; INPUTS:
;  longitude: Galactic Longitude. Currently must be in the range [-180,
;  180] in degrees.
;
;  velocity: Radial velocity in km/s
;
; KEYWORD PARAMETERS:
;  RADIANS: If set, input Longitude is in radians
;  DEBUG: If set, produce debugging plots / information
;
; OUTPUTS:
;  The two element vector [near_distance, far_distance] in kpc.
;
; RESTRICTIONS:
;  Functionality for outer galaxy is not well-characterized.
;
; MODIFICATION HISTORY:
;  Written by: Chris Beaumont, June 2008.
;  June 23, 2008: Changed name from kinematic_distance to kdist. cnb
;  June 23, 2008: Fixed bug in modding l with 2 pi. cnb
;  July 17, 2008: Removed degrees keyword. Added radians. cnb
;  July 17, 2008: Changed things so that, if l is a variable,
;  it isn't modified
;  March 18, 2009: Changed theta_0 and v_sol to reflect arxiv
;  0902.3913
;  March 20, 2009: Changed minor typo in value of a1. Added /DEBUG
;  keyword
;  March 11, 2012: Changed name from kdist to kidst2. Tom Rice,
;  t.rice90@gmail.com (tsr)
;  March 11, 2012: Replaced all references to "latitude" with
;  "longitude". tsr
;  March 11, 2012: Added support for the outer galaxy (l obtuse);
;  whether it's doing the math correctly is unconfirmed, but
;  the code works. tsr
;  March 15, 2012: Now if we get a "zero solution" case the code
;  returns two NaNs rather than crashing with an error.
;  March 15, 2012: Also, if we get a "cannot determine galactocentric
;  distance" error the code returns two NaNs rather than an error.
;-
FUNCTION kdist2, longitude, velocity, radians=radians, debug = debug
compile_opt idl2
on_error, 2

if n_params() != 2 then begin
    message,'Calling Sequence: dist=kdist2(l,v,[/radians])'
endif

if ~keyword_set(radians) then l=longitude/!radeg else l=longitude
l=(l mod (2*!dpi))
if(l < 0) then l += 2 * !dpi

; We used to only allow acute angles...
;if (l >= np.pi/2) && (l < 3*np.pi/2) then message, $ 
;   'Error -- Longitude must be acute'

;- notation:
; v : rotation velocity
; vo: Solar rotation velocity
; vr: Radial velocity wrt LSR
; r : Galactocentric radius
; ro: Solar galactocentric radius

;rotation curve parameters from Brand and Blitz 1993
; v/vo = a1 * (R/Ro)**a2 + a3
; vr = np.sin(l) * vo * [ a1 * (r/ro)^(a2-1) + a3*(r/ro)^-1 - 1]

a1=1.0077D
a2=.0394D
a3=.0077D

;- the following are from arxiv 0902.3913 (VLBI parallax)
;- value taken from section 4, the best fit to this particular roation curve
ro=8.8 
vo=275.
;ro = 8.5 ; - Brand Blitz value
;vo = 225

;determine r
r = (findgen(2000) + 1) / 2000. * 2 * ro
root=sin(l) * vo * (a1 * (r/ro)^(a2-1) + a3 * (r/ro)^(-1) - 1) - velocity
backup = root

;find the zero crossing
root*=shift(root,1)
root[0]=1

root[n_elements(root)-1]=1

hit=where(root < 0, ct)
if keyword_set(debug) then begin
   plot, r, backup, xtitle = 'Galactocentric Distance (kpc)', $
         ytitle = 'Relative radial velocity (galaxy - object)', $
         charsize = 1.8, yra = [-20, 20]
endif

if ct != 1 then begin
   print 'Cannot determine galactocentric distance'
   print, 'returned 0 (NaN) solutions'
   return, [np.nan, np.nan]
endif
r=r[hit[0]]

warnmsg = 'Warning: The galactocentric distance extrapolates the measured galactic roation curve'
if (r < 2 || r > 17) then $
   message, /continue, warnmsg

if keyword_set(debug) then begin
   oplot, r * [1,1], [-20, 20], color = fsc_color('crimson')
endif


; Having determined "r", we now need to solve the triangle and find
; "d". If "l" is acute, there are (usually) two possible solutions; if
; "l" is obtuse, then there should be one unique solution.


; Case I: obtuse l - 1 solution

if (l >= np.pi/2) && (l < 3*np.pi/2) then begin 

   ; "Gamma" and "Beta" refer to the angles opposite "ro" and "d", 
   ; respectively. (I made a triangle to sketch this all out.)
   gamma = np.arcsin( ro / r * sin(l))
   beta = np.pi - ( l + gamma)
   d = r * sin(beta) / sin( l )

   print, 'returned one solution'
   return, [d, d]

; Case II: acute l - 2 or 0 solutions

endif else begin

   rmin=ro*np.cos(l)
   dr=np.sqrt(r^2-(ro*sin(l))^2)
   if ~np.isfinite(dr) then begin
      print, 'motion cannot be reproduced via galactic rotation'
      print, 'returned 0 (NaN) solutions'
      return, [np.nan, np.nan]
   endif else begin


      print, 'returned 2 solutions'
      return,[rmin-dr,rmin+dr]
      
   endelse
   
   
endelse


end
'''