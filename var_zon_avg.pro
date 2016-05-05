;----------------------------------------------------------------
;                                                                  
; NAME:
;      VAR_ZON_AVG
;                                                                            
; PURPOSE:
;      plot the vertical profile of a variety of variables from 
;      model runs in order to view the difference between the models
;                    
; CALLING SEQUENCE:
;      VAR_ZON_AVG, REFC1=REFC1, VAIRH=VAIRH, PSC=PSC, CLO=CLO,
;                   H2O=H2O, NO2=NO2, XNO=XNO, POS=POS
;
; INPUTS:
;
;  KEYWORD PARAMETERS:
;      REFC1 - to plot output from Kane's REFC1 model
;      VAIRH - to plot output from Matt's varih model
;      PSC  - to plot PSC mass mixing ratio
;      CLO  - to plot ClO mass mixing ratio
;      H2O  - to plot water mass mixing ratio
;      NO2  - to plot NO2 mass mixing ratio
;      XNO  - to plot NO mass mixing ratio
;      POS   - to save as a post script file
;                    
; SUBROUTINES:
;      N/A
;                 
; NOTES:
;      None
;                    
; EXAMPLES
;      var_zon_avg,/vairh,/psc
;                                                                         
; CREATED:      
;      03.02.16 by Kaitlyn Lieshcke           
;                                                                          
;-----------------------------------------------------------------  

pro var_zon_avg, refc1=refc1,vairh=vairh,psc=psc,clo=clo,h2o=h2o,$
                 no2=no2,xno=xno,pos=pos

;=========================================================      
; READ IN THE VAIRH MODEL OUTPUT
;=========================================================

if keyword_set(vairh) then begin

; Setup
modelpath = '/g/data1/p66/mxw599/ACCESS/output/vairh/'
nyrs = 2010-2005+1
;nmms = 12 * nyrs
var = fltarr(192,145,85,4)
pr = fltarr(192,145,85,4)
model = 'VAIRH'
if keyword_set(psc) then begin
   field = 'field34158'
   species = 'PSC'
endif
if keyword_set(clo) then begin
   field = 'field34042'
   species = 'ClO'
endif
if keyword_set(h2o) then begin
   field = 'hus_tl'
   species = 'H2O'
endif
if keyword_set(no2) then begin
   field = 'field34153'
   species = 'NO2'
endif
if keyword_set(xno) then begin
   field = 'field34002'
   species = 'NO'
endif

;Create variable to count passing through years
yr=0

for year = 2005,2010 do begin

   ;Create string variable of year
   syyyy = string(year,'(i4.4)')

   print,'Reading '+syyyy+' files'

   for mon = 0,11 do begin

      ;Create string variable of month
      smm = string(mon+1,'(i2.2)')

      ; Read model data
      ncdf_read, struct, file=modelpath+'vairha.pm-'+syyyy+'-'+smm+'.nc',$
                 variables=[field,'lon','lat','pfull']

      ;Save output
      ;Calculate number of days in nyears of each season
      ndsum=(31+31+28.25)*nyrs
      ndaut=(31+30+31)*nyrs
      ndwin=(30+31+31)*nyrs
      ndspr=(30+31+30)*nyrs

      ;Add fraction of seasonal mean to final array
      ;Assign the command to execute based on the variable desired
      jan = 'var[*,*,*,0]=var[*,*,*,0]+(struct.'+field+'*(31./ndsum))'
      feb = 'var[*,*,*,0]=var[*,*,*,0]+(struct.'+field+'*(28.25/ndsum))'
      mar = 'var[*,*,*,1]=var[*,*,*,1]+(struct.'+field+'*(31./ndaut))'
      apr = 'var[*,*,*,1]=var[*,*,*,1]+(struct.'+field+'*(30./ndaut))'
      may = 'var[*,*,*,1]=var[*,*,*,1]+(struct.'+field+'*(31./ndaut))'
      jun = 'var[*,*,*,2]=var[*,*,*,2]+(struct.'+field+'*(30./ndwin))'
      jul = 'var[*,*,*,2]=var[*,*,*,2]+(struct.'+field+'*(31./ndwin))'
      aug = 'var[*,*,*,2]=var[*,*,*,2]+(struct.'+field+'*(31./ndwin))'
      sep = 'var[*,*,*,3]=var[*,*,*,3]+(struct.'+field+'*(30./ndspr))'
      oct = 'var[*,*,*,3]=var[*,*,*,3]+(struct.'+field+'*(31./ndspr))'
      nov = 'var[*,*,*,3]=var[*,*,*,3]+(struct.'+field+'*(30./ndspr))'
      dec = 'var[*,*,*,3]=var[*,*,*,3]+(struct.'+field+'*(31./ndsum))'

      if mon eq 0 then begin
         com = execute(jan)
         pr[*,*,*,0]=pr[*,*,*,0]+(struct.pfull*(31./ndsum))
      endif
      if mon eq 1 then begin
         com = execute(feb)
         pr[*,*,*,0]=pr[*,*,*,0]+(struct.pfull*(28.25/ndsum))
       endif
      if mon eq 2 then begin
         com = execute(mar)
         pr[*,*,*,1]=pr[*,*,*,1]+(struct.pfull*(31./ndaut))
      endif
      if mon eq 3 then begin
         com = execute(apr)
         pr[*,*,*,1]=pr[*,*,*,1]+(struct.pfull*(30./ndaut))
      endif
      if mon eq 4 then begin
         com = execute(may)
         pr[*,*,*,1]=pr[*,*,*,1]+(struct.pfull*(31./ndaut))
      endif
      if mon eq 5 then begin
         com = execute(jun)
         pr[*,*,*,2]=pr[*,*,*,2]+(struct.pfull*(30./ndwin))
      endif
      if mon eq 6 then begin
         com = execute(jul)
         pr[*,*,*,2]=pr[*,*,*,2]+(struct.pfull*(31./ndwin))
      endif
      if mon eq 7 then begin
         com = execute(aug)
         pr[*,*,*,2]=pr[*,*,*,2]+(struct.pfull*(31./ndwin))
      endif
      if mon eq 8 then begin
         com = execute(sep)
         pr[*,*,*,3]=pr[*,*,*,3]+(struct.pfull*(30./ndspr))
      endif
      if mon eq 9 then begin
         com = execute(oct)
         pr[*,*,*,3]=pr[*,*,*,3]+(struct.pfull*(31./ndspr))
      endif
      if mon eq 10 then begin
         com = execute(nov)
         pr[*,*,*,3]=pr[*,*,*,3]+(struct.pfull*(30./ndspr))
      endif
      if mon eq 11 then begin
         com = execute(dec)
         pr[*,*,*,0]=pr[*,*,*,0]+(struct.pfull*(31./ndsum))
      endif

      if mon eq 0 then begin
      lat=struct.lat
      lon=struct.lon
      endif

   endfor

;Add one to variable counting years
yr = yr + 1

endfor

endif

;=========================================================      
; READ IN THE REFC1 MODEL OUTPUT
;=========================================================

if keyword_set(refc1) then begin

; Setup *NB: file is mislabelled and PSCs are 'temp' variable
modelpath = '/g/data1/p66/kjl574/Kanes_output/'
nyrs = 2010-2005+1
;nmms = 12 * nyrs
model = 'REFC1'
if keyword_set(psc) then begin
   varpath = 'ACCESSCCM_REFC1_NATPSC_seas.nc'
   species = 'PSC'
   field = 'temp'
endif
if keyword_set(clo) then begin
   varpath = 'ACCESSCCM_REFC1_ClO_seas.nc'
   species = 'ClO'
   field = 'field542'
endif
if keyword_set(h2o) then begin
   varpath = 'ACCESSCCM_REFC1_SpecHum_seas.nc'
   species = 'H2O'
   field = 'q'
endif
if keyword_set(no2) then begin
   varpath = 'ACCESSCCM_REFC1_NO2_seas.nc'
   species = 'NO2'
   field = 'field1861'
endif
if keyword_set(xno) then begin
   varpath = 'ACCESSCCM_REFC1_NO_seas.nc'
   species = 'NO'
   field = 'tracer2'
endif

; Read model data
ncdf_read, struct, file=modelpath+varpath,$
           variables=[field,'latitude','longitude']
print,'Reading PSC file'

; Save output
lat = struct.latitude
lon = struct.longitude
if keyword_set(psc) then var = struct.temp
if keyword_set(clo) then var = struct.field542
if keyword_set(h2o) then var = struct.q
if keyword_set(no2) then var = struct.field1861
if keyword_set(xno) then var = struct.tracer2

; Read model data
ncdf_read, struct, file=modelpath+'ACCESSCCM_REFC1_Pressure_seas.nc',$
           variables=['p']
print,'Reading pressure file'

; Save output
pr = struct.p

endif

;=========================================================      
; INTERPOLATE PRESSURES ONTO CONSISTENT LEVELS
;=========================================================

print,'Interpolating pressures'

; Setup for pressure interpolation
intpr = reverse(10^((indgen(390)-80.)/100))*1d2

;Create variable to store interpolated values
var_int = fltarr(n_elements(lon),n_elements(lat),n_elements(intpr),4)

; Loop over lat,lon and seasons to regrid pressure levels
for time = 0,3 do begin

   for tmplat = 0,n_elements(lat)-1 do begin

      for tmplon = 0,n_elements(lon)-1 do begin
         
         var_int[tmplon,tmplat,*,time]=interpol(var[tmplon,tmplat,*,time],pr[tmplon,tmplat,*,time],intpr)
      endfor

   endfor

endfor

; Test interpolation code
;; plot, var[73,7,*,2], logpr[73,7,*,2], yrange=[12,-2],xrange=[0,5e-10]
;; oplot, var_int[73,7,*,2], intpr,color=2
;; plot, var[73,7,*,2], pr[73,7,*,2], yrange=[1000,100],xrange=[0,4e-10]
;; oplot, var_int[73,7,*,2], intpr,color=2

; Average longitudinally and change units from MMR if necessary
if keyword_set(psc) then begin
zon_var_int = reform(mean(var_int,1)) * 1e12
units='MMR * 1e12'
endif
if keyword_set(clo) then begin
zon_var_int = reform(mean(var_int,1)) * 1e9
units='MMR * 1e9'
mindata=0.
maxdata=1.5
endif
if keyword_set(h2o) then begin
zon_var_int = reform(mean(var_int,1)) * 1e5
units='MMR * 1d5'
mindata=0
maxdata=10.
endif
if keyword_set(no2) then begin
zon_var_int = reform(mean(var_int,1)) * 1e11
units='MMR * 1e11'
endif
if keyword_set(xno) then begin
zon_var_int = reform(mean(var_int,1)) * 1e8
units='MMR * 1e8'
mindata = 0.
maxdata = 3.
endif

; Convert pressure from Pa to hPa
intpr = intpr / 100

;=========================================================      
; PLOTTING
;=========================================================
if keyword_set(pos) then begin
ps_setup,filename=model+'_'+species+'_zon_avg.ps',/open,/portrait,$
         /color
!p.font=0
!p.charsize=1.
endif else begin
window,0,xsize=1300,ysize=800
!p.charsize=2
endelse
!x.margin=[3,3]
!y.margin=[3,3]
!x.omargin=[2,0]
multipanel,rows=2,cols=2

tvplot,reform(zon_var_int[*,*,0]),lat,intpr,/cbar,cbunit=units,$
       title='Summer '+species,/fcontour,$
       yrange=[1000,1],/ylog,ytitle='Log Pressure (hPa)',$
       xtitle='Latitude',mindata=mindata,maxdata=maxdata

tvplot,reform(zon_var_int[*,*,1]),lat,intpr,/cbar,cbunit=units,$
       title='Autumn '+species,/fcontour,$
       yrange=[1000,1],/ylog,ytitle='Pressure (hPa)',$
       xtitle='Latitude',mindata=mindata,maxdata=maxdata

tvplot,reform(zon_var_int[*,*,2]),lat,intpr,/cbar,cbunit=units,$
       title='Winter '+species,/fcontour,$
       yrange=[1000,1],/ylog,ytitle='Pressure (hPa)',$
       xtitle='Latitude',mindata=mindata,maxdata=maxdata

tvplot,reform(zon_var_int[*,*,3]),lat,intpr,/cbar,cbunit=units,$
       title='Spring '+species,/fcontour,$
       yrange=[1000,1],ytitle='Pressure (hPa)',/ylog,$
       xtitle='Latitude',mindata=mindata,maxdata=maxdata

if keyword_set(pos) then ps_setup,/close,/noview

stop
end
