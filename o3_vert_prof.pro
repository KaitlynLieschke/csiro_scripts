;----------------------------------------------------------------      
;                                                                   
; NAME:                                                          
;      O3_VERT_PROF
;                                                                  
; PURPOSE:                                                          
;      plot the vertical profile of ozone in number density from
;      model runs and ozone sondes in order to view the bias in
;      the models
;                                                                     
; CALLING SEQUENCE:                                             
;      O3_VERT_PROF,VAIRH=VAIRH,VALNZ=VALNZ,REFC1=REFC1,POS=POS
;
; INPUTS:
;
;  KEYWORD PARAMETERS:
;     VAIRH - to plot output from Matt's VAIRH model
;     VALNZ - to plot output from Matt's VALNZ model
;     REFC1 - to plot output from Kane's REFC1 model
;     POS   - to save as a post script file
;                                                              
; SUBROUTINES:                                                 
;      N/A                                                       
;                                                                
; NOTES:                                                         
;      None                                                       
;                                                                   
; EXAMPLES:                                              
;      o3_vert_prof,/vair,/refc,/pos
;                                                                  
; CREATED:                                                         
;      12.01.16 by Kaitlyn Lieshcke                                 
;                                                                          
;-----------------------------------------------------------------      

pro o3_vert_prof,vair=vair,valn=valn,refc=refc,pos=pos

;=========================================================     
; OVERALL SETUP
;========================================================= 

latitudes = [-37.5,-45,-54.6,-68.5,-90]
longitudes = [145,169.7,158.9,79,169]

;=========================================================     
; READ IN THE VAIRH MODEL OUTPUT (NUDGED MODEL)
;========================================================= 

if keyword_set(vairh) then begin

; SETUP                                               
nyrs = 2010 - 2005 + 1
variables = ['lon','lat','field34001', 'pfull','theta']
modelpath = '/g/data1/p66/mxw599/ACCESS/output/vairh/'
vh_o3_tmp_all = fltarr(5,nyrs,12,85)
vh_pr_tmp_all = fltarr(5,nyrs,12,85)
vh_th_tmp_all = fltarr(5,nyrs,12,85)
vh_o3_seas_avg = fltarr(5,4,85)
vh_pr_seas_avg = fltarr(5,4,85)
vh_tp_seas_avg = fltarr(5,4,85)
vh_o3_seas_std = fltarr(5,4,85)
vh_tp_seas_std = fltarr(5,4,85)

; Count the years passed through in the loop
vyear = 0

; Loop over years to average
for yyyy = 2005, 2010 do begin

   ; Count the months passed through in each loop
   vmonth = 0

   ; Setup string for year
   syyyy = string(yyyy, '(i4.4)')
   print, 'Reading files for '+syyyy

   ; Loop over each onth in the year
   for mm = 1,12 do begin

      ; Setup string for month
      smm = string(mm, '(i2.2)')

      ; Read model file
      ncdf_read, UKvair, file=modelpath+'vairha.pm-'+syyyy+'-'+smm+'.nc',$
                 variable = variables

      ; Save the mass mixing ratio and air pressure
      o3_tmp = UKvair.field34001
      pr_tmp = UKvair.pfull
      th_tmp = UKvair.theta
      if yyyy eq 2005 and mm eq 1 then begin
         UKlat = UKvair.lat
         UKlon = UKvair.lon
      endif

      ; Loop through coordinates for each gridbox
      for ll = 0,4 do begin

         ; Find grid box containing location
         near = min(abs(UKlat - latitudes[ll]),latindex)
         near = min(abs(UKlon - longitudes[ll]),lonindex)

         ; Add ozone mixing ratio to calculate average over season
         vh_o3_tmp_all[ll,vyear,vmonth,*] = reform(o3_tmp[lonindex,latindex,*])
         vh_pr_tmp_all[ll,vyear,vmonth,*] = reform(pr_tmp[lonindex,latindex,*])
         vh_th_tmp_all[ll,vyear,vmonth,*] = reform(th_tmp[lonindex,latindex,*])

      endfor

      ; Add one to variable counting months
      vmonth = vmonth + 1

   endfor

   ; Add one to variable counting years
   vyear = vyear + 1

endfor

; CONVERSIONS FOR PLOTTING
; Convert ozone mass mixing ratio to number density (molecules*d12 /cm-3)
; Calculate air temperature from potential tempersture
vh_tp_tmp_all = vh_th_tmp_all / ((100000/vh_pr_tmp_all)^(0.286))
;Calculate air density from temperature and pressure (in kg/cm3)
vh_dn_tmp_all = vh_pr_tmp_all / (vh_tp_tmp_all * 287.058) * 1e-6
;MMR / molar mass * avogadro's number * density
vh_o3_tmp_all = (vh_o3_tmp_all / 0.048) * 6.022d23 * vh_dn_tmp_all


; Seasonally average the vairh ozone vertical profile over all 10 years
; Set months of seasons over which to average
ssmm = [[0,1,11],[2,3,4],[5,6,7],[8,9,10]]
for isite = 0,4 do begin
   for ilev = 0,84 do begin
      for ii = 0,3 do begin
         vh_o3_seas_avg[isite,ii,ilev] = mean(vh_o3_tmp_all[isite,*,ssmm[*,ii],ilev])
         vh_o3_seas_std[isite,ii,ilev] = stddev(vh_o3_tmp_all[isite,*,ssmm[*,ii],ilev])
         vh_pr_seas_avg[isite,ii,ilev] = mean(vh_pr_tmp_all[isite,*,ssmm[*,ii],ilev])
         vh_tp_seas_avg[isite,ii,ilev] = mean(vh_tp_tmp_all[isite,*,ssmm[*,ii],ilev])
         vh_tp_seas_std[isite,ii,ilev] = stddev(vh_tp_tmp_all[isite,*,ssmm[*,ii],ilev])
      endfor
   endfor
endfor

;convert Pa to hPa
vh_pr_seas_avg = vh_pr_seas_avg / 1d2
;Convert molecs/cm3 to molecs*d12/cm3
vh_o3_seas_avg = vh_o3_seas_avg / 1d12
vh_o3_seas_std = vh_o3_seas_std / 1d12

endif

;=========================================================     
; READ IN THE VALNZ MODEL OUTPUT (UN-NUDGED MODEL)
;========================================================= 

if keyword_set(valnz) then begin

; SETUP                                               
nyrs = 2010 - 2005 + 1
variables = ['lon','lat','field34001', 'pfull','theta']
modelpath = '/g/data1/p66/mxw599/ACCESS/output/valnz/'
vl_o3_tmp_all = fltarr(5,nyrs,12,85)
vl_pr_tmp_all = fltarr(5,nyrs,12,85)
vl_th_tmp_all = fltarr(5,nyrs,12,85)
vl_o3_seas_avg = fltarr(5,4,85)
vl_pr_seas_avg = fltarr(5,4,85)
vl_tp_seas_avg = fltarr(5,4,85)
vl_o3_seas_std = fltarr(5,4,85)
vl_tp_seas_std = fltarr(5,4,85)

; Count the years passed through in the loop
vyear = 0

; Loop over years to average
for yyyy = 2005, 2010 do begin

   ; Count the months passed through in each loop
   vmonth = 0

   ; Setup string for year
   syyyy = string(yyyy, '(i4.4)')
   print, 'Reading files for '+syyyy

   ; Loop over each onth in the year
   for mm = 1,12 do begin

      ; Setup string for month
      smm = string(mm, '(i2.2)')

      ; Read model file
      ncdf_read, UKvair, file=modelpath+'valnza.pm-'+syyyy+'-'+smm+'.nc',$
                 variable = variables

      ; Save the mass mixing ratio and air pressure
      o3_tmp = UKvair.field34001
      pr_tmp = UKvair.pfull
      th_tmp = UKvair.theta
      if yyyy eq 2005 and mm eq 1 then begin
         UKlat = UKvair.lat
         UKlon = UKvair.lon
      endif

      ; Loop through coordinates for each gridbox
      for ll = 0,4 do begin

         ; Find grid box containing location
         near = min(abs(UKlat - latitudes[ll]),latindex)
         near = min(abs(UKlon - longitudes[ll]),lonindex)

         ; Add ozone mixing ratio to calculate average over season
         vl_o3_tmp_all[ll,vyear,vmonth,*] = reform(o3_tmp[lonindex,latindex,*])
         vl_pr_tmp_all[ll,vyear,vmonth,*] = reform(pr_tmp[lonindex,latindex,*])
         vl_th_tmp_all[ll,vyear,vmonth,*] = reform(th_tmp[lonindex,latindex,*])

      endfor

      ; Add one to variable counting months
      vmonth = vmonth + 1

   endfor

   ; Add one to variable counting years
   vyear = vyear + 1

endfor

; CONVERSIONS FOR PLOTTING
; Convert ozone mass mixing ratio to number density (molecules*d12 /cm-3)
; Calculate air temperature from potential tempersture
vl_tp_tmp_all = vl_th_tmp_all / ((100000/vl_pr_tmp_all)^(0.286))
;Calculate air density from temperature and pressure (in kg/cm3)
vl_dn_tmp_all = vl_pr_tmp_all / (vl_tp_tmp_all * 287.058) * 1e-6
;MMR / molar mass * avogadro's number * density
vl_o3_tmp_all = (vl_o3_tmp_all / 0.048) * 6.022d23 * vl_dn_tmp_all

; Seasonally average the valnz ozone vertical profile over all 10 years
; Set months of seasons over which to average
ssmm = [[0,1,11],[2,3,4],[5,6,7],[8,9,10]]
for isite = 0,4 do begin
   for ilev = 0,84 do begin
      for ii = 0,3 do begin
         vl_o3_seas_avg[isite,ii,ilev] = mean(vl_o3_tmp_all[isite,*,ssmm[*,ii],ilev])
         vl_o3_seas_std[isite,ii,ilev] = stddev(vl_o3_tmp_all[isite,*,ssmm[*,ii],ilev])
         vl_pr_seas_avg[isite,ii,ilev] = mean(vl_pr_tmp_all[isite,*,ssmm[*,ii],ilev])
         vl_tp_seas_avg[isite,ii,ilev] = mean(vl_tp_tmp_all[isite,*,ssmm[*,ii],ilev])
         vl_tp_seas_std[isite,ii,ilev] = stddev(vl_tp_tmp_all[isite,*,ssmm[*,ii],ilev])
      endfor
   endfor
endfor

;convert Pa to hPa
vl_pr_seas_avg = vl_pr_seas_avg / 1d2
;Convert molecs/cm3 to molecs*d12/cm3
vl_o3_seas_avg = vl_o3_seas_avg / 1d12
vl_o3_seas_std = vl_o3_seas_std / 1d12

endif

;=========================================================     
; READ IN KANES MODEL OUTPUT
;========================================================= 

if keyword_set(refc1) then begin

; SETUP
files = ['ozone','Temperature','Pressure']
variables = ['tracer1','temp','p']
pathname = '/g/data1/p66/kjl574/Kanes_output/'
kn_o3_tmp_all = fltarr(5,60,6,12)
kn_tp_tmp_all = fltarr(5,60,6,12)
kn_pr_tmp_all = fltarr(5,60,6,12)
kn_o3_seas_avg = fltarr(5,60,4)
kn_tp_seas_avg = fltarr(5,60,4)
kn_pr_seas_avg = fltarr(5,60,4)
kn_o3_seas_std = fltarr(5,60,4)
kn_tp_seas_std = fltarr(5,60,4)

; NB: time is in days since 1940-12-01 00:00:00 (360 day calendar)

;Loop over each file to be read
for vv = 0,2 do begin

   ; Read model file
   ncdf_read, knstruct,file=pathname+'ACCESSCCM_REFC1_'+files[vv]+'_short_lonlat.nc', $
              variable = ['longitude','latitude',variables[vv]]
   print, 'Reading '+files[vv]+' file'

   ; Save latitudes and longitudes from file
   if vv eq 0 then begin
      knlat = knstruct.latitude
      knlon = knstruct.longitude
   endif

   ; Loop through coordinates for each gridbox
   for ll = 0,4 do begin

      ; Find grid box containing location
      near = min(abs(knlat - latitudes[ll]),latindex)
      near = min(abs(knlon - longitudes[ll]),lonindex)

      ; Save and reform each variable to only include spring mean
      if vv eq 0 then begin
         o3_tmp = knstruct.tracer1

         ; Count total months in file
         totmm = 0

         ; Loop through years
         for years=0,5 do begin
            ; Loop through months in a year
            for yrmm = 0,11 do begin
               kn_o3_tmp_all[ll,*,years,yrmm]=o3_tmp[lonindex,latindex,*,totmm]

               totmm = totmm + 1
            endfor
         endfor
      endif
      if vv eq 1 then begin
         tp_tmp = knstruct.temp

         ; Count total months in file
         totmm = 0

         ; Loop through years
         for years=0,5 do begin
            ; Loop through months in a year
            for yrmm = 0,11 do begin
               kn_tp_tmp_all[ll,*,years,yrmm]=tp_tmp[lonindex,latindex,*,totmm]

               totmm = totmm + 1
            endfor
         endfor
      endif
      if vv eq 2 then begin
         pr_tmp = knstruct.p

         ; Count total months in file
         totmm = 0

         ; Loop through years
         for years=0,5 do begin
            ; Loop through months in a year
            for yrmm = 0,11 do begin
               kn_pr_tmp_all[ll,*,years,yrmm]=pr_tmp[lonindex,latindex,*,totmm]

               totmm = totmm + 1
            endfor
         endfor
      endif

   endfor

endfor

; Convert ozone mass mixing ratio to number density (molecules*d12 /cm-3)
;Calculate air density from temperature and pressure (in kg/cm3)
kn_dn_tmp_all = kn_pr_tmp_all / (kn_tp_tmp_all * 287.058) * 1e-6
;MMR / molar mass * avogadro's number * density
kn_o3_tmp_all = (kn_o3_tmp_all / 0.048) * 6.022d23 * kn_dn_tmp_all

; Seasonally average the REFC1 ozone vertical profile over all 10 years
; Set months of seasons over which to average
ssmm = [[0,1,11],[2,3,4],[5,6,7],[8,9,10]]
for isite = 0,4 do begin
   for ilev = 0,59 do begin
      for ii = 0,3 do begin
         kn_o3_seas_avg[isite,ilev,ii] = mean(kn_o3_tmp_all[isite,ilev,*,ssmm[*,ii]])
         kn_o3_seas_std[isite,ilev,ii] = stddev(kn_o3_tmp_all[isite,ilev,*,ssmm[*,ii]])
         kn_pr_seas_avg[isite,ilev,ii] = mean(kn_pr_tmp_all[isite,ilev,*,ssmm[*,ii]])
         kn_tp_seas_avg[isite,ilev,ii] = mean(kn_tp_tmp_all[isite,ilev,*,ssmm[*,ii]])
         kn_tp_seas_std[isite,ilev,ii] = stddev(kn_tp_tmp_all[isite,ilev,*,ssmm[*,ii]])
      endfor
   endfor
endfor

;convert Pa to hPa
kn_pr_seas_avg = kn_pr_seas_avg / 1d2
;Convert molecs/cm3 to molecs*d12/cm3
kn_o3_seas_avg = kn_o3_seas_avg / 1d12
kn_o3_seas_std = kn_o3_seas_std / 1d12

endif

;=========================================================     
; READ IN SONDES DATA
;========================================================= 

; Assign path for ozone sondes data
ospath = '/g/data1/p66/mxw599/obs/ozonesondes_from_Kane/Seasonal_means/'

; Set up variables
os_pr_seas_avg = fltarr(5,64,4)
os_o3_seas_avg = fltarr(5,64,4)
os_o3_seas_std = fltarr(5,64,4)
os_tp_seas_avg = fltarr(5,64,4)
os_tp_seas_std = fltarr(5,64,4)

; Loop over locations and seasons
sites = ['Melbourne','Lauder','Macquarie','Davis','South_Pole']
seasons=['sum','aut','win','spr']

for pp = 0,4 do begin

   print,'Reading '+sites[pp]+' ozone sondes files'

   for jj = 0,3 do begin
   
      ; Open file
      openr,lun,ospath+sites[pp]+'_'+seasons[jj]+'_2005-2011.txt',/get_lun

      count=0
      ; Set up variables
      tmp_pr=0.0D
      os_pr=0.0D
      tmp_o3=0.0D
      os_o3=0.0D
      tmpostd=0.0D
      os_ostd=0.0D
      tmp_tp=0.0D
      os_tp=0.0D
      tmp_tstd=0.0D
      os_tstd=0.0D

      ; Read file
;      print,'Reading '+sites[pp]+' '+seasons[jj]+' ozone sondes file'
      while ~eof(lun) do begin
         readf,lun,tmp_pr,tmp_o3,tmp_ostd,tmp_tp,tmp_tstd,$
               format='(3x,e13.1,3x,d13.7,3x,d13.7,3x,d13.5,3x,d13.5)'
         os_pr=[os_pr,tmp_pr]
         os_o3=[os_o3,tmp_o3]
         os_ostd=[os_ostd,tmp_ostd]
         os_tp=[os_tp,tmp_tp]
         os_tstd=[os_tstd,tmp_tstd]
         count=count+1
      endwhile

      ; Close file and free memory
      close,lun
      free_lun,lun

      ; Remove zero values from start of array
      os_pr_tmp=os_pr[1:*]
      os_o3_tmp=os_o3[1:*]/1d12 ;change to molecs x 10^12 /cm3 
      os_ostd_tmp=os_o3[1:*]/1d12 ;change to molecs x 10^12 /cm
      os_tp_tmp=os_tp[1:*]
      os_tstd_tmp=os_tstd[1:*]

      ; Add data to variable containing all data
      os_pr_seas_avg[pp,*,jj] = os_pr_tmp
      os_o3_seas_avg[pp,*,jj] = os_o3_tmp
      os_o3_seas_std[pp,*,jj] = os_ostd_tmp
      os_tp_seas_avg[pp,*,jj] = os_tp_tmp
      os_tp_seas_std[pp,*,jj] = os_tstd_tmp

   endfor
endfor

;=========================================================     
; PLOTTING
;========================================================= 
if keyword_set(pos) then begin
ps_setup,filename='o3_vert_prof.ps',/portrait,/color,/open
!p.font=0
!p.charsize=1.
endif else begin
; Set window size
window,0,xsize=1400,ysize=800
;Adjust font size
!p.charsize=2
endelse
!x.margin=[3,3]
!y.margin=[3,3]
!x.omargin=[2,0]
; Plot twenty graphs
multipanel,rows=5,cols=5

; Set up legend
legend, line=intarr(4),thick=intarr(4)+2,label=['vairh','valnz','REFC1','Sondes'],$
        lcolor=[1,3,4,2],position=[0.13,0.84,0.18,0.9]
   multipanel,/advance
   multipanel,/advance
   multipanel,/advance
   multipanel,/advance
   multipanel,/advance

; Assign axes
seasons = ['Summer','Autumn','Winter','Spring']
sites = ['Melbourne','Lauder','Macquarie','Davis','South Pole']

; Loop over season
for jj = 0,3 do begin

   ; Loop over location
   for pp = 0,4 do begin

      if keyword_set(vair) then begin
      ; Plot ozone vertical profile (number density vs pressure) for vairh
      plot, reform(vh_o3_seas_avg[pp,jj,*]),reform(vh_pr_seas_avg[pp,jj,*]),$
            /ylog,yrange=[1000,1],thick=2,xrange=[0,7]
      endif

      if keyword_set(valn) then begin
      plot,reform(vl_o3_seas_avg[pp,*,jj]),reform(vl_pr_seas_avg[pp,*,jj]),$
            color='3',thick=2,/ylog,yrange=[1000,1]
      endif

      if keyword_set(refc) then begin
      oplot,reform(kn_o3_seas_avg[pp,*,jj]),reform(kn_pr_seas_avg[pp,*,jj]),$
            color='4',thick=2
      endif

      oplot,reform(os_o3_seas_avg[pp,*,jj]),reform(os_pr_seas_avg[pp,*,jj]),$
            color='2',thick=2

   endfor
endfor

if keyword_set(pos) then ps_setup,/close,/noview

stop
end
