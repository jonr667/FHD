FUNCTION fhd_struct_init_meta,file_path_vis,hdr,params,lon=lon,lat=lat,alt=alt,$
    zenra=zenra,zendec=zendec,obsra=obsra,obsdec=obsdec,phasera=phasera,phasedec=phasedec,$
    rephase_to_zenith=rephase_to_zenith,precess=precess,degpix=degpix,dimension=dimension,elements=elements,$
    obsx=obsx,obsy=obsy,instrument=instrument,mirror_X=mirror_X,mirror_Y=mirror_Y,no_rephase=no_rephase,$
    meta_data=meta_data,meta_hdr=meta_hdr,_Extra=extra

IF N_Elements(instrument) EQ 0 THEN instrument=''
IF N_Elements(lon) EQ 0 THEN lon=116.67081524 & lon=Float(lon);degrees
IF N_Elements(lat) EQ 0 THEN lat=-26.7033194 & lat=Float(lat);degrees
IF N_Elements(alt) EQ 0 THEN alt=377.83 & alt=Float(alt);altitude (meters)
metafits_ext='.metafits'
metafits_dir=file_dirname(file_path_vis)
metafits_name=file_basename(file_path_vis,'.uvfits',/fold_case)
metafits_name=file_basename(metafits_name,'_cal',/fold_case) ;sometimes "_cal" is present, sometimes not.
metafits_path=metafits_dir+path_sep()+metafits_name+metafits_ext

time=params.time
b0i=Uniq(time)
jdate=double(hdr.jd0)+time[b0i]

IF N_Elements(dimension) EQ 0 THEN dimension=1024.
IF N_Elements(elements) EQ 0 THEN elements=dimension
IF N_Elements(obsx) EQ 0 THEN obsx=dimension/2.
IF N_Elements(obsy) EQ 0 THEN obsy=elements/2.

degpix2=[degpix,degpix]
IF Keyword_Set(mirror_X) THEN degpix2[0]*=-1
IF Keyword_Set(mirror_Y) THEN degpix2[1]*=-1
n_pol=hdr.n_pol

IF file_test(metafits_path) THEN BEGIN
    meta_hdr=headfits(metafits_path,exten=0,/silent)
    
    meta_data=mrdfits(metafits_path,1,hdr1,/silent)
    tile_nums=meta_data.antenna
    tile_nums=radix_sort(tile_nums,index=tile_order)
    meta_data=meta_data[tile_order]
    pol_names=meta_data.pol
    single_i=where(pol_names EQ pol_names[0],n_single)
    tile_names=meta_data.tile
    tile_names=tile_names[single_i]
    tile_height=meta_data.height
    tile_height=tile_height[single_i]-alt
    tile_flag=Ptrarr(n_pol) & FOR pol_i=0,n_pol-1 DO tile_flag[pol_i]=Ptr_new(meta_data(single_i+pol_i).flag)
    
    obsra=sxpar(meta_hdr,'RA')
    obsdec=sxpar(meta_hdr,'Dec')
    phasera=sxpar(meta_hdr,'RAPHASE')
    phasedec=sxpar(meta_hdr,'DECPHASE')
    
;    LST=sxpar(meta_hdr,'LST')
;    HA=sxpar(meta_hdr,'HA')
;    HA=ten([Fix(Strmid(HA,0,2)),Fix(Strmid(HA,3,2)),Fix(Strmid(HA,6,2))])*15.
    date_obs=sxpar(meta_hdr,'DATE-OBS')
    JD0=date_conv(date_obs,'JULIAN')
    
;    IF instrument EQ 'mwa' THEN BEGIN
;        ;if the MWA, check if the uvfits files were created during a period when there was a one-day offset in the time
;        cotter_date=date_conv(sxpar(meta_hdr,'DATE'),'REAL')
;        IF (cotter_date LT 2014115.) AND (cotter_date GT 2013115.) THEN BEGIN
;            JD0+=1.
;            date_obs=date_conv(jd0,'FITS')
;        ENDIF
;    ENDIF
    ct2lst,LST_hr,lon,0,JD0
    LST=LST_hr*360./24.
    
    zenra=LST
    zendec=lat
    epoch=date_conv(date_obs,'REAL')/1000.
    Precess,zenra,zendec,epoch,2000.    
    
    beamformer_delays=sxpar(meta_hdr,'DELAYS')
    beamformer_delays=Ptr_new(Float(Strsplit(beamformer_delays,',',/extract)))
ENDIF ELSE BEGIN
    ;use hdr and params to guess metadata
    print,'### NOTE ###'
    print,'Metafits file not found! Calculating obs settings from the uvfits header instead'
    
;    print,metafits_path+' not found. Calculating obs settings from the uvfits header instead'
    ;256 tile upper limit is hard-coded in CASA format
    ;these tile numbers have been verified to be correct
    tile_A1=Long(Floor(params.baseline_arr/256)) ;tile numbers start from 1
    tile_B1=Long(Fix(params.baseline_arr mod 256))
    hist_A1=histogram(tile_A1,min=0,max=256,/binsize,reverse_ind=ria)
    hist_B1=histogram(tile_B1,min=0,max=256,/binsize,reverse_ind=rib)
    hist_AB=hist_A1+hist_B1
    tile_names=where(hist_AB,n_tile)
    tile_height=Fltarr(n_tile)
    tile_flag=Ptrarr(n_pol) & FOR pol_i=0,n_pol-1 DO tile_flag[pol_i]=Ptr_new(intarr(n_tile))
    date_obs=hdr.date
    JD0=date_string_to_julian(date_obs)
    epoch=date_conv(hdr.date)/1000.
    
    IF ~Keyword_Set(time_offset) THEN time_offset=0d
    time_offset/=(24.*3600.)
    JD0=Min(Jdate)+time_offset
    
    obsra=hdr.obsra
    obsdec=hdr.obsdec
    IF Keyword_Set(precess) THEN Precess,obsra,obsdec,epoch,2000.
    IF N_Elements(phasera) EQ 0 THEN phasera=obsra
    IF N_Elements(phasedec) EQ 0 THEN phasedec=obsdec
;    Precess,obsra,obsdec,2000.,epoch

    IF Keyword_Set(zenra) THEN BEGIN
        IF Keyword_Set(precess) THEN BEGIN
            IF N_Elements(zendec) EQ 0 THEN zendec=lat
            Precess,zenra,zendec,epoch,2000.
        ENDIF ELSE BEGIN
            IF N_Elements(zendec) EQ 0 THEN BEGIN
                zendec=lat
                zenra0=zenra
                Precess,zenra0,zendec,epoch,2000. ;slight error, since zenra0 is NOT in J2000, but assume the effect on zendec is small
            ENDIF
        ENDELSE
    ENDIF ELSE zenpos2,JD0,zenra,zendec, lat=lat, lng=lon,/degree,/J2000
    beamformer_delays=Ptr_new()
ENDELSE

orig_phasera=phasera
orig_phasedec=phasedec
IF Keyword_Set(no_rephase) THEN BEGIN
    phasera=obsra
    phasedec=obsdec
ENDIF

;IF Abs(obsra-zenra) LT degpix THEN zenra=obsra
;IF Abs(obsdec-zendec) LT degpix THEN zendec=obsdec

IF Keyword_Set(rephase_to_zenith) THEN BEGIN
    phasera=obsra
    phasedec=obsdec
    obsra=zenra
    obsdec=zendec
ENDIF
projection_slant_orthographic,astr=astr,degpix=degpix2,obsra=obsra,obsdec=obsdec,zenra=zenra,zendec=zendec,$
    dimension=dimension,elements=elements,obsx=obsx,obsy=obsy,zenx=zenx,zeny=zeny,phasera=phasera,phasedec=phasedec,$
    epoch=2000.,JDate=JD0,date_obs=date_obs

Eq2Hor,obsra,obsdec,JD0,obsalt,obsaz,lat=lat,lon=lon,alt=Mean(alt)
meta={obsra:Float(obsra),obsdec:Float(obsdec),zenra:Float(zenra),zendec:Float(zendec),phasera:Float(phasera),phasedec:Float(phasedec),$
    epoch:Float(epoch),tile_names:tile_names,lon:Float(lon),lat:Float(lat),alt:Float(alt),JD0:Double(JD0),Jdate:Double(Jdate),astr:astr,$
    obsx:Float(obsx),obsy:Float(obsy),zenx:Float(zenx),zeny:Float(zeny),obsaz:Float(obsaz),obsalt:Float(obsalt),$
    delays:beamformer_delays,tile_height:Float(tile_height),tile_flag:tile_flag,orig_phasera:Float(orig_phasera),orig_phasedec:Float(orig_phasedec)}

RETURN,meta
END