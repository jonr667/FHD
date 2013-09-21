PRO fhd_quickview,obs,psf,cal,image_uv_arr=image_uv_arr,weights_arr=weights_arr,source_array=source_array,$
    model_uv_arr=model_uv_arr,file_path_fhd=file_path_fhd,silent=silent,$
    gridline_image_show=gridline_image_show,pad_uv_image=pad_uv_image,$
    zoom_radius=zoom_radius,zoom_low=zoom_low,zoom_high=zoom_high,grid_spacing=grid_spacing,reverse_image=reverse_image,$
    _Extra=extra
n_pol=obs.npol
dimension_uv=obs.dimension
astr=obs.astr

IF Keyword_Set(pad_uv_image) THEN BEGIN
    pad_uv_image=pad_uv_image>1.

    obs_out=obs
    astr_out=astr
    astr_out.cdelt/=pad_uv_image
    astr_out.crpix*=pad_uv_image
    
    obs_out.astr=astr_out
    obs_out.dimension*=pad_uv_image
    obs_out.elements*=pad_uv_image
    obs_out.obsx*=pad_uv_image
    obs_out.obsy*=pad_uv_image
    obs_out.zenx*=pad_uv_image
    obs_out.zeny*=pad_uv_image
    obs_out.degpix/=pad_uv_image
ENDIF ELSE obs_out=obs
dimension=obs_out.dimension
elements=obs_out.elements
degpix=obs_out.degpix
astr_out=obs_out.astr

zoom_radius=Round(18./(degpix)/16.)*16.
zoom_low=dimension/2.-zoom_radius
zoom_high=dimension/2.+zoom_radius-1
;stats_radius=10. ;degrees
pol_names=['xx','yy','xy','yx','I','Q','U','V']

grid_spacing=10.
offset_lat=5.;15. paper 10 memo
offset_lon=5.;15. paper 10 memo
reverse_image=0   ;1: reverse x axis, 2: y-axis, 3: reverse both x and y axes
map_reverse=0;1 paper 3 memo
label_spacing=1.




;Build a fits header
mkhdr,fits_header,*residual_array[0]
putast, fits_header, astr_out;, cd_type=1




IF N_Elements(source_array) GT 0 THEN BEGIN
    si_use=where(source_array.ston GE fhd.sigma_cut,ns_use)
    source_arr=source_array[si_use]
    source_arr_out=source_array
    sx=(source_array.x-obs.dimension/2.)*2.+obs_out.dimension/2.
    sy=(source_array.y-obs.elements/2.)*2.+obs_out.elements/2.
    source_arr_out.x=sx & source_arr_out.y=sy
    
    extend_test=where(Ptr_valid(source_arr_out.extend),n_extend)
    IF n_extend GT 0 THEN BEGIN
        FOR ext_i=0L,n_extend-1 DO BEGIN
            comp_arr_out=*source_array[extend_test[ext_i]].extend
            ad2xy,comp_arr_out.ra,comp_arr_out.dec,astr_out,cx,cy
            comp_arr_out.x=cx & comp_arr_out.y=cy
            source_arr_out[extend_test[ext_i]].extend=Ptr_new(/allocate)
            *source_arr_out[extend_test[ext_i]].extend=comp_arr_out
        ENDFOR
    ENDIF
ENDIF
END
;
;
;;,fhd,obs,image_uv_arr,model_uv_holo,source_array,comp_arr,beam_base,$
;;    file_path_fhd=file_path_fhd,image_filter_fn=image_filter_fn,_Extra=extra
;
;IF N_Elements(show_grid) EQ 0 THEN show_grid=1
;
;np=N_Params() 
;IF N_Elements(beam_base) GT 0 THEN np=3 ELSE np=np<2
;SWITCH np OF
;    0:restore,file_path_fhd+'_fhd_params.sav'
;    1:restore,file_path_fhd+'_obs.sav'
;    2:restore,file_path_fhd+'_fhd.sav'
;    3:
;ENDSWITCH
;
;basename=file_basename(file_path_fhd)
;dirpath=file_dirname(file_path_fhd)
;export_path=filepath(basename,root=dirpath,sub='export')
;export_dir=file_dirname(export_path)
;image_path=filepath(basename,root=dirpath,sub='images')
;image_dir=file_dirname(image_path)
;IF file_test(image_dir) EQ 0 THEN file_mkdir,image_dir
;IF file_test(export_dir) EQ 0 THEN file_mkdir,export_dir
;
;IF Keyword_Set(image_filter_fn) THEN BEGIN
;    dummy_img=Call_function(image_filter_fn,fltarr(2,2),name=filter_name)
;    IF Keyword_Set(filter_name) THEN filter_name='_'+filter_name ELSE filter_name=''
;ENDIF ELSE filter_name=''
;;filter_name=''
;
;;IF N_Elements(normalization_arr) GT 0 THEN normalization=Mean(normalization_arr)/2.
;npol=fhd.npol
;dimension=obs.dimension
;elements=obs.elements
;degpix=obs.degpix
;
;stats_radius=10. ;degrees
;pol_names=['xx','yy','xy','yx','I','Q','U','V']
;
;grid_spacing=10.
;offset_lat=5.;15. paper 10 memo
;offset_lon=5.;15. paper 10 memo
;reverse_image=0   ;1: reverse x axis, 2: y-axis, 3: reverse both x and y axes
;map_reverse=0;1 paper 3 memo
;label_spacing=1.
;
;instr_images_filtered=Ptrarr(npol,/allocate)
;instr_images=Ptrarr(npol,/allocate)
;instr_sources=Ptrarr(npol,/allocate)
;restored_beam_width=(!RaDeg/(obs.MAX_BASELINE/obs.KPIX)/obs.degpix)/(2.*Sqrt(2.*Alog(2.)))
;restored_beam_width=restored_beam_width>0.75
;FOR pol_i=0,npol-1 DO BEGIN
;;    *instr_images[pol_i]=dirty_image_generate(*residual_array[pol_i],image_filter_fn=image_filter_fn,degpix=degpix,$
;;        _Extra=extra)*weight_invert(*beam_base[pol_i])
;    *instr_images[pol_i]=dirty_image_generate(*image_uv_arr[pol_i]-*model_uv_holo[pol_i],$
;        degpix=degpix,image_filter_fn=image_filter_fn)*weight_invert(*beam_base[pol_i])
;    *instr_images_filtered[pol_i]=*instr_images[pol_i];-Median(*instr_images[pol_i],fhd.smooth_width,/even)
;    *instr_sources[pol_i]=source_image_generate(comp_arr,obs,pol_i=pol_i,resolution=16,$
;        dimension=dimension,width=restored_beam_width)
;ENDFOR
;
;stokes_images=stokes_cnv(instr_images,beam=beam_base)
;stokes_images_filtered=stokes_cnv(instr_images_filtered,beam=beam_base)
;stokes_sources=stokes_cnv(instr_sources,beam=beam_base)
;
;beam_mask=fltarr(dimension,elements)+1
;beam_avg=fltarr(dimension,elements)
;alias_mask=fltarr(dimension,elements) 
;alias_mask[dimension/4:3.*dimension/4.,elements/4:3.*elements/4.]=1
;FOR pol_i=0,(npol<2)-1 DO BEGIN
;    beam_mask_test=fltarr(dimension,elements)
;    beam_i=where(*beam_base[pol_i]*alias_mask GE fhd.beam_threshold/2.)
;;    beam_i=region_grow(*beam_base[pol_i],dimension/2.+dimension*elements/2.,threshold=[fhd.beam_threshold/2.,Max(*beam_base[pol_i])])
;    beam_avg+=*beam_base[pol_i]/(2.<npol)
;    beam_mask_test[beam_i]=1.
;    beam_mask*=beam_mask_test
;ENDFOR
;x_inc=where(beam_mask) mod dimension
;y_inc=Floor(where(beam_mask)/dimension)
;zoom_low=min(x_inc)<min(y_inc)
;zoom_high=max(x_inc)>max(y_inc)
;
;mkhdr,fits_header,(*stokes_images[0])[zoom_low:zoom_high,zoom_low:zoom_high]
;astr=obs.astr
;astr.crpix-=zoom_low
;putast, fits_header, astr
;beam_avg_use=beam_avg[zoom_low:zoom_high,zoom_low:zoom_high]
;beam_mask_use=beam_mask[zoom_low:zoom_high,zoom_low:zoom_high]
;
;FOR pol_i=0,npol-1 DO BEGIN
;    stokes_residual=((*stokes_images[pol_i])*beam_mask)[zoom_low:zoom_high,zoom_low:zoom_high]
;    stokes_residual_filtered=((*stokes_images_filtered[pol_i])*beam_mask)[zoom_low:zoom_high,zoom_low:zoom_high]
;    stokes_source=((*stokes_sources[pol_i])*beam_mask)[zoom_low:zoom_high,zoom_low:zoom_high]
;    stokes_restored=stokes_residual_filtered+stokes_source
;    
;    stokes_low=Min((stokes_residual_filtered*Sqrt(beam_avg_use>0))[where(beam_mask_use)])
;    stokes_high=Max((stokes_residual_filtered*Sqrt(beam_avg_use>0))[where(beam_mask_use)])
;    stokesS_high=Max(stokes_restored[where(beam_mask_use)])
;    IF pol_i EQ 0 THEN log=1 ELSE log=0
;    IF pol_i EQ 1 THEN Imagefast,stokes_residual_filtered,file_path=image_path+filter_name+'_Residual_'+pol_names[pol_i+4],$
;        /right,sig=2,color_table=0,back='white',reverse_image=reverse_image,low=stokes_low,high=stokes_high,$
;        lat_center=obs.obsdec,lon_center=obs.obsra,rotation=0,grid_spacing=grid_spacing,degpix=obs.degpix,$
;        offset_lat=offset_lat,offset_lon=offset_lon,label_spacing=label_spacing,map_reverse=map_reverse,show_grid=show_grid,/sphere,/no_ps
;    IF pol_i EQ 0 THEN Imagefast,stokes_restored,file_path=image_path+filter_name+'_Restored_'+pol_names[pol_i+4],$
;        /right,sig=2,color_table=0,back='white',reverse_image=reverse_image,log=log,low=stokes_low,high=stokesS_high,$
;        lat_center=obs.obsdec,lon_center=obs.obsra,rotation=0,grid_spacing=grid_spacing,degpix=obs.degpix,$
;        offset_lat=offset_lat,offset_lon=offset_lon,label_spacing=label_spacing,map_reverse=map_reverse,show_grid=show_grid,/sphere,/no_ps
;    FitsFast,stokes_residual,fits_header,/write,file_path=export_path+'_Residual_'+pol_names[pol_i+4]
;ENDFOR
;
;;write sources to a text file
;Ires=(*stokes_images_filtered[0])[source_array.x,source_array.y]
;IF npol GT 1 THEN Qres=(*stokes_images_filtered[1])[source_array.x,source_array.y]
;radius=angle_difference(obs.obsdec,obs.obsra,source_array.dec,source_array.ra,/degree)
;source_array_export,source_array,beam_avg,radius=radius,Ires=Ires,Qres=Qres,file_path=export_path+'_source_list'
;
;Ires=(*stokes_images_filtered[0])[comp_arr.x,comp_arr.y]
;IF npol GT 1 THEN Qres=(*stokes_images_filtered[1])[comp_arr.x,comp_arr.y]
;radius=angle_difference(obs.obsdec,obs.obsra,comp_arr.dec,comp_arr.ra,/degree)
;source_array_export,comp_arr,beam_avg,radius=radius,Ires=Ires,Qres=Qres,file_path=export_path+'_component_list'
;
;END