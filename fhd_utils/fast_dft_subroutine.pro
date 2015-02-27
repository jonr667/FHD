FUNCTION fast_dft_subroutine,x_vec,y_vec,amp_vec,dft_kernel_threshold=dft_kernel_threshold,dimension=dimension,$
    elements=elements,dft_approximation_resolution=dft_approximation_resolution,conserve_memory=conserve_memory,return_kernel=return_kernel

IF N_Elements(elements) EQ 0 THEN elements=dimension
IF N_Elements(dft_approximation_resolution) EQ 0 THEN resolution=32. ELSE resolution=Float(Round(dft_approximation_resolution))
IF N_Elements(dft_kernel_threshold) EQ 0 THEN dft_kernel_threshold=2./(!Pi*dimension) ;value of kernel_test along either axis at the edge of the image. 
IF resolution LE 1 THEN resolution=32.

xv_test=Abs(meshgrid(dimension,elements,1)-dimension/2.)
yv_test=Abs(meshgrid(dimension,elements,2)-elements/2.)

kernel_test=1./(((!Pi*xv_test)>1.)*((!Pi*yv_test)>1.)) 
kernel_i=where(kernel_test GE dft_kernel_threshold,n_k)

IF Keyword_Set(conserve_memory) THEN BEGIN
    IF conserve_memory GE 1E7 THEN mem_threshold=conserve_memory ELSE mem_threshold=1E8
    WHILE n_k*resolution^2. GT mem_threshold DO BEGIN
        resolution=Float(Round(resolution/Sqrt(2.)))
        dft_kernel_threshold*=Sqrt(2.)
        kernel_test=1./(((!Pi*xv_test)>1.)*((!Pi*yv_test)>1.)) 
        kernel_i=where(kernel_test GE dft_kernel_threshold,n_k)
    ENDWHILE
ENDIF

xv_k=(kernel_i mod dimension)-dimension/2.
yv_k=Floor(kernel_i/dimension)-elements/2.

kernel_recalc=1
IF Keyword_Set(return_kernel) THEN BEGIN
    IF Min(Ptr_valid(return_kernel)) EQ 1 THEN BEGIN
        kernel_arr=return_kernel
        kernel_free=0
        kernel_recalc=0
    ENDIF
ENDIF   
IF Keyword_Set(kernel_recalc) THEN BEGIN
    kernel_arr=Ptrarr(resolution,resolution)
;    kernel_norm=0.
    FOR i=0.,resolution-1 DO BEGIN
        IF i EQ 0 THEN BEGIN
            kernel_x=Sin(!DPi*(xv_k+i/resolution))/((!DPi*(xv_k+i/resolution))>1.)
            kernel_x[where(xv_k EQ 0)]=1.
        ENDIF ELSE kernel_x=Sin(!DPi*(xv_k+i/resolution))/(!DPi*(xv_k+i/resolution))
        FOR j=0.,resolution-1 DO BEGIN
            IF j EQ 0 THEN BEGIN
                kernel_y=Sin(!DPi*(yv_k+j/resolution))/((!DPi*(yv_k+j/resolution))>1.)
                kernel_y[where(yv_k EQ 0)]=1.
            ENDIF ELSE kernel_y=Sin(!DPi*(yv_k+j/resolution))/(!DPi*(yv_k+j/resolution))
            kernel_single=kernel_x*kernel_y
            kernel_norm=Total(kernel_single,/double)
            kernel_arr[i,j]=Ptr_new(Float(kernel_single/kernel_norm))
        ENDFOR
    ENDFOR
    
    IF Keyword_Set(return_kernel) THEN BEGIN
        return_kernel=kernel_arr
        kernel_free=0
    ENDIF ELSE kernel_free=1
ENDIF

x_offset=Round((Ceil(x_vec)-x_vec)*resolution) mod resolution    
y_offset=Round((Ceil(y_vec)-y_vec)*resolution) mod resolution
xcen0=Round(x_vec+x_offset/resolution) ;do this after offset, in case it has rounded to the next grid point
ycen0=Round(y_vec+y_offset/resolution)

si1=where((xcen0 GE 0) AND (ycen0 GE 0) AND (xcen0 LE dimension-1) AND (ycen0 LE elements-1),ns)

;test if any gridding kernels would extend beyond image boudaries
xv_test=Minmax(xcen0[si1])+Minmax(xv_k)
yv_test=Minmax(ycen0[si1])+Minmax(yv_k)

IF xv_test[0] LT 0 OR xv_test[1] GT dimension-1 OR yv_test[0] LT 0 OR yv_test[1] GT elements-1 THEN BEGIN
    mod_flag=1
    dimension_use=xv_test[1]-xv_test[0]
    elements_use=yv_test[1]-yv_test[0]
    xcen0-=xv_test[0]
    ycen0-=yv_test[0]
ENDIF ELSE BEGIN
    mod_flag=0 
    dimension_use=dimension
    elements_use=elements
ENDELSE

model_img_use=fltarr(dimension_use,elements_use)
FOR si=0L,ns-1L DO BEGIN
    model_img_use[xcen0[si1[si]]+xv_k,ycen0[si1[si]]+yv_k]+=amp_vec[si1[si]]*(*kernel_arr[x_offset[si1[si]],y_offset[si1[si]]])
ENDFOR

IF Keyword_Set(mod_flag) THEN BEGIN
    model_img=Fltarr(dimension,elements)
    x_low0=xv_test[0]>0
    y_low0=yv_test[0]>0
    x_high0=(xv_test[0]+dimension_use-1)<(dimension-1)
    y_high0=(yv_test[0]+elements_use-1)<(elements-1)
    x_low1=-xv_test[0]>0
    y_low1=-yv_test[0]>0
    x_high1=x_high0-x_low0+x_low1
    y_high1=y_high0-y_low0+y_low1
    model_img[x_low0:x_high0,y_low0:y_high0]=model_img_use[x_low1:x_high1,y_low1:y_high1]
    
    ;add in aliasing!
    IF x_low1 GT 0 THEN BEGIN
        model_img[dimension-x_low1:dimension-1,y_low0:y_high0]=model_img_use[0:x_low1-1,y_low1:y_high1]
    ENDIF
    IF y_low1 GT 0 THEN BEGIN
        model_img[x_low0:x_high0,elements-y_low1:elements-1]=model_img_use[x_low1:x_high1,0:y_low1-1]
    ENDIF
    IF x_high1 LT dimension_use-1 THEN BEGIN
        model_img[0:dimension_use-x_high1-2,y_low0:y_high0]=model_img_use[x_high1+1:dimension_use-1,y_low1:y_high1]
    ENDIF
    IF y_high1 LT elements_use-1 THEN BEGIN
        model_img[x_low0:x_high0,0:elements_use-y_high1-2]=model_img_use[x_low1:x_high1,y_high1+1:elements_use-1]
    ENDIF
ENDIF ELSE model_img=model_img_use
IF Keyword_Set(kernel_free) THEN Ptr_free,kernel_arr

RETURN,model_img
END