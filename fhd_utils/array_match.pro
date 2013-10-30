FUNCTION array_match,array1,array2,value_match=value_match
;IF array2 is supplied, it should have the SAME number of elements as array1

min_use=Min(array1)<Min(array2)
max_use=Max(array1)>Max(array2)

hist1=histogram(array1,binsize=1,min=min_use,max=max_use,reverse_ind=ri1)
hist2=histogram(array2,binsize=1,min=min_use,max=max_use,reverse_ind=ri2)

hist12=hist1+hist2
bins=where(hist12 GT 0, nb)

;select values to be used
hist_v1=histogram(bins+min_use,binsize=1,omin=omin,omax=omax)
hist_v2=histogram(value_match,binsize=1,min=omin,max=omax)
vals=where(hist_v1 AND hist_v2,nv)+omin-min_use

IF nv EQ 0 THEN RETURN,-1
ind_arr=intarr(size(array1,/dimension))
FOR vi=0L,nv-1 DO BEGIN
    i=vals[vi]
    IF hist1[i] GT 0 THEN ind_arr[ri1[ri1[i]]:ri1[ri1[i+1]-1]]=1
    IF hist2[i] GT 0 THEN ind_arr[ri2[ri2[i]]:ri2[ri2[i+1]-1]]=1
ENDFOR

match_indices=where(ind_arr)
RETURN,match_indices
END