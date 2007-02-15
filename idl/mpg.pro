PRO  make_mpg, prefix=prefix, suffix=suffix, n_start=n_start, n_end=n_end, digits=digits, $
	dims=dims, format=format, frame_rate=frame_rate, mpeg_file=mpeg_file, tmp_dir=tmp_dir
;-----------------------------------------------------------
; the image filenames are assumed of the format: image#.ext, 
; where # is the index of the sequence
; and ext is one of the suported image types. 
;-----------------------------------------------------------

;------- Check the arguments 
if keyword_set(mpeg_file) eq 0  then begin
	mpeg_file='outfile.mpg' 
end 


if keyword_set(prefix) eq 0  then begin
	print, 'prefix is  missing'
  	return
end
if keyword_set(suffix) eq 0  then begin
  print, 'suffix is  missing'
  return
end

if keyword_set(n_start) eq 0 then begin
	n_start=0
	print, 'n_start=', n_start
end

if keyword_set(n_end)  eq 0 then begin
   n_end=0
	print, 'n_end=', n_end
end

if keyword_set(format) eq 0 then begin
	format = 0
end
if keyword_set(frame_rate) eq 0 then begin
	frame_rate = 5
end

if keyword_set(tmp_dir)  eq 0  then begin
       	tmp_dir='.'
end

if (n_start > n_end) then begin
	print, 'n_start, and n_end do not make sense'
  return
end
;------- Create the MPEG
if keyword_set(dims)  ne 1  then begin
	mympeg = obj_new('IDLgrMPEG', filename = mpeg_file, format=format, frame_rate=frame_rate, temp_directory=tmp_dir)
endif else begin
	mympeg = obj_new('IDLgrMPEG', filename = mpeg_file, dimensions=dims, format=format, frame_rate=frame_rate, temp_directory=tmp_dir)
endelse
;------- Read the images
for j = n_start, n_end do begin
  n = string(digits)
  format_str = '(I' + n+ '.' +n +')'
  index = string(format=format_str, j)
  image_name = prefix +strcompress(index, /remove_all) + '.' + suffix
  image = READ_IMAGE (image_name)
  image_size = SIZE(image)

  if (image_size[0] eq  0) then begin
    print, 'Cannot read image file: ', image_name
    return
  end
;------- Add the image to the sequence 
mympeg -> Put, reverse(image,3) , j
print, '.'
endfor

;------- Generate the Mpeg
mympeg -> Save
obj_destroy, mympeg
print, 'done'
return
END


