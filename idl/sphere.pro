;PROCEDURE SPHERE
;A procedure to render a sphere

;set color parameters for iMac
device,decomposed=0
loadct,1 ;blue-white color table

; Create an empty, 3D array:
SPHERE = FLTARR(20, 20, 20)

; Create the spherical dataset:
FOR X=0,19 DO FOR Y=0,19 DO FOR Z=0,19 DO $
   SPHERE(X, Y, Z) = SQRT((X-10)2 + (Y-10)2 + (Z-10)2)

   ; Find the vertices and polygons for a density level of 8:
   SHADE_VOLUME, SPHERE, 8, V, P

   ; Set up an appropriate 3D transformation so we can see the
   ; sphere. This step is very important:
   SCALE3, XRANGE=[0,20], YRANGE=[0,20], ZRANGE=[0,20]

   ; Render the image. Note that the T3D keyword has been set so that
   ; the previously-established 3D transformation is used:
   image = POLYSHADE(V, P, /T3D)

   ; Display the image:
   TV, image

   ;end sphere.pro
   end

