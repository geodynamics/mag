
;PROCEDURE TUBE
;A procedure to render a tube along a trajectory
;For a single tube, a monochromatic color table is
;recommended, such as loadct=0(b&w) or 1(blue-white) 

;set color parameters for iMac
device,decomposed=0
loadct,1


;DEFINE COORDINATE AND VARIABLE ARRAYS
x=fltarr(10,10,10)
y=fltarr(10,10,10)
z=fltarr(10,10,10)
vx=fltarr(10,10,10)
vy=fltarr(10,10,10)
vz=fltarr(10,10,10)

;SPECIFY COORDINATE GRID
for k=0,9 do begin
x(k,*,*)=k & y(*,k,*)=k & z(*,*,k)=k
endfor

;SPECIFY FIELD VARIABLE - INCOMPRESSIBLE!

;case#1: spiral
vz=replicate(0.1,10,10,10)
vx=2.*sin(0.5*x)*cos(0.5*y)
vy=-2.*cos(0.5*x)*sin(0.5*y)

;LOAD FIELD VARIABLE INTO DATA ARRAYS
data=fltarr(3,10,10,10)
data(0,*,*,*)=vx
data(1,*,*,*)=vy
data(2,*,*,*)=vz

;SPECIFY THE STARTING POINT FOR THE TRAJECTORY
seeds=[1.,1.,0] 

;create the trajectory
particle_trace,data,seeds,verts,conn,normals, $
	max_stepsize=0.25, max_iterations=500	
print,'particle trace done'

;define a circular tube around the trajectory
tubrad=0.1
pi=3.1416
xp=fltarr(2,11)
np=findgen(11)/5 & np(10)=0 
xp(0,*)=tubrad*cos(np*pi) & xp(1,*)=tubrad*sin(np*pi)

;create the tube
streamline,verts,conn,normals,outverts,outconn,profile=xp
print,'streamline done'

;render the tube
SCALE3,XRANGE=[0,9],YRANGE=[0,9],ZRANGE=[0,9],Ax=15
set_shading,light=[0,0,1]
tvscl,polyshade(outverts,outconn,/t3d)

end
