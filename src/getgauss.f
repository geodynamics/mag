	subroutine getgauss(l,m,alm,blm,glm,hlm,id)
	
c	this subroutine converts between Gauss coefficients (glm, hlm)
c      	and the spherical harmonic coefficients (alm,blm)
c	of the magnetic potential at harmonic degree l 
c	and order m used by the numerical dynamo model MAG. 
c 
c	Inputs:
c		l,m: integer spherical harmonic degree and order;
c
c		id: transformation direction:
c			id = 1 converts (alm,blm) to (glm,hlm)
c			id =-1 converts (glm,hlm) to (alm,blm)             
c
c		alm,blm: dimensionless spherical harmonic coefficients 
c			of the radial magnetic field potential on the
c			outer boundary of the dynamo model. alm is the
c			coefficient of the real (cosine) part and blm 
c			is for the imaginary (sine) part.
c
c		glm,hlm: Gauss coefficients in nanotessla	(nT).	
c                        glm multiplies with cos(m*phi)
c                        hlm multiplies with sin(m*phi)
c		
c	The calculation uses SI units and nominal Earth values for
c	surface and core radii. The 
c	conversion from dimensionless (alm,blm) to (glm,hlm) in
c	nT uses Earth values for density, rotation, electrical 
c	conductivity, and uses the depth of the outer core as the
c	length scale.

		if (l .le. 0 .or. m .lt. 0 .or. m .gt. l) then
			write(6,'(''bad l or m in getgauss'')')
			return
		endif	 
c
c	constants
c
		pi=4.*atan(1.)
		pi4inv=1./(4.*pi)
		anano=1.0e+09
c
c	nominal Earth parameters
c
		sigma=6.0e+05     ! Core electrical conductivity in S/m
		rho=11.e+03       ! Core mean density in kg m^{-3}
		omega=7.2924e-05  ! Rotation angular frequency
		escale=sqrt(rho*omega/sigma)  ! "Elsasser number" scaling	
		re=6.371e+06      ! Earth radius in m
		rc=3.485e+06      ! Core radius in m
		ri=1.215e+06      ! Inner core radius in m
		cdepth=rc-ri      ! Depth of outer core
		rcre=rc/re        ! core/surface radius ratio
		rcsqinv=(cdepth/rc)**2  ! dim-less inverse squared core radius
c     				
c	harmonic factors	
c
		fl=float(l)
		flp2=fl + 2.
		fl2p1=(2.*fl) + 1.
		fm=float(m)
		fact=sqrt(fl2p1*pi4inv)
		fact1=((-1.)**fm)*fl*fact
                if(m.gt.0) fact1=fact1*sqrt(2.)
		fact2=rcsqinv*(rcre**flp2)
c
c	define MAG harmonic-Gauss harmonic conversion factors 
c	conalm,conblm
c
		conalm=fact1
		conblm=-fact1
c
		if (id .gt. 0) then  ! form Gauss coeffs in nT	
			glm=anano*escale*fact2*conalm*alm
			hlm=anano*escale*fact2*conblm*blm
		        return
		else ! form dimensionless fully normalized potential coeffs
			alm=glm/(anano*escale*fact2*conalm)
			blm=hlm/(anano*escale*fact2*conblm)
		return
		endif
	end
