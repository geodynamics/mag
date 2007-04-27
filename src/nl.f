      program nl
c
c     A dynamic dynamo model driven by thermal convection
c     in a rotating spherical fluid shell.
c     This version is restricted to Boussinesq fluids and
c     non-dimensional variables are used throughout.
c
c     The following set of equations is solved:
c
c     E {dv/dt + v.grad(v)} = -grad(p) - 2e_z x v
c         +1/Pm rot(B) x B + RaE/Pr g/g_o T
c
c     dB/dt = rot(v x B) + 1/Pm Lapl(B)
c
c     dT/dt + v.grad(T) = 1/Pr Lapl(T) + epsc0
c
c     div(v)=0          div(B)=0
c
c       subject to the following boundary conditions
c       at the inner and outer radii:
c
c       v_r=0, and either no slip or stress free
c       T=0 / T=1  or fixed heat flux (the latter not tested!)
c       B fitted to exterior potential fields, or parts of B
c       specified on the boundaries
c
c     List of symbols:
c     
c     v: velocity          p: pressure        B: magnetic field
c     g: gravity           g_o: reference value at outer radius
c     T: temperature       epsc0: rate of internal heating
c     e_z: unit vector parallel to the rotation axis
c     d/dt: partial time derivative  Lapl: Laplace operator
c     
c     Scaling properties:
c
c     nu: kinematic viscosity         d: shell width
c     omega: angular frequency        alpha: thermal expansion coeff
c     delta_T: temperature contrast   kappa: thermal diffusivity
c     eta: magnetic diffusivity       rho: density
c     mu_o: magnetic permeability
c
c     Scaling:
c
c     Length:   d              time:      d^2/nu
c     Velocity: nu/d           pressure:  rho*nu*omega
c     Temperature: delta_T     mag.field: sqrt(rho*mu_o*eta*omega)
c     
c
c     Non-dimensional numbers:
c
c     E: Ekman number     E= nu*d^2/omega
c     Ra: Rayleigh number Ra = alpha*g_o*delta_T*d^3/(kappa*nu)
c     Pr: Prandtl number  Pr = nu/kappa
c     Pm: Magnetic Prandtl number    Pm=nu/eta      
c
c
c
c
c     Numerical simulations via a nonlinear, multimode,
c     initial-boundary value problem.
c
c *** entropy boundary condtions (tops and bots on input)
c     if ktops = 1, entropy specified on outer boundary
c     if ktops = 2, radial heat flux specified on outer boundary
c     if kbots = 1, entropy specified on inner boundary
c     if kbots = 2, radial heat flux specified on inner boundary
c     for example: ktops=1,
c           the spherically-symmetric temperature
c           on the outer boundary relative to the reference state
c
c *** velocity boundary condtions
c     if ktopv = 1, stress-free outer boundary
c     if ktopv = 2, non-slip outer boundary
c     if kbotv = 1, stress-free inner boundary
c     if kbotv = 2, non-slip inner boundary
c
c *** magnetic boundary condtions
c     if ktopb = 1, insulating outer boundary (mag coupling if cmb.gt.0)
c     if kbotb = 1, perfectly insulating inner boundary
c     if kbotb = 2, perfectly conducting inner boundary
c
c *** magneto-convection
c     bpeak = max amplitude of imposed magnetic field
c     if imagcon .eq. 1, imposed toroidal field via inner bc on J(l=2,m=0)
c     if imagcon .eq.10, imposed tor. field on both icb and cmb J(l=2,m=0)
c     if imagcon .eq.11, imposed tor. field on both icb and cmb J(l=2,m=0)
c                        opposite sign
c     if imagcon .eq.12, imposed tor. field on both icb and cmb J(l=1,m=0)
c     if imagcon .lt. 0, imposed poloidal field via inner bc on B(l=1,m=0)
c
c
c     if init .eq. 0, initial conditions are read from "in".
c     if init .gt. 0, random initial entropy (and magnetic) conditions.
c     if init .lt. 0, initial hydro conditions are read from "in"
c                     with random initial magnetic conditions.
c
c     since nj .ge. (3*mmax+1)
c       and ni .ge. (3*mmax+1)/2 for triangular truncation,
c       horizontal transforms are alias free.
c     if nnaf .lt. nn, aliasing in radial transform is reduced.
c
c     if symmetry is forced in longitude (minc .gt. 1)
c        (i.e., longitudinal periodicity of order minc)
c        then jc = 1 to nja=nj/minc.
c
c     nstep  = number of timesteps per printout (even)
c     nprnt  = number of printouts per data storage
c     nstor  = number of data storages per run
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com3.f'
      include 'com8.f'
c
      character chd*3,chg*3
      character*72 movefile,movmfile,movafile,logfile,grafile,
     $ lsfile,lpfile,ccfile, cgfile
c---------------------------------------------------------------
c
      ifirst=1
      kc0=0
      ispc1=1
c     open(5,file='par',form='formatted',status='old')
c     open(6,file='po',form='formatted',status='new')
c
      call prep
c
      if(nstep .lt. 1) then
         stop '00'
      endif
c
c  open various output files
c
      logfile='l.'//outfile
      lsfile='ls.'//outfile
      movefile='me.'//outfile
      movafile='ma.'//outfile
      movmfile='mm.'//outfile
      ccfile='cc.'//outfile
      cgfile='cg.'//outfile
      lpfile='lp.'//outfile
      open(15,file=logfile,status='new',form='formatted')
      open(16,file=lsfile,status='new',form='formatted')
      if(nplog.gt.0) 
     $ open(17,file=lpfile,status='new',form='formatted')
      if(mod(imovopt,10).ge.1)
     $ open(18,file=movefile,status='new',form='formatted')
      if(mod(imovopt,100).ge.10)
     $ open(19,file=movafile,status='new',form='formatted')
      if(mod(imovopt,1000).ge.100)
     $ open(20,file=movmfile,status='new',form='formatted')
      if(imovopt.ge.1000)
     $ open(21,file=ccfile,status='new',form='formatted')
       open(22,file=cgfile,status='new',form='formatted')
c
c  write header for movie file
c
      if(imovopt.lt.1) iframes=0
      imovct=1
      tmovnext=tmovstart
      kvp=mod(imovopt,1000)/100
      if(mod(imovopt,10).ge.1) call moveout(kc0)
      if(mod(imovopt,100).ge.10) call movaout(kc0)
      if(mod(imovopt,1000).ge.100) call movmout(kc0,kvp)
c
c ***  start of iteration loop
c
      do 1 istor=1,nstor
c
         if(istor .gt. 1) then
            close(14)
         endif
         if(nstor.gt.1) then
           write(chd,'(''d'',i1,''.'')') istor-1
           write(chg,'(''g'',i1,''.'')') istor-1
           rstfile=chd//outfile
           grafile=chg//outfile
         else
           rstfile='d.'//outfile
           grafile='g.'//outfile
         endif
         open(10,file=rstfile,status='new',form='unformatted')
         if(abs(ngform).eq.1)
     $    open(14,file=grafile,status='new',form='formatted')
         if(ngform.ge.2)
     $    open(14,file=grafile,status='new',form='unformatted')
c
         do 2 iprnt=1,nprnt
            call amhd
            if(istop.gt.0) go to 1
    2    continue
c
    1 continue
c
      if(istop.gt.0)
     &   write(6,'(/'' Terminated by stop signal'')')
      stop 'Regular end of program mag'
      end
