      subroutine prep
c
c     nj must be multiple of 4 
c     ni should be nj/2
c     nn must be of form 4*i+1, where i is an integer
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com1.f'
      include 'com2.f'
      include 'com3.f'
      include 'com4.f'
      include 'com5.f'
      include 'com6.f'
      include 'com7.f'
      include 'com8.f'
c
      dimension colat(ni),clm(nlafp1,nmafa),gauss(ni),
     $bleg1(nlafp1),bleg2(nlafp1),bleg3(nlaf),
     $rv1(nnp2),rv2(nnp2),rv3(nnp2)
c
      character*72 infile
c
      logical treset
c
      namelist /contrl/ outfile,infile,runid,
     $init,nstep,nprnt,nstor,iframes,tmovstart,tmovstep,imovopt,
     $alpha,courfac,alffac,tipdipole,ra,pr,prmag,ek,radratio,cmb,
     $epsc0,dtmax,dtstart,bpeak,difamp,ldif,ldifexp,
     $icour,samp,nlogstep,treset,iscale,enscale,
     $amps,ampb,ampj,ampw,ampz,imagcon,
     $ktops,kbots,ktopv,kbotv,ktopb,kbotb,ifvfrz,ifbfrz,ifsfrz,
     $ngform,ngcolat,ngrad,nglon,nplog,nfilt,alfilt,ivfilt,dipfilt
c
      namelist /bounds/ tops,bots
c
      runid="Dummy run name for magb"
      outfile="dummy"
      infile="in"
      init=0
      alpha=0.6
      courfac=2.5
      alffac=1.0
      dtmax=1.e-3
      radratio=0.35
      ra=1.e5
      ek=1.e-3
      pr=1.
      prmag=1.
      epsc0=0.0
      cmb=0.0 
      bpeak=0.
      difamp=0.
      ldif=1
      ldifexp=-1
      icour=2
      samp=1.
      amps=1.
      ampb=1.
      ampj=1.
      ampw=1.
      ampz=1.
      imagcon=0
      ktops=1
      kbots=1
      ktopv=2
      kbotv=2
      ktopb=1
      kbotb=1
      ifvfrz=.false.
      ifbfrz=.false.
      ifsfrz=.false.
      ngform=0
      ngrad=1
      ngcolat=1
      nglon=1
      imovopt=0
      iframes=0
      tmovstart=0.
      tmovstep=0.5e-3
      treset=.false.
      nlogstep=25
      nplog=0
      nfilt=0
      alfilt=9999.
      ivfilt=-1
      dipfilt=1.0
      iscale=1
      enscale=1.
      dtstart=0.0
      tipdipole=0.0
c
      i01=1
c---------------------------------------------------------------
c
      write(6,*) 'magb'
c
      if(minc.lt.1) stop '04'
      if((minc.eq.1) .and. (nlma.ne.nlm)) stop '05'
      if(mod(nj,minc) .ne. 0) stop '06'
      if(mod(nja,2) .ne. 0) stop '07'
      njd2=nja/2
      if((mod(njd2,2) .ne. 0) .or. (njd2 .le. 2)) stop '08'
      if((mod(nn1,4) .ne. 0) .or. (nn .le. 3)) stop '09'
      if(ni .le. 2) stop '09'
      if(ni .lt. nj/2) stop '10'
      if((nnaf .gt. nn) .or. (nnaf .lt. 1)) stop '11'
c     note if (nja+1)*ni .gt. lot*nnp2, then dimension work(lot,nwk)
c     where nwk=(nja+1)*ni/lot+1
c     since work is used by both fourtf and chebtf
      if((nja+1)*ni .gt. lot*nnp2) stop '12'
c
      write(6,267)
  267 format(/,2x,"triangular truncation of spherical harmonics")
      nmax=nnaf-1
      write(6,266) nn,ni,nj,nmax,lmax,mmax
  266 format(/,2x,"  nn =",i4,"     ni =",i4,"     nj =",i4,
     $         2x,"nmax =",i4,"   lmax =",i4,"   mmax =",i4,/)
      if(minc.gt.1) write(6,268) minc
  268 format(2x,i3," - fold symmetry in longitude",/)
c
c *** control parameters
c
      do m=0,mmax
         do l=0,lmax
            tops(l,m)=0.
            bots(l,m)=0.
         enddo
      enddo
c
      read(5,contrl)
      read(5,bounds)
c
      opr=1./pr
      opm=1./prmag
      oek=1./ek
      oekpm=oek*opm
      rapr=ra/pr
      dtmax=min(dtmax,0.25*ek)
      dtmin=dtmax/10000.
c
      do l=0,lmax
         tops(l,0)=real(tops(l,0))
         bots(l,0)=real(bots(l,0))
      enddo
c **************************************************
c     for this magneto-convection version:
c     the bpeak bc on aj(3,..), aj(2,..) or b(2,nn) requires
      kbotb=1
c **************************************************
      ldif=max(1,ldif)
      ldif=min(lmax,ldif)
      if(mod(icour,2) .ne. 0) icour=icour+1
      write(6,contrl)
      if(mod(nstep,2) .eq. 0) nstep=nstep-1
      do m=0,mmax
         do l=m,lmax
            if(tops(l,m) .ne. cmplx(0.,0.))
     $         write(6,901) m,l,tops(l,m)
  901       format(" m =",i3,"  l =",i3,
     $         "  tops =",2(1x,1pe12.5))
            if(bots(l,m) .ne. cmplx(0.,0.))
     $         write(6,902) m,l,bots(l,m)
  902       format(" m =",i3,"  l =",i3,
     $         "  bots =",2(1x,1pe12.5))
         enddo
      enddo
c
c *** parameters
c
      ai=cmplx(0.,1.)
      pi=4.*atan(1.)
      pik=pi/float(nn1)
      sqrt2pi=sqrt(2.*pi)
      y00=1./sqrt(4.*pi)
      p00co=4./sqrt(3.)
      anorm=sqrt(2./float(nn1))
c
c  radial grid structure
c
      radtop=1./(1.-radratio)
      radbot=radtop-1.
      ocorevol =4.*pi/3. *(radtop**3 - radbot**3)
c
      do kc=1,nn
         r(kc)=radbot+0.5*(1.+cos(float(kc-1)*pik))
      enddo
c
c *** input data
c
c -----------------------------------------------------------------
      if(init .le. 0) then    ! initial condition from restart file
c -----------------------------------------------------------------
         open(8,file=infile,status='old',form='unformatted')
         read(8) time,dt,raold,prold,pmold,ekold,radratiold,
     $           kstep,nnold,niold,njold,minco
         write(6,260) infile
  260    format(/,2x,"input file: ",a64,/)
c
         if(ra.ne.raold.or.pr.ne.prold.or.prmag.ne.pmold.or.
     $      ek.ne.ekold.or.radratio.ne.radratiold) then
          write(6,262)  
  262     format("Parameter values changed  New / Old : ")
          write(6,263) radratio,ra,pr,prmag,ek
          write(6,263) radratiold,raold,prold,pmold,ekold
  263     format('radratio= ',f8.5,'  ra= ',f9.0,'  pr= ',f7.3,
     $    '  prmag= ',f7.3,'  ek= ',e12.4)
         endif
c
         if(dtstart.gt.0.0) dt=dtstart
         if(ni.eq.niold.and.nj.eq.njold.and.minc.eq.minco
     $     .and.nn.eq.nnold) then
          read(8) w,z,p,s
          read(8) dsdt1,dwdt1,dzdt1,dpdt1
          if(init .gt. -10) read(8) b,aj,dbdt1,djdt1
         else
c
c  mapping from  different grid
c
ccc
             write(6,'(4i4)') nnold,niold,njold,minco
          call mapdata(nnold,niold,njold,minco)
         endif
c
         if(init.lt.-6.and.ifvfrz) call zerorot(init)
c
         if(treset) then
           time=0.
           kstep=0
         endif
c -----------------------------------------------------------------
      else       !  initial condition from scratch
c -----------------------------------------------------------------
         kstep=0
         time=0.
         if(dtstart.gt.0.0) then
           dt=dtstart
         else
           dt=dtmax
         endif
c -----------------------------------------------------------------
      endif
c -----------------------------------------------------------------
c
      dtold=dt
      oodt=1./dt
c
c *** chebyshev polynomials and derivatives
c
      do 10 kc=1,nn
         do 10 nc=1,nn
            cheb(nc,kc)=cos(float((nc-1)*(kc-1))*pik)
   10 continue
c
      do 70 kc=1,nn
      dcheb(1,kc)=0.
      do 71 nc=2,nn
      if(mod(nc,2) .eq. 0) go to 72
      dcheb(nc,kc)=0.
      n1=2
      go to 73
   72 dcheb(nc,kc)=0.5
      n1=3
   73 n2=nc-1
      if(n2 .lt. n1) go to 69
      do 74 ncc=n1,n2,2
      dcheb(nc,kc)=dcheb(nc,kc)+cheb(ncc,kc)
   74 continue
   69 dcheb(nc,kc)=float(2*(nc-1))*dcheb(nc,kc)
   71 continue
      d2cheb(1,kc)=0.
      d2cheb(2,kc)=0.
      do 75 nc=3,nn
      if(mod(nc,2) .eq. 0) go to 76
      d2cheb(nc,kc)=0.5*float((nc-1)**2)
      n1=3
      go to 77
   76 d2cheb(nc,kc)=0.
      n1=2
   77 n2=nc-2
      if(n2 .lt. n1) go to 79
      do 78 ncc=n1,n2,2
      d2cheb(nc,kc)=d2cheb(nc,kc)+float((nc-ncc)*(nc+ncc-2))
     $*cheb(ncc,kc)
   78 continue
   79 d2cheb(nc,kc)=float(nc-1)*d2cheb(nc,kc)
   75 continue
   70 continue
      do 90 ncc=1,nn
         do 91 kc=1,nn
            rv1(kc)=d2cheb(ncc,kc)
   91    continue
         call rderiv(1.,anorm,rv1,rv2)
         do 92 kc=1,nn
            d3cheb(ncc,kc)=0.5*rv2(kc)
   92    continue
   90 continue
      do 94 kc=1,nn
         d3cheb(1,kc)=0.
         d3cheb(2,kc)=0.
         d3cheb(3,kc)=0.
   94 continue
c
      call chebi(nn,wsave,trigsc,ifaxc,k2k)
c
c *** complex fourier polynomials
c
      call fax(ifaxf,nja,3)
      k=ifaxf(1)
      if((k .lt. 1) .or. (ifaxf(k+1) .gt. 5)) stop '17'
ctest
      write(6,'(/'' Fourier factors= '',10I2/)')
     $ (ifaxf(kk),kk=2,k+1)
      call fftrig(trigsf,nja,3)
c
c *** radially dependent variables
c
      do kc=1,nn
         grav(kc)=r(kc)/radtop
         qk(kc,1)=1./r(kc)**2   
         qk(kc,2)=4./r(kc)
         qk(kc,3)=1./r(kc)
         qk(kc,4)=1./r(kc)**4
         qk(kc,5)=2./r(kc)**3
         qk(kc,6)=2./r(kc)
      enddo    
c
c *** scaling parameters
c
      zscl=dt*radtop**2
      wscl=zscl
      pscl=radtop**2
      sscl=dt
      oosscl=1./sscl
      bscl=dt*radtop**2
      ampnu=2.*oek*radtop*radtop/real(lmax*(lmax+1))
      ampnu=max(1.d0,ampnu)
c
      lm0=3
      if(imagcon .gt. 0) then
         if(imagcon.eq.12) lm0=2
         b0norm=4./float(lm0) * sqrt(pi/(2.*float(lm0)-1.))
         bpeakbot=b0norm*radbot/bscl*bpeak
         if(imagcon.ge.10) then
           bpeaktop=b0norm*radtop/bscl*bpeak
           if(imagcon.eq.11) bpeaktop=-bpeaktop
         else
           bpeaktop=.0
         endif
      else if(imagcon.lt.0) then
         bpeakbot=-sqrt(pi/3.)*radbot*radbot/bscl*bpeak
         bpeaktop=.0
      else
         bpeaktop=.0
         bpeakbot=.0
      endif
c
      if((ktops.eq.1) .and. (kbots.eq.1)) then
         tops(0,0)=-radbot*radbot/(radtop*radtop+radbot*radbot)/y00
         bots(0,0)= radtop*radtop/(radtop*radtop+radbot*radbot)/y00
      endif
c
      do l=0,lmax
         do m=0,min(l,mmax)
           tops(l,m)=tops(l,m)*oosscl
           bots(l,m)=bots(l,m)*oosscl
         enddo
      enddo
c
c *** courant lengths
c
      c1=1./float(lmax*(lmax+1))
      delxh2(1)=c1*r(1)**2
      delxh2(nn)=c1*r(nn)**2
      delxr(1)=r(1)-r(2)
      delxr(nn)=r(nn1)-r(nn)
      do 120 kc=2,nn1
      delxh2(kc)=c1*r(kc)**2
      delxr(kc)=min((r(kc-1)-r(kc)),(r(kc)-r(kc+1)))
  120 continue
c
c *** chebyshev integrals
c
      c1=0.0625
      c2=0.5*(r(nn)+0.5)
      c3=(r(nn)+0.5)**2+0.125
      do 18 nc=1,nnp2
      rv1(nc)=0.
   18 continue
      do 19 nc=1,nnp2,2
      n=nc-1
      rv1(nc)=2./(1.-float(n**2))
   19 continue
      do 20 nc=1,nn
      n=nc-1
      qn(nc,1)=(c1*rv1(iabs(n-2)+1)+c2*rv1(iabs(n-1)+1)+
     $c3*rv1(nc)+c2*rv1(nc+1)+c1*rv1(nc+2))
      qn(nc,3)=0.
      qn(nc,4)=0.
   20 continue
      do 28 nc=1,nn,2
      qn(nc,3)=2./(1.-float((nc-1)**2))
   28 continue
c
      do 22 ncc=1,nn
      do 23 kc=1,nn
      rv1(kc)=0.0
      rv3(kc)=cheb(ncc,kc)*r(kc)**2
   23 continue
      call chebtf(i01,ns2,i01,nn1,ns2,rv1,rv1,rv1,wsave,
     $work,work(2,1),work(2,1),work(2,1),trigsc,ifaxc,k2k)
      call chebtf(i01,ns2,i01,nn1,ns2,rv3,rv3,rv3,wsave,
     $work,work(2,1),work(2,1),work(2,1),trigsc,ifaxc,k2k)
      do 26 nc=2,nn3,2
      rv2(nc)=(rv1(nc-1)-rv1(nc+1))/float(2*(nc-1))
      rva(nc)=(rv3(nc-1)-rv3(nc+1))/float(2*(nc-1))
   26 continue
      qn(ncc,2)=rv1(nn)/float(4*nn)
      qn(ncc,5)=rv3(nn)/float(4*nn)
      qn(ncc,6)=0.0                  
      rv2(nn1)=(rv1(nn2)-0.5*rv1(nn))/float(2*nn2)
      rva(nn1)=(rv3(nn2)-0.5*rv3(nn))/float(2*nn2)
      do 27 ncb=1,nn2,2
      nc=nn-ncb
      qn(ncc,2)=qn(ncc,2)+rv2(nc)
      qn(ncc,5)=qn(ncc,5)+rva(nc)
   27 continue
      qn(ncc,2)=2.*anorm*qn(ncc,2)
      qn(ncc,5)=2.*anorm*qn(ncc,5)
   22 continue
c
c *** horizontally dependent coefficients
c
      call gquad(ni,colat,gauss)
c
      do 32 ic=1,ni
         qi(ic,1)=1./(sin(colat(ic)))**2
         qi(ic,2)=cos(colat(ic))*qi(ic,1)
         qi(ic,3)=sin(colat(ic))
         qi(ic,4)=cos(colat(ic))
         qi(ic,5)=    colat(ic) 
   32 continue
c
c
      do 33 mc=1,nmaf,minc
      m=mc-1
      mca=m/minc+1
      do 34 lc=mc,nlafp1
      l=lc-1
      clm(lc,mca)=sqrt(float((l+m)*(l-m))/float((2*l-1)*(2*l+1)))
   34 continue
   33 continue
c
      lm=0
      lmp=0
      do 35 mc=1,nmaf,minc
      m=mc-1
      mca=m/minc+1
      do 31 lc=mc,nlaf
      l=lc-1
      lm=lm+1
      lmp=lmp+1
      mclm(lm)=mc
      mclma(lm)=mca
      mcalmp(lmp)=mca
      ql(lm,1)=float(l+2)*clm(lc+1,mca)
         if(lc .eq. nlaf) ql(lm,1)=0.
      ql(lm,2)=float(l-1)*clm(lc,mca)
      ql(lm,3)=float(l*lc)
      ql(lm,4)=float(l)
      ql(lm,5)=float(lc)
      ql(lm,6)=ql(lm,3)
         if(mc .eq. 1) ql(lm,6)=0.5*ql(lm,6)
      ql(lm,7)=float(lc)*clm(lc,mca)
      ql(lm,8)=float(l)*clm(lc+1,mca)
c        note: if(lc .eq. nlaf) ql(lm,8) .ne. 0, use only for nl terms
      ql(lm,9)=1.
         if(l .eq. lmax) ql(lm,9)=0.
      ql(lm,10)=float(m)
      if((ldifexp .gt. 0) .and. (l .ge. ldif)) then
         ql(lm,11)=ql(lm,3)*(1.+
     $      difamp*(float(l+1-ldif)/float(lmax+1-ldif))**ldifexp)
         ql(lm,12)=ql(lm,3)*(1.+
     $      difamp*(float(l+1-ldif)/float(lmax+1-ldif))**ldifexp)
      else
         ql(lm,11)=ql(lm,3)
         ql(lm,12)=ql(lm,3)
      endif
         if(l .eq. lmax) ql(lm,12)=ampnu*ql(lm,12)
      ql(lm,13)=ql(lm,10)*ql(lm,9)
      ql(lm,14)=ql(lm,1)*ql(lm,9)
      ql(lm,15)=ql(lm,2)*ql(lm,9)
      ql(lm,16)=ql(lm,4)*ql(lm,1)*ql(lm,9)
      ql(lm,17)=ql(lm,5)*ql(lm,9)
      ql(lm,18)=ql(lm,5)*ql(lm,2)*ql(lm,9)
      ql(lm,19)=ql(lm,4)*ql(lm,9)
      ql(lm,20)=ql(lm,3)*ql(lm,9)
   31 continue
      lmp=lmp+1
      mcalmp(lmp)=mca
   35 continue
      if(lm .ne. nlma) then
         write(6,*) lm,' .ne. nlma .eq. ',nlma
         stop '18'
      endif
      if(lmp .ne. nlmpa) then
         write(6,*) lmp,' .ne. nlmpa .eq. ',nlmpa
         stop '19'
      endif
c
c *** legendre functions
c
      do 36 ic=1,ni
      lm=0
      lmp=0
      do 36 mc=1,nmaf,minc
         m=mc-1
         mca=m/minc+1
         xm=(-1.)**m
         do 37 lc=mc,nlafp1
            l=lc-1
            call pbar(colat(ic),l,m,plm)
            bleg1(lc)=xm/sqrt2pi*plm
            bleg2(lc)=xm*sqrt2pi*gauss(ic)*plm
   37    continue
         bleg3(mc)=float(m)*clm(mc+1,mca)*bleg1(mc+1)
         if(mc .lt. nmaf) then
            do 39 lc=mc+1,nlaf
               l=lc-1
               bleg3(lc)=float(l)*clm(lc+1,mca)*bleg1(lc+1)-
     $         float(lc)*clm(lc,mca)*bleg1(lc-1)
   39       continue
         endif
         do 52 lc=mc,nlaf
            lm=lm+1
            lmp=lmp+1
            aleg1(lm,ic)=bleg1(lc)
            aleg2(lmp,ic)=bleg2(lc)
            aleg3(lm,ic)=bleg3(lc)
   52    continue
         lmp=lmp+1
         aleg2(lmp,ic)=bleg2(nlaf+1)
   36 continue
c
c *** initial conditions
c
      do 644 kc=1,nnp1
         w(1,kc)=0.
         dw(1,kc)=0.
         ddw(1,kc)=0.
         z(1,kc)=0.
         dz(1,kc)=0.
         b(1,kc)=0.
         db(1,kc)=0.
         ddb(1,kc)=0.
         aj(1,kc)=0.
         dj(1,kc)=0.
  644 continue
      do 645 kc=1,nn
         dpdt(1,kc,1)=0.
         dpdt(1,kc,2)=0.
         dzdt(1,kc,1)=0.
         dzdt(1,kc,2)=0.
         dbdt(1,kc,1)=0.
         dbdt(1,kc,2)=0.
         djdt(1,kc,1)=0.
         djdt(1,kc,2)=0.
  645 continue
c
c    -radial derivatives
c
      if(init .le. 0) then
         if(amps .ne. 1.) then
            do kc=1,nn
               do lm=1,nlma
                  s(lm,kc)=amps*s(lm,kc)
                  dsdt(lm,kc,1)=amps*dsdt(lm,kc,1)
               enddo
            enddo
         endif
         if((ampw .ne. 1.) .or. (ampz .ne. 1.)) then
            do kc=1,nn
               do lm=1,nlma
                  w(lm,kc)=ampw*w(lm,kc)
                  dwdt(lm,kc,1)=ampw*dwdt(lm,kc,1)
                  z(lm,kc)=ampz*z(lm,kc)
                  dzdt(lm,kc,1)=ampz*dzdt(lm,kc,1)
               enddo
            enddo
         endif
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   w,w,w,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   z,z,z,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         do lm=1,nlma
            dw(lm,nn)=0.
            dw(lm,nn1)=float(nn1)*w(lm,nn)
            ddw(lm,nn)=0.
            ddw(lm,nn1)=0.
            dz(lm,nn)=0.
            dz(lm,nn1)=float(nn1)*z(lm,nn)
         enddo
         do n=nn2,1,-1
            do lm=1,nlma
               dw(lm,n)=dw(lm,n+2)+float(2*n)*w(lm,n+1)
               ddw(lm,n)=ddw(lm,n+2)+float(2*n)*dw(lm,n+1)
               dz(lm,n)=dz(lm,n+2)+float(2*n)*z(lm,n+1)
            enddo
         enddo
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   w,w,w,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   dw,dw,dw,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   ddw,ddw,ddw,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   z,z,z,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   dz,dz,dz,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         do kc=1,nnp1
            w(1,kc)=0.
            dw(1,kc)=0.
            ddw(1,kc)=0.
            z(1,kc)=0.
            dz(1,kc)=0.
         enddo
      endif
c
      if(init .gt. -10) then
         if((ampb .ne. 1.) .or. (ampj .ne. 1.)) then
            do kc=1,nn
               do lm=1,nlma
                  b(lm,kc)=ampb*b(lm,kc)
                  dbdt(lm,kc,1)=ampb*dbdt(lm,kc,1)
                  aj(lm,kc)=ampj*aj(lm,kc)
                  djdt(lm,kc,1)=ampj*djdt(lm,kc,1)
               enddo
            enddo
         endif
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   b,b,b,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   aj,aj,aj,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         do lm=1,nlma
            db(lm,nn)=0.
            db(lm,nn1)=float(nn1)*b(lm,nn)
            ddb(lm,nn)=0.
            ddb(lm,nn1)=0.
            dj(lm,nn)=0.
            dj(lm,nn1)=float(nn1)*aj(lm,nn)
         enddo
         do n=nn2,1,-1
            do lm=1,nlma
               db(lm,n)=db(lm,n+2)+float(2*n)*b(lm,n+1)
               ddb(lm,n)=ddb(lm,n+2)+float(2*n)*db(lm,n+1)
               dj(lm,n)=dj(lm,n+2)+float(2*n)*aj(lm,n+1)
            enddo
         enddo
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   b,b,b,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   db,db,db,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   ddb,ddb,ddb,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   aj,aj,aj,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         call chebtf(lot,ns2,lot,nnp1,nps2,
     $   dj,dj,dj,wsave,work,
     $   work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
         do kc=1,nnp1
            b(1,kc)=0.
            db(1,kc)=0.
            ddb(1,kc)=0.
            aj(1,kc)=0.
            dj(1,kc)=0.
         enddo
      endif
c
c    -initial entropy field and boundary conditions
c
      do 40 nc=1,nn
         do 40 kc=1,nn
            s0mat(kc,nc)=sscl*anorm*opr*
     $         (4.*d2cheb(nc,kc)+qk(kc,2)*dcheb(nc,kc))
   40 continue
      do 41 nc=1,nnaf
         if(ktops .eq. 1) then
            s0mat(1,nc)=anorm
         else
            s0mat(1,nc)=2.*dcheb(nc,1)*anorm
         endif
         if(kbots .eq. 1) then
            s0mat(nn,nc)=cheb(nc,nn)*anorm
         else
            s0mat(nn,nc)=2.*dcheb(nc,nn)*anorm
         endif
   41 continue
      if(nnaf .lt. nn) then
         do 911 nc=nnaf+1,nn
            s0mat(1,nc)=0.
            s0mat(nn,nc)=0.
  911    continue
      endif
      do 42 kc=1,nn
         s0mat(kc,1)=0.5*s0mat(kc,1)
         s0mat(kc,nn)=0.5*s0mat(kc,nn)
   42 continue
      call sgefa(s0mat,nn,nn,is0,info)
      if(info .ne. 0) stop '20'
      do 44 kc=2,nn1
         rva(kc)=-epsc0
   44 continue
      rva(1)=real(tops(0,0))
      rva(nn)=real(bots(0,0))
      call sgesl(s0mat,nn,nn,is0,rva,0)
      do nc=1,nnaf
         rva(nc)=rva(nc)*sscl
      enddo
      if(nnaf .lt. nn) then
         do 912 nc=nnaf+1,nn
            rva(nc)=0.
  912    continue
      endif
      call chebtf(1,ns2,1,nn1,ns2,rva,rva,rva,wsave,
     $work,work(2,1),work(2,1),work(2,1),trigsc,ifaxc,k2k)
c
      do 58 lm=1,nlma
         s(lm,nnp1)=0.
         w(lm,nnp1)=0.
         z(lm,nnp1)=0.
         p(lm,nnp1)=0.
         dw(lm,nnp1)=0.
         ddw(lm,nnp1)=0.
         dz(lm,nnp1)=0.
         b(lm,nnp1)=0.
         aj(lm,nnp1)=0.
         db(lm,nnp1)=0.
         ddb(lm,nnp1)=0.
         dj(lm,nnp1)=0.
   58 continue
c
      if(init .gt. 0 .or. init .le. -10) then
c
c    -initial toroidal magnetic field
c
       if(imagcon .ge. 0) then
       do nc=1,nn
         do kc=1,nn
            s0mat(kc,nc)=ql(lm0,3)*opm*qk(kc,1)*
     $         (4.*d2cheb(nc,kc)-
     $         ql(lm0,3)*qk(kc,1)*cheb(nc,kc))*anorm*bscl
         enddo
       enddo
       do nc=1,nnaf
         s0mat(1,nc)=anorm
         s0mat(nn,nc)=cheb(nc,nn)*anorm
       enddo
       if(nnaf .lt. nn) then
         do nc=nnaf+1,nn
            s0mat(1,nc)=0.
            s0mat(nn,nc)=0.
         enddo
       endif
       do kc=1,nn
         s0mat(kc,1)=0.5*s0mat(kc,1)
         s0mat(kc,nn)=0.5*s0mat(kc,nn)
       enddo
       call sgefa(s0mat,nn,nn,is0,info)
       if(info .ne. 0) stop '21'
       do kc=1,nn1
         rvb(kc)=0.
       enddo
       rvb(nn)=bpeakbot                   ! Inner boundary
       rvb(1)= bpeaktop                   ! Outer boundary
       call sgesl(s0mat,nn,nn,is0,rvb,0)
       do nc=1,nnaf
         rvb(nc)=rvb(nc)*bscl
       enddo
       if(nnaf .lt. nn) then
         do nc=nnaf+1,nn
            rvb(nc)=0.
         enddo
       endif
       call chebtf(1,ns2,1,nn1,ns2,rvb,rvb,rvb,wsave,
     $ work,work(2,1),work(2,1),work(2,1),trigsc,ifaxc,k2k)
       endif
c
c  initial entropy and velocity
c
       if(init .gt. 0) then
         do 400 kc=1,nn
            do 402 lm=1,nlma
               s(lm,kc)=0.
               w(lm,kc)=0.
               z(lm,kc)=0.
               p(lm,kc)=0.
               dw(lm,kc)=0.
               ddw(lm,kc)=0.
               dz(lm,kc)=0.
               dwdt(lm,kc,1)=0.
               dwdt(lm,kc,2)=0.
               dzdt(lm,kc,1)=0.
               dzdt(lm,kc,2)=0.
               dpdt(lm,kc,1)=0.
               dpdt(lm,kc,2)=0.
               dsdt(lm,kc,1)=0.
               dsdt(lm,kc,2)=0.
  402       continue
  400    continue
       endif
c
       do 410 kc=1,nn
         do 412 lm=1,nlma
            b(lm,kc)=0.
            aj(lm,kc)=0.
            db(lm,kc)=0.
            ddb(lm,kc)=0.
            dj(lm,kc)=0.
            dbdt(lm,kc,1)=0.
            dbdt(lm,kc,2)=0.
            djdt(lm,kc,1)=0.
            djdt(lm,kc,2)=0.
  412    continue
  410  continue
c
c radial dependence of initial perturbation in rv1
c
       do 49 kc=1,nn
         x=2.*r(kc)-r(1)-r(nn)
         rv1(kc)=1.-3.*x**2+3.*x**4-x**6
   49  continue
       qllm4max=float(lmax)+0.1
c
c  random noise initialization
c
       if(init .gt. 0 .and. init .lt. 100) then
         do 85 lm=2,nlma
            if(ql(lm,4) .gt. qllm4max) go to 85
            ra1=(-1.+2.*random(0.))*samp
            ra2=(-1.+2.*random(0.))*samp
            do 86 kc=1,nn
               c1=ra1*rv1(kc)
               c2=ra2*rv1(kc)
               if(mclm(lm) .gt. 1) then
                  s(lm,kc)=cmplx(c1,c2)
               else
                  s(lm,kc)=c1
               endif
   86       continue
   85    continue
         do kc=1,nn
c            s(1,kc)=0.
             s(1,kc)=rva(kc)
         enddo
       endif
c
c  initialize one mode specifically
c
       if(init .ge. 100) then
        m=mod(init,100)
        if(mod(m,minc).ne.0) stop 'm_init incompatible with minc'
        l=init/100
        if(l.gt.lmax.or.l.lt.m) stop 'l_init > lmax or < m_init'
        lm=m*(lmax+1)/minc-m*(m-minc)/(2*minc)+l-m+1
            do 88 kc=1,nn
               c1=rv1(kc)*samp
               c2=0.0
               if(mclm(lm) .gt. 1) then
                  s(lm,kc)=cmplx(c1,c2)
               else
                  s(lm,kc)=c1
               endif
   88       continue
            do kc=1,nn
c              s(1,kc)=0.
               s(1,kc)=rva(kc)
            enddo
         write(6,'(/'' Initialized at mode l= '',i3,''  m= '',i3/)')
     $   nint(ql(lm,4)),nint(ql(lm,10))
       endif
c
       if(imagcon .ge. 0) then
         if(bpeak.gt.0.0) then
         do kc=1,nn
            aj(lm0,kc)=rvb(kc)
         enddo
         else if(bpeak.lt.0.0) then
         bpk=2.*bpeak*sqrt(pi/3.)
         bpp=-4./3.*bpeak*sqrt(pi/5.)
         do kc=1,nn
            b(2,kc)=bpk*(.375*r(kc)**3-0.5*r(1)*r(kc)**2
     $                    +r(nn)**4/(8.*r(kc)))
            aj(3,kc)=bpp*r(kc)*sin(pi*(r(kc)-r(nn)))
         enddo
         endif
       else
         concof=bpeakbot*bscl*r(nn)
         do kc=1,nn
            b(2,kc)=concof*qk(kc,3)
         enddo
       endif
      endif
c
c  urc: determine scaling factors
c
      if(iscale.eq.1) then
        scdiff=1.       
      else if(iscale.eq.2) then
        scdiff=opr        
      else
        scdiff=opm         
      endif
      tscale= 1./scdiff 
      vscale=scdiff
      pscale=scdiff*oek
      escale=scdiff**2 * ocorevol  / enscale
      alum0=opr*radtop*radbot/y00**2
c
c    -energies
c
      if(init .le. 0) then
         call kei(envp,envt,adrke,amcke)
         env=envp+envt
         enbp=0.
         enbt=0.
         enb=0.
         if(init .gt. -10) then
            call mei(enbp,enbt,apome,atome)
            enb=enbp+enbt
         endif
      else
         env=0.
         enb=0.
      endif
      write(6,'('' Ekin & Emag= '',2f9.2)') env/escale,enb/escale
c
      if(nstep.lt.1) then
        write(6,'(/'' R   T  ''/)')
        do kc=1,nn
          write(6,'(f8.5,1x,f8.5)')
     &    r(kc)/radtop,real(s(1,kc)-s(1,1))*y00
        enddo
      endif
c
      write(6,900) ra,ek,pr,prmag
  900 format(3x,"rayleigh =",1pe10.3,3x,"ekman =",1pe10.3,/,
     $3x,"prandtl =",1pe9.2,3x,"mag prandtl =",1pe9.2,/)
      write(6,904) vscale,pscale,escale
  904 format(3x,"scales, V= ",1pe10.3,"  P= ",1pe10.3,
     $ "  E= ",1pe10.3)
      if(enscale.ne.1.)
     $ write(6,'(3x,'' energy multiplied by '',f8.4)') enscale
      write(6,6) dt/tscale,kstep,time/tscale
    6 format(3x,"dt =",f10.8,3x,"kstep =",i7,3x,
     $"time =",f10.6/)
c
c *** construct lu decomposed matrices for w, z, s, b, aj, and p equations
c
      call ludc
c
      return
      end
c
c*********************************************************************
c
      subroutine zerorot(init)
c
c  For init=-7 or -9, reset time-derivatives of magnetic field to zero
c  For init<=-8, set m=0 terms of toroidal velocity to zero
c  To be used in conjunction with kinematic dynamo runs
c
      include 'param.f'
      include 'com5.f'
c
      if(init.eq.-8) go to 10
      do kc=1,nn
        do lm=1,nlma
           dbdt1(lm,kc)=(0.0,0.0)
           djdt1(lm,kc)=(0.0,0.0)
        enddo
      enddo
c
      if(init.gt.-8) return
c
   10 do kc=1,nnp1
        do lm=1,lmax+1
          z(lm,kc)=(0.0,0.0)
        enddo
      enddo
c
      return
      end
