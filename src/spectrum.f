      subroutine spectrum(imode)  
c
c---------------------------------------------------------------
c
c  calculation of power spectrum of kinetic and magnetic energy as
c  function of harmonic degree  (urc)
c
c  imode=0: sum all modes of given l and various m
c  imode=1: sum all modes of given m and various l
c
c   called in amhd
c
c---------------------------------------------------------------
c
      include 'param.f'
      include 'com1.f'
      include 'com3.f'
      include 'com4.f'
      include 'com5.f'
      include 'com8.f'
c
      complex c
c
      dimension ekin(0:nlaf)
      dimension emag(0:nlaf)
      dimension emcmb(0:9)
      dimension vtcmb(0:nlaf)
      dimension vpcmb(0:nlaf)
      dimension rvl(nn,0:nlaf)
      dimension rvm(nn,0:nlaf)
c
c  rvl(nn,0) assembles the axisymmetric portion of kinetic energy
c  rvm(nn,0) assembles the axisymmetric portion of kinetic energy
c
      cabssq(c)=real(c)**2+aimag(c)**2
c---------------------------------------------------------------
c
      isort=5
      lblock=(nlaf-1)/10+1
      lstep=1
      lc0=0
      if(imode.gt.0) then
       isort=10
       lblock=(nlaf-1)/(10*minc)+1
       lstep=minc
       lc0=1
      endif
c
      do 20 lc=0,nlaf
      do 20 kc=1,nn
         rvl(kc,lc)=0.
         rvm(kc,lc)=0.
   20 continue
      do 22 lc=0,9 
        emcmb(lc)=0.0
   22 continue
      polmag=0.0
      polaxi=0.0
c
      do 30 lm=nlma,2,-1
         lc=ql(lm,isort)+lc0
         do 31 kc=1,nn
            rvl(kc,lc)=rvl(kc,lc)+
     $         ql(lm,6)*(ql(lm,3)*qk(kc,1)*cabssq(w(lm,kc))+
     $         cabssq(z(lm,kc))+4.*cabssq(dw(lm,kc)))
            rvm(kc,lc)=rvm(kc,lc)+
     $         ql(lm,6)*(ql(lm,3)*qk(kc,1)*cabssq(b(lm,kc))+
     $         cabssq(aj(lm,kc))+4.*cabssq(db(lm,kc)))
   31    continue
         if(imode.lt.1.and.ql(lm,10).lt.0.1) then
           do 32 kc=1,nn
              rvl(kc,0)=rvl(kc,0)+
     $           ql(lm,6)*(ql(lm,3)*qk(kc,1)*cabssq(w(lm,kc))+
     $           cabssq(z(lm,kc))+4.*cabssq(dw(lm,kc)))
              rvm(kc,0)=rvm(kc,0)+
     $           ql(lm,6)*(ql(lm,3)*qk(kc,1)*cabssq(b(lm,kc))+
     $           cabssq(aj(lm,kc))+4.*cabssq(db(lm,kc)))
   32      continue
           polaxi=polaxi +
     $           ql(lm,6)*(ql(lm,3)*qk(1,1)*cabssq(b(lm,1))
     $           +4.*cabssq(db(lm,1)))
         endif
         if(imode.lt.1) then
           polmag=polmag +
     $           ql(lm,6)*(ql(lm,3)*qk(1,1)*cabssq(b(lm,1))
     $           +4.*cabssq(db(lm,1)))
           lc=nint(ql(lm,4))
           vpcmb(lc)=vtcmb(lc)+
     $         ql(lm,6)*(ql(lm,3)*qk(1,1)*cabssq(w(lm,1))
     $         +4.*cabssq(dw(lm,1)))
           vtcmb(lc)=vtcmb(lc)+
     $         ql(lm,6)*cabssq(z(lm,1))
         endif
         if(imode.lt.1 .and. ql(lm,4).lt.10) then 
           lc=nint(ql(lm,4))
           emcmb(lc)=emcmb(lc) +
     $           ql(lm,6)*(ql(lm,3)*qk(1,1)*cabssq(b(lm,1))
     $           +4.*cabssq(db(lm,1)))
         endif
c
   30 continue
c
      surface=4.*pi*radtop*radtop
      dipolax=sqrt(2.*ql(2,6)*(ql(2,3)*qk(1,1)*cabssq(b(2,1))
     $  +4.*cabssq(db(2,1))) /surface)
      dipole=dipolax
      if(minc.eq.1) then 
        lm=lmax+2
        dipole=sqrt(dipolax**2 +
     $   2.*ql(lm,6)*(ql(lm,3)*qk(1,1)*cabssq(b(lm,1))
     $  +4.*cabssq(db(lm,1))) /surface)
      endif
      polmag=sqrt(2.*polmag/surface)
      polaxi=sqrt(2.*polaxi/surface)
c
      do 33 lc=0,9
        emcmb(lc)=2.*emcmb(lc)/surface
   33 continue
      do 34 lc=1,nlaf
        vpcmb(lc)=2.*vpcmb(lc)/surface
        vtcmb(lc)=2.*vtcmb(lc)/surface
   34 continue
c
      ekinmax=0.0
      ekinaxi=0.0
      ekinsum=0.0
      emagmax=0.0
      emagaxi=0.0
      emagsum=0.0
c
      do 50 lc=0,nlaf
c
        call chebtf(1,ns2,1,nn1,ns2,rvl(1,lc),rvl(1,lc),rvl(1,lc),
     $   wsave,work,work(2,1),work(2,1),work(2,1),trigsc,ifaxc,k2k)
        call chebtf(1,ns2,1,nn1,ns2,rvm(1,lc),rvm(1,lc),rvm(1,lc),
     $   wsave,work,work(2,1),work(2,1),work(2,1),trigsc,ifaxc,k2k)
c
        ekin(lc)=0.0
        emag(lc)=0.0
        rvl(1,lc)=0.5*rvl(1,lc)
        rvl(nn,lc)=0.5*rvl(nn,lc)
        rvm(1,lc)=0.5*rvm(1,lc)
        rvm(nn,lc)=0.5*rvm(nn,lc)
        do 45 ncb=1,nn,2
         nc=nnp1-ncb
         ekin(lc)=ekin(lc)+rvl(nc,lc)*qn(nc,3)
         emag(lc)=emag(lc)+rvm(nc,lc)*qn(nc,3)
   45   continue
        ekin(lc)=0.5*anorm*ekin(lc)
        if(ekin(lc).gt.ekinmax) ekinmax=ekin(lc)
        if(lc.gt.0) ekinsum=ekinsum+ekin(lc)
        emag(lc)=0.5*oekpm*anorm*emag(lc)
        if(emag(lc).gt.emagmax) emagmax=emag(lc)
        if(lc.gt.0) emagsum=emagsum+emag(lc)
   50 continue
c
c  print
c
      ekinaxi=ekin(0)/escale
      ekinsum=ekinsum/escale
      emagaxi=emag(0)/escale
      emagsum=emagsum/escale
      dipolax=sign(dipolax,real(b(2,1)))
      dipole=dipole
      polmag=polmag
      polaxi=polaxi
      tiltdipole=0.0
      if(minc.eq.1) 
     $  tiltdipole=atan2(abs(b(lmax+2,1)),real(b(2,1)))
c
      if(imode.lt.1) then
        write(6,'('' Spectrum Ekin vs l    Total/Axisym.= '',
     $   2f9.2 )') ekinsum,ekinaxi
        write(6,'(10(2x,i3,3x))') (l,l=0,9)
      else
        write(6,'('' Spectrum Ekin vs m    Total  = '',
     $   1f9.2 )') ekinsum
        write(6,'(10(2x,i3,3x))') (l,l=0,9*lstep,lstep)
      endif
      do 60 lb=1,lblock
      l0=(lb-1)*10*lstep
      l1=min(l0+9*lstep,nlaf-1)
      write(6,'(10(f8.2))') (ekin(l+1)/escale,l=l0,l1,lstep)
   60 continue
c
      if(imode.lt.1) then
        write(6,'(/'' Poloidal/Toroidal Ekin (x10) at CMB'')')
        write(6,'(10(2x,i3,3x))') (l,l=0,9)
      do 61 lb=1,lblock
      l0=(lb-1)*10*lstep
      l1=min(l0+9*lstep,nlaf-1)
      write(6,'(10(f8.2))') (10.*vpcmb(l+1)/escale,l=l0,l1,lstep)
      write(6,'(10(f8.2)/)') (10.*vtcmb(l+1)/escale,l=l0,l1,lstep)
   61 continue
        write(6,'(/'' CMB field: Total='',f7.4,
     $  ''  Axi='',f7.4,''  Dipole='',f7.4,''  Ax.Dip.='',
     $  f7.4,''  Tilt='',f5.1)') 
     $  polmag,polaxi,dipole,dipolax,tiltdipole*180./pi
        write(6,'('' Spectrum Emag vs l    Total/Axisym.= '',
     $   2f9.2 )') emagsum,emagaxi
        write(6,'(10(2x,i3,3x))') (l,l=0,9)
        write(6,'('' CMB:   '',9f8.5)') (emcmb(l),l=1,9)
      else
        write(6,'('' Spectrum Emag vs m    Total  = '',
     $   1f9.2 )') emagsum
        write(6,'(10(2x,i3,3x))') (l,l=0,9*lstep,lstep)
      endif
      do 70 lb=1,lblock
      l0=(lb-1)*10*lstep
      l1=min(l0+9*lstep,nlaf-1)
      write(6,'(10(f8.2))') (emag(l+1)/escale,l=l0,l1,lstep)
   70 continue
c
c  output on logs-file
c
      if(imode.lt.1) then
        write(16,'(f10.6,1x,65f10.4)') time/tscale,(ekin(l)/escale,
     &   l=1,nlaf)
        write(16,'(f10.6,1x,65f10.4)') time/tscale,(emag(l)/escale,
     &   l=1,nlaf,lstep)
      else
        write(16,'(f10.6,1x,65f10.4)') time/tscale,(ekin(l)/escale,
     &   l=1,nlaf,lstep)
        write(16,'(f10.6,1x,65f10.4)') time/tscale,(emag(l)/escale,
     &   l=1,nlaf,lstep)
      endif
c
      return
      end
