      subroutine moveout(kc)
c
c  supplied by urc
c
c  called in nl and amhd
c
c  output of entropy, z-vorticity, and z-field in equatorial plane
c  for procuding movie           
c  for kc=0 a header is written
c  for kc>0 values at radial level kc are written
c
      include 'param.f'
      include 'com1.f'
      include 'com3.f'
      include 'com4.f'
      include 'com6.f'  
      include 'com8.f'
c
      dimension vr(nrp,ni),vp(nrp,ni)
      dimension bt(nrp,ni)
      dimension sr(nrp,ni)
c
      dimension ur(0:nja+1),wz(nja)
c
      equivalence (vr,vrc),(vp,vpc)
      equivalence (bt,btc)
      equivalence (sr,sc)
c
      if(kc.gt.0) go to 300
c
c  write header
c
      write(18,'(A64)') runid
c
      write(18,'(2e13.5)') tscale,vscale
      write(18,'(5i6)') nn,ni,nj,minc,iframes
c
c write radii           
c
      write(18,'(65(1x,f8.5))') (r(ir)/radtop,ir=1,nn)
      return
c
c  output for radial levels
c
  300 iequat=ni/2
      dphi2=16.*atan(1.)/real(nj)
c
c  calculate z-vorticity
c
      kup=mod(kc+1,3)+1
      kmd=mod(kc+2,3)+1
      klw=mod(kc+3,3)+1 
      do jc=1,nja
        up(jc,klw)=vp(jc,iequat)*qk(kc,3)
        ur(jc)=    vr(jc,iequat)*qk(kc,1)
      enddo
      ur(0)=ur(nja)
      ur(nja+1)=ur(1)
c
      if(kc.eq.1) go to 320
      if(kc.eq.2) then
        if(ktopv.eq.1) then
          do jc=1,nja
            wz(jc)=0.0
          enddo
        else
          do jc=1,nja
            wz(jc) =
     &      + (up(jc,kmd)-up(jc,klw))/(r(1)-r(2))
          enddo
        endif
        go to 320
      endif
c
c  kc > 2
c
      drup=r(kc-2)-r(kc-1)
      drlw=r(kc-1)-r(kc)
      frlw=-drup/(drlw*(drlw+drup))
      frmd=(drup-drlw)/(drup*drlw) + 1./r(kc-1)
      frup= drlw/(drup*(drlw+drup))
      do jc=1,nja
        wz(jc)= -urdp(jc)
     &  +up(jc,klw)*frlw+up(jc,kmd)*frmd+up(jc,kup)*frup
      enddo
c
  320 do jc=1,nja
        urdp(jc)=(ur(jc+1)-ur(jc-1))/(dphi2*r(kc))
      enddo
c
      if(kc.eq.1) then        
        write(18,'(I5,1x,f8.5)') imovct,time/tscale
      else
c
c  write vorticity for level kc-1
c
        write(18,900) (wz(jc)/vscale,jc=1,nja)
      endif
c
c write entropy 
c
      write(18,901) (sr(jc,iequat),jc=1,nja)
c
c  calculate and write z-field           
c
      write(18,902) (bt(jc,iequat)*qk(kc,3),jc=1,nja)
c
      if(kc.lt.nn) return
c
c  calculate and write z-vorticity for final radius kc=nn
c
      if(kbotv.eq.1) then 
        do jc=1,nja
          wz(jc)=0.0
        enddo
      else
        do jc=1,nja
          wz(jc)=
     &     + (up(jc,kmd)-up(jc,klw))/(r(nn-1)-r(nn))
        enddo
      endif
      write(18,900) (wz(jc)/vscale,jc=1,nja)
c
  900 format(256(1X,f7.1))
  901 format(256(1X,f7.4))
  902 format(256(1X,f7.3))
c
      return
      end
