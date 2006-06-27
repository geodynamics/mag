      subroutine movaout(kc)
c
c  supplied by urc
c
c  called in nl and amhd
c
c  output of longitudinally averaged B_phi, j_phi, and v_phi      
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
      dimension vp(nrp,ni)
      dimension bp(nrp,ni),bt(nrp,ni),br(nrp,ni)
c
      dimension bpa(ni),bta(ni),bra(ni),vpa(ni)
      dimension ajp(ni)
c
      equivalence (vp,vpc)
      equivalence (bp,bpc),(bt,btc),(br,brc)
c
      if(kc.gt.0) go to 300
c
c  write header
c
      write(19,'(A64)') runid
c
      write(19,'(2e13.5)') tscale,vscale
      write(19,'(5i6)') nn,ni,nj,minc,iframes
c
c write radii and colatitudes           
c
      write(19,'(65(1x,f8.5))') (r(ir)/radtop,ir=1,nn)
      write(19,'(65(1x,f8.5))') (qi(ic,5),ic=1,ni)     
      return
c
c  calculate averages
c
  300 do ic=1,ni
        bpa(ic)=0.0
        bta(ic)=0.0
        bra(ic)=0.0
        vpa(ic)=0.0
        do jc=1,nja
          bpa(ic)=bpa(ic)+bp(ic,jc)
          bta(ic)=bta(ic)+bt(ic,jc)
          bra(ic)=bra(ic)+br(ic,jc)
          vpa(ic)=vpa(ic)+vp(ic,jc)
        enddo
        bpa(ic)=bpa(ic)/real(nja)  *qk(kc,3)/qi(ic,3)      
        bta(ic)=bta(ic)/real(nja)  *qk(kc,3)/qi(ic,3)
        bra(ic)=bra(ic)/real(nja)  *qk(kc,1)
        vpa(ic)=vpa(ic)/real(nja)  *qk(kc,3)/qi(ic,3)
      enddo
c
c  calculate phi-component of curl(B) for radius kc-1
c
      kup=mod(kc+1,3)+1
      kmd=mod(kc+2,3)+1
      klw=mod(kc+3,3)+1 
      do ic=1,ni
        bts(ic,klw)=bta(ic)
      enddo
c
      if(kc.eq.1) go to 320
      if(kc.eq.2) then
          do ic=1,ni   
            ajp(ic) = -brdt(ic)
     &      + (bts(ic,kmd)-bts(ic,klw))/(r(1)-r(2))
     &      +  bts(ic,kmd)/r(1)
          enddo
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
      do ic=1,ni
        ajp(ic)= -brdt(ic)
     &  +bts(ic,klw)*frlw+bts(ic,kmd)*frmd+bts(ic,kup)*frup
      enddo
c
  320 do ic=2,ni-1
        brdt(ic)=(bra(ic+1)-bra(ic-1))/
     &   (r(kc)*(qi(ic+1,5)-qi(ic-1,5)))
      enddo
        brdt(1)=(bra(2)-bra(1))/
     &   (r(kc)*(qi(2,5)-qi(1,5)))
        brdt(ni)=(bra(ni)-bra(ni-1))/
     &   (r(kc)*(qi(ni,5)-qi(ni-1,5)))
c
      if(kc.eq.1) then        
        write(19,'(I5,1x,f8.5)') imovct,time/tscale
      else
c
c  write j_phi for level kc-1
c
        write(19,900) (ajp(jc),ic=1,ni)
      endif
c
c write v_phi     
c
      write(19,902) (vpa(ic)/vscale,ic=1,ni)
c
c write b_phi     
c
      write(19,901) (bpa(ic),ic=1,ni)
c
      if(kc.lt.nn) return
c
c  calculate and write j_phi for final radius kc=nn
c
        do ic=1,ni
          ajp(ic)= -brdt(ic)
     &     + (bts(ic,kmd)-bts(ic,klw))/(r(nn-1)-r(nn))
     &     + bts(ic,klw)/r(nn)
        enddo
      write(19,900) (ajp(ic),ic=1,ni)
c
  900 format(256(1X,f7.2))
  901 format(256(1X,f7.4))
  902 format(256(1X,f7.3))
c
      return
      end
