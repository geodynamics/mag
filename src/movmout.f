      subroutine movmout(kc,kc0)
c
c  called in nl and amhd
c
c  supplied by urc
c
c  output of v_r at mid-depth, and B_r at mid-depth and the top    
c  surface for producing movie           
c  for kc=0 a header is written
c  for kc>0 values at radial level kc are written
c
      include 'param.f'
      include 'com1.f'
      include 'com3.f'
      include 'com4.f'
      include 'com5.f'
      include 'com6.f'  
      include 'com8.f'
c
      dimension vr(nrp,ni)
      dimension vp(nrp,ni)
      dimension vt(nrp,ni)
      dimension br(nrp,ni)
      dimension brf(nrp,ni)
c
      double complex brfc(ncp,ni),blm(nlma)
c
      equivalence (vr,vrc)
      equivalence (vp,vpc)
      equivalence (vt,vtc)
      equivalence (br,brc)
      equivalence (brf,brfc)
c
      if(kc.eq.0) then           
c
c  write header
c
        ndat=3
        if(kc0.eq.9) ndat=1
        write(20,'(A64)') runid
c
        write(20,'(2e13.5)') tscale,vscale
        write(20,'(8i6)') nn,ni,nj,minc,iframes,ndat,nglon,ngcolat
c
c write colatitudes           
c
        write(20,'(128(1x,f8.5))') (qi(ic,5),ic=1,ni,ngcolat)     
        return
c
c  output of B_r at outer surface  
c
      else if (kc.eq.1) then 
        write(20,'(I5,1x,f10.6)') imovct,time/tscale
        if(nfilt.le.0) then
         do ic=1,ni,ngcolat
          write(20,901) (br(jc,ic)*qk(1,1),jc=1,nja,nglon)
         enddo
        else
         call filter(b,blm,1,alfilt,nfilt,dipfilt)
         call spherictf(blm,brfc)
         do ic=1,ni,nglon
          write(20,901) (brf(jc,ic)*qk(1,1),jc=1,nja,nglon)
         enddo
        endif
c
c  output of phi-component of velocity at level kc0
c
      endif
      if(kc0.eq.9) return
      if (kc.eq.kc0) then 
        do ic=1,ni,ngcolat
          write(20,902) (vp(jc,ic)*qk(kc0,3)/(qi(ic,3)*vscale),
     &      jc=1,nja,nglon)
        enddo
        do ic=1,ni,ngcolat
          write(20,902) (vt(jc,ic)*qk(kc0,3)/(qi(ic,3)*vscale),
     &      jc=1,nja,nglon)
        enddo
      endif
c
      return
c     
  900 format(256(1X,f7.1))
  901 format(256(1X,f7.4))
  902 format(256(1X,f7.2))
c
      end
