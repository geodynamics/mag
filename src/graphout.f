      subroutine graphout(kc,iwrit)
c
c  supplied by urc
c
c  called in amhd
c
c  output of components of velocity, magnetic field vector and
c  entropy for graphics          
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
      dimension vr(nrp,ni),vt(nrp,ni),vp(nrp,ni)
      dimension br(nrp,ni),bt(nrp,ni),bp(nrp,ni)
      dimension cbr(nrp,ni),cbt(nrp,ni),cbp(nrp,ni)
      dimension sr(nrp,ni)
      dimension brf(nrp,ni)
      complex brfc(ncp,ni)
      complex blm(nlma)
c
      real*4 dumm0,dumm1(ni),dummy(nj,ni)
c
      equivalence (vr,vrc),(vt,vtc),(vp,vpc)
      equivalence (br,brc),(bt,btc),(bp,bpc)
      equivalence (cbr,cbrc),(cbt,cbtc),(cbp,cbpc)
      equivalence (sr,sc)
      equivalence (brf,brfc)
c
      save smin
c
      if(iwrit.eq.0) return
c
      nim=ni
      nis=ngcolat
      njm=nja
      njs=nglon
      if(kc.eq.1) smin=sr(1,1)
c
      if(kc.gt.0) go to 300
c
c  write header & colatitudes
c
      if(iwrit.lt.2) then
        write(14,'(A64)') runid
        if(iwrit.lt.0) 
     $  write(14,'(/'' Time,nn,ni,nj,ngrad,nglon,minc:'')')
        write(14,'(e12.5,1x,7i6)') time,nn,ni,nj,ngrad,ngcolat,
     $    nglon,minc
        write(14,'(5e14.6)') ra,ek,pr,prmag,radratio

        if(iwrit.lt.0) write(14,'(/'' Colatitudes '')')
        write(14,'(129(1x,f8.5))') (qi(ic,5),ic=1,nim,nis)
      else
        write(14) runid
        dumm0=time
        write(14) dumm0,nn,ni,nj,ngrad,ngcolat,nglon,minc
        dumm1(1)=ra
        dumm1(2)=ek
        dumm1(3)=pr
        dumm1(4)=prmag
        dumm1(5)=radratio
        write(14) (dumm1(i),i=1,5)
        do ic=1,nim,nis
          dumm1(ic)=qi(ic,5)
        enddo
        write(14) (dumm1(ic),ic=1,nim,nis)
      endif
c
      return
c
c  output for radial levels
c
  300 if(iwrit.lt.0) write(14,'(/'' Radial level, radius        '')')
      if(iwrit.lt.2) then
        write(14,'(i4,1x,f9.5)') kc,r(kc)/radtop
      else
        dumm0=r(kc)/radtop
        write(14) kc,dumm0          
      endif
c
c write entropy 
c
      if(iwrit.lt.0) write(14,'(/'' S: '')')
      if(iwrit.lt.2) then
        do ic=1,nim,nis
         write(14,902) ((sr(jc,ic)-smin),jc=1,njm,njs)
        enddo
      else
        do ic=1,nim,nis
         do jc=1,njm,njs
           dummy(jc,ic)=sr(jc,ic)-smin
         enddo
        enddo
        write(14) ((dummy(jc,ic),jc=1,njm,njs),ic=1,nim,nis)
       endif
  350 continue
c
c  calculate and write radial velocity
c
      if(iwrit.lt.0) write(14,'(/'' Vr: '')')
      if(kc.eq.1.and.nfilt.gt.0.and.ivfilt.gt.0) then
       call filter(w,blm,ivfilt,alfilt,nfilt,1.0)
       call spherictf(blm,brfc)
       do ic=1,nim,nis    
        do jc=1,njm,njs  
          dummy(jc,ic)=brf(jc,ic)*qk(ivfilt,1)/vscale
        enddo
       enddo
      else  
      do ic=1,nim,nis    
        do jc=1,njm,njs  
          dummy(jc,ic)=vr(jc,ic)*qk(kc,1)/vscale
        enddo
      enddo
      endif
      if(iwrit.lt.2) then
        do ic=1,nim,nis
         write(14,900) (dummy(jc,ic),jc=1,njm,njs)
        enddo
      else
        write(14) ((dummy(jc,ic),jc=1,njm,njs),ic=1,nim,nis)
      endif
c
c  calculate and write latitudinal velocity
c
      if(iwrit.lt.0) write(14,'(/'' Vt: '')')
      do ic=1,nim,nis
       do jc=1,njm,njs
         dummy(jc,ic)=vt(jc,ic)*qk(kc,3)/(vscale*qi(ic,3))
       enddo      
      enddo      
      if(iwrit.lt.2) then
        do ic=1,nim,nis
         write(14,900) (dummy(jc,ic),jc=1,njm,njs)
        enddo
      else
        write(14) ((dummy(jc,ic),jc=1,njm,njs),ic=1,nim,nis)
      endif
c
c  calculate and write longitudinal velocity
c
      if(iwrit.lt.0) write(14,'(/'' Vp: '')')
      do ic=1,nim,nis
       do jc=1,njm,njs
        dummy(jc,ic)=vp(jc,ic)*qk(kc,3)/(vscale*qi(ic,3))
       enddo   
      enddo   
      if(iwrit.lt.2) then
        do ic=1,nim,nis
         write(14,900) (dummy(jc,ic),jc=1,njm,njs)
        enddo
      else
        write(14) ((dummy(jc,ic),jc=1,njm,njs),ic=1,nim,nis)
      endif
c
c  calculate and write radial magnetic field
c
      if(iwrit.lt.0) write(14,'(/'' Br: '')')
      if(kc.eq.1.and.nfilt.gt.0) then
       call filter(b,blm,1,alfilt,nfilt,dipfilt)
       call spherictf(blm,brfc)
       do ic=1,nim,nis    
        do jc=1,njm,njs  
          dummy(jc,ic)=brf(jc,ic)*qk(1,1)
        enddo
       enddo
      else  
      do ic=1,nim,nis
       do jc=1,njm,njs
        dummy(jc,ic)=br(jc,ic)*qk(kc,1)
       enddo      
      enddo      
      endif
      if(iwrit.lt.2) then
        do ic=1,nim,nis
         write(14,901) (dummy(jc,ic),jc=1,njm,njs)
        enddo
      else
        write(14) ((dummy(jc,ic),jc=1,njm,njs),ic=1,nim,nis)
      endif
c
c  calculate and write latitudinal magnetic field
c
      if(iwrit.lt.0) write(14,'(/'' Bt: '')')
      do ic=1,nim,nis
       do jc=1,njm,njs
        dummy(jc,ic)=bt(jc,ic)*qk(kc,3)/qi(ic,3)
       enddo
      enddo
      if(iwrit.lt.2) then
        do ic=1,nim,nis
         write(14,901) (dummy(jc,ic),jc=1,njm,njs)
        enddo
      else
        write(14) ((dummy(jc,ic),jc=1,njm,njs),ic=1,nim,nis)
      endif
c
c  calculate and write longitudinal magnetic field
c
      if(iwrit.lt.0) write(14,'(/'' Bp: '')')
      do ic=1,nim,nis
       do jc=1,njm,njs
        dummy(jc,ic)=bp(jc,ic)*qk(kc,3)/qi(ic,3)
       enddo
      enddo
      if(iwrit.lt.2) then
        do ic=1,nim,nis
         write(14,901) (dummy(jc,ic),jc=1,njm,njs)
        enddo
      else
        write(14) ((dummy(jc,ic),jc=1,njm,njs),ic=1,nim,nis)
      endif
c
c  if iwrit > 2 write also components of curl(B)
c
      if(iwrit.le.2) return
c
c  r-component
      do ic=1,nim,nis
       do jc=1,njm,njs
        dummy(jc,ic)=cbr(jc,ic)*qk(kc,1)
       enddo
      enddo
      write(14) ((dummy(jc,ic),jc=1,njm,njs),ic=1,nim,nis)
c  theta-component
      do ic=1,nim,nis
       do jc=1,njm,njs
        dummy(jc,ic)=cbt(jc,ic)*qk(kc,3)/qi(ic,3)
       enddo
      enddo
      write(14) ((dummy(jc,ic),jc=1,njm,njs),ic=1,nim,nis)
c  phi-component
      do ic=1,nim,nis
       do jc=1,njm,njs
        dummy(jc,ic)=cbp(jc,ic)*qk(kc,3)/qi(ic,3)
       enddo
      enddo
      write(14) ((dummy(jc,ic),jc=1,njm,njs),ic=1,nim,nis)

  900 format(256(1X,f7.2))
  901 format(256(1X,f7.3))
  902 format(256(1X,f7.5))
c
      return
      end
