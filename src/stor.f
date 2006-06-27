      subroutine stor
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com1.f'
      include 'com3.f'
      include 'com4.f'
      include 'com5.f'
c---------------------------------------------------------------
c
c *** store output
c
      rewind(10)
      write(10) time,dt,ra,pr,prmag,ek,radratio,
     $          kstep,nn,ni,nj,minc
      write(10) w,z,p,s
      write(10) dsdt1,dwdt1,dzdt1,dpdt1
      write(10) b,aj,dbdt1,djdt1
c
      if(iprnt .eq. nprnt) close(10)
      write(6,1) rstfile
    1 format(/,8x,"**  restart data stored in ",a32,"  **",/)
c
      return
      end
