      subroutine dtchck(kstep,newdt,dt,dtnew,
     $dtmin,dtmax,dtr,dth,ifirst,kcour)
c
c *** Check if Courant criterion based on combined
c *** fluid and Alfven velocity is satisfied
c *** Returns new value of time step dtnew
c
      kcour=max(kcour-1,0)
c
      if(ifirst .eq. 1) then
         ifirst=0
         if(abs((dt-dtmax)/dt) .lt. 1.e-7) dt=dtmax
      endif
c
      dtnew=dt
c
      dtlo=0.5
      dtrh=min(dtr,dth)
      dt2=min(0.5*(dtlo+1.0)*dtrh,dtmax)
      if(dt .le. dtmax) go to 40
      write(6,20) kstep,dtmax
      go to 50
c
   40 dtlim=     dtrh
      if(dt .le. dtlim) go to 10
      write(6,20) kstep,dtr,dth
   20 format(/,1x,"step=",i6,2x,"dt> dtr or dth=",2e12.4)
      go to 50
c
   10 dtlim=dtlo*dtrh
      if((dt .ge. dtlim) .or. (dt .ge. dtmax)) return
      write(6,30) kstep,dtr,dth
   30 format(/,1x,"step=",i6,2x,"2*dt< dtr and dth=",2e12.4)
c
   50 dtnew=dt2
      if(dtnew .lt. dtmin) dtmin=0.  ! Signal to stop
      newdt=1
      write(6,25) dtnew
   25 format(12x,"dt changed to ",e12.4)
c
      kcour=2
c
      return
      end
