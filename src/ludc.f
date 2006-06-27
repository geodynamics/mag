      subroutine ludc
c
c  contruct matrices for chebychev collocation of linear terms in
c  the governing equations
c  perform LU-decomposition of matrices
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com1.f'
      include 'com2.f'
      include 'com3.f'
      include 'com4.f'
      include 'com7.f'
c---------------------------------------------------------------
c
c *** chebyshev collocation, l=0 terms for s and p
c
      do 342 nc=1,nn
         do 343 kc=1,nn
            s0mat(kc,nc)=(cheb(nc,kc)*oodt-
     $         alpha*opr*(4.*d2cheb(nc,kc)+
     $         qk(kc,2)*dcheb(nc,kc)))*anorm*sscl
            p0mat(kc,nc)=2.*dcheb(nc,kc)*anorm
  343    continue
  342 continue
c
c *** boundary conditions
c
      do 351 nc=1,nnaf
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
         p0mat(nps2,nc)=anorm
  351 continue
      if(nnaf .lt. nn) then
         do 353 nc=nnaf+1,nn
            s0mat(1,nc)=0.
            s0mat(nn,nc)=0.
            p0mat(nps2,nc)=0.
  353    continue
      endif
c
c *** normalization
c
      do 361 kc=1,nn
         s0mat(kc,1)=0.5*s0mat(kc,1)
         s0mat(kc,nn)=0.5*s0mat(kc,nn)
         p0mat(kc,1)=0.5*p0mat(kc,1)
         p0mat(kc,nn)=0.5*p0mat(kc,nn)
  361 continue
c
c *** construct lu decomposed matrices
c
      call sgefa(s0mat,nn,nn,is0,info)
      if(info .ne. 0) stop '28'
      call sgefa(p0mat,nn,nn,ip0,info)
      if(info .ne. 0) stop '29'
c
c *** chebyshev collocation, l>0 terms for s
c  
      do 200 l=1,lmax
c
      lm=l+1
c
      do 242 nc=1,nn
         do 243 kc=1,nn
            smat(kc,nc,l)=(cheb(nc,kc)*oodt-
     $         alpha*opr*(4.*d2cheb(nc,kc)+
     $         qk(kc,2)*dcheb(nc,kc)-
     $         ql(lm,3)*qk(kc,1)*cheb(nc,kc)))*anorm*sscl
  243    continue
  242 continue
c
c *** boundary conditions
c
      do 251 nc=1,nnaf
         if(ktops .eq. 1) then
            smat(1,nc,l)=anorm
         else
            smat(1,nc,l)=2.*dcheb(nc,1)*anorm
         endif
         if(kbots .eq. 1) then
            smat(nn,nc,l)=cheb(nc,nn)*anorm
         else
            smat(nn,nc,l)=2.*dcheb(nc,nn)*anorm
         endif
  251 continue
      if(nnaf .lt. nn) then
         do 253 nc=nnaf+1,nn
            smat(1,nc,l)=0.
            smat(nn,nc,l)=0.
  253    continue
      endif
c
c *** normalization
c
      do 261 kc=1,nn
         smat(kc,1,l)=0.5*smat(kc,1,l)
         smat(kc,nn,l)=0.5*smat(kc,nn,l)
  261 continue
c
c *** construct lu decomposed matrix   
c
      call sgefa(smat(1,1,l),nn,nn,is(1,l),info)
      if(info .ne. 0) stop '31'
c
  200 continue
c
c *** chebyshev collocation, l>0 terms for magnetic potentials
c          
      do 400 l=1,lmax
c
      lm=l+1
c
      do 442 nc=1,nn
         do 443 kc=1,nn
            bmat(kc,nc,l)=(oodt*ql(lm,3)*qk(kc,1)*cheb(nc,kc)-
     $         alpha*ql(lm,11)*opm*qk(kc,1)*
     $         (4.*d2cheb(nc,kc)-
     $         ql(lm,3)*qk(kc,1)*cheb(nc,kc)))*anorm*bscl
            ajmat(kc,nc,l)=(oodt*ql(lm,3)*qk(kc,1)*cheb(nc,kc)-
     $         alpha*ql(lm,11)*opm*qk(kc,1)*
     $         (4.*d2cheb(nc,kc)-
     $         ql(lm,3)*qk(kc,1)*cheb(nc,kc)))*anorm*bscl
  443    continue
  442 continue
c
c  *** boundary conditions
c
      do 451 nc=1,nnaf
         bmat(1,nc,l)=((2.*dcheb(nc,1)+float(l)/r(1)*cheb(nc,1))+
     $      cmb*(4.*d2cheb(nc,1)-float(l*(l+1))/(r(1)*r(1))*
     $      cheb(nc,1)))*anorm
         ajmat(1,nc,l)=(2.*cmb*dcheb(nc,1)+cheb(nc,1))*anorm
         if(kbotb .eq. 2) then
            bmat(nn1,nc,l)=4.*d2cheb(nc,nn)*anorm
            ajmat(nn,nc,l)=2.*dcheb(nc,nn)*anorm
         else
            bmat(nn,nc,l)=(2.*dcheb(nc,nn)-float(l+1)/r(nn)*
     $         cheb(nc,nn))*anorm
            ajmat(nn,nc,l)=cheb(nc,nn)*anorm
         endif
         if(l.eq.1.and.imagcon.lt.0) bmat(nn,nc,1)=cheb(nc,nn)*anorm
  451 continue
c
      if(nnaf .lt. nn) then
         do 453 nc=nnaf+1,nn
            bmat(1,nc,l)=0.
            ajmat(1,nc,l)=0.
            if(kbotb .eq. 2) then
               bmat(nn1,nc,l)=0.
               ajmat(nn,nc,l)=0.
            else
               bmat(nn,nc,l)=0.
               ajmat(nn,nc,l)=0.
            endif
  453    continue
      endif
c
c  *** normalization
c
      do 461 kc=1,nn
         bmat(kc,1,l)=0.5*bmat(kc,1,l)
         bmat(kc,nn,l)=0.5*bmat(kc,nn,l)
         ajmat(kc,1,l)=0.5*ajmat(kc,1,l)
         ajmat(kc,nn,l)=0.5*ajmat(kc,nn,l)
  461 continue
c
c  *** lu-decomposition
c
      call sgefa(bmat(1,1,l),nn,nn,ib(1,l),info)
      if(info .ne. 0) stop '32'
      call sgefa(ajmat(1,1,l),nn,nn,ij(1,l),info)
      if(info .ne. 0) stop '33'
  400 continue
c
c  *** chebycheff collocation, l>0 terms for velocity potentials
c
      do 500 l=1,lmax
c
      lm=l+1
c
      do 542 nc=1,nn
         nd=nc+nn
         do 543 kc=1,nn
            kd=kc+nn
            zmat(kc,nc,l)=(oodt*ql(lm,3)*qk(kc,1)*cheb(nc,kc)-
     $         alpha*ql(lm,12)*qk(kc,1)*
     $         (4.*d2cheb(nc,kc)-
     $         (ql(lm,3)*qk(kc,1))*cheb(nc,kc)))*anorm*zscl
            wpmat(kc,nc,l)=(oodt*ql(lm,3)*qk(kc,1)*cheb(nc,kc)-
     $         alpha*ql(lm,12)*qk(kc,1)*
     $         (4.*d2cheb(nc,kc)-
     $         (ql(lm,3)*qk(kc,1))*cheb(nc,kc)))*
     $         anorm*wscl
            wpmat(kc,nd,l)=alpha*(2.*dcheb(nc,kc))*anorm*pscl
            wpmat(kd,nc,l)=(-2.*oodt*ql(lm,3)*qk(kc,1)*dcheb(nc,kc)-
     $         alpha*ql(lm,12)*qk(kc,1)*
     $         (-8.*d3cheb(nc,kc)+
     $         2.*(ql(lm,3)*qk(kc,1))*dcheb(nc,kc)-
     $         ql(lm,3)*qk(kc,5)*cheb(nc,kc)))*anorm*wscl
            wpmat(kd,nd,l)=-alpha*ql(lm,3)*qk(kc,1)*cheb(nc,kc)*anorm*
     $         pscl
  543    continue
  542 continue
c
c  *** boundary conditions
c
      do 551 nc=1,nnaf
         nd=nc+nn
         wpmat(1,nc,l)=cheb(nc,1)*anorm
         wpmat(1,nd,l)=0.
         wpmat(nn,nc,l)=cheb(nc,nn)*anorm
         wpmat(nn,nd,l)=0.
         if(ktopv.eq.1) then
            zmat(1,nc,l)=(2.*dcheb(nc,1)-
     $         qk(1,6)*cheb(nc,1))*anorm
            wpmat(nnp1,nc,l)=(4.*d2cheb(nc,1)-
     $         2.*dcheb(nc,1)*qk(1,6))*anorm
         else
            zmat(1,nc,l)=cheb(nc,1)*anorm
            wpmat(nnp1,nc,l)=2.*dcheb(nc,1)*anorm
         endif
         wpmat(nnp1,nd,l)=0.
         if(kbotv.eq.1) then
            zmat(nn,nc,l)=(2.*dcheb(nc,nn)-
     $         qk(nn,6)*cheb(nc,nn))*anorm
            wpmat(nnx2,nc,l)=(4.*d2cheb(nc,nn)-
     $         2.*qk(nn,6)*dcheb(nc,nn))*anorm
         else
            zmat(nn,nc,l)=cheb(nc,nn)*anorm
            wpmat(nnx2,nc,l)=2.*dcheb(nc,nn)*anorm
         endif
         wpmat(nnx2,nd,l)=0.
  551 continue
c
      if(nnaf .lt. nn) then
         do 553 nc=nnaf+1,nn
         nd=nc+nn
         wpmat(1,nc,l)=0.
         wpmat(1,nd,l)=0.
         wpmat(nn,nc,l)=0.
         wpmat(nn,nd,l)=0.
         zmat(1,nc,l)=0.
         wpmat(nnp1,nc,l)=0.
         wpmat(nnp1,nd,l)=0.
         zmat(nn,nc,l)=0.
         wpmat(nnx2,nc,l)=0.
         wpmat(nnx2,nd,l)=0.
  553    continue
      endif
c
c  *** normalization
c
      do 561 kc=1,nn
         kd=kc+nn
         zmat(kc,1,l)=0.5*zmat(kc,1,l)
         zmat(kc,nn,l)=0.5*zmat(kc,nn,l)
         wpmat(kc,1,l)=0.5*wpmat(kc,1,l)
         wpmat(kc,nn,l)=0.5*wpmat(kc,nn,l)
         wpmat(kc,nnp1,l)=0.5*wpmat(kc,nnp1,l)
         wpmat(kc,nnx2,l)=0.5*wpmat(kc,nnx2,l)
         wpmat(kd,1,l)=0.5*wpmat(kd,1,l)
         wpmat(kd,nn,l)=0.5*wpmat(kd,nn,l)
         wpmat(kd,nnp1,l)=0.5*wpmat(kd,nnp1,l)
         wpmat(kd,nnx2,l)=0.5*wpmat(kd,nnx2,l)
  561 continue
c
c  *** lu-decomposition of matrices
c
      call sgefa(zmat(1,1,l),nn,nn,iz(1,l),info)
      if(info .ne. 0) stop '34'
      call sgefa(wpmat(1,1,l),nnx2,nnx2,iwp(1,l),info)
      if(info .ne. 0) stop '35'
c
  500 continue
c
      return
      end
