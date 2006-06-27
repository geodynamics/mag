      subroutine chebtf(lot,ns2,jmp1,jmp2,js2,x,x3,xh,wsave,
     $work,x2,xh1,xh2,trigs,ifax,k2k)
c
c     performs multiple fast chebyshev transforms
c
c     called in amhd, prep, keimei, spectrum
c
c
c     x(l,k) = sqrt(2/nm1) * sum(j=0 to nm1) of x(l,j) * cheb(j,k)
c              with the j=0 and nm1 terms multiplied by 1/2,
c              where cheb(j,k)=cos(pi*j*k/nm1)
c
c     chebtf is the normalized inverse of itself
c
c     sub chebi must be called once to initialize wsave, trigs, ifax, k2k
c              these three and n must not be changed
c
c     n = nm1+1 = 2*ns2+1 = length of the vectors to be tranformed
c              n must be odd and .gt. 3
c              nm1 must be divisible by 4
c
c     lot is the number of vectors to be tranformed
c
c     jmp1 is the first dim of x and is .ge. lot
c
c     jmp2 = (2*js2) is the first dim of x2 and x3 and is .ge. (2*ns2)
c              and should not be a multiple of 8 on the cray
c              but must be a multiple of 2
c
c     real  x(jmp1,jmp2/n),wsave(n),work(jmp1,jmp2+1),trigs(n),ifax(13),k2k(nm1)
c              note x2,xh1,xh2 are equivalenced with w(1,2)
c              and x3 and xh are equivalenced with x
c
c     if array x is complex x(nlm,n) then jmp1=2*nlm
c
      dimension x(jmp1,*),wsave(*),work(jmp1,*),trigs(*),ifax(*),k2k(*),
     $x2(jmp2,jmp1),x3(jmp2,jmp1),xh(jmp1,2,*),xh1(2,js2,jmp1),
     $xh2(jmp1,2,js2)
c
      if(2*js2 .ne. jmp2) stop '40'
      nm1=2*ns2
      n=nm1+1
      np1=n+1
      ns2m1=ns2-1
      ns2p1=ns2+1
      ns2p2=ns2+2
      nq=ns2/2
      nqp1=nq+1
      nfax=ifax(1)
      la=1
      igo=110
      i02=2
      con=sqrt(1./(8.*float(nm1)))
c
      do 10 l=1,lot
      work(l,1)=x(l,1)-x(l,n)
   10 continue
c
      do 102 k=2,ns2
      k1=np1-k
      k2=ns2p1-k
      do 102 l=1,lot
      work(l,1)=work(l,1)+wsave(k2)*(x(l,k)-x(l,k1))
  102 continue
c
      k3=k2k(1)
      k4=k2k(ns2p1)
c
      do 11 l=1,lot
      work(l,1)=work(l,1)+work(l,1)
      x2(k3,l)=x(l,1)+x(l,n)
      x2(k4,l)=x(l,ns2p1)+x(l,ns2p1)
   11 continue
c
      do 103 k=2,ns2
      k1=np1-k
      k2=k-1
      k3=k2k(k)
      k4=k2k(k1)
      do 103 l=1,lot
      x2(k3,l)=x(l,k)+x(l,k1)-wsave(k2)*(x(l,k)-x(l,k1))
      x2(k4,l)=x(l,k)+x(l,k1)+wsave(k2)*(x(l,k)-x(l,k1))
  103 continue
c
      if(mod(nfax,2) .eq. 0) go to 15
      igo=120
c
      do 16 k=1,nm1
      do 16 l=1,lot
      x3(k,l)=x2(k,l)
   16 continue
c
   15 continue
c
      do 140 k=1,nfax
      if(igo .eq. 120) go to 120
c
      call vpassm(x2(1,1),x2(2,1),x3(1,1),x3(2,1),trigs,i02,i02,
     $jmp2,jmp2,lot,ns2,ifax(k+1),la)
c
      igo=120
      go to 130
  120 continue
c
      call vpassm(x3(1,1),x3(2,1),x2(1,1),x2(2,1),trigs,i02,i02,
     $jmp2,jmp2,lot,ns2,ifax(k+1),la)
c
      igo=110
  130 continue
      la=la*ifax(k+1)
  140 continue
c
      do 12 l=1,lot
      xh(l,1,nqp1)=xh1(1,nqp1,l)+xh1(1,nqp1,l)
      xh(l,2,nqp1)=xh1(2,nqp1,l)+xh1(2,nqp1,l)
      xh(l,1,1)=xh1(1,1,l)+xh1(1,1,l)+xh1(2,1,l)+xh1(2,1,l)
      xh(l,2,1)=xh1(1,1,l)+xh1(1,1,l)-xh1(2,1,l)-xh1(2,1,l)
   12 continue
c
      if(nq .lt. 2) go to 107
c
      do 101 k=2,nq
      k1=ns2p2-k
      do 101 l=1,lot
      xh(l,1,k)=xh1(1,k,l)+xh1(1,k1,l)
      xh(l,2,k)=xh1(2,k,l)-xh1(2,k1,l)
      xh(l,1,k1)=xh1(2,k,l)+xh1(2,k1,l)
      xh(l,2,k1)=xh1(1,k1,l)-xh1(1,k,l)
  101 continue
c
      do 100 k=2,nq
      k1=ns2p2-k
      k2=k+ns2m1
      k3=k1+ns2m1
      do 100 l=1,lot
      xh2(l,1,k1)=wsave(k2)*xh(l,1,k1)+wsave(k3)*xh(l,2,k1)
      xh2(l,2,k1)=wsave(k2)*xh(l,2,k1)-wsave(k3)*xh(l,1,k1)
  100 continue
c
      do 104 k=2,nq
      k1=ns2p2-k
      do 104 l=1,lot
      xh(l,1,k1)=xh2(l,1,k1)
      xh(l,2,k1)=xh2(l,2,k1)
  104 continue
c
      do 105 k=2,nq
      k1=ns2p2-k
      do 105 l=1,lot
      xh2(l,1,k1)=xh(l,1,k)-xh(l,1,k1)
      xh2(l,2,k1)=xh(l,2,k)-xh(l,2,k1)
      xh2(l,1,k)=xh(l,1,k)+xh(l,1,k1)
      xh2(l,2,k)=-xh(l,2,k)-xh(l,2,k1)
  105 continue
c
  107 continue
c
      do 13 l=1,lot
      xh2(l,1,1)=xh(l,1,1)
      xh2(l,2,1)=xh(l,2,1)
      xh2(l,1,nqp1)=xh(l,1,nqp1)
      xh2(l,2,nqp1)=xh(l,2,nqp1)
   13 continue
c
      do 14 l=1,lot
cc    x(l,1)=work(l,2)*con
cc    x(l,2)=work(l,1)*con
cc    x(l,n)=work(l,3)*con
c
c  message from Gary: if problem with optimizer (e.g. on T3E),
c  replace the previous 3 lines of code by the following 3 lines
c
      xh(l,1,1)=work(l,2)*con
      xh(l,2,1)=work(l,1)*con
      xh(l,n,1)=work(l,3)*con
   14 continue
c
      do 106 i=4,nm1,2
      i1=i+1
      i2=i-2
      i3=i-1
      do 106 l=1,lot
      x(l,i3)=work(l,i)*con
      x(l,i)=work(l,i1)*con+x(l,i2)
  106 continue
c
      return
      end
