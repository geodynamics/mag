      subroutine wpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
c
c     called in chebtf and fourtf
c
      parameter(mmx=96)
      dimension a(*),b(*),c(*),d(*),trigs(*)
      dimension c1(mmx),c2(mmx),c3(mmx),c4(mmx)
      dimension s1(mmx),s2(mmx),s3(mmx),s4(mmx)
      dimension iindex(mmx),jindex(mmx)
      data sin36/0.587785252292473/,cos36/0.809016994374947/,
     *     sin72/0.951056516295154/,cos72/0.309016994374947/,
     *     sin60/0.866025403784437/
c
      m=n/ifac
      if(m.gt.mmx) stop 'wpassm: m>mmx'
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.gt.4) return
      go to (10,50,90,130),igo
c
c     coding for factor 2
c
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
c
      ims=1
      if(la.lt.m.and.la.lt.16) go to 25
      do 20 ijk=1,lot
      iadd=(ijk-1)*inc3
      jadd=(ijk-1)*inc4
      do 15 l=1,la      
      i=(l-1)*inc1+iadd
      j=(l-1)*inc2+jadd
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
   15 continue
   20 continue
c
      if (la.eq.m) return
      ims=la+1  
c
   25 do 26 im=ims,m     
      kc=(im-1)/la
      kk=kc*la
      c1(im)=trigs(2*kk+1)
      s1(im)=trigs(2*kk+2)
      iindex(im)=im*inc1
      jindex(im)=im*inc2+jump*kc
   26 continue
c
      do 30 ijk=1,lot
      iadd=(ijk-1)*inc3 - inc1
      jadd=(ijk-1)*inc4 - inc2
      do 28 im=ims,m
      i= iindex(im)+iadd
      j= jindex(im)+jadd
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1(im)*(a(ia+i)-a(ib+i))-s1(im)*(b(ia+i)-b(ib+i))
      d(jb+j)=s1(im)*(a(ia+i)-a(ib+i))+c1(im)*(b(ia+i)-b(ib+i))
   28 continue
   30 continue
      return
c
c     coding for factor 3
c
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
c
      ims=1
      if(la.lt.m.and.la.lt.16) go to 65
      do 60 ijk=1,lot
      iadd=(ijk-1)*inc3
      jadd=(ijk-1)*inc4
      do 55 l=1,la      
      i=(l-1)*inc1+iadd
      j=(l-1)*inc2+jadd
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
   55 continue
   60 continue
c
      if (la.eq.m) return
      ims=la+1  
c
   65 do 66 im=ims,m     
      kc=(im-1)/la
      kk=kc*la
      c1(im)=trigs(2*kk+1)
      s1(im)=trigs(2*kk+2)
      c2(im)=trigs(4*kk+1)
      s2(im)=trigs(4*kk+2)
      iindex(im)=im*inc1
      jindex(im)=im*inc2+jump*kc
   66 continue
c
      do 70 ijk=1,lot
      iadd=(ijk-1)*inc3 - inc1
      jadd=(ijk-1)*inc4 - inc2
      do 68 im=ims,m
      i= iindex(im)+iadd
      j= jindex(im)+jadd
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=
     *    c1(im)*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
     *           -(sin60*(b(ib+i)-b(ic+i))))
     *   -s1(im)*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))
     *           +(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=
     *    s1(im)*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
     *          -(sin60*(b(ib+i)-b(ic+i))))
     *   +c1(im)*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))
     *          +(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=
     *    c2(im)*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
     *          +(sin60*(b(ib+i)-b(ic+i))))
     *   -s2(im)*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))
     *          -(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=
     *    s2(im)*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))
     *          +(sin60*(b(ib+i)-b(ic+i))))
     *   +c2(im)*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))
     *          -(sin60*(a(ib+i)-a(ic+i))))
   68 continue
   70 continue
      return
c
c     coding for factor 4
c
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
c
      ims=1
      if(la.lt.m.and.la.lt.16) go to 105
      do 100 ijk=1,lot
      iadd=(ijk-1)*inc3
      jadd=(ijk-1)*inc4
      do  95 l=1,la
      i=(l-1)*inc1+iadd
      j=(l-1)*inc2+jadd
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
   95 continue
  100 continue
 
      if (la.eq.m) return
      ims=la+1
c      
  105 do 106 im=ims,m
      kc=(im-1)/la
      kk=kc*la
      c1(im)=trigs(2*kk+1)
      s1(im)=trigs(2*kk+2)
      c2(im)=trigs(4*kk+1)
      s2(im)=trigs(4*kk+2)
      c3(im)=trigs(6*kk+1)
      s3(im)=trigs(6*kk+2)
      iindex(im)=im*inc1
      jindex(im)=im*inc2+jump*kc
  106 continue
c
      do 110 ijk=1,lot
      iadd=(ijk-1)*inc3 - inc1
      jadd=(ijk-1)*inc4 - inc2
      do 108 im=ims,m     
      i= iindex(im)+iadd
      j= jindex(im)+jadd
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=
     *    c2(im)*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   -s2(im)*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=
     *    s2(im)*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   +c2(im)*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=
     *    c1(im)*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   -s1(im)*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=
     *    s1(im)*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   +c1(im)*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=
     *    c3(im)*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   -s3(im)*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=
     *    s3(im)*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   +c3(im)*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
  108 continue
  110 continue
      return
c
c     coding for factor 5
c
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
c
      ims=1
      if(la.lt.m.and.la.lt.16) go to 145
      do 140 ijk=1,lot
      iadd=(ijk-1)*inc3
      jadd=(ijk-1)*inc4
      do 135 l=1,la
      i=(l-1)*inc1+iadd
      j=(l-1)*inc2+jadd
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
  135 continue
  140 continue
c
      if (la.eq.m) return
      ims=la+1
c      
  145 do 146 im=ims,m
      kc=(im-1)/la
      kk=kc*la
      c1(im)=trigs(2*kk+1)
      s1(im)=trigs(2*kk+2)
      c2(im)=trigs(4*kk+1)
      s2(im)=trigs(4*kk+2)
      c3(im)=trigs(6*kk+1)
      s3(im)=trigs(6*kk+2)
      c4(im)=trigs(8*kk+1)
      s4(im)=trigs(8*kk+2)
      iindex(im)=im*inc1
      jindex(im)=im*inc2+jump*kc
  146 continue
c
      do 150 ijk=1,lot
      iadd=(ijk-1)*inc3 - inc1
      jadd=(ijk-1)*inc4 - inc2
      do 148 im=ims,m     
      i= iindex(im)+iadd
      j= jindex(im)+jadd
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=
     *    c1(im)*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))
     *      -cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s1(im)*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))
     *      -cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)=
     *    s1(im)*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))
     *      -cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c1(im)*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))
     *      -cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)=
     *    c4(im)*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))
     *      -cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s4(im)*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))
     *      -cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)=
     *    s4(im)*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))
     *      -cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c4(im)*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))
     *      -cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)=
     *    c2(im)*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))
     *      +cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s2(im)*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))
     *      +cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)=
     *    s2(im)*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))
     *      +cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c2(im)*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))
     *      +cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)=
     *    c3(im)*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))
     *      +cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s3(im)*((b(ia+i)-cos36*(b(ib+i)
     *      +b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)=
     *    s3(im)*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))
     *      +cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c3(im)*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))
     *      +cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
  148 continue
  150 continue
      return
      end
