      subroutine fourtf(a,work,trigs,ifax,inc,jump,n,lot,isign)
c
c     same as fft991
c
c     called in amhd
c
c
c purpose      perform a number of simultaneous real/half-complex
c              periodic fast fourier transforms or corresponding inverse
c              transforms, using ordinary spatial order of
c              gridpoint values.  given a set
c              of real data vectors, the package returns a set of
c              "half-complex" fourier coefficient vectors, or vice
c              versa.  the length of the transforms must be an even
c              number that has no other factors except possibly powers
c              of 2, 3, and 5.  this version of fft991 is
c              optimized for use on the cray-1.
c
c argument     a(lot*(n+2)), work(lot*(n+1)), trigs(3*n/2+1), ifax(13)
c dimensions
c
c arguments
c
c on input     a
c               an array of length lot*(n+2) containing the input data
c               or coefficient vectors.  this array is overwritten by
c               the results.
c
c              work
c               a work array of dimension lot*(n+1)
c
c              trigs
c               an array set up by fftfax, which must be called first.
c
c              ifax
c               an array set up by fftfax, which must be called first.
c
c              inc
c               the increment (in words) between successive elements of
c               each data or coefficient vector (e.g.  inc=1 for
c               consecutively stored data).
c
c              jump
c               the increment (in words) between the first elements of
c               successive data or coefficient vectors.  on the cray-1,
c               try to arrange data so that jump is not a multiple of 8
c               (to avoid memory bank conflicts).
c
c              n
c               the length of each transform (see definition of
c               transforms, below).
c
c              lot
c               the number of transforms to be done simultaneously.
c
c              isign
c               = +1 for a transform from fourier coefficients to
c                    gridpoint values.
c               = -1 for a transform from gridpoint values to fourier
c                    coefficients.
c
c on output    a
c               if isign = +1, and lot coefficient vectors are supplied
c               each containing the sequence
c
c               a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)
c
c               then the result consists of lot data vectors each
c               containing the corresponding n+2 gridpoint values
c
c               for fft991, x(0), x(1), x(2),...,x(n-1),0,0.
c                    (n+2) real values with x(n)=x(n+1)=0
c
c               when isign = +1, the transform is defined by
c                 x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
c                 where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c                 and i=sqrt (-1)
c                    for k=0,...,n/2    i.e., (n/2+1) complex values
c                    with c(0) = c(n) = a(0) and c(n/2)=a(n/2)=0
c
c               if isign = -1, and lot data vectors are supplied each
c               containing a sequence of gridpoint values x(j) as
c               defined above, then the result consists of lot vectors
c               each containing the corresponding fourier cofficients
c               a(k), b(k), 0 .le. k .le n/2.
c
c               when isign = -1, the inverse transform is defined by
c                 c(k)=(1/n)*sum(j=0,...,n-1)(x(j)*exp(-2*i*j*k*pi/n))
c                 where c(k)=a(k)+i*b(k) and i=sqrt(-1)
c                 for k=0,...,n/2
c
c               a call with isign=+1 followed by a call with isign=-1
c               (or vice versa) returns the original data.
c
c               note the fact that the gridpoint values x(j) are real
c               implies that b(0)=b(n/2)=0.  for a call with isign=+1,
c               it is not actually necessary to supply these zeros.
c               note starting from grid with x(n)=x(n+1)=0
c               then transforming to spectral (sign=-1)
c               then c(n/2)=a(n/2) is not necessarily 0
c               unless there is no aliasing.
c
      dimension a(*),work(*),trigs(*),ifax(*)
c
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 10
c
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=1
      jbase=1
      do 21 l=1,lot
      i=ibase
      j=jbase
      do 11 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   11 continue
      ibase=ibase+jump
      jbase=jbase+nx
   21 continue
      igo=60
      go to 40
c
c     preprocessing (isign=+1)
c
   10 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     complex transform
c
   40 continue
      ia=1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call wpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call wpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
c
      if (isign.eq.-1) go to 130
c
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=1
      do 100 l=1,lot
      i=ibase
      j=jbase
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
c
c     fill in zeros
c
  110 continue
      ib=n*inc+1
cdir$ ivdep
      do 120 l=1,lot
      a(ib)=0.0
      a(ib+inc)=0.0
      ib=ib+jump
  120 continue
      go to 140
c
c     postprocessing (isign=-1)
c
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 continue
      return
      end
