      function random(r)
c
c     random number generator
c
c     if(r .eq. 0.) then
c        random(r) = next random number (between 0. and 1.)
c     else if(r .lt. 0.) then
c        random(r) = previous random number
c     else if(r .gt. 0.) then
c        random(r) = a new sequence of random numbers is started
c                  with seed r mod 1
c                  note: r must be a non-integer to get a different seq
c     endif
c
c     called in prep
c
      save ia1, ia0, ia1ma0, ic, ix1, ix0
      data ia1, ia0, ia1ma0 /1536, 1029, 507/
      data ic /1731/
      data ix1, ix0 /0, 0/
c
      if (r.lt.0.) go to 10
      if (r.gt.0.) go to 20
c
      iy0 = ia0*ix0
      iy1 = ia1*ix1 + ia1ma0*(ix0-ix1) + iy0
      iy0 = iy0 + ic
      ix0 = mod (iy0, 2048)
      iy1 = iy1 + (iy0-ix0)/2048
      ix1 = mod (iy1, 2048)
c
 10   random = ix1*2048 + ix0
      random = random / 4194304.
      return
c
 20   ix1 = mod(r,1.)*4194304. + 0.5
      ix0 = mod (ix1, 2048)
      ix1 = (ix1-ix0)/2048
      go to 10
c
      end
