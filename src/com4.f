      common /com4/ ai,tops(0:lmax,0:mmax),bots(0:lmax,0:mmax),
     $qi(ni,5),qk(nn,6),ql(nlma,20),qn(nn,6),
     $aleg1(nlma,ni),aleg2(nlmpa,ni),aleg3(nlma,ni),
     $trigsc(nn),trigsf(ntf),wsave(nn),delxr(nn),delxh2(nn),
     $work(lot,nnp2),bpeak,bpeakbot,bpeaktop,
     $wscl,zscl,pscl,sscl,oosscl,bscl,
     $cmb,pi,y00,p00co,anorm,
     $mclm(nlma),mcalmp(nlmpa),mclma(nlma),
     $ifaxc(13),ifaxf(13),k2k(nn1)
c
      complex ai,tops,bots
