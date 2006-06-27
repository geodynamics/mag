      common /com51/ w(nlma,nnp1),z(nlma,nnp1)
      common /com52/ s(nlma,nnp1),p(nlma,nnp1)
      common /com53/ dw(nlma,nnp1),ddw(nlma,nnp1),dddw(nlma,nnp1)
     $               ,dz(nlma,nnp1),ddz(nlma,nnp1),dp(nlma,nnp1)
     $               ,db(nlma,nnp1),ddb(nlma,nnp1),dj(nlma,nnp1)
      common /com55/ dwdt(nlma,nn,2)
      common /com56/ dzdt(nlma,nn,2)
      common /com57/ dpdt(nlma,nn,2)
      common /com58/ dsdt(nlma,nn,2)
      common /com59/ b(nlma,nnp1),aj(nlma,nnp1)
      common /com512/ dbdt(nlma,nn,2)
      common /com513/ djdt(nlma,nn,2)
c
      complex w,z,s,p,dw,ddw,dddw,dz,ddz,dp,dwdt,dzdt,
     $ dpdt,dsdt
      complex b,aj,db,ddb,dj,dbdt,djdt
c
      complex dwdt1(nlma,nn),dzdt1(nlma,nn),
     $dsdt1(nlma,nn),dpdt1(nlma,nn)
      equivalence (dwdt,dwdt1),(dzdt,dzdt1),(dsdt,dsdt1),(dpdt,dpdt1)
      complex dbdt1(nlma,nn),djdt1(nlma,nn)
      equivalence (dbdt,dbdt1),(djdt,djdt1)
