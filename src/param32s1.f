      parameter (nn=25,ni=48,nj=096,nnaf=23,minc=1)
c
c This file is an example of a input grid parameter file 
c to be linked to 'param.f' using command 'ln -sf param32s1.f param.f'
c after linking, run the makefile with command 'make'.
c
      parameter (nnp2=nn+2,nnp1=nn+1,nn1=nn-1,nn2=nn-2,nn3=nn-3,
     $nps2=nnp1/2,ns2=nn1/2,nnx2=2*nn,nja=nj/minc,
     $nrp=nja+2,ncp=nrp/2,ntf=3*nja/2+1,njp1=nj+1,nip1=ni+1,
     $lmax=(nj)/3,mmax=(lmax/minc)*minc,nmaf=mmax+1,
     $nmafa=mmax/minc+1,nlm=(nmaf*(nmaf+1))/2,nlaf=lmax+1,
     $nlma=mmax*nlaf/minc-mmax*(mmax-minc)/(2*minc)+nlaf-mmax,
     $lot=2*nlma,nlafp1=nlaf+1,nlmpa=nlma+nmafa)
