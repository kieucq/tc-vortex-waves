      program STREAM
c                                                   
c  This program computes the streamfunctions from 
c  the wind field using the relaxation and the 
c  Fourier methods .
c
      parameter (l=21,m=13,np=7,l1=l-1,m1=m-1,l2=l-2,m2=m-2)                                               
      real    datau (l,m,np),datav(l,m,np),zinv(100)                                      
      real    u(l,m),v(l,m),dx(m),z(100)                                               
      real    psi(l,m),a(l,m),work(2*l)                                         
      complex uu(l,m),vv(l,m)                                                  
      open (20,file='uv21.dat',status='old')                      
      open (30,file='psirf.dat',status='unknown')                                 
c
c  read the wind components (1000 to 100 mb) from unit 20.                          
c                                                                               
  878 format(10f8.2)                                                           
      do 4100 ip = 1, np                                                           
         read (20,878) ((datau(i,j,ip),i=1,l),j=1,m)                                
         read (20,878) ((datav(i,j,ip),i=1,l),j=1,m)                                
 4100 continue
c                                                                  
c  select wind field at 500 mb                                                  
c                                                                               
      do 4102 i = 1, l                                                             
      do 4102 j = 1, m                                                             
         u   (i,j) = datau (i,j,4)                                             
         v   (i,j) = datav (i,j,4)                                             
 4102 continue
c                                                                 
c  define the grid spacing and the invariant                                    
c  constants for the domain.                                                   
c                                                                               
      slat      = -15.                                                              
      grid      = 2.5                                                             
      pi        = 4.0*atan(1.0)                                                     
      rad       = pi/180.                                                           
      dy        = 111.1 * 1000. * grid                                              
      do 4104 j = 1, m                                                             
         alat   = (slat + (j-1)*grid)*rad                                        
         dx(j)  = dy * cos(alat)                                                 
4104  continue                                                                 
      do 4106 j = 1, m                                                              
         z(j)   = dx(j)*dx(j)                                                          
         zinv(j)= 1./z(j)                                                           
 4106 continue                                                                  
      zz        = dy*dy                                                                 
      zzinv     = 1./zz
c                                                               
c define the forcing function (relative vorticity)                             
c                                                                               
      do 4108 j = 2, m1                                                           
      do 4110 i = 2, l1                                                         
 4110    a(i,j) = (v(i+1,j)-v(i-1,j))/(2.*dx(j))                              
     &            -(u(i,j+1)-u(i,j-1))/(2.*dy)                                 
         a(1,j) = (v(2,j)-v(l1,j))/(2.*dx(j))                                 
     &            -(u(1,j+1)-u(1,j-1))/(2.*dy)                                 
         a(l,j) = a(1,j)                                                      
 4108 continue                                                                  
      do 4112 i = 1, l                                                              
         a(i,1) = 2.*a(i,2)-a(i,3)                                               
4112     a(i,m) = 2.*a(i,m1)-a(i,m2)                                             
c                                                                               
c  compute the net mass out-flux.the outward velocity                           
c  is corrected to yield a net outward mass flux.                               
c                                                                               
c  vno  is the integral mass flux.                                              
c  uno  is the integral of the magnitude of mass flux.                          
c                                                                               
      vno       = v(1,m)*dx(m)/2.+v(l,m)*dx(m)/2.                                 
      uno       = abs(v(1,m))*dx(m)/2.+abs(v(l,m))*dx(m)/2.                       
c                                                                               
      do 4114 i = 2, l1                                                         
         uno    = uno + abs (v(i,m))*dx(m)                                        
 4114    vno    = vno + v(i,m)*dx(m)                                              
      vno       = vno + u(l,m)*dy/2. + u(l,1)*dy/2.                               
      uno       = uno + abs(u(l,m))*dy/2.+abs(u(l,1))*dy/2.                       
c                                                                               
      do 4116 j = 2, m1                                                           
         uno    = uno + abs(u(l,j))*dy                                            
 4116    vno    = vno + u(l,j)*dy                                                 
                                                                                
      vno       = vno - v(l,1)*dx(1)/2.-v(1,1)*dx(1)/2.                           
      uno       = uno + abs(v(l,1))*dx(1)/2.+abs(v(1,1))*dx(1)/2.                 
c                                                                               
      do 4118 i = 2, l1                                                           
         uno    = uno + abs(v(i,1))*dx(1)                                         
 4118    vno    = vno - v(i,1)*dx(1)                                              
         vno    = vno - u(1,1)*dy/2.-u(1,m)*dy/2.                                 
         uno    = uno + abs(u(1,1))*dy/2.+abs(u(1,m))*dy/2.                       
c                                                                               
c  computation of the correction factor epsilon
c                                                       
      do 4120 j = 2, m1                                                           
         uno    = uno + abs(u(1,j))*dy                                            
 4120    vno    = vno - u(1,j)*dy                                                 
         eps    = vno/uno                                                         
      write(6,798)                                                              
      write(6,799) uno,vno,eps                                                  
  798 format(2x,'uno ,vno , eps')                                               
  799 format(6x,3e14.5)                                                         
c                                                                               
c  correction of the outward normal velocity.
c                                   
      do 4122 i = 1, l                                                             
         v(i,1) = v(i,1) + eps*abs(v(i,1))                                      
 4122    v(i,m) = v(i,m) - eps*abs(v(i,m))                                      
c                                                                               
      do 4124 j = 1, m                                                             
         u(1,j) = u(1,j) + eps*abs(u(1,j))                                      
 4124    u(l,j) = u(l,j) - eps*abs(u(l,j))                                      
c                                                                               
c assume psi(1,m) is known and compute the remaining                           
c boundary values using the corrected outward normal                           
c velocitiy.                                                                   
c                                                                               
      psi(1,m)  = 0.                                                               
      do 4126 i = 2, l                                                             
 4126  psi(i,m) = psi(i-1,m) + (v(i,m)+v(i-1,m))*dx(m)/2.                      
      do 4128 jj= 1, m1                                                           
         j      = m-jj                                                                    
 4128  psi(l,j) = psi(l,j+1) + (u(l,j)+u(l,j+1))*dy/2.                         
      do 4130 ii= 1, l1                                                           
         i      = l-ii                                                                    
 4130  psi(i,1) = psi(i+1,1) - (v(i,1)+v(i+1,1))*dx(1)/2.                      
      do 4132 j = 2, m1                                                            
 4132  psi(1,j) = psi(1,j-1) - (u(1,j)+u(1,j-1))*dy/2.                         
c                                                                               
c  solve the poisson equation using the relaxation technique.               
c  the tolerance factor is set to 1000.                                         
c                                                                               
      call RELAX (psi,zzinv,z,a,zinv,l,l1,m1,m)                                 
c                                                                               
c normalize the streamfunction ,and write output to tape 4.                    
c
      do 4134 j = 1, m                                                       
      do 4134 i = 1, l
 4134  psi(i,j) = psi(i,j)                                   
      write (30,222)psi   
 222  format(10e13.6)                                                           
c                                                                               
c  compute the steamfunction via fourier expansions.                            
c                                                                               
      do 4136 j = 1, m                                                            
      do 4136 i = 1, l                                                            
        uu(i,j) = cmplx (u(i,j),0.0)                                     
        vv(i,j) = cmplx (v(i,j),0.0)                                     
 4136 continue
c                                                                                                                                                 
      call PSICHI (uu,vv,l,m,dx,dy,work,1,1)                                    
c                                                                               
      do 4138 j = 1, m                                                            
      do 4138 i = 1, l                                                            
       psi(i,j) = real (uu(i,j))                                       
 4138 continue
c
c  write output .
c
              write (30,222)psi                                                 
c                                                                                
      stop                                                                      
      end
                                                                       
      subroutine FOURT (data,nn,ndim,isign,iform,work)
c
c  This subroutines performs the fast Fourier transform (FFT)
c                           
      dimension data(1),nn(1),ifact(32),work(1)                                 
      real *8 twopi,rthlf
      data twopi/6.2831853071796/,rthlf/0.70710678118655/                       
      if(ndim-1)920,1,1                                                         
1     ntot      = 2                                                                    
      do 2 idim = 1, ndim                                                          
      if(nn(idim))920,920,2                                                     
2     ntot=ntot*nn(idim)                                                        
c                                                                               
c     main loop for each dimension                                              
c                                                                               
      np1       = 2                                                                     
      do 910 idim = 1, ndim                                                        
         n        = nn(idim)                                                                
         np2      = np1*n                                                                 
         if (n-1) 920,900,5                                                          
c                                                                               
c     is n a power of two and if not, what are its factors                      
c                                                                               
5     m           = n                                                                       
      ntwo        = np1                                                                  
      if          = 1                                                                      
      idiv        = 2                                                                    
10    iquot       = m/idiv                                                              
      irem        = m-idiv*iquot                                                         
      if (iquot-idiv) 50,11,11                                                    
11    if (irem) 20,12,20                                                          
12    ntwo        = ntwo+ntwo                                                            
      ifact(if)   = idiv                                                            
      if          = if+1                                                                   
      m           = iquot                                                                   
      go to 10                                                                  
20    idiv        = 3                                                                    
      inon2       = if                                                                   
30    iquot       = m/idiv                                                               
      irem        = m-idiv*iquot                                                         
      if (iquot-idiv) 60,31,31                                                    
31    if (irem) 40,32,40                                                          
32    ifact(if)   = idiv                                                            
      if          = if+1                                                                   
      m           = iquot                                                                   
      go to 30                                                                  
40    idiv        = idiv+2                                                               
      go to 30                                                                  
50    inon2       = if                                                                  
      if (irem) 60,51,60                                                          
51    ntwo        = ntwo+ntwo                                                            
      go to 70                                                                  
60    ifact(if)   = m                                                               
c                                                                               
c     separate four cases--                                                     
c        1. complex transform or real transform for the 4th, 9th,etc.           
c           dimensions.                                                               
c        2. real transform for the 2nd or 3rd dimension.  method--              
c           transform half the data, supplying the other half by con-           
c           jugate symmetry.                                                    
c        3. real transform for the 1st dimension, n odd.  method--              
c           set the imaginary parts to zero.                                    
c        4. real transform for the 1st dimension, n even.  method--             
c           transform a complex array of length n/2 whose real parts            
c           are the even numbered real values and whose imaginary parts         
c           are the odd numbered real values.  separate and supply              
c           the second half by conjugate symmetry.                              
c                                                                               
70    icase     = 1                                                                   
      ifmin     = 1                                                                   
      i1rng     = np1                                                                 
      if (idim-4) 7 1,100,100                                                      
71    if (iform) 72,72,100                                                        
72    icase     = 2                                                                   
      i1rng     = npo*(1+nprev/2)                                                     
      if (idim-1) 73,73,100                                                       
73    icase     = 3                                                                   
      i1rng     = np1                                                                 
      if (ntwo-np1) 100,100,74                                                    
74    icase     = 4                                                                   
      ifmin     = 2                                                                   
      ntwo      = ntwo/2                                                               
      n         = n/2                                                                      
      np2       = np2/2                                                                 
      ntot      = ntot/2                                                               
      i         = 1                                                                       
      do 80 j   = 1, ntot                                                            
      data(j)   = data(i)                                                           
80    i         = i+2                                                                     
c                                                                               
c     shuffle data by bit reversal, since n=2**k.  as the shuffling             
c     can be done by simple interchange, no working array is needed             
c                                                                               
100   if (ntwo-np2) 200,110,110                                                   
110   np2hf     = np2/2                                                               
      j         = 1                                                                       
      do 150 i2 = 1, np2, np1                                                       
         if (j-i2) 120,130,130                                                       
120   i1max     = i2+np1-2                                                            
      do 125 i1 = i2, i1max, 2                                                      
      do 125 i3 = i1, ntot, np2                                                     
         j3     = j+i3-i2                                                                
          tempr = data(i3)                                                            
          tempi = data(i3+1)                                                          
       data(i3) = data(j3)                                                         
      data(i3+1)= data(j3+1)                                                     
       data(j3) = tempr                                                            
125   data(j3+1)= tempi                                                          
130   m         = np2hf                                                                   
140   if (j-m) 150,150,145                                                        
145   j         = j-m                                                                     
      m         = m/2                                                                     
      if (m-np1) 150,140,140                                                      
150   j         = j+m                                                                     
      go to 300                                                                 
c                                                                               
c     shuffle data by digit reversal for general n                              
c                                                                               
200   nwork     = 2*n                                                                 
      do 270 i1 = 1, np1, 2                                                         
      do 270 i3 = i1,ntot, np2                                                     
         j      = i3                                                                      
      do 260 i  = 1, nwork, 2                                                        
         if (icase-3) 210,220,210                                                    
210      work(i)= data(j)                                                           
      work(i+1) = data(j+1)                                                       
      go to 230                                                                 
220   work(i)   = data(j)                                                           
      work(i+1) = 0.                                                              
230   ifp2      = np2                                                                   
      if        = ifmin                                                                  
240   ifp1      = ifp2/ifact(if)                                                       
      j         = j+ifp1                                                                  
      if (j-i3-ifp2) 260,250,250                                                  
250   j         = j-ifp2                                                                  
      ifp2      = ifp1                                                                 
      if        = if+1                                                                   
      if (ifp2-np1) 260,260,240                                                   
260   continue                                                                  
      i2max     = i3+np2-np1                                                          
      i         = 1                                                                       
      do 270 i2 = i3,i2max,np1                                                    
      data(i2)  = work(i)                                                          
      data(i2+1)= work(i+1)                                                      
270   i         = i+2                                                                     
c                                                                               
c     main loop for factors of two.  perform fourier transforms of              
c     length four, with one of length two if needed.  the twiddle factor        
c     w = exp(isign*2*pi*sqrt(-1)*m/(4*mmax)).  check for w=isign*sqrt(-1)        
c     and repeat for w=w*(1+isign*sqrt(-1))/sqrt(2).                            
c                                                                               
300   if (ntwo-np1) 600,600,305                                                   
305   np1tw     = np1+np1                                                             
      ipar      = ntwo/np1                                                             
310   if (ipar-2) 350,330,320                                                     
320   ipar      = ipar/4                                                               
      go to 310                                                                 
330   do 340 i1 = 1,i1rng,2                                                       
      do 340 k1 = i1,ntot,np1tw                                                   
         k2     = k1+np1                                                                 
      tempr     = data(k2)                                                            
      tempi     = data(k2+1)                                                          
      data(k2)  = data(k1)-tempr                                                   
      data(k2+1)= data(k1+1)-tempi                                               
      data(k1)  = data(k1)+tempr                                                   
340   data(k1+1)= data(k1+1)+tempi                                               
350   mmax      = np1                                                                  
360   if (mmax-ntwo/2) 370,600,600                                                
370   lmax      = max0(np1tw,mmax/2)                                                   
      do 570 l  = np1, lmax, np1tw                                                   
      m         = l                                                                       
      if (mmax-np1) 420,420,380                                                   
380   theta     = -twopi*float(l)/float(4*mmax)                                       
      if (isign) 400,390,390                                                      
390   theta     = -theta                                                              
400   wr        = cos(theta)                                                             
      wi        = sin(theta)                                                             
410   w2r       = wr*wr-wi*wi                                                           
      w2i       = 2.*wr*wi                                                              
      w3r       = w2r*wr-w2i*wi                                                         
      w3i       = w2r*wi+w2i*wr                                                         
420   do 530 i1 = 1,i1rng,2                                                       
      kmin      = i1+ipar*m                                                            
      if (mmax-np1) 430,430,440                                                   
430   kmin      = i1                                                                   
440   kdif      = ipar*mmax                                                            
450   kstep     = 4*kdif                                                              
      if (kstep-ntwo) 460,460,530                                                 
460   do 520 k1 = kmin, ntot, kstep                                                 
      k2        = k1+kdif                                                                
      k3        = k2+kdif                                                                
      k4        = k3+kdif                                                                
      if (mmax-np1) 470,470,480                                                   
470   u1r       = data(k1)+data(k2)                                                     
      u1i       = data(k1+1)+data(k2+1)                                                 
      u2r       = data(k3)+data(k4)                                                     
      u2i       = data(k3+1)+data(k4+1)                                                 
      u3r       = data(k1)-data(k2)                                                     
      u3i       = data(k1+1)-data(k2+1)                                                 
      if (isign) 471,472,472                                                      
471   u4r       = data(k3+1)-data(k4+1)                                                 
      u4i       = data(k4)-data(k3)                                                     
      go to 510                                                                 
472   u4r       = data(k4+1)-data(k3+1)                                                 
      u4i       = data(k3)-data(k4)                                                     
      go to 510                                                                 
480   t2r       = w2r*data(k2)-w2i*data(k2+1)                                           
      t2i       = w2r*data(k2+1)+w2i*data(k2)                                           
      t3r       = wr*data(k3)-wi*data(k3+1)                                             
      t3i       = wr*data(k3+1)+wi*data(k3)                                             
      t4r       = w3r*data(k4)-w3i*data(k4+1)                                           
      t4i       = w3r*data(k4+1)+w3i*data(k4)                                           
      u1r       = data(k1)+t2r                                                          
      u1i       = data(k1+1)+t2i                                                        
      u2r       = t3r+t4r                                                               
      u2i       = t3i+t4i                                                               
      u3r       = data(k1)-t2r                                                          
      u3i       = data(k1+1)-t2i                                                        
      if (isign) 490,500,500                                                      
490   u4r       = t3i-t4i                                                               
      u4i       = t4r-t3r                                                               
      go to 510                                                                 
500   u4r       = t4i-t3i                                                               
      u4i       = t3r-t4r                                                               
510   data(k1)  = u1r+u2r                                                          
      data(k1+1)= u1i+u2i                                                        
      data(k2)  = u3r+u4r                                                          
      data(k2+1)= u3i+u4i                                                        
      data(k3)  = u1r-u2r                                                          
      data(k3+1)= u1i-u2i                                                        
      data(k4)  = u3r-u4r                                                          
520   data(k4+1)= u3i-u4i                                                        
      kdif      = kstep                                                                
      kmin      = 4*(kmin-i1)+i1                                                       
      go to 450                                                                 
530   continue                                                                  
      m         = m+lmax                                                                  
      if (m-mmax) 540,540,570                                                     
540   if(isign)550,560,560                                                      
550   tempr     = wr                                                                  
      wr        = (wr+wi)*rthlf                                                          
      wi        = (wi-tempr)*rthlf                                                       
      go to 410                                                                 
560   tempr     = wr                                                                  
      wr        = (wr-wi)*rthlf                                                          
      wi        = (tempr+wi)*rthlf                                                       
      go to 410                                                                 
570   continue                                                                  
      ipar      = 3-ipar                                                               
      mmax      = mmax+mmax                                                            
      go to 360                                                                 
c                                                                               
c     main loop for factors not equal to two.  apply the twiddle factor         
c     w=exp(isign*2*pi*sqrt(-1)*(j1-1)*(j2-j1)/(ifp1+ifp2)),then                
c     perform a fourier transform of length ifact(if), making use of            
c     conjugate symmetries.                                                     
c                                                                               
600   if (ntwo-np2) 605,700,700                                                   
605   ifp1      = ntwo                                                                 
      if        = inon2                                                                  
      np1hf     = np1/2                                                               
610   ifp2      = ifact(if)*ifp1                                                       
      j1min     = np1+1                                                               
      if (j1min-ifp1) 615,615,640                                                 
615   do 635 j1 = j1min, ifp1, np1                                                  
      theta     = -twopi*float(j1-1)/float(ifp2)                                      
      if (isign) 625,620,620                                                      
620   theta     = -theta                                                              
625   wstpr     = cos(theta)                                                          
      wstpi     = sin(theta)                                                          
      wr        = wstpr                                                                  
      wi        = wstpi                                                                    
      j2min     = j1+ifp1                                                             
      j2max     = j1+ifp2-ifp1                                                        
      do 635 j2 = j2min, j2max, ifp1                                                
      i1max     = j2+i1rng-2                                                          
      do 630 i1 = j2, i1max, 2                                                      
      do 630 j3 = i1, ntot, ifp2                                                    
      tempr     = data(j3)                                                            
      data(j3)  = data(j3)*wr-data(j3+1)*wi                                        
630   data(j3+1)= tempr*wi+data(j3+1)*wr                                         
      tempr     = wr                                                                  
      wr        = wr*wstpr-wi*wstpi                                                      
635   wi        = tempr*wstpi+wi*wstpr                                                   
640   theta     = -twopi/float(ifact(if))                                             
      if (isign) 650,645,645                                                      
645   theta     = -theta                                                              
650   wstpr     = cos(theta)                                                          
      wstpi     = sin(theta)                                                          
      j2rng     = ifp1*(1+ifact(if)/2)                                                
      do 695 i1 = 1, i1rng, 2                                                       
      do 695 i3 = i1, ntot, np2                                                     
      j2max     = i3+j2rng-ifp1                                                       
      do 690 j2 = i3,j2max,ifp1                                                   
      j1max     = j2+ifp1-np1                                                         
      do 680 j1 = j2, j1max, np1                                                    
      j3max     = j1+np2-ifp2                                                         
      do 680 j3 = j1,j3max,ifp2                                                   
      jmin      = j3-j2+i3                                                             
      jmax      = jmin+ifp2-ifp1                                                       
      i         = 1+(j3-i3)/np1hf                                                         
      if (j2-i3) 655,655,665                                                      
655   sumr      = 0.                                                                   
      sumi      = 0.                                                                   
      do 660 j  = jmin, jmax, ifp1                                                   
      sumr      = sumr+data(j)                                                         
660   sumi      = sumi+data(j+1)                                                       
      work(i)   = sumr                                                              
      work(i+1) = sumi                                                            
      go to 680                                                                 
665   iconj     = 1+(ifp2-2*j2+i3+j3)/np1hf                                           
      j         = jmax                                                                    
      sumr      = data(j)                                                              
      sumi      = data(j+1)                                                            
      oldsr     = 0.                                                                  
      oldsi     = 0.                                                                  
      j         = j-ifp1                                                                  
670   tempr     = sumr                                                                
      tempi     = sumi                                                                
      sumr      = twowr*sumr-oldsr+data(j)                                             
      sumi      = twowr*sumi-oldsi+data(j+1)                                           
      oldsr     = tempr                                                               
      oldsi     = tempi                                                               
      j         = j-ifp1                                                                  
      if (j-jmin) 675,675,670                                                     
675   tempr     = wr*sumr-oldsr+data(j)                                               
      tempi     = wi*sumi                                                             
      work(i)   = tempr-tempi                                                       
      work(iconj)=tempr+tempi                                                   
      tempr     = wr*sumi-oldsi+data(j+1)                                             
      tempi     = wi*sumr                                                             
      work(i+1) = tempr+tempi                                                     
      work(iconj+1)=tempr-tempi                                                 
680   continue                                                                  
      if (j2-i3) 685,685,686                                                      
685   wr        = wstpr                                                                  
      wi        = wstpi                                                                  
      go to 690                                                                 
686   tempr     = wr                                                                  
      wr        = wr*wstpr-wi*wstpi                                                      
      wi        = tempr*wstpi+wi*wstpr                                                   
690   twowr     = wr+wr                                                               
      i         = 1                                                                       
      i2max     = i3+np2-np1                                                          
      do 695 i2 = i3,i2max,np1                                                    
      data(i2)  = work(i)                                                          
      data(i2+1)= work(i+1)                                                      
695   i         = i+2                                                                     
      if        = if+1                                                                   
      ifp1      = ifp2                                                                 
      if (ifp1-np2) 610,700,700                                                   
c                                                                               
c     complete a real transform in the 1st dimension, n even, by con-           
c     jugate symmetries.                                                        
c                                                                               
700   go to (900,800,900,701),icase                                             
701   nhalf     = n                                                                   
      n         = n+n                                                                     
      theta     = -twopi/float(n)                                                     
      if (isign) 703,702,702                                                      
702   theta     = -theta                                                              
703   wstpr     = cos(theta)                                                          
      wstpi     = sin(theta)                                                          
      wr        = wstpr                                                                  
      wi        = wstpi                                                                  
      imin      = 3                                                                    
      jmin      = 2*nhalf-1                                                            
      go to 725                                                                 
710   j         = jmin                                                                    
      do 720 i  = imin, ntot, np2                                                    
      sumr      = (data(i)+data(j))/2.                                                 
      sumi      = (data(i+1)+data(j+1))/2.                                             
      difr      = (data(i)-data(j))/2.                                                 
      difi      = (data(i+1)-data(j+1))/2.                                             
      tempr     = wr*sumi+wi*difr                                                     
      tempi     = wi*sumi-wr*difr                                                     
      data(i)   = sumr+tempr                                                        
      data(i+1) = difi+tempi                                                      
      data(j)   = sumr-tempr                                                        
      data(j+1) = -difi+tempi                                                     
720   j         = j+np2                                                                   
      imin      = imin+2                                                               
      jmin      = jmin-2                                                               
      tempr     = wr                                                                  
      wr        = wr*wstpr-wi*wstpi                                                      
      wi        = tempr*wstpi+wi*wstpr                                                   
725   if (imin-jmin) 710,730,740                                                  
730   if (isign) 731,740,740                                                      
731   do 735 i  = imin, ntot, np2                                                    
735   data(i+1) = -data(i+1)                                                      
740   np2       = np2+np2                                                               
      ntot      = ntot+ntot                                                            
      j         = ntot+1                                                                  
      imax      = ntot/2+1                                                             
745   imin      = imax-2*nhalf                                                         
      i         = imin                                                                    
      go to 755                                                                 
750   data(j)   = data(i)                                                           
      data(j+1) = -data(i+1)                                                      
755   i         = i+2                                                                     
      j         = j-2                                                                     
      if (i-imax) 750,760,760                                                     
760   data(j)   = data(imin)-data(imin+1)                                           
      data(j+1) = 0.                                                              
      if (i-j) 770,780,780                                                        
765   data(j)   = data(i)                                                           
      data(j+1) = data(i+1)                                                       
770   i         = i-2                                                                     
      j         = j-2                                                                     
      if (i-imin) 775,775,765                                                     
775   data(j)   = data(imin)+data(imin+1)                                           
      data(j+1) = 0.                                                              
      imax      = imin                                                                 
      go to 745                                                                 
780   data(1)   = data(1)+data(2)                                                   
      data(2)   = 0.                                                                
      go to 900                                                                 
c                                                                               
c     complete a real transform for the 2nd or 3rd dimension by                 
c     conjugate symmetries.                                                     
c                                                                               
800   if (i1rng-np1)805,900,900                                                 
805   do 860 i3 = 1, ntot, np2                                                      
      i2max=i3+np2-np1                                                          
      do 860 i2 = i3, i2max,np1                                                    
      imin      = i2+i1rng                                                             
      imax      = i2+np1-2                                                             
      jmax      = 2*i3+np1-imin                                                        
      if (i2-i3) 820,820,810                                                      
810   jmax      = jmax+np2                                                             
820   if (idim-2) 850,850,830                                                     
830   j         = jmax+npo                                                                
      do 840 i  = imin,imax,2                                                      
      data(i)   = data(j)                                                           
      data(i+1) = -data(j+1)                                                      
840   j         = j-2                                                                     
850   j         = jmax                                                                    
      do 860 i  = imin,imax,npo                                                    
      data(i)   = data(j)                                                           
      data(i+1) = -data(j+1)                                                      
860   j         = j-npo                                                                   
c                                                                               
c     end of loop on each dimension                                             
c                                                                               
900   npo       = np1                                                                   
      np1       = np2                                                                   
910   nprev     = n                                                                   
920   return                                                                    
      end
      subroutine RELAX (x,zzinv,z,y,zinv,l,l1,m1,m)                             
c                                                                               
c  This subroutine solve the Poisson equation using 
c  the overrelaxation method .
c                                                   
c  definition of variables :                                                     
c                                                                               
c x         contains first guess and eventually                                 
c           the output of the field to be relaxed                                   
c y         forcing function of l by m matrix                                   
c z         latitudinal increment square                                        
c zinv      inverse of dx**2                                                    
c zzinv     inverse of dy**2                                                    
c l         east-west dimension                                                 
c m         north-south dimension                                               
c                                                                               
c  definition of constants :                                                                   
c npts      number of points to be relaxed                            
c nrel      number of points relaxed                                  
c alfa      relaxation factor                                         
c ia        maximum number of iterations allowed                      
c eps       tolerance error                                           
c nsc       count of number of scan for convergence                   
c lsc       last scan after convergence                               
c                                                                               
      real x(l,m),y(l,m),z(m),zinv(m)                                                   
      npts      = l*(m-2)                                                            
      nlax      = 1                                                                  
      mm        = 2                                                                    
      mmm       = m1                                                                  
      alfa      = .46                                                                
      ia        = 1000                                                                 
      eps       = 1.                                                                  
      nsc       = 0                                                                   
      lsc       = -1                                                                  
   15 nrel      = 0                                                                  
      do 4110 j = mm, mmm                                                             
      do 4110 i = 1, l                                                                
         im1    = i-1                                                                 
         ip1    = i+1                                                                 
         if (im1.lt.1) im1 = l1                                                     
         if (ip1.gt.l) ip1 = 2                                                      
         r      = (x(ip1,j)+x(im1,j)-2.*x(i,j))*zinv(j)+
     &            (x(i,j+1)+x(i,j-1)-2.*x(i,j))*zzinv                                                            
         r      = (r-y(i,j))*z(j)                                                       
         if (lsc-nsc) 29,29,30                                                     
   29    x(i,j) = x(i,j) + alfa*r                                                  
   30    if (abs(r).le.eps) nrel = nrel+1                                           
 4110 continue                                                                  
c                                                                               
c nrel gives the number of points that have converged .
c nsc gives the number of scan made over the domain and
c if it is less then the maximum number of iteration 
c allowed for complete convergence it keeps going to                                        
c statement number 15. lsc is the final check once 
c convergence has been achieved it does one more loop 
c before finally jumping out to statement 300 .                                                         
c                                                                               
      nsc       = nsc+1                                                               
      if (nrel-npts) 13,14,14                                                   
   14 if (lsc .ge. nsc) go to 300                                               
   18 lsc       = nsc+1                                                               
   13 if (nsc.lt.ia) go to 15                                                   
  201 format(50h   progress of relaxation npts,nrel,nsc,ia        )             
  300 continue                                                                  
      write(6,201)                                                              
  200 format(6x,4i9)                                                            
      write (6,200)npts,nrel,nsc,ia                                             
      return                                                                    
      end                                                                       
      subroutine PSICHI (u,v,l,m,dx,dy,work,if1,if2)                             
c                                                                               
c  This subroutine is composed by two parts,                                   
c  psichi and uuvv .uuvv does not stand alone and 
c  can be used only after psichi has been called .                                                   
c  one external routine needed is the fft,subroutine fourt.                 
c                                                                               
c  subroutine psichi..computes (psi or zta) and (chi or div)                            
c  psi    : streamfunction                                                    
c  zta    : vorticity                                                         
c  chi    : velocity potential                                                
c  div    : divergence                                                        
c                                                                               
c  calling sequence :                                                        
c                                                                               
c  call PSICHI (u,v,l,m,dx,dy,work,if1,if2)                                   
c                                                                               
c  u      : complex array dimensioned (l,m) on input, 
c           the real part contains zonal wind component            
c           on output, the real part contains psi or zta                     
c  v      : complex array dimensioned (l,m) on input, 
c           the real part contains meridional wind component                                        
c           on output, the real part contains chi or div                     
c
c  variables definition :
c  l      : number of grid point in zonal direction                          
c  m      : number of grid point in meridional direction                     
c  dx     : grid distance in zonal direction, dimensioned m                  
c  dy     : grid distance in meridional direction                            
c  work   : working array needed for fourt, dimensioned 2*                   
c           2*(larger value of l and m that is not a power of 2)               
c  if1 = 1: compute psi. 0 : compute zta.                                
c  if2 = 1: compute chi. 0 : compute div.                                
c                                                                               
c  this routine is consistent with the following                             
c  finite difference scheme                                                  
c                                                                               
c  psi to zta         : five point scheme with one dx and one dy               
c  chi to div         : same as above                                          
c  u,v to zta         : forward difference scheme                              
c  u,v to div         : backward scheme                                        
c  psi to upsi,vpsi   : forward scheme                                         
c  chi to uchi,vchi   : backward scheme                                        
c                                                                               
c                                                                               
c  this routine works only with non-cyclic domain.for cyclic domain
c  user should create the proper array by removing one of the two 
c  cyclic boundaries before calling this routine.                                                             
c  one assumption is made in the computation of psi and chi:                 
c   ubar(area averaged u) is contained solely in psi, and                  
c   vbar(area averaged v) is contained solely in chi.                      
c                                                                               
c  subroutine UUVV computes u and v from (psi or zta) and (chi or div)                 
c  which has been computed through psichi. when psichi has been called
c  more than once,two important data are lost, namely ubar and vbar, 
c  which are needed for the reconstruction of the wind field.
c  so if user needed to retrieve wind field at a later time,he/she should 
c  use the labeled common /uvbar/ ubar, vbar in the main program and save
c  the ubar and vbar in anotherarray and set the proper ubar and vbar 
c  before calling UUVV.               
c                                                                               
c  calling sequence :                                                        
c                                                                               
c  call UUVV (u,v,l,m,dx,dy,work,if1,if2)                                     
c                                                                               
c  u      : complex array dimensioned (l,m)                                  
c           on input, the real part contains psi or zta                      
c           on output, the real part contains u                              
c  v      : complex array dimensioned (l,m)                                  
c           on input, the real part contains chi or div                      
c           on output, the real part contains v                              
c  l      : same as in psichi                                                
c  m      : same as in psichi                                                
c  dx     : same as in in psichi                                             
c  dy     : same as in psichi                                                
c  work   : same as in psichi                                                
c  if1 = 1: input is psi. 0 : input is zta.                              
c  if2 = 1: input is chi. 0 : input is div.                              
c                                                                               
      dimension dx(1), work(1), nn(2)                                           
      complex u(l,m), v(l,m)                                                    
      complex temp, xi, xj                                                      
      common /uvbar/ ubar, vbar                                                 
      data pi/3.1415926535898/                                                  
      nn(1)     = l                                                                 
      nn(2)     = m                                                                 
      fliv      = 1. / float(l)                                                      
      fmiv      = 1. / float(m)                                                      
      flmiv     = fliv * fmiv                                                       
      twopibl   = 2. * pi * fliv                                                  
      twopibm   = 2. * pi * fmiv                                                  
      dyiv      = 1. / dy                                                            
      dysqiv    = dyiv * dyiv                                                      
      call FOURT (u,nn,2,-1,0,work)                                              
      call FOURT (v,nn,2,-1,0,work)                                              
      do 4120 j = 1, m                                                            
      do 4120 i = 1, l                                                            
         u(i,j) = u(i,j) * flmiv                                                   
         v(i,j) = v(i,j) * flmiv                                                   
 4120 continue                                                                  
      ubar      = real(u(1,1))                                                       
      vbar      = real(v(1,1))                                                       
      u(1,1)    = (0.,0.)                                                          
      v(1,1)    = (0.,0.)                                                          
      do 4121 j = 1, m                                                            
         xcj    = cos(twopibm * (j - 1))                                              
         xsj    = sin(twopibm * (j - 1))                                              
         xcj1   = 1. - xcj                                                           
         xj     = cmplx(xcj1,xsj)                                                      
         dxiv   = 1. / dx(j)                                                         
         dxsqiv = dxiv * dxiv                                                      
      do 4121 i = 1, l                                                            
      if (i.eq.1.and.j.eq.1) go to 4121                                            
          xci   = cos(twopibl * (i - 1))                                              
          xsi   = sin(twopibl * (i - 1))                                              
          xci1  = 1. - xci                                                           
          xi    = cmplx(xci1,xsi)                                                      
          term  = -2. * xci1 * dxsqiv - 2. * xcj1 * dysqiv                           
      termivx   = 1. / term                                                        
      termivy   = termivx                                                         
      if (if1.eq.0) termivx = 1.                                                 
      if (if2.eq.0) termivy = 1.                                                 
         temp   = (v(i,j)*xi*dxiv-u(i,j)*xj*dyiv)*termivx                            
         v(i,j) = (u(i,j)*conjg(xi)*dxiv+v(i,j)*conjg(xj)*dyiv)*termivy            
         u(i,j) = temp                                                             
      if (if2.eq.0) v(i,j) = -v(i,j)                                             
 4121 continue                                                                  
      call FOURT (u,nn,2,1,1,work)                                               
      call FOURT (v,nn,2,1,1,work)                                               
   30 if (if1.eq.0) go to 50                                                     
      do 4122 j = 1, m                                                            
      do 4122 i = 1, l                                                            
         u(i,j) = u(i,j) - (j - 1) * ubar * dy                                     
 4122 continue                                                                  
   50 if (if2.eq.0) go to 70                                                     
      do 4123 j = 1, m                                                            
      do 4123 i = 1, l                                                            
         v(i,j) = v(i,j) - (j - 1) * vbar * dy                                     
 4123 continue                                                                  
   70 continue                                                                  
      return                                                                    
      entry uuvv                                                                
      if (if1.eq.1) go to 100                                                    
      call FOURT (u,nn,2,-1,0,work)                                              
      do 4130 j = 1, m                                                            
         xcj1   = 1. - cos(twopibm * (j - 1))                                        
         dxiv   = 1. / dx(j)                                                         
         dxsqiv = dxiv * dxiv                                                      
      do 4130 i = 1, l                                                            
         if (i.eq.1.and.j.eq.1) go to 4130                                            
         xci1   = 1. - cos(twopibl * (i - 1))                                        
         term   = -2. * xci1 * dxsqiv - 2. * xcj1 * dysqiv                           
         termiv = 1. / term                                                        
         u(i,j) = u(i,j) * termiv * flmiv                                          
 4130 continue                                                                  
      call FOURT (u,nn,2,1,1,work)                                               
      do 4131 j = 1, m                                                            
      do 4131 i = 1, l                                                            
         u(i,j) = u(i,j) - (j - 1) * ubar * dy                                     
 4131 continue                                                                  
  100 if (if2.eq.1) go to 130                                                    
      call FOURT (v,nn,2,-1,0,work)                                              
      do 4132 j = 1, m                                                           
         xcj1   = 1. - cos(twopibm * (j-1))                                          
         dxiv   = 1. / dx(j)                                                         
         dxsqiv = dxiv * dxiv                                                      
      do 4132 i = 1, l                                                           
         if (i.eq.1.and.j.eq.1) go to 4132                                           
         xci1   = 1. - cos(twopibl * (i - 1))                                        
         term   = -2. * xci1 * dxsqiv - 2. * xcj1 * dysqiv                           
         termiv = 1. / term                                                        
         v(i,j) = -v(i,j) * termiv * flmiv                                         
 4132 continue                                                                  
      call FOURT (v,nn,2,1,1,work)                                               
      do 4133 j = 1, m                                                           
      do 4133 i = 1, l                                                           
         v(i,j) = v(i,j) - (j - 1) * vbar * dy                                     
 4133 continue                                                                  
  130 continue                                                                  
      do 4134 j = 1, m                                                           
         jp1    = j + 1                                                               
         jm1    = j - 1                                                               
         if (j.eq.1) jm1 = m                                                        
         if (j.eq.m) jp1 = 1                                                        
         dxiv   = 1. / dx(j)                                                         
      do 4134 i = 1, l                                                           
         ip1    = i + 1                                                               
         im1    = i - 1                                                               
         if (i.eq.1) im1 = l                                                        
         if (i.eq.l) ip1 = 1                                                        
         tt1    = -(real(u(i,jp1)) - real(u(i,j))) * dyiv                             
     &            -(real(v(i,j)) - real(v(im1,j))) * dxiv                             
         tt2    = (real(u(ip1,j)) - real(u(i,j))) * dxiv                              
     &           -(real(v(i,j)) - real(v(i,jm1))) * dyiv                              
         if (j.eq.m) tt1 = tt1 + m * ubar                                           
         if (j.eq.1) tt2 = tt2 + m * vbar                                           
         tt3    = real(u(i,j))                                                        
         tt4    = real(v(i,j))                                                        
         u(i,j) = cmplx(tt3,tt1)                                                   
         v(i,j) = cmplx(tt4,tt2)                                                   
 4134 continue                                                                  
      temp      = cmplx(0.,-1.)                                                      
      do 4135 j = 1, m                                                           
      do 4135 i = 1, l                                                           
         u(i,j) = temp * u(i,j)                                                    
         v(i,j) = temp * v(i,j)                                                    
 4135 continue                                                                  
      return                                                                    
      end                                                                       
