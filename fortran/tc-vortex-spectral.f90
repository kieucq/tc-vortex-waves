      program tc_binary
      implicit none
      integer, parameter :: NN=30,MM=30,nx=201,ny=201
      real               :: c(MM,NN),lambda,delta,ct(MM,NN)
      integer            :: i,j,m,n,CM,CN,k,l,loop,Nmax   
      real               :: pi, LL, dt, irec, psi(nx,ny) 
      real               :: xx, yy, dx, dy
!
! initialize spectrum tendency matrix and some coefficients
!      
      pi           = 4.*atan(1.)
      LL           = 1000.e3
      dt           = 50.
      dx           = LL/(nx-1)
      dy           = LL/(ny-1)
      c            = 0
      c(10,10)     = 1.e+5
      c(12,12)     = -1.e+5
      c(14,14)     = -1.e+5
      Nmax         = 10000
!
! re-construct the initial streamline field from the spectrum coefficient
! and print out
!
      open(99,file='tc.dat',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=nx*ny)
      irec         = 1
      do j         = 1,ny
       yy          = pi*(j-1)*dy/LL
       do i        = 1,nx
        xx         = pi*(i-1)*dx/LL
        psi(i,j)   = 0.
        do m       = 1,MM
         do n      = 1,NN
          psi(i,j) = psi(i,j) + c(m,n)*sin(m*xx)*sin(n*yy)
         enddo
        enddo
       enddo
      enddo
      write(99,rec=irec)((psi(i,j),i=1,nx),j=1,ny)
!
! loop thru all the cycles
!
      loop         = 1
18    continue     
      print*,'Looping at: ',loop
      ct           = 0
      do CM        = 1,MM
       do CN       = 1,NN
!         
! ... compute tendency of the spectrum coefficient using Eq. (14)
!
         ct(CM,CN) = 0
         do m      = 1,MM
          do n     = 1,NN
           if (c(m,n).ne.0) then
             do k  = 1,MM
              do l = 1,NN
               if (c(k,l).ne.0) then
                 ct(CM,CN) = n*k*(lambda(k,l,LL)-lambda(m,n,LL))*c(m,n)*c(k,l)*      &
                             (delta(CM,k+m)+delta(CM,m-k)-delta(CM,k-m))*            &
                             (delta(CN,n+l)+delta(CN,l-n)-delta(CN,n-l)) + ct(CM,CN)
               endif
              enddo
             enddo
           endif 
          enddo
         enddo
         ct(CM,CN) = ct(CM,CN)*pi*pi/(2*LL*LL*lambda(CM,CN,LL))
       enddo
      enddo 
!
! ... update c with ct, using a simple fwd scheme.
! 
      c = c + ct*dt
!
! ... re-construct the streamline field from the spectrum coefficient
!     and print out every 100 time steps
!
      if (mod(loop,100).eq.0) then
       do j        = 1,ny
        yy         = pi*(j-1)*dy/LL
        do i       = 1,nx
         xx        = pi*(i-1)*dx/LL
         psi(i,j)  = 0.
         do m      = 1,MM
          do n     = 1,NN
           psi(i,j)= psi(i,j) + c(m,n)*sin(m*xx)*sin(n*yy)
          enddo
         enddo 
        enddo
       enddo
!
! ... write out the psi field
!
       irec        = irec + 1
       write(99,rec=irec)((psi(i,j),i=1,nx),j=1,ny)
      endif
!
! ... next loop
!
      loop = loop + 1
      if (loop.le.Nmax) goto 18
      print*,'Done'
      end

      real function lambda(m,n,LL)
      integer n,m
      real LL, pi
      pi     = 4.*atan(1.)
      lambda = pi*pi*(m*m+n*n)/(LL*LL)
      return 
      end function

      real function delta(m,n)
      integer n,m
      if (m.eq.n) then
       delta = 1
      else
       delta = 0
      endif
      return
      end function


