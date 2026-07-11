      program wave_decomposition
      implicit none
      integer            :: nx,ny              ! domain dimension
      real, allocatable  :: uc(:,:),vc(:,:)    ! cartesian raw u,v input 
      real, allocatable  :: ur(:,:),vr(:,:)    ! cylindrical raw u,v input
      real, allocatable  :: ur0(:,:),vr0(:,:)  ! wave 0 component in cylindral coordinate
      real, allocatable  :: u0(:,:),v0(:,:)    ! wave 0 component in Cartesian coordinate
      real, allocatable  :: tm(:,:),t0(:,:)    ! temperature
      real, allocatable  :: st(:,:),s0(:,:)    ! stream fucntion
      real, allocatable  :: st0(:,:),vo0(:,:)  ! axisymmtric decomposition
      real, allocatable  :: vo(:,:)            ! vorticity
      real, allocatable  :: dx2i(:)            ! inverted squared of dx
      integer            :: i,j,k              ! running index
      integer            :: ic,jc              ! storm center
      integer            :: irec               ! record index for plotting 
      character*300      :: atcfline           ! line of atcf format 
      real               :: vmax               ! wave decom maximum wind 
      integer            :: imax0,imax1,imax2  ! rounded maximum wind (kt)
      integer            :: imaxr              ! residual maximum wind (kt)
      integer            :: vmax0              ! total maximum wind (kt)
      real               :: pi,sta,sta1,xc,yc,rmw
      real               :: dx,dy,dx2,dy2i,radius
!
! initialized vorticity field and model parameters
!
      nx          = 251
      ny          = 251
      dx          = 1000.
      dy          = 1000.
      rmw         = 30000
      allocate(uc(nx,ny),vc(nx,ny),tm(nx,ny),st(nx,ny),vo(nx,ny))
      allocate(vr(nx,ny),ur(nx,ny))
      allocate(vr0(nx,ny),ur0(nx,ny))
      allocate(v0(nx,ny),u0(nx,ny))
      allocate(t0(nx,ny),s0(nx,ny))
      allocate(vo0(nx,ny),st0(nx,ny))
      allocate(dx2i(ny))
      ic          = (nx+1)/2
      jc          = (ny+1)/2
      do j        = 1,ny
       do i       = 1,nx
         radius   = sqrt(((i-ic)*dx)**2 + ((j-jc)*dy)**2)
         if (radius < rmw) then
          vo(i,j) = 1.0e-3
         else
          vo(i,j) = 0.
         endif
       enddo
       dx2i(j)    = 1/(dx*dx)
      enddo
      dx2         = dx*dx
      dy2i        = 1/(dy*dy) 
!     call RELAX(s0,dy2i,dx2,vo,dx2i,nx,nx-1,ny-1,ny)
      call p2d_dd_relax(s0,nx,ny,dx,dy,vo,2500,1.e-3,0)
      do j        = 2,ny-1
       do i       = 2,nx-1
        uc(i,j)   = -(s0(i,j+1)-s0(i,j-1))/(2*dy)
        vc(i,j)   =  (s0(i+1,j)-s0(i-1,j))/(2*dx)
       enddo
      enddo
      uc(1,1:ny)  = uc(2,1:ny)
      uc(1:nx,2)  = uc(1:nx,1)
      uc(nx,1:ny) = uc(nx-1,1:ny)
      uc(1:nx,ny) = uc(1:nx,ny-1)
      vc(1,1:ny)  = vc(2,1:ny)
      vc(1:nx,2)  = vc(1:nx,1)
      vc(nx,1:ny) = vc(nx-1,1:ny)
      vc(1:nx,ny) = vc(1:nx,ny-1)
!
! plot the axisymmetric initialization for checking
!
      open(91,file='wave.dat',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=nx*ny)
      irec        = 1
      write(91,rec=irec)((uc(i,j),i=1,nx),j=1,ny)
      irec        = irec + 1
      write(91,rec=irec)((vc(i,j),i=1,nx),j=1,ny)
      irec        = irec + 1
      write(91,rec=irec)((s0(i,j),i=1,nx),j=1,ny)
      irec        = irec + 1
      write(91,rec=irec)((vo(i,j),i=1,nx),j=1,ny)
!
! Adding a few small spot for vorticity next to see the impacts of small
! eddies around RMW = 30km
!
      pi          = 4.*atan(1.)
      do j        = 1,ny
       do i       = 1,nx
          xc=float(i)-ic
          yc=float(j)-jc
          radius   = sqrt(((i-ic)*dx)**2 + ((j-jc)*dy)**2)
          if(xc.ne.0.0) sta1=abs(atan(yc/xc))
          if(xc.ge.0.0.and.yc.eq.0.0) sta=0.0
          if(xc.eq.0.0.and.yc.ge.0.0) sta=pi/2.0
          if(xc.eq.0.0.and.yc.lt.0.0) sta=pi*1.5
          if(xc.lt.0.0.and.yc.eq.0.0) sta=pi
          if(xc.gt.0.0.and.yc.ge.0.0) sta=sta1
          if(xc.lt.0.0.and.yc.gt.0.0) sta=pi-sta1
          if(xc.lt.0.0.and.yc.lt.0.0) sta=pi+sta1
          if(xc.gt.0.0.and.yc.lt.0.0) sta=pi*2.0-sta1
          if (0.lt.sta.and.sta.lt.pi/2.and.        &
              rmw.lt.radius.and.radius.lt.1.2*rmw) then
              vo(i,j) = vo(i,j) + 5.e-3*exp(-(radius-rmw)**2/1e5)
          endif
       enddo
      enddo
      call p2d_dd_relax(st,nx,ny,dx,dy,vo,2500,1.e-3,0)
      do j        = 2,ny-1
       do i       = 2,nx-1
        uc(i,j)   = -(s0(i,j+1)-s0(i,j-1))/(2*dy)
        vc(i,j)   =  (s0(i+1,j)-s0(i-1,j))/(2*dx)
       enddo
      enddo
      uc(1,1:ny)  = uc(2,1:ny)
      uc(1:nx,2)  = uc(1:nx,1)
      uc(nx,1:ny) = uc(nx-1,1:ny)
      uc(1:nx,ny) = uc(1:nx,ny-1)
      vc(1,1:ny)  = vc(2,1:ny)
      vc(1:nx,2)  = vc(1:nx,1)
      vc(nx,1:ny) = vc(nx-1,1:ny)
      vc(1:nx,ny) = vc(1:nx,ny-1)
      irec        = irec + 1
      write(91,rec=irec)((uc(i,j),i=1,nx),j=1,ny)
      irec        = irec + 1
      write(91,rec=irec)((vc(i,j),i=1,nx),j=1,ny)
      irec        = irec + 1
      write(91,rec=irec)((st(i,j),i=1,nx),j=1,ny)
      irec        = irec + 1
      write(91,rec=irec)((vo(i,j),i=1,nx),j=1,ny)
!
! projecting wind field from Cartesian to cylindrical coordinate first before decompositing
!
      pi          = 4.*atan(1.)
      do j        = 1,ny
       do i       = 1,nx
          xc=float(i)-ic
          yc=float(j)-jc
          if(xc.ne.0.0) sta1=abs(atan(yc/xc))
          if(xc.ge.0.0.and.yc.eq.0.0) sta=0.0
          if(xc.eq.0.0.and.yc.ge.0.0) sta=pi/2.0
          if(xc.eq.0.0.and.yc.lt.0.0) sta=pi*1.5
          if(xc.lt.0.0.and.yc.eq.0.0) sta=pi
          if(xc.gt.0.0.and.yc.ge.0.0) sta=sta1
          if(xc.lt.0.0.and.yc.gt.0.0) sta=pi-sta1
          if(xc.lt.0.0.and.yc.lt.0.0) sta=pi+sta1
          if(xc.gt.0.0.and.yc.lt.0.0) sta=pi*2.0-sta1
          vr(i,j)=vc(i,j)*cos(sta) - uc(i,j)*sin(sta)
          ur(i,j)=vc(i,j)*sin(sta) + uc(i,j)*cos(sta)
       enddo
      enddo
!
! Now decomposing and computing the azimuhtal average
!
      call decom(vr,vr0,nx,ny,0,0)
      call decom(ur,ur0,nx,ny,0,0)
      call decom(vo,vo0,nx,ny,0,0)
      call decom(st,st0,nx,ny,0,0)
!
! projecting back the wind field from cylindrical to Cartesian to
! coordinate for plotting
!
      ic          = (nx+1)/2
      jc          = (ny+1)/2
      do j        = 1,ny
       do i       = 1,nx
          xc=float(i)-ic
          yc=float(j)-jc
          if(xc.ne.0.0) sta1=abs(atan(yc/xc))
          if(xc.ge.0.0.and.yc.eq.0.0) sta=0.0
          if(xc.eq.0.0.and.yc.ge.0.0) sta=pi/2.0
          if(xc.eq.0.0.and.yc.lt.0.0) sta=pi*1.5
          if(xc.lt.0.0.and.yc.eq.0.0) sta=pi
          if(xc.gt.0.0.and.yc.ge.0.0) sta=sta1
          if(xc.lt.0.0.and.yc.gt.0.0) sta=pi-sta1
          if(xc.lt.0.0.and.yc.lt.0.0) sta=pi+sta1
          if(xc.gt.0.0.and.yc.lt.0.0) sta=pi*2.0-sta1
          v0(i,j)=ur0(i,j)*sin(sta) + vr0(i,j)*cos(sta)
          u0(i,j)=ur0(i,j)*cos(sta) - vr0(i,j)*sin(sta)
       enddo
      enddo
!
! plot the decomposition
!
      irec        = irec + 1
      write(91,rec=irec)((ur0(i,j),i=1,nx),j=1,ny)
      irec        = irec + 1
      write(91,rec=irec)((vr0(i,j),i=1,nx),j=1,ny)
      irec        = irec + 1
      write(91,rec=irec)((ur(i,j),i=1,nx),j=1,ny)
      irec        = irec + 1
      write(91,rec=irec)((vr(i,j),i=1,nx),j=1,ny)
      end
!================================================================
! SUBROUTINES
!================================================================
     SUBROUTINE p2d_dd_relax(phi,nx,ny,dx,dy,rhs,nmax,epsi,debug)
!
! Solver for 2D poisson with DD bnd type, using Lipman iteration
!
     IMPLICIT NONE
     INTEGER nx,ny
     REAL phi(nx,ny),rhs(nx,ny),cx(nx,ny),cy(nx,ny),cxodx2(nx,ny)
     REAL phio(nx,ny),dx,dy,dx2,dy2,denom(nx,ny),numer,cyody2(nx,ny)
     REAL err,err0,epsi,res,relaxfac
     INTEGER i,j,k,loop,nmax,debug
     dx2         = dx*dx
     dy2         = dy*dy
     cx          = 1.
     cy          = 1. 
     relaxfac    = 1.8
     DO i        = 1,nx
      DO j       = 1,ny
       denom(i,j)= 2*(cx(i,j)/dx2+cy(i,j)/dy2)
       cxodx2(i,j) = cx(i,j)/(dx2*denom(i,j))
       cyody2(i,j) = cy(i,j)/(dy2*denom(i,j))
      ENDDO
     ENDDO
     iteration_loop : DO loop = 1,nmax
      phio      = phi
      err       = 0
      DO i      = 2,nx-1
       DO j     = 2,ny-1
         res     = cxodx2(i,j)*(phi(i+1,j) + phi(i-1,j))     &
                 + cyody2(i,j)*(phi(i,j+1) + phi(i,j-1))     &
                 - rhs(i,j)/denom(i,j) - phi(i,j)
   
      phi(i,j)=phi(i,j) + relaxfac*res 
        err     = err + abs(phi(i,j)-phio(i,j))
       ENDDO
      ENDDO
   
      IF (debug.eq.1) PRINT*,loop,err0,err
      IF (loop.eq.1) err0 = err
      IF (err0.ne.0.and.err/err0.le.epsi) goto 55 
     ENDDO iteration_loop
   55 CONTINUE
     IF (loop.lt.nmax) THEN
       write(*,*)"Iteration converged at loop:",loop
     ELSE
       write(*,*)"Reaching max iteration. Stop"
     ENDIF
     RETURN
     END SUBROUTINE p2d_dd_relax

      subroutine RELAX (x,zzinv,z,y,zinv,l,l1,m1,m)
!
!  This subroutine solve the Poisson equation using
!  the overrelaxation method .
!
!  definition of variables :
!
! x         contains first guess and eventually
!           the output of the field to be relaxed
! y         forcing function of l by m matrix
! z         latitudinal increment square
! zinv      inverse of dx**2
! zzinv     inverse of dy**2
! l         east-west dimension
! m         north-south dimension
!
!  definition of constants :
! npts      number of points to be relaxed
! nrel      number of points relaxed
! alfa      relaxation factor
! ia        maximum number of iterations allowed
! eps       tolerance error
! nsc       count of number of scan for convergence
! lsc       last scan after convergence
!
      real x(l,m),y(l,m),z(m),zinv(m)
      npts      = l*(m-2)
      nlax      = 1
      mm        = 2
      mmm       = m1
      alfa      = .46
      ia        = 1000
      eps       = 1.e-5
      nsc       = 0
      lsc       = -1
   15 nrel      = 0
      do 4110 j = mm, mmm
      do 4110 i = 1, l
         im1    = i-1
         ip1    = i+1
         if (im1.lt.1) im1 = l1
         if (ip1.gt.l) ip1 = 2
         r      = (x(ip1,j)+x(im1,j)-2.*x(i,j))*zinv(j)+ &
                  (x(i,j+1)+x(i,j-1)-2.*x(i,j))*zzinv
         r      = (r-y(i,j))*z(j)
         if (lsc-nsc) 29,29,30
   29    x(i,j) = x(i,j) + alfa*r
   30    if (abs(r).le.eps) nrel = nrel+1
 4110 continue
!
! nrel gives the number of points that have converged .
! nsc gives the number of scan made over the domain and
! if it is less then the maximum number of iteration
! allowed for complete convergence it keeps going to
! statement number 15. lsc is the final check once
! convergence has been achieved it does one more loop
! before finally jumping out to statement 300 .
!
      nsc       = nsc+1
      if (nrel-npts) 13,14,14
   14 if (lsc .ge. nsc) go to 300
   18 lsc       = nsc+1
      write(*,*)"Iteration at step",nsc,r,eps
   13 if (nsc.lt.ia) go to 15
  201 format(50h   progress of relaxation npts,nrel,nsc,ia        )
  300 continue
      write(6,201)
  200 format(6x,4i9)
      write (6,200)npts,nrel,nsc,ia
      return
      end

      subroutine decom(fi,fir,mm,nn,mc,uvt)
      integer mm,nn,mc
      parameter(lm=124,ln=180)
      real fi(mm,nn),fir(mm,nn)
      real fd(2*mm-1,2*nn-1),fdr(2*mm-1,2*nn-1)
      real fi0(2*mm-1,2*nn-1),fi1(2*mm-1,2*nn-1)
      real fi2(2*mm-1,2*nn-1),fi3(2*mm-1,2*nn-1)
      real fi4(2*mm-1,2*nn-1),fi5(2*mm-1,2*nn-1)
      real fi6(2*mm-1,2*nn-1)
      real fs(lm)
      real ft0(lm,ln),ft(lm,ln),fts(lm,ln),fta(lm,ln)
      real rst1(lm,ln),rst2(lm,ln),rst3(lm,ln)
      real rst4(lm,ln),rst5(lm,ln),rst6(lm,ln)
      real x1(2*mm-1),x2(2*nn-1)
      integer uvt,mm2,nn2

      mm2=2*mm-1
      nn2=2*nn-1

      j0=mm
      i0=nn

      dd=360.0/float(ln)

      do i=1,mm
      do j=1,nn
       fd(2*i-1,2*j-1)=fi(i,j)
      end do
      end do

      do i=1,mm-1
      do j=1,nn-1
       fd(2*i,  2*j  )=(fi(i,j)+fi(i+1,j)+fi(i,j+1)+fi(i+1,j+1))*0.25
      end do
      end do

      do i=1,mm-1
      do j=1,nn
       fd(2*i  ,2*j-1)=(fi(i,j)+fi(i+1,j))*0.5
      end do
      end do

      do i=1,mm
      do j=1,nn-1
       fd(2*i-1,2*j  )=(fi(i,j)+fi(i,j+1))*0.5
      end do
      end do
      fd(i0+1,j0+1)= fd(i0+1,j0+1)*0.4+fd(i0+2,j0+2)*0.6
      fd(i0+1,j0)=  fd(i0+1,j0)*0.4+fd(i0+2,j0)*0.6
      fd(i0+1,j0-1)=  fd(i0+1,j0-1)*0.4+fd(i0+2,j0-2)*0.6
      fd(i0-1,j0+1)=  fd(i0-1,j0+1)*0.4+fd(i0-2,j0+2)*0.6
      fd(i0-1,j0)=  fd(i0-1,j0)*0.4+fd(i0-2,j0)*0.6
      fd(i0-1,j0-1)=  fd(i0-1,j0-1)*0.4+fd(i0-2,j0-2)*0.6
      fd(i0,j0+1)=  fd(i0,j0+1)*0.4+fd(i0,j0+2)*0.6
      fd(i0,j0-1)=  fd(i0,j0-1)*0.4+fd(i0,j0-2)*0.6
!------------------------------------------------------------
! angle is measued from due north (j0,i0) the center of frame
!************************************************************
      pi1=atan(1.0)/45.0
      do 10 i=1,lm
      r1=float(i-1)
      do 10 j=1,ln
      phi=float(j-1)*dd

      x=float(j0)-r1*sin(phi*pi1)
      y=float(i0)+r1*cos(phi*pi1)
      call SCINEX(x,y,fd,scint0,mm2,nn2)

      ft0(i,j)=scint0
      ft(i,j)=scint0

  10  continue

      do 11 j=1,ln

          ft0(1,j)=fd(j0,i0)
          ft(1,j)=ft0(1,j)

  11  continue

      do i=1,lm
      do j=1,ln

       ftsmooth=0.0

! this few lines is for rotating angle of fields

       do jr=-2,3
        j00=j+jr
        if(j00.le.0)j00=ln+j00
        if(j00.gt.ln)j00=j00-ln
        ftsmooth=ftsmooth+ft0(i,j00)
       end do

       ft(i,j)=ftsmooth/6.0

      end do
      end do

      if(mc.eq.0) then

      call transf000(fts,fs,ft,lm,ln)    ! get the symmetric component
      call transi(fdr,fts,lm,ln,mm2,nn2)  ! tansform back to Cartesian

      else if(mc.eq.1) then

      call transf00(fts,fs,ft,lm,ln)    ! get the symmetric component
      call transi(fi0,fts,lm,ln,mm2,nn2)  ! tansform back to Cartesian
      call transf0(fta,fs,ft,lm,ln)     ! get the asymmetric component
      call transf(rst1,fta,1.0,lm,ln)   ! take the  wave number 1 comp
      call transi(fi1,rst1,lm,ln,mm2,nn2)
      do j=1,mm2
      do i=1,nn2
         fdr(j,i)=fi0(j,i)+fi1(j,i)           ! the sum of wn0 and wn1 component
      end do
      end do

      else if(mc.eq.2) then

      call transf00(fts,fs,ft,lm,ln)    ! get the symmetric component
      call transi(fi0,fts,lm,ln,mm2,nn2)  ! tansform back to Cartesian
      call transf0(fta,fs,ft,lm,ln)     ! get the asymmetric component
      call transf(rst2,fta,2.0,lm,ln)
      call transi(fi2,rst2,lm,ln,mm2,nn2)
      do j=1,mm2
      do i=1,nn2
         fdr(j,i)=fi0(j,i)+fi2(j,i)           ! the sum of wn0 and wn2
      end do
      end do

      else if(mc.eq.3) then

      call transf00(fts,fs,ft,lm,ln)    ! get the symmetric component
      call transi(fi0,fts,lm,ln,mm2,nn2)  ! tansform back to Cartesian
      call transf0(fta,fs,ft,lm,ln)     ! get the asymmetric component
      call transf(rst3,fta,3.0,lm,ln)
      call transi(fi3,rst3,lm,ln,mm2,nn2)

      do j=1,mm2
      do i=1,nn2
         fdr(j,i)=fi0(j,i)+fi3(j,i)           ! the sum of wn0 and wn2
      end do
      end do

      else if(mc.eq.4) then

      call transf00(fts,fs,ft,lm,ln)    ! get the symmetric component
      call transi(fi0,fts,lm,ln,mm2,nn2)  ! tansform back to Cartesian
      call transf0(fta,fs,ft,lm,ln)     ! get the asymmetric component
      call transf(rst4,fta,4.0,lm,ln)
      call transi(fi4,rst4,lm,ln,mm2,nn2)

      do j=1,mm2
      do i=1,nn2
         fdr(j,i)=fi0(j,i)+fi4(j,i)           ! the sum of wn0 and wn4
      end do
      end do

      else if(mc.eq.5) then

      call transf00(fts,fs,ft,lm,ln)    ! get the symmetric component
      call transi(fi0,fts,lm,ln,mm2,nn2)  ! tansform back to Cartesian
      call transf0(fta,fs,ft,lm,ln)     ! get the asymmetric component
      call transf(rst5,fta,5.0,lm,ln)
      call transi(fi5,rst5,lm,ln,mm2,nn2)

      do j=1,mm2
      do i=1,nn2
         fdr(j,i)=fi0(j,i)+fi5(j,i)           ! the sum of wn0 and wn5
      end do
      end do


      else if(mc.eq.6) then

      call transf00(fts,fs,ft,lm,ln)    ! get the symmetric component
      call transi(fi0,fts,lm,ln,mm2,nn2)  ! tansform back to Cartesian
      call transf0(fta,fs,ft,lm,ln)     ! get the asymmetric component
      call transf(rst1,fta,1.0,lm,ln)   ! take the  wave number 1
      call transi(fi1,rst1,lm,ln,mm2,nn2)
      call transf(rst2,fta,2.0,lm,ln)
      call transi(fi2,rst2,lm,ln,mm2,nn2)
      call transf(rst3,fta,3.0,lm,ln)
      call transi(fi3,rst3,lm,ln,mm2,nn2)
      call transf(rst4,fta,4.0,lm,ln)
      call transi(fi4,rst4,lm,ln,mm2,nn2)
      call transf(rst5,fta,5.0,lm,ln)
      call transi(fi5,rst5,lm,ln,mm2,nn2)
      call transf(rst6,fta,6.0,lm,ln)
      call transi(fi6,rst6,lm,ln,mm2,nn2)

      do 20 j=1,mm2
      do 20 i=1,nn2
      fdr(j,i)=fi0(j,i)+fi1(j,i)+fi2(j,i)+fi3(j,i)+fi4(j,i)+fi5(j,i)+fi6(j,i)  ! the sum of 0,1,2,3,4,5,6 components
  20  continue

      else if (mc.eq.7) then

      call transi(fdr,ft,lm,ln,mm2,nn2)  ! tansform back to Cartesian

      end if

      do i=1,mm
      do j=1,nn
       fir(i,j)=fdr(2*i-1,2*j-1)
      end do
      end do

      return
      end subroutine decom

        subroutine transf000(out,amean,u,lm,ln)
        real u(lm,ln),out(lm,ln),amean(lm)
        do 20 j=1,lm
        amean(j)=0.
        do i=1,ln
        amean(j)=amean(j)+u(j,i)
        enddo
        amean(j)=amean(j)/float(ln)
  20    continue

        do 26 i=1,ln
        do 25 j=1,10
        out(j,i)=amean(1)
  25    continue
        out(lm,i)=amean(lm)
  26    continue

       do 35 ii=1,50
        do 30 i=1,ln
        do 30 j=2,lm-1
        out(j,i)=0.25*amean(j-1)+0.5*amean(j)+0.25*amean(j+1)
  30    continue
  35   continue
       return
       end subroutine transf000


        subroutine transf00(out,amean,u,lm,ln)
        real u(lm,ln),out(lm,ln),amean(lm)
        do 20 j=1,lm
        amean(j)=0.
        do i=1,ln
        amean(j)=amean(j)+u(j,i)
        enddo
        amean(j)=amean(j)/float(ln)
  20    continue

        do 30 i=1,ln
        do 30 j=1,lm
        out(j,i)=amean(j)
  30    continue
        return
        end subroutine transf00

        subroutine transf0(out,amean,u,lm,ln)
        real u(lm,ln),out(lm,ln),amean(lm)
        do 20 j=1,lm
        amean(j)=0.
        do i=1,ln
        amean(j)=amean(j)+u(j,i)
        enddo
        amean(j)=amean(j)/float(ln)
  20    continue
        do 30 i=1,ln
        do 30 j=1,lm
        out(j,i)=u(j,i)-amean(j)
  30    continue
        return
        end subroutine transf0


        subroutine transf(out,psia,ai,l,m)
        real out(l,m),ps1a(m),psia(l,m)
        pi=4.0*atan(1.0)
        dl=2*pi/float(m)
!
! to get the component you want
!
        do 30 j=1,l
        do i=1,m
        ps1a(i)=psia(j,i)
        enddo
        call coscoeff(ps1a,cc1a,m,dl,pi,ai)
        call sincoeff(ps1a,cs1a,m,dl,pi,ai)
        do i=1,m
        ang=float(i-1)*dl
        out(j,i)=cc1a*cos(ang*ai)+cs1a*sin(ang*ai)
        enddo 
  30    continue
        return
        end subroutine transf


        subroutine coscoeff(psi,psim,m,dl,pi,ai)
        dimension psi(m)
        psim=0.0
        do 20 j=1,m
        ang=float(j-1)*dl
        psim=psim+psi(j)*cos(ang*ai)*dl
  20    continue
        psim=psim/pi
        return
        end subroutine coscoeff


        subroutine sincoeff(psi,psim,m,dl,pi,ai)
        dimension psi(m)
        psim=0.0
        do 20 j=1,m
        ang=float(j-1)*dl
        psim=psim+psi(j)*sin(ang*ai)*dl
  20    continue
        psim=psim/pi
        return
        end subroutine sincoeff


      subroutine transi0(u1,uu,im,in,mm,nn)
      real uu(im),y2(im),rd(im),u1(mm,nn)
!
! j is x-direction (zonal), i is y-direction (meridional)
! angle is measued from due north,(j0,i0) the center of frame
!
       pi=4.0*atan(1.0)
       pi1=180.0/pi
       m0=(mm+1)/2
       n0=(nn+1)/2

       do i=1,im
        rd(i)=float(i)
       end do

!     call SPLINE to get second derivatives

        call spline(rd,uu,im,yp1,ypn,y2)

       do 110 j=1,mm
       do 110 i=1,nn

       x=float(j-m0)
       y=float(i-n0)
       r1=sqrt(x**2+y**2)

!     call SPLINT for interpolations

       call splint(rd,uu,y2,im,r1,sz1)

       u1(j,i)=sz1
110    continue

       u1(m0,n0)=uu(1)

       return
       end subroutine transi0

!*****************************************************************
!  This code is used to transfer vector in Cylindrical to Cartesian 
!  cooridinates, originally coded by Liguang Wu (1996.8) and modified
!  by Yuqing Wang (1997.1)
!  lm,ln are grid number in cyl. co., respectively. the thansfering
!  domain is 2000km*2000km for 10km*10km grid. 
!*****************************************************************
      subroutine transi(u1,uu,im,in,mm,nn)
      real uu(im,in),u1(mm,nn)
      dd=360.0/float(in)
!************************************************************* 
! j is x-direction (zonal), i is y-direction (meridional)
! angle is measued from due north,(j0,i0) the center of frame
!*************************************************************
       pi=4.0*atan(1.0)
       pi1=180.0/pi
       m0=(mm+1)/2
       n0=(nn+1)/2
       do 110 j=1,mm
       do 110 i=1,nn
       if(j.eq.m0.and.i.eq.n0) then
       u1(j,i)=uu(1,1)
       else
       x=float(j-m0)
       y=float(i-n0)
       r1=sqrt(x**2+y**2)
       phi=asin(x/r1)
       if(x.le.0..and.y.ge.0.) phi=-phi
       if(x.lt.0..and.y.lt.0.) phi=pi+phi
       if(x.ge.0..and.y.le.0.) phi=pi+phi
       if(x.ge.0..and.y.gt.0.) phi=2*pi-phi
       phi=phi*pi1
!************************************************
! (lm,ln) the position in cylindrical coordinates
!************************************************
       ln=int(phi/dd+0.5)+1
       lm=int(r1+0.5)+1
       dr=r1-lm
       da=phi-ln*dd
        mp1=lm+1
        np1=ln+1
        mm1=lm-1
        nm1=ln-1
        if(mp1.gt.im) mp1=im
        if(np1.gt.in) np1=np1-in
        if(mm1.lt.1) mm1=1 
        if(nm1.lt.1) nm1=in-nm1 
        if(ln.gt.in) ln=ln-in
        if(lm.gt.im) lm=im
!.......................................................
      d1=0.5*(uu(mp1,ln)-uu(mm1,ln))
      d2=0.5*(uu(lm,np1)-uu(lm,nm1))/dd
      d3=0.5*(uu(mp1,ln)-2.*uu(lm,ln)+uu(mm1,ln))
      d4=0.5*(uu(lm,np1)-2.*uu(lm,ln)+uu(lm,nm1))/dd/dd
      d5=0.125*(uu(mp1,np1)+uu(mm1,nm1)-uu(mm1,np1)-uu(mp1,nm1))/dd
      u1(j,i)=uu(lm,ln)+d1*dr+d2*da+d3*dr**2+d4*da**2+d5*da*dr
!.......................................................
      endif
      phi=phi/pi1
110   continue
      return
      end

      subroutine d2spline(x1,x2,xx1,xx2,y,f,m,n)
      INTEGER M,N,i,j
      REAL xx1,xx2,x1(M),x2(N),y(M,N),y2(M,N),f
      call splie2(x1,x2,y,M,N,y2)
      call splin2(x1,x2,y,y2,M,N,xx1,xx2,f)
      return
      END subroutine d2spline


      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      INTEGER m,n,NN
      REAL x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100)
      INTEGER j,k
      REAL y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      call spline(x1a,yytmp,m,1.e30,1.e30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      END SUBROUTINE splin2

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END SUBROUTINE splint 

      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      INTEGER m,n,NN
      REAL x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100)
      INTEGER j,k
      REAL y2tmp(NN),ytmp(NN)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.e30,1.e30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return
      END SUBROUTINE splie2

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END  SUBROUTINE spline


      SUBROUTINE SCINEX(GM,GN,SCALA,SCINTO,lq,lp)
!
! THIS SUBROUTINE PRODUCES THE VALUE SCINTO OF A SCALAR FIELD AT A POINT
! GM,GN BY INTERPOLATION OR EXTRAPOLATION OF THE FIELD SCALA  (2-DIRECTI
! BESSEL INTERPOLATION FORMULA). MMIN,MMAX AND NMIN,NMAX ARE THE BOUNDAR
! OF THE GRID ARRAY.
!
      REAL SCALA(lq,lp)
      MMIN=1
      NMIN=1
      MMAX=lq
      NMAX=lp
      IGM=int(GM)
      JGN=int(GN)
      FM=GM-IGM
      FN=GN-JGN
      IF(FM.LT.1.E-06)FM=0.
      IF(FN.LT.1.E-06)FN=0.
      MS=MMAX-1
      NS=NMAX-1
      MR=MMIN+1
      NR=NMIN+1
      IF(GM.LT.MMAX)GO TO 60
      IF(GN.LT.NMAX)GO TO 20
      E=GM-MMAX
      T1=E*(SCALA(MMAX,NMAX)-SCALA(MS,NMAX))
      E=GN-NMAX
      T2=E*(SCALA(MMAX,NMAX)-SCALA(MMAX,NS))
      SCINTO=SCALA(MMAX,NMAX)+T1+T2
      RETURN
   20 IF(GN.GE.NMIN)GO TO 40
      E=GM-MMAX
      T1=E*(SCALA(MMAX,NMIN)-SCALA(MS,NMIN))
      E=NMIN-GN
      T2=E*(SCALA(MMAX,NMIN)-SCALA(MMAX,NR))
      SCINTO=SCALA(MMAX,NMIN)+T1+T2
      RETURN
   40 P=SCALA(MMAX,JGN)+FN*(SCALA(MMAX,JGN+1)-SCALA(MMAX,JGN))
      H=SCALA(MS,JGN)+FN*(SCALA(MS,JGN+1)-SCALA(MS,JGN))
      E=GM-MMAX
      SCINTO=P+E*(P-H)
      RETURN
   60 IF(GM.GE.MMIN)GO TO 140
      IF(GN.LT.NMAX)GO TO 80
      E=GN-NMAX
      T2=E*(SCALA(MMIN,NMAX)-SCALA(MMIN,NS))
      E=MMIN-GM
      T1=E*(SCALA(MMIN,NMAX)-SCALA(MR,NMAX))
      SCINTO=SCALA(MMIN,NMAX)+T1+T2
      RETURN
   80 IF(GN.GE.NMIN)GO TO 100
      E=NMIN-GN
      T2=E*(SCALA(MMIN,NMIN)-SCALA(MMIN,NR))
      E=MMIN-GM
      T1=E*(SCALA(MMIN,NMIN)-SCALA(MR,NMIN))
      SCINTO=SCALA(MMIN,NMIN)+T1+T2
      RETURN
  100 E=MMIN-GM
      P=SCALA(MMIN,JGN)+FN*(SCALA(MMIN,JGN+1)-SCALA(MMIN,JGN))
      H=SCALA(MR,JGN)+FN*(SCALA(MR,JGN+1)-SCALA(MR,JGN))
      SCINTO=P+E*(P-H)
      RETURN
  120 E=GN-NMAX
      P=SCALA(IGM,NMAX)+FM*(SCALA(IGM+1,NMAX)-SCALA(IGM,NMAX))
      H=SCALA(IGM,NS)+FM*(SCALA(IGM+1,NS)-SCALA(IGM,NS))
      SCINTO=P+E*(P-H)
      RETURN
  140 IF(GN.GE.NMAX)GO TO 120
      IF(GN.GE.NMIN)GO TO 160
      E=NMIN-GN
      P=SCALA(IGM,NMIN)+FM*(SCALA(IGM+1,NMIN)-SCALA(IGM,NMIN))
      H=SCALA(IGM,NR)+FM*(SCALA(IGM+1,NR)-SCALA(IGM,NR))
      SCINTO=P+E*(P-H)
      RETURN
  160 IF(GM.LT.MS.AND.GM.GE.MR.AND.GN.LT.NS.AND.GN.GE.NR)GO TO 180
      P=SCALA(IGM+1,JGN)+FN*(SCALA(IGM+1,JGN+1)-SCALA(IGM+1,JGN))
         H=SCALA(IGM,JGN)+FN*(SCALA(IGM,JGN+1)-SCALA(IGM,JGN))
      SCINTO=H+FM*(P-H)
      RETURN
  180    FQ=0.25*(FM*FM-FM)
      A=SCALA(IGM,JGN-1)+FM*(SCALA(IGM+1,JGN-1)-SCALA(IGM,JGN-1))      &
      +FQ*(SCALA(IGM+2,JGN-1)+SCALA(IGM-1,JGN-1)-SCALA(IGM+1,JGN-1)-   &
      SCALA(IGM,JGN-1))
      B=SCALA(IGM,JGN)+FM*(SCALA(IGM+1,JGN)-SCALA(IGM,JGN))            &
      +FQ*(SCALA(IGM+2,JGN)+SCALA(IGM-1,JGN)-SCALA(IGM+1,JGN)-         &
      SCALA(IGM,JGN))
      C=SCALA(IGM,JGN+1)+FM*(SCALA(IGM+1,JGN+1)-SCALA(IGM,JGN+1))      &
      +FQ*(SCALA(IGM+2,JGN+1)+SCALA(IGM-1,JGN+1)-SCALA(IGM+1,JGN+1)    &
      -SCALA(IGM,JGN+1))
      D=SCALA(IGM,JGN+2)+FM*(SCALA(IGM+1,JGN+2)-SCALA(IGM,JGN+2))      &
      +FQ*(SCALA(IGM+2,JGN+2)+SCALA(IGM-1,JGN+2)-SCALA(IGM+1,JGN+2)    &
      -SCALA(IGM,JGN+2))
      SCINTO=B+FN*(C-B)+0.25*(FN*FN-FN)*(A+D-B-C)
      RETURN
      END
       
      subroutine car2cyl(nx,ny,ua,va,upr,vpr)
      integer           :: nx,ny
      real, intent(in)  :: ua(nx,ny),va(nx,ny)
      real, intent(out) :: upr(nx,ny),vpr(nx,ny)
      real              :: vt(nx,ny),vr(nx,ny)
      real              :: vtr(nx,ny),ur0(nx,ny)
      integer           :: i,j,k
      real              :: xc,yc,txc0,tyc0
      real              :: sta,pil,a,b
      txc0         = float(nx)/2.
      tyc0         = float(ny)/2.
      pil          = 4.*atan(1.)
      do i=1,ny
       do j=1,nx
        xc=float(j)-txc0
        yc=float(i)-tyc0
        if(xc.ne.0.0) sta1=abs(atan(yc/xc))
        if(xc.ge.0.0.and.yc.eq.0.0) sta=0.0
        if(xc.eq.0.0.and.yc.ge.0.0) sta=pil/2.0
        if(xc.eq.0.0.and.yc.lt.0.0) sta=pil*1.5
        if(xc.lt.0.0.and.yc.eq.0.0) sta=pil
        if(xc.gt.0.0.and.yc.ge.0.0) sta=sta1
        if(xc.lt.0.0.and.yc.gt.0.0) sta=pil-sta1
        if(xc.lt.0.0.and.yc.lt.0.0) sta=pil+sta1
        if(xc.gt.0.0.and.yc.lt.0.0) sta=pil*2.0-sta1
        a=sin(sta)
        b=cos(sta)
        vt(j,i)=va(j,i)*b-ua(j,i)*a
        vr(j,i)=va(j,i)*a+ua(j,i)*b
       end do
       end do

       call decom(vt,vtr,nx,ny,0,1)
       call decom(vr,ur0,nx,ny,0,0)

       do i=1,ny
       do j=1,nx
        xc=float(j)-txc0
        yc=float(i)-tyc0
        if(xc.ne.0.0) sta1=abs(atan(yc/xc))
        if(xc.ge.0.0.and.yc.eq.0.0) sta=0.0
        if(xc.eq.0.0.and.yc.ge.0.0) sta=pil/2.0
        if(xc.eq.0.0.and.yc.lt.0.0) sta=pil*1.5
        if(xc.lt.0.0.and.yc.eq.0.0) sta=pil
        if(xc.gt.0.0.and.yc.ge.0.0) sta=sta1
        if(xc.lt.0.0.and.yc.gt.0.0) sta=pil-sta1
        if(xc.lt.0.0.and.yc.lt.0.0) sta=pil+sta1
        if(xc.gt.0.0.and.yc.lt.0.0) sta=pil*2.0-sta1
        a=sin(sta)
        b=cos(sta)
        upr(j,i)=ur0(j,i)*b-vtr(j,i)*a
        vpr(j,i)=ur0(j,i)*a+vtr(j,i)*b
       end do
       end do
       return
       end subroutine car2cyl


