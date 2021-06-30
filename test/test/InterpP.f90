!!  InterP.f90 
!!
!!  FUNCTIONS:
!!  InterP - Entry point of console application.
!!
!
!!****************************************************************************
!!
!!  PROGRAM: InterP
!!
!!  PURPOSE:  Entry point for the console application.
!!
!!****************************************************************************
!module double_precision
!    integer, parameter :: dp_t = selected_real_kind(18,100)
!    private
!    public :: dp_t
!end module
!    
!program InterP
!    use double_precision
!    implicit none
!    integer  ::  m,n,l
!    real(dp_t), allocatable :: X(:,:,:), Y(:,:,:), Z(:,:,:), C0(:,:,:), &
!                               WindX(:,:,:), WindY(:,:,:), WindZ(:,:,:)
!    
!    m=1000
!    n=100
!    l=100
!    
!    !allocate X(:,:,:)
!    allocate ( X (m,n,l) )
!    allocate ( Y (m,n,l) )
!    allocate ( Z (m,n,l) )
!    allocate ( WindX (m,n,l) )
!    allocate ( WindY (m,n,l) )
!    allocate ( WindZ (m,n,l) )
!    allocate ( C0 (m,n,l) )
!    
!   
!    X = 0._dp_t
!    Y(m,n,l) = 0._dp_t
!    Z(m,n,l) = 0._dp_t
!    C0(m,n,l) = 340-log(Z+1)
!    WindX(:,:,:) = log(Z)
!    WindY(m,n,l) = 0._dp_t
!    WindZ(m,n,l) = 0._dp_t
!    
!    ! Variables
!
!    ! Body of InterP
!    print *, X(100,10,10)
!    
!
!    end program InterP


    
    !*****************************************************************************************
!>
!  Units test for 1d-6d tensor product b-spline interpolation (with extrapolation).

    subroutine interpP(x,y,z,xval,yval,zval,fcn_3d,val)

    !*****************************************************************************************


    use bspline_module
    use bspline_kinds_module, only: wp

    implicit none

    integer,parameter :: nx = 4     !! number of points in x
    integer,parameter :: ny = 4     !! number of points in y
    integer,parameter :: nz = 4     !! number of points in z
    integer,parameter :: kx = 3     !! order in x
    integer,parameter :: ky = 3     !! order in y
    integer,parameter :: kz = 3     !! order in z
    integer,parameter :: iknot = 0  !! automatically select the knots
    
    real(wp),parameter :: tol = 1.0e-2_wp !! tolerance for interp/extrap failure
    real(wp) :: x(nx),y(ny),z(nz)
    real(wp) :: xval,yval,zval
    real(wp) :: tx(nx+kx),ty(ny+ky),tz(nz+kz)
    real(wp) :: fcn_3d(nx,ny,nz)
    real(wp) :: val,tru,err,errmax
    
    logical :: fail
    integer :: i,j,k,l,m,n,idx,idy,idz
    integer :: iflag
    integer :: inbvx,inbvy,inbvz
    integer :: iloy,iloz

    fail = .false.
    idx = 0
    idy = 0
    idz = 0

!**************************************************************************
!change my serching domain and function value here
     !do i=1,nx
     !   x(i) = dble(i-1)/dble(nx-1)
     !end do
     !do j=1,ny
     !   y(j) = dble(j-1)/dble(ny-1)
     !end do
     !do k=1,nz
     !   z(k) = dble(k-1)/dble(nz-1)
     !end do
!**************************************************************************     
     
     !do i=1,nx
     !   do j=1,ny
     !      do k=1,nz
     !         fcn_3d(i,j,k) = f3(x(i),y(j),z(k))
     !      end do
     !   end do
     !end do
     
!**************************************************************************



    !have to set these before the first evaluate call:
    inbvx = 1
    inbvy = 1
    inbvz = 1
    iloy  = 1
    iloz  = 1


    ! initialize
    iflag = 0

    call db3ink(x,nx,y,ny,z,nz,fcn_3d,kx,ky,kz,iknot,tx,ty,tz,fcn_3d,iflag)

    ! compute max error at interpolation points

     val = 0.0_wp
     errmax = 0.0_wp
     err = -99999.9_wp

    call db3val(xval,yval,zval,idx,idy,idz,&
                tx,ty,tz,nx,ny,nz,kx,ky,kz,fcn_3d,val,iflag,&
                inbvx,inbvy,inbvz,iloy,iloz,extrap=.true.)
    !tru = f3(xval, yval, zval)
    !err = tru - val
    !write(*,*) err
    ! check max error against tolerance


    !contains
    !
    !
    !    real(wp) function f3 (x,y,z) !! 3d test function
    !    implicit none
    !    real(wp) x,y,z,piov2
    !    piov2 = 1.2_wp*atan(1.0_wp)
    !    f3 = 0.5_wp*( y*exp(-x) + z*sin(piov2*y) )
    !    end function f3


    end subroutine interpP
