!This is a 3D ray tracing code for nolinear acoustic propagation
!developed by Tianshu Zhang, TFDTG, University of Florida

program test
    use bspline_kinds_module, only: wp
    implicit none
    
    integer     :: nx, ny, nz, nt, i, j , k, yd
    integer     :: ii, jj, kk, NTT, ReflectionTimeIndex
    integer, allocatable   ::  ReflectionPointMax(:)
    real        :: ZT, XT, YT
    real, dimension(2)     ::  w
    real(wp)    :: d, d1, xver, yver, zver, dt
    real(wp)    :: ValC, ValU, ValV, Valduz, Valdvz, Valdwz, Valduy, Valdvy, Valdwy
    real(wp)    :: Valdcx, Valdcy, Valdcz, Valdux, Valdvx, Valdwx
    real(wp)    :: omega, Kx, Ky, Kz, Knormal
    real(wp), dimension(4) ::  Xnx,Yny,Znz
    real(wp), allocatable  ::  x(:,:,:),y(:,:,:),z(:,:,:), dudz(:,:,:), dcdz(:,:,:)
    real(wp), allocatable  ::  xc(:,:,:),yc(:,:,:),zc(:,:,:)
    real(wp), allocatable  ::  C0(:,:,:),Uwind(:,:,:),Vwind(:,:,:)
    real(wp), allocatable  ::  XX(:,:,:),YY(:,:,:),ZZ(:,:,:),UU(:,:,:),CC(:,:,:), VV(:,:,:)
    real(wp), allocatable  ::  Xrecord(:), Yrecord(:), Zrecord(:)
    real(wp), allocatable  ::  Vout(:,:,:,:,:), Location(:,:)
    real(wp), allocatable  ::  C0out(:,:,:,:), XgridLoc(:,:,:,:,:), NNormalN(:,:)  
    real(wp), dimension(4,4,4) :: Duz, Dcz, Dvz, Duy, Dcy, Dvy, Dux, Dcx, Dvx, Dwx, Dwy, Dwz
    real(wp), dimension(3) ::  n, Nnew, Term1, Term2    
    real(wp)  ::  soundspeed, Uwindspeed, Vwindspeed
    real(wp)  ::  xt1,yt1,zt1,angle_theta,angle_phi
    real(wp)  ::  xp,yp,zp,xm,ym,zm, Cdown,  Udown,  Vdown, Cup,    Vup,    Uup
    real(wp)  ::  xl,yl,zl,xr,yr,zr, Cleft,  Uleft,  Vleft, Cright, Uright, Vright
    real(wp)  ::  xf,yf,zf,xb,yb,zb, Cfront, Ufront, Vfront,Cback,  Uback , Vback
    character(100)  ::  str2theta, str2phi

    call get_command_argument(1, str2theta)
    call get_command_argument(2, str2phi) 
    read(str2phi,*)angle_phi
    read(str2theta,*)angle_theta
    angle_phi = angle_phi*4.0 * ATAN(1.0)/180.0
    angle_theta = angle_theta*4.0 * ATAN(1.0)/180.0

    nx = 1000
    ny = 50
    nz = 1000
    nt = 10
    dt = 0.01_wp
    omega = 150.
    allocate (ReflectionPointMax(1000))
    allocate (Location(3,nt))
    allocate (x(nx,ny,nz))
    allocate (y(nx,ny,nz))
    allocate (z(nx,ny,nz))
    allocate (xc(nx-1,ny-1,nz-1))
    allocate (yc(nx-1,ny-1,nz-1))
    allocate (zc(nx-1,ny-1,nz-1))
    allocate (C0(nx,ny,nz))
    allocate (Uwind(nx,ny,nz))
    allocate (Vwind(nx,ny,nz))
    allocate (XX(4,4,4))
    allocate (YY(4,4,4))
    allocate (ZZ(4,4,4))
    allocate (UU(4,4,4))
    allocate (VV(4,4,4))
    allocate (CC(4,4,4))
    allocate (dudz(nx,ny,(nz-1)))
    allocate (dcdz(nx,ny,(nz-1)))
    allocate (Xrecord(nt))
    allocate (Yrecord(nt))
    allocate (Zrecord(nt))
    allocate (Vout(4,4,4,3,nt))
    allocate (C0out(4,4,4,nt))
    allocate (XgridLoc(4,4,4,3,nt))
    allocate (NNormalN(3,nt))

    do i=1,nx
      x(i,:,:) = 100000._wp*(dble(i-1)/dble(nx-1)) 
    end do
    do i=1,ny
      y(:,i,:) = 50000._wp*(dble(i-1)/dble(ny-1))
    end do
    do i=1,nz
      z(:,:,i) = 100000._wp*(dble(i-1)/dble(nz-1))
    end do

    !do i=1,nz
    !  C0(:,:,i) = 343._wp-dlog(1._wp+z(:,:,i))
    !end do
    !
    !do i=1,nz
    !    Uwind(:,:,i) = (z(:,:,i)-2000)**2./100000.
    !end do
    !do i=1,nz
    !    Vwind(:,:,i) = log(1.+abs(z(:,:,i)-2000))/100000.
    !end do
    
    !z_step = zrange/nz
    xver = 100000._wp*.13_wp
    yver = 50000._wp*.08_wp
    zver = 20000.+10.

    n(1) = cos(angle_theta)!1._wp/(2._wp)**.5_wp
    n(2) = 0!sin(angle_phi)
    n(3) = sin(angle_theta)!1._wp/(2._wp)**.5_wp
    Nnew = n
    print*, n(1), n(2), n(3)
     
    call cellcenter(x,y,z,xc,yc,zc,nx,ny,nz)
    
    ReflectionTimeIndex = 1
    !open(91,file = '/home/tianshu/Documents/DATATRANSFER/ReflectionTime.txt')

    do NTT=1,nt
        Location(1,NTT) = xver
        Location(2,NTT) = yver
        Location(3,NTT) = zver
        
        omega=30.
        d1 = 1000000._wp
        ii = 0
        jj = 0
        kk = 0

        NNormalN(1,NTT) = n(1) 
        NNormalN(2,NTT) = n(2)
        NNormalN(3,NTT) = n(3)
        !boundary condition judgement (perfect reflecction)
        If (Zver<(200._wp)) then
            Zver=400.-abs(Zver)
            Nnew(3)=-Nnew(3)
            !write(91, *) NTT
            ReflectionPointMax(ReflectionTimeIndex) = NTT
            ReflectionTimeIndex = ReflectionTimeIndex + 1
        endif
        !write(3,'(3e16.8)') xver, yver, zver
        do i = 1,nx-1
            do j = 1,ny-1
                do,k = 1,nz-1
                    d = dabs((xver-xc(i,j,k))**2+(yver-yc(i,j,k))**2+(zver-zc(i,j,k))**2)
                    if (d<d1) then
                        d1 = d
                        ii = i
                        jj = j
                        kk = k
                        !print*, KK
                    endif
                end do
            end do
        end do


        !print*,xver
        !print*, kk
        do i=1,4
            do j=1,4
                do k=1,4
                    !print*, ii-2+i,jj-2+j,kk-2+k
                    XX(i,j,k) = x(ii-2+i,jj-2+j,kk-2+k)
                    YY(i,j,k) = y(ii-2+i,jj-2+j,kk-2+k)
                    ZZ(i,j,k) = z(ii-2+i,jj-2+j,kk-2+k)
                    !UU(i,j,k) = Uwind(ii-2+i,jj-2+j,kk-2+k)
                    !CC(i,j,k) = C0(ii-2+i,jj-2+j,kk-2+k)
                    !print*, ii-2+i,jj-2+j,kk-2+k
                end do
            end do
        end do
        !Output Gird Coordinates



        !
        !calculate Uwind and C0 for 4X4X4 do    ian
        do i=1,4
            do j=1,4
                do k=1,4
                    zt1 = zz(i,j,k) 
                    yt1 = YY(i,j,k)
                    xt1 = XX(i,j,k)
                    xt = real(xt1)
                    yt = real(yt1)
                    zt = real(zt1)
                    UU(i,j,k)  = Uwindspeed(xt1,yt1,zt1)!(zt-2000)**2./100000!w(2)
                    VV(i,j,k)  = Vwindspeed(xt1,yt1,zt1)!log(1.+abs(zt-2000))/100000
                    CC(i,j,k)  = Soundspeed(xt1,yt1,zt1)!343._wp-dlog(1._wp+abs(zt-2000))!340-zt/1000!C0(ii-2+i,jj-2+j,kk-2+k)
                end do
            end do
        end do  
        Vout(:,:,:,1,NTT) = UU
        Vout(:,:,:,2,NTT) = VV
        Vout(:,:,:,3,NTT) = 0.0
        C0out(:,:,:,NTT) = CC
        !XgridLoc(i,j,k,Cor_(1=x,2=y,3=z),Time)
        XgridLoc(:,:,:,1,NTT) = XX
        XgridLoc(:,:,:,2,NTT) = YY
        XgridLoc(:,:,:,3,NTT) = ZZ
        !write(*,*) UU, w(1), w(2), xt
        
        do i=1,4
            do j=1,4
                do k=1,4
                    !For D()/Dz Consideration, get related coordinates
                    xp = x(ii-1+i,jj-1+j,kk+k)      !get coordinate of up point
                    yp = y(ii-1+i,jj-1+j,kk+k)
                    zp = z(ii-1+i,jj-1+j,kk+k)
                    xm = x(ii-1+i,jj-1+j,kk+k-2)    !get coordinate of bottom point
                    ym = y(ii-1+i,jj-1+j,kk+k-2)
                    zm = z(ii-1+i,jj-1+j,kk+k-2)
                    !For D()/Dy consideration, get related coordinates
                    xl = x(ii-1+i, jj+j, kk+k-1)    !get coordinates of left point
                    yl = y(ii-1+1, jj+j, kk+k-1)    
                    zl = z(ii-1+i, jj+j, kk+k-1)
                    xr = x(ii-1+i, jj-2+j, kk+k-1)  !get coordinate of right point
                    yr = y(ii-1+i, jj-2+j, kk+k-1)
                    zr = z(ii-1+i, jj-2+j, kk+k-1)
                    !For D()/Dx consideration, get related coordinates
                    xf = x(ii+i, jj-1+j, kk+k-1)
                    yf = y(ii+i, jj-1+j, kk+k-1)
                    zf = z(ii+i, jj-1+j, kk+k-1)
                    xb = x(ii-2+i, jj-1+j, kk+k-1)
                    yb = y(ii-2+i, jj-1+j, kk+k-1)
                    zb = z(ii-2+i, jj-1+j, kk+k-1)
                    !get Uwind for different points
                    Uup   = Uwindspeed(xp,yp,zp)
                    Udown = Uwindspeed(xm,ym,zm)
                    Uleft = Uwindspeed(xl,yl,zl)
                    Uright= Uwindspeed(xr,yr,zr)
                    Ufront= Uwindspeed(xf,yf,zf)
                    Uback = Uwindspeed(xb,yb,zb)
                    !Get Vwind for different points
                    Vup   = Vwindspeed(xp,yp,zp)
                    Vdown = Vwindspeed(xm,ym,zm)
                    Vleft = Vwindspeed(xl,yl,zl)
                    Vright= Vwindspeed(xr,yr,zr)
                    Vfront= Vwindspeed(xf,yf,zf)
                    Vback = Vwindspeed(xb,yb,zb)
                    !Get Speed of Sound for different points
                    Cup   = Soundspeed(xp,yp,zp)
                    Cdown = Soundspeed(xm,ym,zm)
                    Cleft = Soundspeed(xl,yl,zl)
                    Cright= Soundspeed(xr,yr,zr)
                    Cfront= Soundspeed(xf,yf,zf)
                    Cback = Soundspeed(xb,yb,zb)

                    Duz(i,j,k) = (Uup-Udown)    / 200.
                    Dcz(i,j,k) = (Cup-Cdown)    / 200.
                    Dvz(i,j,k) = (Vup-Vdown)    / 200.
                    Duy(i,j,k) = (Uleft-Uright) / 200.
                    Dcy(i,j,k) = (Cleft-Cright) / 200.
                    Dvy(i,j,k) = (Vleft-Vright) / 200.
                    Dux(i,j,k) = (Ufront-Uback) / 200.
                    Dcx(i,j,k) = (Cfront-Cback) / 200.
                    Dvx(i,j,k) = (Vfront-Vback) / 200.

                end do
            end do
        end do

        Xnx = XX(:,1,1)
        Yny = YY(1,:,1)
        Znz = ZZ(1,1,:)


        !call interpD(Xnx(2:3),Yny(2:3),Znz(2:3),xver,yver,zver, UU(2:3,2:3,2:3),ValU)
        call interpP(Xnx,Yny,Znz,xver,yver,zver,CC,ValC)
    
        call interpP(Xnx,Yny,Znz,xver,yver,zver,UU,ValU)

        call interpP(Xnx,Yny,Znz,xver,yver,zver,VV,ValV)

        !call interpP(Xnx,Yny,Znz,xver,yver,zver,WW,valW)
        call interpP(Xnx,Yny,Znz,xver,yver,zver,Dwz,Valdwz)
        call interpP(Xnx,Yny,Znz,xver,yver,zver,Dvz,Valdvz)
        call interpP(Xnx,Yny,Znz,xver,yver,zver,Duz,Valduz)

        !call interpP(Xnx,Yny,Znz,xver,yver,zver,WW,Valdwz)

        call interpP(Xnx,Yny,Znz,xver,yver,zver,Dwy,Valdwy)
        call interpP(Xnx,Yny,Znz,xver,yver,zver,Dvy,Valdvy)
        call interpP(Xnx,Yny,Znz,xver,yver,zver,Duy,Valduy)

        call interpP(Xnx,Yny,Znz,xver,yver,zver,Dwx,Valdwx)
        call interpP(Xnx,Yny,Znz,xver,yver,zver,Dvx,Valdvx)
        call interpP(Xnx,Yny,Znz,xver,yver,zver,Dux,Valdux)

        call interpP(Xnx,Yny,Znz,xver,yver,zver,Dcx,Valdcx)
        call interpP(Xnx,Yny,Znz,xver,yver,zver,Dcy,Valdcy)
        call interpP(Xnx,Yny,Znz,xver,yver,zver,Dcz,Valdcz)
        
        !print*, Valduz
        !**********************************************************************
        ! start here
        ! For Gainesville Area
        ! 1 degree latitude 111.32 km
        ! 1 degree longtitude 96.7421km

        !you have all the parameters to find n^bar
        !for one time step, your finite difference method for n^bar:

        Knormal = omega/ValC
        Kx = Knormal*Nnew(1)
        Ky = Knormal*Nnew(2)
        Kz = Knormal*Nnew(3)

        Term2(1) = Valdux*kx+Valdvx*ky+Valdwx*kz!dudx*kx+dvdx*ky+0._wp
        Term2(2) = Valduy*kx+Valdvy*ky+Valdwy*kz!dudy*kx+dvdy*ky+0._wp
        Term2(3) = Valduz*kx+Valdvz*ky+Valdvz*kz+0._wp

        Term1(1) = Valdcx*Knormal! Valdcx * Knormmal
        Term1(2) = Valdcy*Knormal! Valdcy * Knormal
        Term1(3) = Valdcz*Knormal! 

        !Knew = K(nt)+dt*(-Term1-Term2)
        !print*,valduz
        Nnew = Nnew+dt/Knormal*(-Term1-Term2-sum(Nnew*(-Term1-Term2))*Nnew)
        Nnew(1) = Nnew(1)/(sqrt((Nnew(1)**2+Nnew(2)**2+Nnew(3)**2)))
        Nnew(2) = Nnew(2)/(sqrt((Nnew(1)**2+Nnew(2)**2+Nnew(3)**2)))
        Nnew(3) = Nnew(3)/(sqrt((Nnew(1)**2+Nnew(2)**2+Nnew(3)**2)))
        Xrecord(NTT) = xver
        Yrecord(NTT) = yver
        Zrecord(NTT) = zver

        xver = xver+dt*(ValC*Nnew(1)+ValU)
        yver = yver+dt*(ValC*Nnew(2)+ValV)+0.!ValV
        zver = zver+dt*(ValC*Nnew(3))+0.!ValW
        print*,  zver,ValU, ValC

    enddo

    ! write(91,*) ReflectionTimeIndex
    ! if (ReflectionTimeIndex > 1) then
    !     write(91,*) ReflectionPointMax(1:ReflectionTimeIndex)
    ! endif
    ! close(91)
    call system ( "mkdir -p " // 'Data' )
    open(unit = 3, file ='./Data/DataXYZ'//trim(str2phi)//'-'//trim(str2theta)//'.txt',form='formatted')
    write(3, '(3e16.8)') Location
    close(3)
    
    !open(51,file = '/home/tianshu/Documents/DATATRANSFER/VelocityGeoOut.txt')
    !write(51, '(64e16.8)') Vout
    !close(51)
    !open(52,file = '/home/tianshu/Documents/DATATRANSFER/SpeedofSoundGeoOut.txt')
    !write(52, '(64e16.8)') C0out
    !close(52)
    !open(53,file = '/home/tianshu/Documents/DATATRANSFER/XgridLoc.txt')
    !write(53, '(32e16.8)') XgridLoc
    !close(53)
    !open(54, file = '/home/tianshu/Documents/DATATRANSFER/Normal.txt')
    !write(54, '(3e16.8)') NNormalN
    !close(54)
    ! open(54, file = 'XverLoc.txt')
    ! write(54,'3e16.4')
    ! close(54)

    !write(*,*) w(1), w(2), ZT , ZZ(4,4,4), (XT/1000.)
    
    
    !do i=1,4
    !    do j=1,4
    !        do k=1,4
                !write(*,*) Duz !ValC, ValU, Valduz, Valdcz, w(2)
    !        end do
    !    end do
    !end do

end program test

    
    
!***************************************************************************************
function  Soundspeed(x,y,z)
  use bspline_kinds_module, only: wp
  implicit none
  real               ::  xt, yt, zt, TT
  real(wp)           ::  Soundspeed, x,y,z
  real, dimension(2) ::  T
  real, dimension(7) ::  AP
  real, dimension(9) ::  DDD
  
  AP = 1.0
  xt = real(x)
  yt = real(y)
  zt = real(z)
  call GTD7(1990,5.0,zt/1000.,(29.6516344+yt/1000./111.32),&
       &(-82.3248262+xt/1000./96.74),12.0,150.0,150.0,AP,48,DDD,T)
  TT = 273+T(2)
  Soundspeed = sqrt(TT*1.4*287)!.3_dp_t*Z
  !Soundspeed = 323._dp_t+Z/10._dp_t!.2_dp_t*z
  !print*, TT, zt
end function

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function  Uwindspeed(x,y,z)
  use bspline_kinds_module, only: wp
  implicit none
  real               ::  xt, yt, zt
  real(wp)           ::  Uwindspeed, x,y,z
  real, dimension(2) ::  w
  
  xt = real(x)
  yt = real(y)
  zt = real(z)
  
  call GWS5(1990,5.0,zt/1000.,(29.6516344+yt/1000./111.32),&
       &(-82.3248262+xt/1000./96.74),12.0,150.0,150.0,(/1.0,1.0/),w)
  
  Uwindspeed = w(2)
  !Uwindspeed = (z-2000)**2./100000
end function

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function  Vwindspeed(x,y,z)
  use bspline_kinds_module, only: wp
  implicit none
  real               ::  xt, yt, zt
  real(wp)           ::  Vwindspeed, x,y,z
  real, dimension(2) ::  w
  
  xt = real(x)
  yt = real(y)
  zt = real(z)
  
  call GWS5(1990,5.0,zt/1000.,(29.6516344+yt/1000./111.32),&
       &(-82.3248262+xt/1000./96.74),12.0,150.0,150.0,(/1.0,1.0/),w)
  
  Vwindspeed = w(1)
  !Vwindspeed = 0.!(z-2000)**2./1000000

end function

!***************************************************************************************
    
    
!program test
!    use bspline_kinds_module, only: wp
!    implicit none
!    
!    integer     :: nx, ny, nz, nt, i, j , k, yd
!    integer     :: ii, jj, kk, NTT
!    real        :: ZT, XT, YT
!    real(wp)    :: d, d1, xver, yver, zver, dudx, dvdx, dcdx, dudy, dvdy, dcdy, dvdz, dt
!    real(wp)    :: ValC, ValU, ValV, Valduz, Valdcz, Valdvz
!    real(wp)    :: omega, Kx, Ky, Kz, Knormal
!    real(wp), dimension(4) ::  Xnx,Yny,Znz
!    real(wp), allocatable  ::  x(:,:,:),y(:,:,:),z(:,:,:), dudz(:,:,:), dcdz(:,:,:)
!    real(wp), allocatable  ::  xc(:,:,:),yc(:,:,:),zc(:,:,:)
!    real(wp), allocatable  ::  C0(:,:,:),Uwind(:,:,:),Vwind(:,:,:)
!    real(wp), allocatable  ::  XX(:,:,:),YY(:,:,:),ZZ(:,:,:),UU(:,:,:),CC(:,:,:), VV(:,:,:)
!    real(wp), allocatable  ::  Xrecord(:), Yrecord(:), Zrecord(:)
!    real, dimension(2)     ::  w ,T
!    real(wp), dimension(3) ::  n, Nnew, Term1, Term2
!    real(wp), dimension(4,4,4) :: Duz, Dcz, Dvz
!    real(wp)  ::  soundspeed, Uwindspeed, Vwindspeed
!    real(wp)  ::  xp,yp,zp,xm,ym,zm,Cup,Cdown,Uup,Udown,Vup,Vdown,xt1,yt1,zt1
!    real,  dimension(9)  ::DDD
!    real,  dimension(7)  ::APP
!    
!    App=1
!    DDD=0
!    
!    call GTD7(1999,1800.,13.,22.2,82.2,12.,150.,150.,APP,48,DDD,T)
!    call GWS5(1999,1800.,13.,22.2,82.2,12.,150.,150.,APP,w)
!    print*, W,T,DDD
!    
!    end program

    
