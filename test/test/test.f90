program test
    use bspline_kinds_module, only: wp
    implicit none
    
    integer     :: nx, ny, nz, nt, i, j , k, yd
    integer     :: ii, jj, kk, NTT
    real        :: ZT, XT, YT
    real(wp)    :: d, d1, xver, yver, zver, dudx, dvdx, dcdx, dudy, dvdy, dcdy, dvdz, dt
    real(wp)    :: ValC, ValU, ValV, Valduz, Valdcz, Valdvz
    real(wp)    :: omega, Kx, Ky, Kz, Knormal
    real(wp), dimension(4) ::  Xnx,Yny,Znz
    real(wp), allocatable  ::  x(:,:,:),y(:,:,:),z(:,:,:), dudz(:,:,:), dcdz(:,:,:)
    real(wp), allocatable  ::  xc(:,:,:),yc(:,:,:),zc(:,:,:)
    real(wp), allocatable  ::  C0(:,:,:),Uwind(:,:,:),Vwind(:,:,:)
    real(wp), allocatable  ::  XX(:,:,:),YY(:,:,:),ZZ(:,:,:),UU(:,:,:),CC(:,:,:), VV(:,:,:)
    real(wp), allocatable  ::  Xrecord(:), Yrecord(:), Zrecord(:)
    real, dimension(2)     ::  w
    real(wp), dimension(3) ::  n, Nnew, Term1, Term2      
    real(wp), dimension(4,4,4) :: Duz, Dcz, Dvz
    real(wp)  ::  soundspeed, Uwindspeed, Vwindspeed
    real(wp)  ::  xp,yp,zp,xm,ym,zm,Cup,Cdown,Uup,Udown,Vup,Vdown,xt1,yt1,zt1

    nx=1000
    ny=20
    nz=100
    nt=50
    dt=0.8_wp
    omega=15.
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
    

    
    dudx = 0._wp
    dudy = 0._wp
    !dudz = Uwind(:,:,2:nz) - Uwind(:,:,1:nz-1)
    dvdx = 0._wp
    dvdy = 0._wp
    !dvdz = 0._wp
    dcdx = 0._wp
    dcdy = 0._wp
    
    xver = 100000._wp*.13_wp
    yver = 50000._wp*.08_wp
    zver = 20000.+10.
    n(1) = 1._wp/(2._wp)**.5_wp
    n(2) = 0._wp
    n(3) = 1._wp/(2._wp)**.5_wp
    Nnew = n
    
     
    call cellcenter(x,y,z,xc,yc,zc,nx,ny,nz)
    open(unit = 3, file ='DataXYZ.txt',form='formatted')
    
do NTT=1,nt
    omega=30.
    d1 = 1000000._wp
    ii = 0
    jj = 0
    kk = 0
    
    !boundary condition judgement
    If (Zver<(2000._wp)) then
        Zver=4000.-abs(Zver)
        Nnew(3)=-Nnew(3)
    endif
    write(3,*) xver, yver, zver
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
   
    !calculate Uwind and C0 for 4X4X4 domian
     do i=1,4
        do j=1,4
            do k=1,4
                zt1 = zz(i,j,k) 
                yt1 = yy(i,j,k)
                xt1 = xx(i,j,k)
                xt = real(xt1)
                yt = real(yt1)
                zt = real(zt1)
                !call gws5(yd,-5.0,zt/1000.,(29.6516344+xt/1000./111.32),...
                !...(-82.3248262+yt/1000./96.74),12.0,150.0,150.0,(/4.0,4.0/),w(1:2))
                !call GTD7(yd,-5.0,zt/1000.,(29.6516344+xt/1000./111.32),...
                !...(-82.3248262+yt/1000./96.74),12,F107A,F107,AP,MASS,D,T)
                !call gws5(yd,-5.0,10.,29.6516344,-82.3248262,12.0,150.0,150.0,(/4.0,4.0/),w(1:2))
                UU(i,j,k)  = Uwindspeed(xt1,yt1,zt1)!(zt-2000)**2./100000!w(2)
                VV(i,j,k)  = Vwindspeed(xt1,yt1,zt1)!log(1.+abs(zt-2000))/100000
                CC(i,j,k)  = Soundspeed(xt1,yt1,zt1)!343._wp-dlog(1._wp+abs(zt-2000))!340-zt/1000!C0(ii-2+i,jj-2+j,kk-2+k)
            end do
        end do
     end do  
    !write(*,*) UU, w(1), w(2), xt
     
 
    do i=1,4
        do j=1,4
            do k=1,4
                xp=x(ii-1+i,jj-1+j,kk+k)
                yp=y(ii-1+i,jj-1+j,kk+k)
                zp=z(ii-1+i,jj-1+j,kk+k)
                xm=x(ii-1+i,jj-1+j,kk+k-2)
                ym=y(ii-1+i,jj-1+j,kk+k-2)
                zm=z(ii-1+i,jj-1+j,kk+k-2)
                
                Uup   = Uwindspeed(xp,yp,zp)
                Udown = Uwindspeed(xm,ym,zm)
                
                Vup   = Vwindspeed(xp,yp,zp)
                Vdown = Vwindspeed(xm,ym,zm)
                
                Cup   = Soundspeed(xp,yp,zp)
                Cdown = Soundspeed(xm,ym,zm)
                
                Duz(i,j,k) = (Uup-Udown)/200.
                Dcz(i,j,k) = (Cup-Cdown)/200.
                Dvz(i,j,k) = (Vup-Vdown)/200.
            end do
        end do
    end do
    
    Xnx = XX(:,1,1)
    Yny = YY(1,:,1)

    Znz = ZZ(1,1,:)
    
    call interpP(Xnx,Yny,Znz,xver,yver,zver,CC,ValC)
   
    call interpP(Xnx,Yny,Znz,xver,yver,zver,UU,ValU)
    
    call interpP(Xnx,Yny,Znz,xver,yver,zver,VV,ValV)
    
    call interpP(Xnx,Yny,Znz,xver,yver,zver,Duz,Valduz)
    
    call interpP(Xnx,Yny,Znz,xver,yver,zver,Dvz,Valdvz)
    
    call interpP(Xnx,Yny,Znz,xver,yver,zver,Dcz,Valdcz)
    !print*, Valduz
    !**********************************************************************
    ! start here
    ! For Gainesville Area
    ! 1 degree latitude 111.32 km
    ! 1 degree longtitude 96.7421km
    
    !you have all the parameters to find n^bar
    !for one time step, your finite difference method for n^bar:
    
    Knormal=omega/ValC
    Kx=Knormal*Nnew(1)
    Ky=Knormal*Nnew(2)
    Kz=Knormal*Nnew(3)
    
    Term2(1) = dudx*kx+dvdx*ky+0._wp
    Term2(2) = dudy*kx+dvdy*ky+0._wp
    Term2(3) = Valduz*kx+Valdvz*ky+0._wp
    
    Term1(1) = dcdx*Knormal
    Term1(2) = dcdy*Knormal
    Term1(3) = Valdcz*Knormal
    
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
