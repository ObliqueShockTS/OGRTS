program test
    use bspline_kinds_module, only: wp
    implicit none
    
    integer     :: nx, ny, nz, nt, i, j , k,yd
    integer     :: ii, jj, kk, NTT
    real        :: ZT, XT, YT
    real(wp)    :: d, d1, xver, yver, zver, dudx, dvdx, dcdx, dudy, dvdy, dcdy, dvdz, dt
    real(wp)    :: ValC, ValU, Valduz, Valdcz
    real(wp)    :: omega, Kx, Ky, Kz, Knormal
    real(wp), dimension(4) ::  Xnx,Yny,Znz
    real(wp), allocatable  ::  x(:,:,:),y(:,:,:),z(:,:,:), dudz(:,:,:), dcdz(:,:,:)
    real(wp), allocatable  ::  xc(:,:,:),yc(:,:,:),zc(:,:,:)
    real(wp), allocatable  ::  C0(:,:,:),Uwind(:,:,:),Vwind(:,:,:)
    real(wp), allocatable  ::  XX(:,:,:),YY(:,:,:),ZZ(:,:,:),UU(:,:,:),CC(:,:,:)
    real(wp), allocatable  ::  Xrecord(:), Yrecord(:), Zrecord(:)
    real, dimension(2)     ::  w
    real(wp), dimension(3) ::  n, Nnew, Term1, Term2
    real(wp), dimension(4,4,4) :: Duz, Dcz

    nx=1000
    ny=50
    nz=100
    nt=100
    dt=1._wp
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
    allocate (CC(4,4,4))
    allocate (dudz(nx,ny,(nz-1)))
    allocate (dcdz(nx,ny,(nz-1)))
    allocate (Xrecord(nt))
    allocate (Yrecord(nt))
    allocate (Zrecord(nt))

    
    do i=1,nx
      x(i,:,:) = 1000000._wp*(dble(i-1)/dble(nx-1))
    end do
    do i=1,ny
      y(:,i,:) = 5000._wp*(dble(i-1)/dble(ny-1))
    end do
    do i=1,nz
      z(:,:,i) = 10000._wp*(dble(i-1)/dble(nz-1))
    end do
    do i=1,nz
      C0(:,:,i) = 343._wp-dlog(1._wp+z(:,:,i))
    end do
    
    do i=1,nz
        Uwind(:,:,i) = z(:,:,i)**2./1000000.
    end do
    
    Vwind = 0._wp
    
    dudx = 0._wp
    dudy = 0._wp
    !dudz = Uwind(:,:,2:nz) - Uwind(:,:,1:nz-1)
    dvdx = 0._wp
    dvdy = 0._wp
    dvdz = 0._wp
    dcdx = 0._wp
    dcdy = 0._wp
    
    xver = 100000._wp*.13_wp
    yver = 5000._wp*.78_wp
    zver = 10000._wp*.3_wp
    n(1) = 1._wp/(2._wp)**.5_wp
    n(2) = 0._wp
    n(3) = 1._wp/(2._wp)**.5_wp
    Nnew = n
    
     
    call cellcenter(x,y,z,xc,yc,zc,nx,ny,nz)
    open(unit = 3, file ='DataXYZ.txt',form='formatted')
    
do NTT=1,nt
    
    d1 = 1000000._wp
    ii = 0
    jj = 0
    kk = 0
    
    
    do i = 1,nx-1
        do j = 1,ny-1
            do,k = 1,nz-1
                d = dabs((xver-xc(i,j,k))**2+(yver-yc(i,j,k))**2+(zver-zc(i,j,k))**2)
                if (d<d1) then
                    d1 = d
                    ii = i
                    jj = j
                    kk = k
                endif
            end do
        end do
    end do
    !print*,xver
    !print*, ii,jj,kk
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
                zt = zz(i,j,k) 
                yt = yy(i,j,k)
                xt = xx(i,j,k)
                call gws5(yd,-5.0,zt/1000.,(29.6516344+xt/1000./111.32),(-82.3248262+yt/1000./96.74),12.0,150.0,150.0,(/4.0,4.0/),w(1:2))
                !call gws5(yd,-5.0,10.,29.6516344,-82.3248262,12.0,150.0,150.0,(/4.0,4.0/),w(1:2))
                UU(i,j,k)  = w(2)!zt**2./1000000!w(2)
                CC(i,j,k)  = 343._wp-dlog(1._wp+zt)!340-zt/1000!C0(ii-2+i,jj-2+j,kk-2+k)
            end do
        end do
     end do     
    !write(*,*) UU, w(1), w(2), xt
     
 
    do i=1,4
        do j=1,4
            do k=1,4              
                Duz(i,j,k) = (UWind(ii-1+i,jj-1+j,kk+k)-UWind(ii-1+i,jj-1+j,kk+k-2))/200.
                Dcz(i,j,k) = (C0(ii-1+i,jj-1+j,kk+k)-C0(ii-1+i,jj-1+j,kk+k-2))/200.
            end do
        end do
    end do
    
    Xnx = XX(:,1,1)
    Yny = YY(1,:,1)

    Znz = ZZ(1,1,:)
    
    call interpP(Xnx,Yny,Znz,xver,yver,zver,CC,ValC)
   
    call interpP(Xnx,Yny,Znz,xver,yver,zver,UU,ValU)
    
    call interpP(Xnx,Yny,Znz,xver,yver,zver,Duz,Valduz)
    
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
    Term2(3) = Valduz*kx+dvdz*ky+0._wp
    
    Term1(1) = dcdx*Knormal
    Term1(2) = dcdy*Knormal
    Term1(3) = Valduz*Knormal
    
    !Knew = K(nt)+dt*(-Term1-Term2)
    !print*,valduz
    Nnew = Nnew+dt/Knormal*(-Term1-Term2-sum(Nnew*(-Term1-Term2))*Nnew)
    Xrecord(NTT) = xver
    Yrecord(NTT) = yver
    Zrecord(NTT) = zver
    
    xver = xver+dt*((ValC*Nnew(1))+ValU)
    yver = yver+dt*(ValC*Nnew(2))+0.!ValV
    zver = zver+dt*(ValC*Nnew(3))+0.!ValW
    print*,  Valdcz, Nnew(3)
    write(3,*) xver, yver, zver
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