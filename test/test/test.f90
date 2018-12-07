program test
    use bspline_kinds_module, only: wp
    implicit none
    
    integer     :: nx, ny, nz, i, j , k,yd
    integer     :: ii, jj, kk
    real        :: ZT, XT, YT
    real(wp)    :: d, d1, xver, yver, zver
    real(wp), allocatable  ::  x(:,:,:),y(:,:,:),z(:,:,:)
    real(wp), allocatable  ::  xc(:,:,:),yc(:,:,:),zc(:,:,:)
    real(wp), allocatable  ::  C0(:,:,:),Uwind(:,:,:),Vwind(:,:,:)
    real(wp), allocatable  ::  XX(:,:,:),YY(:,:,:),ZZ(:,:,:),UU(:,:,:),CC(:,:,:)
    real, dimension(2)  ::  w
    
    nx=1000
    ny=50
    nz=100
    dt=.001_wp
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
    
    do i=1,nx
      x(i,:,:) = 200000._wp*(dble(i-1)/dble(nx-1))
    end do
    do i=1,ny
      y(:,i,:) = 5000._wp*(dble(i-1)/dble(ny-1))
    end do
    do i=1,nz
      z(:,:,i) = 100000._wp*(dble(i-1)/dble(nz-1))
    end do
    do i=1,nz
      C0(:,:,i) = 343._wp-dlog(1._wp+z(:,:,i))
    end do
    do i=1,nz
      Uwind(:,:,i) = .005_wp*(z(:,:,i))
    end do
    Vwind = 0._wp
    
    dudx = 0._wp
    dudy = 0._wp
    dudz = Uwind(:,:,2:nz) - Uwind(:,:,1:nz-1)
    dvdx = 0._wp
    dvdy = 0._wp
    dvdz = 0._wp
    dcdx = 0._wp
    dcdy = 0._wp
    dcdz = C0(:,:,2:nz) - C0(:,:,1:nz-1)
    
    xver = 100000._wp*.33_wp
    yver = 5000._wp*.78_wp
    zver = 10000._wp*.9_wp
    n(1) = 1._wp/(2._wp)**.5_wp
    n(2) = 1._wp/(2._wp)**.5_wp
    n(3) = 0._wp
    
    call cellcenter(x,y,z,xc,yc,zc,nx,ny,nz)
    
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
    !write(*,*) Uwind(ii,jj,kk)
    
    do i=1,4
        do j=1,4
            do k=1,4
                
                XX(i,j,k) = x(ii-2+i,jj-2+j,kk-2+k)
                YY(i,j,k) = y(ii-2+i,jj-2+j,kk-2+k)
                ZZ(i,j,k) = z(ii-2+i,jj-2+j,kk-2+k)
                UU(i,j,k) = Uwind(ii-2+i,jj-2+j,kk-2+k)
                CC(i,j,k) = C0(ii-2+i,jj-2+j,kk-2+k)
            end do
        end do
    end do
    do i=1,2
        do j=1,2
            do k=1,2              
                Duz(i,j,k) = Uwind(ii+i-1,jj+j-1,kk+k)-&
                             &Uwind(ii+i-1,jj+j-1,kk-2+k)
                Dcz(i,j,k) = C0i+i-1,jj+j-1,kk+k)-&
                             &C0i+i-1,jj+j-1,kk-2+k)
            end do
        end do
    end do
    Xnx = XX(:,1,1)
    Yny = YY(1,:,1)
    Znz = ZZ(1,1,:)
    
    call interpP(Xnx,Yny,Znz,xver,yver,zver,CC,ValC)
    call interpP(Xnx,Yny,Znz,xver,yver,zver,UU,ValU)
    call interpD(Xnx(2:3),Yny(2:3),Znz(2:3),xver,yver,zver,Duz,Valduz)
    call interpD(Xnx(2:3),Yny(2:3),Znz(2:3),xver,yver,zver,Dcz,Valdcz)
    
    
    ! For Gainesville Area
    ! 1 degree latitude 111.32 km
    ! 1 degree longtitude 96.7421km
    
    !do i=1,4
    !    do j=1,4
    !        do k=1,4
    !            ZT = ZZ(i,j,k) 
    !            YT = YY(i,j,k)
    !            XT = XX(i,j,k)
    !            call GWS5(YD,-5.0,ZT,(29.6516344+XT/1000./111.32),-82.3248262,12.0,150.0,150.0,(/4.0,4.0/),w(1:2))            
    !        end do
    !    end do
    !end do
    
    write(*,*) w(1), w(2), ZT , ZZ(4,4,4), (XT/1000.)
    
    !do i=1,4
    !    do j=1,4
    !        do k=1,4
    !            write(*,*) XX(i,j,k)
    !        end do
    !    end do
    !end do
end program test