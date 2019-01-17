subroutine CellCenter(x,y,z,xc,yc,zc,nx,ny,nz)
    
    use bspline_kinds_module, only: wp
    implicit none
    
    integer,  intent(in)  ::  nx, ny, nz
    real(wp), intent(in)  ::  x(nx,ny,nz), y(nx,ny,nz), z(nx,ny,nz)
    real(wp), intent(out) ::  xc(nx-1,ny-1,nz-1), &
                              yc(nx-1,ny-1,nz-1), zc(nx-1,ny-1,nz-1)
    
    xc=.5_wp*(x(1:nx-1,1:ny-1,1:nz-1)+x(2:nx,2:ny,2:nz))
    yc=.5_wp*(y(1:nx-1,1:ny-1,1:nz-1)+y(2:nx,2:ny,2:nz))
    zc=.5_wp*(z(1:nx-1,1:ny-1,1:nz-1)+z(2:nx,2:ny,2:nz))
    
end subroutine CellCenter