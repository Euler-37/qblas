subroutine qrot(n,dx,incx,dy,incy,c,s)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:31

!     .. scalar arguments ..

integer, intent(in)                      :: n
real(kind=16), intent(in out)         :: dx(*)
integer, intent(in)                      :: incx
real(kind=16), intent(in out)         :: dy(*)
integer, intent(in)                      :: incy
real(kind=16), intent(in)             :: c
real(kind=16), intent(in)             :: s


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     qrot applies a plane rotation.

!  further details
!  ===============

!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
real(kind=16) :: dtemp
integer :: i,ix,iy
!     ..
if (n <= 0) return
if (incx == 1 .and. incy == 1) then
  
!       code for both increments equal to 1
  
  do i = 1,n
    dtemp = c*dx(i) + s*dy(i)
    dy(i) = c*dy(i) - s*dx(i)
    dx(i) = dtemp
  end do
else
  
!       code for unequal increments or equal increments not equal
!         to 1
  
  ix = 1
  iy = 1
  if (incx < 0) ix = (-n+1)*incx + 1
  if (incy < 0) iy = (-n+1)*incy + 1
  do i = 1,n
    dtemp = c*dx(ix) + s*dy(iy)
    dy(iy) = c*dy(iy) - s*dx(ix)
    dx(ix) = dtemp
    ix = ix + incx
    iy = iy + incy
  end do
end if
return
end subroutine qrot
