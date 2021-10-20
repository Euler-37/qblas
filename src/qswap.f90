subroutine qswap(n,dx,incx,dy,incy)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:31

!     .. scalar arguments ..

integer, intent(in)                      :: n
real(kind=16), intent(in out)         :: dx(*)
integer, intent(in)                      :: incx
real(kind=16), intent(in out)         :: dy(*)
integer, intent(in)                      :: incy

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     interchanges two vectors.
!     uses unrolled loops for increments equal one.

!  further details
!  ===============

!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
real(kind=16) :: dtemp
integer :: i,ix,iy,m,mp1
!     ..
!     .. intrinsic functions ..
intrinsic mod
!     ..
if (n <= 0) return
if (incx == 1 .and. incy == 1) then
  
!       code for both increments equal to 1
  
  
!       clean-up loop
  
  m = mod(n,3)
  if (m /= 0) then
    do i = 1,m
      dtemp = dx(i)
      dx(i) = dy(i)
      dy(i) = dtemp
    end do
    if (n < 3) return
  end if
  mp1 = m + 1
  do i = mp1,n,3
    dtemp = dx(i)
    dx(i) = dy(i)
    dy(i) = dtemp
    dtemp = dx(i+1)
    dx(i+1) = dy(i+1)
    dy(i+1) = dtemp
    dtemp = dx(i+2)
    dx(i+2) = dy(i+2)
    dy(i+2) = dtemp
  end do
else
  
!       code for unequal increments or equal increments not equal
!         to 1
  
  ix = 1
  iy = 1
  if (incx < 0) ix = (-n+1)*incx + 1
  if (incy < 0) iy = (-n+1)*incy + 1
  do i = 1,n
    dtemp = dx(ix)
    dx(ix) = dy(iy)
    dy(iy) = dtemp
    ix = ix + incx
    iy = iy + incy
  end do
end if
return
end subroutine qswap
