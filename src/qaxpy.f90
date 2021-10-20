subroutine qaxpy(n,da,dx,incx,dy,incy)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:30

!     .. scalar arguments ..

integer, intent(in)                      :: n
real(kind=16), intent(in)             :: da
real(kind=16), intent(in)             :: dx(*)
integer, intent(in)                      :: incx
real(kind=16), intent(out)            :: dy(*)
integer, intent(in)                      :: incy


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     qaxpy constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.

!  further details
!  ===============

!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
integer :: i,ix,iy,m,mp1
!     ..
!     .. intrinsic functions ..
intrinsic mod
!     ..
if (n <= 0) return
if (da == 0.0_16) return
if (incx == 1 .and. incy == 1) then
  
!        code for both increments equal to 1
  
  
!        clean-up loop
  
  m = mod(n,4)
  if (m /= 0) then
    do i = 1,m
      dy(i) = dy(i) + da*dx(i)
    end do
  end if
  if (n < 4) return
  mp1 = m + 1
  do i = mp1,n,4
    dy(i) = dy(i) + da*dx(i)
    dy(i+1) = dy(i+1) + da*dx(i+1)
    dy(i+2) = dy(i+2) + da*dx(i+2)
    dy(i+3) = dy(i+3) + da*dx(i+3)
  end do
else
  
!        code for unequal increments or equal increments
!          not equal to 1
  
  ix = 1
  iy = 1
  if (incx < 0) ix = (-n+1)*incx + 1
  if (incy < 0) iy = (-n+1)*incy + 1
  do i = 1,n
    dy(iy) = dy(iy) + da*dx(ix)
    ix = ix + incx
    iy = iy + incy
  end do
end if
return
end subroutine qaxpy
