real(kind=16) function qdot(n,dx,incx,dy,incy)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:30

!     .. scalar arguments ..

integer, intent(in)                      :: n
real(kind=16), intent(in)             :: dx(*)
integer, intent(in)                      :: incx
real(kind=16), intent(in)             :: dy(*)
integer, intent(in)                      :: incy

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     qdot forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.

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
qdot = 0.0_16
dtemp = 0.0_16
if (n <= 0) return
if (incx == 1 .and. incy == 1) then
  
!        code for both increments equal to 1
  
  
!        clean-up loop
  
  m = mod(n,5)
  if (m /= 0) then
    do i = 1,m
      dtemp = dtemp + dx(i)*dy(i)
    end do
    if (n < 5) then
      qdot=dtemp
      return
    end if
  end if
  mp1 = m + 1
  do i = mp1,n,5
    dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) +  &
        dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
  end do
else
  
!        code for unequal increments or equal increments
!          not equal to 1
  
  ix = 1
  iy = 1
  if (incx < 0) ix = (-n+1)*incx + 1
  if (incy < 0) iy = (-n+1)*incy + 1
  do i = 1,n
    dtemp = dtemp + dx(ix)*dy(iy)
    ix = ix + incx
    iy = iy + incy
  end do
end if
qdot = dtemp
return
end function qdot
