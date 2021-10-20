complex(kind=16) function odotu(n,zx,incx,zy,incy)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:35

!     .. scalar arguments ..

integer, intent(in)                      :: n
complex(kind=16), intent(in)               :: zx(*)
integer, intent(in)                      :: incx
complex(kind=16), intent(in)               :: zy(*)
integer, intent(in)                      :: incy

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     odotu forms the dot product of two vectors.

!  further details
!  ===============

!     jack dongarra, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
complex(kind=16) :: ztemp
integer :: i,ix,iy
!     ..
ztemp = (0.0_16,0.0_16)
odotu = (0.0_16,0.0_16)
if (n <= 0) return
if (incx == 1 .and. incy == 1) then
  
!        code for both increments equal to 1
  
  do i = 1,n
    ztemp = ztemp + zx(i)*zy(i)
  end do
else
  
!        code for unequal increments or equal increments
!          not equal to 1
  
  ix = 1
  iy = 1
  if (incx < 0) ix = (-n+1)*incx + 1
  if (incy < 0) iy = (-n+1)*incy + 1
  do i = 1,n
    ztemp = ztemp + zx(ix)*zy(iy)
    ix = ix + incx
    iy = iy + incy
  end do
end if
odotu = ztemp
return
end function odotu
