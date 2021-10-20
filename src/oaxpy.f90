subroutine oaxpy(n,za,zx,incx,zy,incy)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:35

!     .. scalar arguments ..

integer, intent(in)                      :: n
complex(kind=16), intent(in)               :: za
complex(kind=16), intent(in)               :: zx(*)
integer, intent(in)                      :: incx
complex(kind=16), intent(out)              :: zy(*)
integer, intent(in)                      :: incy


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     oaxpy constant times a vector plus a vector.

!  further details
!  ===============

!     jack dongarra, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
integer :: i,ix,iy
!     ..
!     .. external functions ..
!     ..
if (n <= 0) return
if (abs(za) == 0.0_16) return
if (incx == 1 .and. incy == 1) then

    !        code for both increments equal to 1

    do i = 1,n
        zy(i) = zy(i) + za*zx(i)
    end do
else

    !        code for unequal increments or equal increments
    !          not equal to 1

    ix = 1
    iy = 1
    if (incx < 0) ix = (-n+1)*incx + 1
    if (incy < 0) iy = (-n+1)*incy + 1
    do i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
    end do
end if

return
end subroutine oaxpy
