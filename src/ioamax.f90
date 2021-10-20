integer function ioamax(n,zx,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:33

!     .. scalar arguments ..

integer, intent(in)                      :: n
complex(kind=16), intent(in out)         :: zx(*)
integer, intent(in)                      :: incx

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     ioamax finds the index of element having max. absolute value.

!  further details
!  ===============

!     jack dongarra, 1/15/85.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
real(kind=16) :: dmax
integer :: i,ix
!     ..
!     .. external functions ..
!     ..
ioamax = 0
if (n < 1 .or. incx <= 0) return
ioamax = 1
if (n == 1) return
if (incx == 1) then
  
!        code for increment equal to 1
  
  dmax = abs(zx(1))
  do i = 2,n
    if (abs(zx(i)) > dmax) then
      ioamax = i
      dmax = abs(zx(i))
    end if
  end do
else
  
!        code for increment not equal to 1
  
  ix = 1
  dmax = abs(zx(1))
  ix = ix + incx
  do i = 2,n
    if (abs(zx(ix)) > dmax) then
      ioamax = i
      dmax = abs(zx(ix))
    end if
    ix = ix + incx
  end do
end if
return
end function ioamax
