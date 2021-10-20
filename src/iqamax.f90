integer function iqamax(n,dx,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:33

!     .. scalar arguments ..

integer, intent(in)                      :: n
real(kind=16), intent(in out)             :: dx(*)
integer, intent(in)                      :: incx

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     iqamax finds the index of element having max. absolute value.

!  further details
!  ===============

!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
real(kind=16) :: dmax
integer :: i,ix
!     ..
!     .. intrinsic functions ..
intrinsic abs
!     ..
iqamax = 0
if (n < 1 .or. incx <= 0) return
iqamax = 1
if (n == 1) return
if (incx == 1) then
  
!        code for increment equal to 1
  
  dmax = abs(dx(1))
  do i = 2,n
    if (abs(dx(i)) > dmax) then
      iqamax = i
      dmax = abs(dx(i))
    end if
  end do
else
  
!        code for increment not equal to 1
  
  ix = 1
  dmax = abs(dx(1))
  ix = ix + incx
  do i = 2,n
    if (abs(dx(ix)) > dmax) then
      iqamax = i
      dmax = abs(dx(ix))
    end if
    ix = ix + incx
  end do
end if
return
end function iqamax
