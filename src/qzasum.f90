real(kind=16) function qzasum(n,zx,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:32

!     .. scalar arguments ..

integer, intent(in)                      :: n
complex(kind=16), intent(in out)           :: zx(*)
integer, intent(in)                      :: incx

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     qzasum takes the sum of the absolute values.

!  further details
!  ===============

!     jack dongarra, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
real(kind=16) :: stemp
integer :: i,nincx
!     ..
!     .. external functions ..
!     ..
qzasum = 0.0_16
stemp = 0.0_16
if (n <= 0 .or. incx <= 0) return
if (incx == 1) then
  
!        code for increment equal to 1
  
  do i = 1,n
    stemp = stemp + abs(zx(i))
  end do
else
  
!        code for increment not equal to 1
  
  nincx = n*incx
  do i = 1,nincx,incx
    stemp = stemp + abs(zx(i))
  end do
end if
qzasum = stemp
return
end function qzasum
