real(kind=16) function qasum(n,dx,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:30

!     .. scalar arguments ..

integer, intent(in)                      :: n
real(kind=16), intent(in out)             :: dx(*)
integer, intent(in)                      :: incx

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     qasum takes the sum of the absolute values.

!  further details
!  ===============

!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
real(kind=16) :: dtemp
integer :: i,m,mp1,nincx
!     ..
!     .. intrinsic functions ..
intrinsic abs,mod
!     ..
qasum = 0.0_16
dtemp = 0.0_16
if (n <= 0 .or. incx <= 0) return
if (incx == 1) then
!        code for increment equal to 1
  
  
!        clean-up loop
  
  m = mod(n,6)
  if (m /= 0) then
    do i = 1,m
      dtemp = dtemp + abs(dx(i))
    end do
    if (n < 6) then
      qasum = dtemp
      return
    end if
  end if
  mp1 = m + 1
  do i = mp1,n,6
    dtemp = dtemp + abs(dx(i)) + abs(dx(i+1)) +  &
        abs(dx(i+2)) + abs(dx(i+3)) + abs(dx(i+4)) + abs(dx(i+5))
  end do
else
  
!        code for increment not equal to 1
  
  nincx = n*incx
  do i = 1,nincx,incx
    dtemp = dtemp + abs(dx(i))
  end do
end if
qasum = dtemp
return
end function qasum
