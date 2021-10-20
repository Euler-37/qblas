subroutine qscal(n,da,dx,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:31

!     .. scalar arguments ..

integer, intent(in)                      :: n
real(kind=16), intent(in)             :: da
real(kind=16), intent(out)            :: dx(*)
integer, intent(in)                      :: incx


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     qscal scales a vector by a constant.
!     uses unrolled loops for increment equal to one.

!  further details
!  ===============

!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
integer :: i,m,mp1,nincx
!     ..
!     .. intrinsic functions ..
intrinsic mod
!     ..
if (n <= 0 .or. incx <= 0) return
if (incx == 1) then
  
!        code for increment equal to 1
  
  
!        clean-up loop
  
  m = mod(n,5)
  if (m /= 0) then
    do i = 1,m
      dx(i) = da*dx(i)
    end do
    if (n < 5) return
  end if
  mp1 = m + 1
  do i = mp1,n,5
    dx(i) = da*dx(i)
    dx(i+1) = da*dx(i+1)
    dx(i+2) = da*dx(i+2)
    dx(i+3) = da*dx(i+3)
    dx(i+4) = da*dx(i+4)
  end do
else
  
!        code for increment not equal to 1
  
  nincx = n*incx
  do i = 1,nincx,incx
    dx(i) = da*dx(i)
  end do
end if
return
end subroutine qscal
