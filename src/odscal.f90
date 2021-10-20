subroutine odscal(n,da,zx,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:36

!     .. scalar arguments ..

integer, intent(in)                      :: n
real(kind=16), intent(in out)         :: da
complex(kind=16), intent(out)              :: zx(*)
integer, intent(in)                      :: incx


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     odscal scales a vector by a constant.

!  further details
!  ===============

!     jack dongarra, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
integer :: i,nincx
!     ..
!     .. intrinsic functions ..
intrinsic dcmplx
!     ..
if (n <= 0 .or. incx <= 0) return
if (incx == 1) then
  
!        code for increment equal to 1
  
  do i = 1,n
    zx(i) = dcmplx(da,0.0_16)*zx(i)
  end do
else
  
!        code for increment not equal to 1
  
  nincx = n*incx
  do i = 1,nincx,incx
    zx(i) = dcmplx(da,0.0_16)*zx(i)
  end do
end if
return
end subroutine odscal
