subroutine oscal(n,za,zx,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:37

!     .. scalar arguments ..

integer, intent(in)                      :: n
complex(kind=16), intent(in)               :: za
complex(kind=16), intent(out)              :: zx(*)
integer, intent(in)                      :: incx


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     oscal scales a vector by a constant.

!  further details
!  ===============

!     jack dongarra, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

!  =====================================================================

!     .. local scalars ..
integer :: i,nincx
!     ..
if (n <= 0 .or. incx <= 0) return
if (incx == 1) then
  
!        code for increment equal to 1
  
  do i = 1,n
    zx(i) = za*zx(i)
  end do
else
  
!        code for increment not equal to 1
  
  nincx = n*incx
  do i = 1,nincx,incx
    zx(i) = za*zx(i)
  end do
end if
return
end subroutine oscal
