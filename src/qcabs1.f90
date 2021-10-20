real(kind=16) function qabs1(z)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:30

!     .. scalar arguments ..

complex(kind=16), intent(in out)           :: z

!     ..
!     ..
!  purpose
!  =======

!  qabs1 computes absolute value of a double complex number

!  =====================================================================

!     .. intrinsic functions ..
intrinsic abs,aimag

qabs1 = abs(real(z,kind=16)) + abs(aimag(z))
return
end function qabs1
