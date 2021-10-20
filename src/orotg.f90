subroutine orotg(ca,cb,c,s)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:37

!     .. scalar arguments ..

complex(kind=16), intent(out)              :: ca
complex(kind=16), intent(in)               :: cb
real(kind=16), intent(out)            :: c
complex(kind=16), intent(out)              :: s


!     ..

!  purpose
!  =======

!     orotg determines a double complex givens rotation.

!  =====================================================================

!     .. local scalars ..
complex(kind=16) :: alpha
real(kind=16) :: norm,scale
!     ..
!     .. intrinsic functions ..
intrinsic abs,dcmplx,conjg,sqrt
!     ..
if (abs(ca) == 0.0_16) then
  c = 0.0_16
  s = (1.0_16,0.0_16)
  ca = cb
else
  scale = abs(ca) + abs(cb)
  norm = scale*sqrt((abs(ca/dcmplx(scale,0.0_16)))**2+  &
      (abs(cb/dcmplx(scale,0.0_16)))**2)
  alpha = ca/abs(ca)
  c = abs(ca)/norm
  s = alpha*conjg(cb)/norm
  ca = alpha*norm
end if
return
end subroutine orotg
