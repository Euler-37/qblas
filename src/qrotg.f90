subroutine qrotg(da,db,c,s)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:31

!     .. scalar arguments ..

real(kind=16), intent(in out)         :: da
real(kind=16), intent(in out)         :: db
real(kind=16), intent(out)            :: c
real(kind=16), intent(out)            :: s

!     ..

!  purpose
!  =======

!     qrotg construct givens plane rotation.

!  further details
!  ===============

!     jack dongarra, linpack, 3/11/78.

!  =====================================================================

!     .. local scalars ..
real(kind=16) :: r,roe,scale,z
!     ..
!     .. intrinsic functions ..
intrinsic abs,sign,sqrt
!     ..
roe = db
if (abs(da) > abs(db)) roe = da
scale = abs(da) + abs(db)
if (scale == 0.0_16) then
  c = 1.0_16
  s = 0.0_16
  r = 0.0_16
  z = 0.0_16
else
  r = scale*sqrt((da/scale)**2+ (db/scale)**2)
  r = sign(1.0_16,roe)*r
  c = da/r
  s = db/r
  z = 1.0_16
  if (abs(da) > abs(db)) z = s
  if (abs(db) >= abs(da) .and. c /= 0.0_16) z = 1.0_16/c
end if
da = r
db = z
return
end subroutine qrotg
