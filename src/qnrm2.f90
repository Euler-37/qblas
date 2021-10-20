real(kind=16) function qnrm2(n,x,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:31

!     .. scalar arguments ..

integer, intent(in out)                  :: n
real(kind=16), intent(in out)         :: x(*)
integer, intent(in)                      :: incx

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  qnrm2 returns the euclidean norm of a vector via the function
!  name, so that

!     qnrm2 := sqrt( x'*x )

!  further details
!  ===============

!  -- this version written on 25-october-1982.
!     modified on 14-october-1993 to inline the call to dlassq.
!     sven hammarling, nag ltd.

!  =====================================================================

!     .. parameters ..

real(kind=16), parameter :: one=1.0d+0
real(kind=16), parameter :: zero=0.0d+0
!     ..
!     .. local scalars ..
real(kind=16) :: absxi,norm,scale,ssq
integer :: ix
!     ..
!     .. intrinsic functions ..
intrinsic abs,sqrt
!     ..
if (n < 1 .or. incx < 1) then
  norm = zero
else if (n == 1) then
  norm = abs(x(1))
else
  scale = zero
  ssq = one
!        the following loop is equivalent to this call to the lapack
!        auxiliary routine:
!        call dlassq( n, x, incx, scale, ssq )
  
  do  ix = 1,1 + (n-1)*incx,incx
    if (x(ix) /= zero) then
      absxi = abs(x(ix))
      if (scale < absxi) then
        ssq = one + ssq* (scale/absxi)**2
        scale = absxi
      else
        ssq = ssq + (absxi/scale)**2
      end if
    end if
  end do
  norm = scale*sqrt(ssq)
end if

qnrm2 = norm
return

!     end of qnrm2.

end function qnrm2
