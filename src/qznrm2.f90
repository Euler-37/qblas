real(kind=16) function qznrm2(n,x,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:32

!     .. scalar arguments ..

integer, intent(in out)                  :: n
complex(kind=16), intent(in out)           :: x(*)
integer, intent(in)                      :: incx

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  qznrm2 returns the euclidean norm of a vector via the function
!  name, so that

!     qznrm2 := sqrt( x**h*x )

!  further details
!  ===============

!  -- this version written on 25-october-1982.
!     modified on 14-october-1993 to inline the call to zlassq.
!     sven hammarling, nag ltd.

!  =====================================================================

!     .. parameters ..

real(kind=16), parameter :: one=1.0d+0
real(kind=16), parameter :: zero=0.0d+0
!     ..
!     .. local scalars ..
real(kind=16) :: norm,scale,ssq,temp
integer :: ix
!     ..
!     .. intrinsic functions ..
intrinsic abs,dble,aimag,sqrt
!     ..
if (n < 1 .or. incx < 1) then
  norm = zero
else
  scale = zero
  ssq = one
!        the following loop is equivalent to this call to the lapack
!        auxiliary routine:
!        call zlassq( n, x, incx, scale, ssq )
  
  do  ix = 1,1 + (n-1)*incx,incx
    if (real(x(ix),kind=16) /= zero) then
      temp = abs(real(x(ix),kind=16))
      if (scale < temp) then
        ssq = one + ssq* (scale/temp)**2
        scale = temp
      else
        ssq = ssq + (temp/scale)**2
      end if
    end if
    if (aimag(x(ix)) /= zero) then
      temp = abs(aimag(x(ix)))
      if (scale < temp) then
        ssq = one + ssq* (scale/temp)**2
        scale = temp
      else
        ssq = ssq + (temp/scale)**2
      end if
    end if
  end do
  norm = scale*sqrt(ssq)
end if

qznrm2 = norm
return

!     end of qznrm2.

end function qznrm2
