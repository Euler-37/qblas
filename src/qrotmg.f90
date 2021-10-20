subroutine qrotmg(dd1,dd2,dx1,dy1,dparam)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:31

!     .. scalar arguments ..

real(kind=16), intent(out)            :: dd1
real(kind=16), intent(out)            :: dd2
real(kind=16), intent(out)            :: dx1
real(kind=16), intent(in)             :: dy1
real(kind=16), intent(out)            :: dparam(5)

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     construct the modified givens transformation matrix h which zeros
!     the second component of the 2-vector  (sqrt(dd1)*dx1,sqrt(dd2)*
!     dy2)**t.
!     with dparam(1)=dflag, h has one of the following forms..

!     dflag=-1.d0     dflag=0.q0        dflag=1.d0     dflag=-2.d0

!       (dh11  dh12)    (1.d0  dh12)    (dh11  1.d0)    (1.d0  0.q0)
!     h=(          )    (          )    (          )    (          )
!       (dh21  dh22),   (dh21  1.d0),   (-1.d0 dh22),   (0.q0  1.d0).
!     locations 2-4 of dparam contain dh11, dh21, dh12, and dh22
!     respectively. (values of 1.d0, -1.d0, or 0.q0 implied by the
!     value of dparam(1) are not stored in dparam.)

!     the values of gamsq and rgamsq set in the data statement may be
!     inexact.  this is ok as they are only used for testing the size
!     of dd1 and dd2.  all actual scaling of data is done using gam.


!  arguments
!  =========

!  dd1    (input/output) real(kind=16)

!  dd2    (input/output) real(kind=16)

!  dx1    (input/output) real(kind=16)

!  dy1    (input) real(kind=16)

!  dparam (input/output)  real(kind=16) array, dimension 5
!     dparam(1)=dflag
!     dparam(2)=dh11
!     dparam(3)=dh21
!     dparam(4)=dh12
!     dparam(5)=dh22

!  =====================================================================

!     .. local scalars ..
real(kind=16) :: dflag,dh11,dh12,dh21,dh22,dp1,dp2,dq1,dq2,dtemp,  &
    du
!     ..
!     .. intrinsic functions ..
intrinsic abs
!     ..
!     .. data statements ..

real(kind=16),parameter::zero=0.0_16,one=1.0_16,two=2.0_16
real(kind=16),parameter::gam=4096.0_16,gamsq=16777216.0_16,rgamsq=one/gamsq
!     ..

if (dd1 < zero) then
    !        go zero-h-d-and-dx1..
    dflag = -one
    dh11 = zero
    dh12 = zero
    dh21 = zero
    dh22 = zero

    dd1 = zero
    dd2 = zero
    dx1 = zero
else
    !        case-dd1-nonnegative
    dp2 = dd2*dy1
    if (dp2 == zero) then
        dflag = -two
        dparam(1) = dflag
        return
    end if
    !        regular-case..
    dp1 = dd1*dx1
    dq2 = dp2*dy1
    dq1 = dp1*dx1

    if (abs(dq1) > abs(dq2)) then
        dh21 = -dy1/dx1
        dh12 = dp2/dp1

        du = one - dh12*dh21

        if (du > zero) then
            dflag = zero
            dd1 = dd1/du
            dd2 = dd2/du
            dx1 = dx1*du
        end if
    else

        if (dq2 < zero) then
            !              go zero-h-d-and-dx1..
            dflag = -one
            dh11 = zero
            dh12 = zero
            dh21 = zero
            dh22 = zero

            dd1 = zero
            dd2 = zero
            dx1 = zero
        else
            dflag = one
            dh11 = dp1/dp2
            dh22 = dx1/dy1
            du = one + dh11*dh22
            dtemp = dd2/du
            dd2 = dd1/du
            dd1 = dtemp
            dx1 = dy1*du
        end if
    end if

    !     procedure..scale-check
    if (dd1 /= zero) then
        do while ((dd1 <= rgamsq) .or. (dd1 >= gamsq))
            if (dflag == zero) then
                dh11 = one
                dh22 = one
                dflag = -one
            else
                dh21 = -one
                dh12 = one
                dflag = -one
            end if
            if (dd1 <= rgamsq) then
                dd1 = dd1*gam**2
                dx1 = dx1/gam
                dh11 = dh11/gam
                dh12 = dh12/gam
            else
                dd1 = dd1/gam**2
                dx1 = dx1*gam
                dh11 = dh11*gam
                dh12 = dh12*gam
            end if
        end do
    end if

    if (dd2 /= zero) then
        do while ( (abs(dd2) <= rgamsq) .or. (abs(dd2) >= gamsq) )
            if (dflag == zero) then
                dh11 = one
                dh22 = one
                dflag = -one
            else
                dh21 = -one
                dh12 = one
                dflag = -one
            end if
            if (abs(dd2) <= rgamsq) then
                dd2 = dd2*gam**2
                dh21 = dh21/gam
                dh22 = dh22/gam
            else
                dd2 = dd2/gam**2
                dh21 = dh21*gam
                dh22 = dh22*gam
            end if
        end do
    end if

end if

if (dflag < zero) then
    dparam(2) = dh11
    dparam(3) = dh21
    dparam(4) = dh12
    dparam(5) = dh22
else if (dflag == zero) then
    dparam(3) = dh21
    dparam(4) = dh12
else
    dparam(2) = dh11
    dparam(5) = dh22
end if

dparam(1) = dflag
return
end subroutine qrotmg




