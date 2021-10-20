subroutine qrotm(n,dx,incx,dy,incy,dparam)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:31

!     .. scalar arguments ..

integer, intent(in)                      :: n
real(kind=16), intent(in out)         :: dx(*)
integer, intent(in)                      :: incx
real(kind=16), intent(in out)         :: dy(*)
integer, intent(in)                      :: incy
real(kind=16), intent(in)             :: dparam(5)

!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!     apply the modified givens transformation, h, to the 2 by n matrix

!     (dx**t) , where **t indicates transpose. the elements of dx are in
!     (dy**t)

!     dx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
!     lx = (-incx)*n, and similarly for sy using ly and incy.
!     with dparam(1)=dflag, h has one of the following forms..

!     dflag=-1.d0     dflag=0.q0        dflag=1.d0     dflag=-2.d0

!       (dh11  dh12)    (1.d0  dh12)    (dh11  1.d0)    (1.d0  0.q0)
!     h=(          )    (          )    (          )    (          )
!       (dh21  dh22),   (dh21  1.d0),   (-1.d0 dh22),   (0.q0  1.d0).
!     see qrotmg for a description of data storage in dparam.

!  arguments
!  =========

!  n      (input) integer
!         number of elements in input vector(s)

!  dx     (input/output) real(kind=16) array, dimension n
!         double precision vector with n elements

!  incx   (input) integer
!         storage spacing between elements of dx

!  dy     (input/output) real(kind=16) array, dimension n
!         double precision vector with n elements

!  incy   (input) integer
!         storage spacing between elements of dy

!  dparam (input/output)  real(kind=16) array, dimension 5
!     dparam(1)=dflag
!     dparam(2)=dh11
!     dparam(3)=dh21
!     dparam(4)=dh12
!     dparam(5)=dh22

!  =====================================================================

!     .. local scalars ..
real(kind=16) :: dflag,dh11,dh12,dh21,dh22,w,z
integer :: i,kx,ky,nsteps
!     ..
!     .. data statements ..
! data zero,two/0.q0,2.d0/
real(kind=16),parameter::zero=0.0_16,two=2.0_16
!     ..

dflag = dparam(1)
if (n <= 0 .or. (dflag+two == zero)) return
if (incx == incy.and.incx > 0) then
  
  nsteps = n*incx
  if (dflag < zero) then
    dh11 = dparam(2)
    dh12 = dparam(4)
    dh21 = dparam(3)
    dh22 = dparam(5)
    do i = 1,nsteps,incx
      w = dx(i)
      z = dy(i)
      dx(i) = w*dh11 + z*dh12
      dy(i) = w*dh21 + z*dh22
    end do
  else if (dflag == zero) then
    dh12 = dparam(4)
    dh21 = dparam(3)
    do i = 1,nsteps,incx
      w = dx(i)
      z = dy(i)
      dx(i) = w + z*dh12
      dy(i) = w*dh21 + z
    end do
  else
    dh11 = dparam(2)
    dh22 = dparam(5)
    do i = 1,nsteps,incx
      w = dx(i)
      z = dy(i)
      dx(i) = w*dh11 + z
      dy(i) = -w + dh22*z
    end do
  end if
else
  kx = 1
  ky = 1
  if (incx < 0) kx = 1 + (1-n)*incx
  if (incy < 0) ky = 1 + (1-n)*incy
  
  if (dflag < zero) then
    dh11 = dparam(2)
    dh12 = dparam(4)
    dh21 = dparam(3)
    dh22 = dparam(5)
    do i = 1,n
      w = dx(kx)
      z = dy(ky)
      dx(kx) = w*dh11 + z*dh12
      dy(ky) = w*dh21 + z*dh22
      kx = kx + incx
      ky = ky + incy
    end do
  else if (dflag == zero) then
    dh12 = dparam(4)
    dh21 = dparam(3)
    do i = 1,n
      w = dx(kx)
      z = dy(ky)
      dx(kx) = w + z*dh12
      dy(ky) = w*dh21 + z
      kx = kx + incx
      ky = ky + incy
    end do
  else
    dh11 = dparam(2)
    dh22 = dparam(5)
    do i = 1,n
      w = dx(kx)
      z = dy(ky)
      dx(kx) = w*dh11 + z
      dy(ky) = -w + dh22*z
      kx = kx + incx
      ky = ky + incy
    end do
  end if
end if
return
end subroutine qrotm
