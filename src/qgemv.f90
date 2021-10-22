subroutine qgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:31

!     .. scalar arguments ..

character(len=*),intent(in)        :: trans
integer, intent(in)                      :: m
integer, intent(in)                      :: n
real(kind=16), intent(in)             :: alpha
real(kind=16), intent(in)             :: a(lda,*)
integer, intent(in out)                  :: lda
real(kind=16), intent(in)             :: x(*)
integer, intent(in)                      :: incx
real(kind=16), intent(in)             :: beta
real(kind=16), intent(out)            :: y(*)
integer, intent(in)                      :: incy



!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  qgemv  performs one of the matrix-vector operations

!     y := alpha*a*x + beta*y,   or   y := alpha*a**t*x + beta*y,

!  where alpha and beta are scalars, x and y are vectors and a is an
!  m by n matrix.

!  arguments
!  ==========

!  trans  - character*1.
!           on entry, trans specifies the operation to be performed as
!           follows:

!              trans = 'n' or 'n'   y := alpha*a*x + beta*y.

!              trans = 't' or 't'   y := alpha*a**t*x + beta*y.

!              trans = 'c' or 'c'   y := alpha*a**t*x + beta*y.

!           unchanged on exit.

!  m      - integer.
!           on entry, m specifies the number of rows of the matrix a.
!           m must be at least zero.
!           unchanged on exit.

!  n      - integer.
!           on entry, n specifies the number of columns of the matrix a.
!           n must be at least zero.
!           unchanged on exit.

!  alpha  - real(kind=16).
!           on entry, alpha specifies the scalar alpha.
!           unchanged on exit.

!  a      - real(kind=16) array of dimension ( lda, n ).
!           before entry, the leading m by n part of the array a must
!           contain the matrix of coefficients.
!           unchanged on exit.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least
!           max( 1, m ).
!           unchanged on exit.

!  x      - real(kind=16) array of dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
!           before entry, the incremented array x must contain the
!           vector x.
!           unchanged on exit.

!  incx   - integer.
!           on entry, incx specifies the increment for the elements of
!           x. incx must not be zero.
!           unchanged on exit.

!  beta   - real(kind=16).
!           on entry, beta specifies the scalar beta. when beta is
!           supplied as zero then y need not be set on input.
!           unchanged on exit.

!  y      - real(kind=16) array of dimension at least
!           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
!           before entry with beta non-zero, the incremented array y
!           must contain the vector y. on exit, y is overwritten by the
!           updated vector y.

!  incy   - integer.
!           on entry, incy specifies the increment for the elements of
!           y. incy must not be zero.
!           unchanged on exit.

!  further details
!  ===============

!  level 2 blas routine.
!  the vector and matrix arguments are not referenced when n = 0, or m = 0

!  -- written on 22-october-1986.
!     jack dongarra, argonne national lab.
!     jeremy du croz, nag central office.
!     sven hammarling, nag central office.
!     richard hanson, sandia national labs.

!  =====================================================================

!     .. parameters ..

real(kind=16), parameter :: one=1.0d+0
real(kind=16), parameter :: zero=0.0d+0
!     ..
!     .. local scalars ..
real(kind=16) :: temp
integer :: i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny
!     ..
!     .. external functions ..
logical :: lsame
external lsame
!     ..
!     .. external subroutines ..
external xerbla
!     ..
!     .. intrinsic functions ..
intrinsic max
!     ..

!     test the input parameters.

info = 0
if (.not.lsame(trans,'n') .and. .not.lsame(trans,'t') .and.  &
      .not.lsame(trans,'c')) then
  info = 1
else if (m < 0) then
  info = 2
else if (n < 0) then
  info = 3
else if (lda < max(1,m)) then
  info = 6
else if (incx == 0) then
  info = 8
else if (incy == 0) then
  info = 11
end if
if (info /= 0) then
  call xerbla('qgemv ',info)
  return
end if

!     quick return if possible.

if ((m == 0) .or. (n == 0) .or. ((alpha == zero).and. (beta == one))) return

!     set  lenx  and  leny, the lengths of the vectors x and y, and set
!     up the start points in  x  and  y.

if (lsame(trans,'n')) then
  lenx = n
  leny = m
else
  lenx = m
  leny = n
end if
if (incx > 0) then
  kx = 1
else
  kx = 1 - (lenx-1)*incx
end if
if (incy > 0) then
  ky = 1
else
  ky = 1 - (leny-1)*incy
end if

!     start the operations. in this version the elements of a are
!     accessed sequentially with one pass through a.

!     first form  y := beta*y.

if (beta /= one) then
  if (incy == 1) then
    if (beta == zero) then
      do  i = 1,leny
        y(i) = zero
      end do
    else
      do  i = 1,leny
        y(i) = beta*y(i)
      end do
    end if
  else
    iy = ky
    if (beta == zero) then
      do  i = 1,leny
        y(iy) = zero
        iy = iy + incy
      end do
    else
      do  i = 1,leny
        y(iy) = beta*y(iy)
        iy = iy + incy
      end do
    end if
  end if
end if
if (alpha == zero) return
if (lsame(trans,'n')) then
  
!        form  y := alpha*a*x + y.
  
  jx = kx
  if (incy == 1) then
    do  j = 1,n
      if (x(jx) /= zero) then
        temp = alpha*x(jx)
        do  i = 1,m
          y(i) = y(i) + temp*a(i,j)
        end do
      end if
      jx = jx + incx
    end do
  else
    do  j = 1,n
      if (x(jx) /= zero) then
        temp = alpha*x(jx)
        iy = ky
        do  i = 1,m
          y(iy) = y(iy) + temp*a(i,j)
          iy = iy + incy
        end do
      end if
      jx = jx + incx
    end do
  end if
else
  
!        form  y := alpha*a**t*x + y.
  
  jy = ky
  if (incx == 1) then
    do  j = 1,n
      temp = zero
      do  i = 1,m
        temp = temp + a(i,j)*x(i)
      end do
      y(jy) = y(jy) + alpha*temp
      jy = jy + incy
    end do
  else
    do  j = 1,n
      temp = zero
      ix = kx
      do  i = 1,m
        temp = temp + a(i,j)*x(ix)
        ix = ix + incx
      end do
      y(jy) = y(jy) + alpha*temp
      jy = jy + incy
    end do
  end if
end if

return

!     end of qgemv .

end subroutine qgemv
