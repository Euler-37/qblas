subroutine qsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:32

!     .. scalar arguments ..

character(len=1),intent(in)        :: uplo
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

!  qsymv  performs the matrix-vector  operation

!     y := alpha*a*x + beta*y,

!  where alpha and beta are scalars, x and y are n element vectors and
!  a is an n by n symmetric matrix.

!  arguments
!  ==========

!  uplo   - character*1.
!           on entry, uplo specifies whether the upper or lower
!           triangular part of the array a is to be referenced as
!           follows:

!              uplo = 'u' or 'u'   only the upper triangular part of a
!                                  is to be referenced.

!              uplo = 'l' or 'l'   only the lower triangular part of a
!                                  is to be referenced.

!           unchanged on exit.

!  n      - integer.
!           on entry, n specifies the order of the matrix a.
!           n must be at least zero.
!           unchanged on exit.

!  alpha  - real(kind=16).
!           on entry, alpha specifies the scalar alpha.
!           unchanged on exit.

!  a      - real(kind=16) array of dimension ( lda, n ).
!           before entry with  uplo = 'u' or 'u', the leading n by n
!           upper triangular part of the array a must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of a is not referenced.
!           before entry with uplo = 'l' or 'l', the leading n by n
!           lower triangular part of the array a must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of a is not referenced.
!           unchanged on exit.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least
!           max( 1, n ).
!           unchanged on exit.

!  x      - real(kind=16) array of dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ).
!           before entry, the incremented array x must contain the n
!           element vector x.
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
!           ( 1 + ( n - 1 )*abs( incy ) ).
!           before entry, the incremented array y must contain the n
!           element vector y. on exit, y is overwritten by the updated
!           vector y.

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
real(kind=16) :: temp1,temp2
integer :: i,info,ix,iy,j,jx,jy,kx,ky
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
if (.not.lsame(uplo,'u') .and. .not.lsame(uplo,'l')) then
  info = 1
else if (n < 0) then
  info = 2
else if (lda < max(1,n)) then
  info = 5
else if (incx == 0) then
  info = 7
else if (incy == 0) then
  info = 10
end if
if (info /= 0) then
  call xerbla('qsymv ',info)
  return
end if

!     quick return if possible.

if ((n == 0) .or. ((alpha == zero).and. (beta == one))) return

!     set up the start points in  x  and  y.

if (incx > 0) then
  kx = 1
else
  kx = 1 - (n-1)*incx
end if
if (incy > 0) then
  ky = 1
else
  ky = 1 - (n-1)*incy
end if

!     start the operations. in this version the elements of a are
!     accessed sequentially with one pass through the triangular part
!     of a.

!     first form  y := beta*y.

if (beta /= one) then
  if (incy == 1) then
    if (beta == zero) then
      do  i = 1,n
        y(i) = zero
      end do
    else
      do  i = 1,n
        y(i) = beta*y(i)
      end do
    end if
  else
    iy = ky
    if (beta == zero) then
      do  i = 1,n
        y(iy) = zero
        iy = iy + incy
      end do
    else
      do  i = 1,n
        y(iy) = beta*y(iy)
        iy = iy + incy
      end do
    end if
  end if
end if
if (alpha == zero) return
if (lsame(uplo,'u')) then
  
!        form  y  when a is stored in upper triangle.
  
  if ((incx == 1) .and. (incy == 1)) then
    do  j = 1,n
      temp1 = alpha*x(j)
      temp2 = zero
      do  i = 1,j - 1
        y(i) = y(i) + temp1*a(i,j)
        temp2 = temp2 + a(i,j)*x(i)
      end do
      y(j) = y(j) + temp1*a(j,j) + alpha*temp2
    end do
  else
    jx = kx
    jy = ky
    do  j = 1,n
      temp1 = alpha*x(jx)
      temp2 = zero
      ix = kx
      iy = ky
      do  i = 1,j - 1
        y(iy) = y(iy) + temp1*a(i,j)
        temp2 = temp2 + a(i,j)*x(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
      jx = jx + incx
      jy = jy + incy
    end do
  end if
else
  
!        form  y  when a is stored in lower triangle.
  
  if ((incx == 1) .and. (incy == 1)) then
    do  j = 1,n
      temp1 = alpha*x(j)
      temp2 = zero
      y(j) = y(j) + temp1*a(j,j)
      do  i = j + 1,n
        y(i) = y(i) + temp1*a(i,j)
        temp2 = temp2 + a(i,j)*x(i)
      end do
      y(j) = y(j) + alpha*temp2
    end do
  else
    jx = kx
    jy = ky
    do  j = 1,n
      temp1 = alpha*x(jx)
      temp2 = zero
      y(jy) = y(jy) + temp1*a(j,j)
      ix = jx
      iy = jy
      do  i = j + 1,n
        ix = ix + incx
        iy = iy + incy
        y(iy) = y(iy) + temp1*a(i,j)
        temp2 = temp2 + a(i,j)*x(ix)
      end do
      y(jy) = y(jy) + alpha*temp2
      jx = jx + incx
      jy = jy + incy
    end do
  end if
end if

return

!     end of qsymv .

end subroutine qsymv
