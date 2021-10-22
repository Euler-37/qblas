subroutine qsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:32

!     .. scalar arguments ..

character(len=*),intent(in)        :: uplo
integer, intent(in)                      :: n
real(kind=16), intent(in)             :: alpha
real(kind=16), intent(in)             :: x(*)
integer, intent(in)                      :: incx
real(kind=16), intent(in)             :: y(*)
integer, intent(in)                      :: incy
real(kind=16), intent(out)            :: a(lda,*)
integer, intent(in out)                  :: lda



!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  qsyr2  performs the symmetric rank 2 operation

!     a := alpha*x*y**t + alpha*y*x**t + a,

!  where alpha is a scalar, x and y are n element vectors and a is an n
!  by n symmetric matrix.

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

!  x      - real(kind=16) array of dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ).
!           before entry, the incremented array x must contain the n
!           element vector x.
!           unchanged on exit.

!  incx   - integer.
!           on entry, incx specifies the increment for the elements of
!           x. incx must not be zero.
!           unchanged on exit.

!  y      - real(kind=16) array of dimension at least
!           ( 1 + ( n - 1 )*abs( incy ) ).
!           before entry, the incremented array y must contain the n
!           element vector y.
!           unchanged on exit.

!  incy   - integer.
!           on entry, incy specifies the increment for the elements of
!           y. incy must not be zero.
!           unchanged on exit.

!  a      - real(kind=16) array of dimension ( lda, n ).
!           before entry with  uplo = 'u' or 'u', the leading n by n
!           upper triangular part of the array a must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of a is not referenced. on exit, the
!           upper triangular part of the array a is overwritten by the
!           upper triangular part of the updated matrix.
!           before entry with uplo = 'l' or 'l', the leading n by n
!           lower triangular part of the array a must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of a is not referenced. on exit, the
!           lower triangular part of the array a is overwritten by the
!           lower triangular part of the updated matrix.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least
!           max( 1, n ).
!           unchanged on exit.

!  further details
!  ===============

!  level 2 blas routine.

!  -- written on 22-october-1986.
!     jack dongarra, argonne national lab.
!     jeremy du croz, nag central office.
!     sven hammarling, nag central office.
!     richard hanson, sandia national labs.

!  =====================================================================

!     .. parameters ..

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
else if (incx == 0) then
  info = 5
else if (incy == 0) then
  info = 7
else if (lda < max(1,n)) then
  info = 9
end if
if (info /= 0) then
  call xerbla('qsyr2 ',info)
  return
end if

!     quick return if possible.

if ((n == 0) .or. (alpha == zero)) return

!     set up the start points in x and y if the increments are not both
!     unity.

if ((incx /= 1) .or. (incy /= 1)) then
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
  jx = kx
  jy = ky
end if

!     start the operations. in this version the elements of a are
!     accessed sequentially with one pass through the triangular part
!     of a.

if (lsame(uplo,'u')) then
  
!        form  a  when a is stored in the upper triangle.
  
  if ((incx == 1) .and. (incy == 1)) then
    do  j = 1,n
      if ((x(j) /= zero) .or. (y(j) /= zero)) then
        temp1 = alpha*y(j)
        temp2 = alpha*x(j)
        do  i = 1,j
          a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
        end do
      end if
    end do
  else
    do  j = 1,n
      if ((x(jx) /= zero) .or. (y(jy) /= zero)) then
        temp1 = alpha*y(jy)
        temp2 = alpha*x(jx)
        ix = kx
        iy = ky
        do  i = 1,j
          a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
          ix = ix + incx
          iy = iy + incy
        end do
      end if
      jx = jx + incx
      jy = jy + incy
    end do
  end if
else
  
!        form  a  when a is stored in the lower triangle.
  
  if ((incx == 1) .and. (incy == 1)) then
    do  j = 1,n
      if ((x(j) /= zero) .or. (y(j) /= zero)) then
        temp1 = alpha*y(j)
        temp2 = alpha*x(j)
        do  i = j,n
          a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
        end do
      end if
    end do
  else
    do  j = 1,n
      if ((x(jx) /= zero) .or. (y(jy) /= zero)) then
        temp1 = alpha*y(jy)
        temp2 = alpha*x(jx)
        ix = jx
        iy = jy
        do  i = j,n
          a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
          ix = ix + incx
          iy = iy + incy
        end do
      end if
      jx = jx + incx
      jy = jy + incy
    end do
  end if
end if

return

!     end of qsyr2 .

end subroutine qsyr2
