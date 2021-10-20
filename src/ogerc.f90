subroutine ogerc(m,n,alpha,x,incx,y,incy,a,lda)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:36

!     .. scalar arguments ..

integer, intent(in)                      :: m
integer, intent(in)                      :: n
complex(kind=16), intent(in)               :: alpha
complex(kind=16), intent(in)               :: x(*)
integer, intent(in)                      :: incx
complex(kind=16), intent(in out)           :: y(*)
integer, intent(in)                      :: incy
complex(kind=16), intent(out)              :: a(lda,*)
integer, intent(in out)                  :: lda


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  ogerc  performs the rank 1 operation

!     a := alpha*x*y**h + a,

!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and a is an m by n matrix.

!  arguments
!  ==========

!  m      - integer.
!           on entry, m specifies the number of rows of the matrix a.
!           m must be at least zero.
!           unchanged on exit.

!  n      - integer.
!           on entry, n specifies the number of columns of the matrix a.
!           n must be at least zero.
!           unchanged on exit.

!  alpha  - complex(kind=16)      .
!           on entry, alpha specifies the scalar alpha.
!           unchanged on exit.

!  x      - complex(kind=16)       array of dimension at least
!           ( 1 + ( m - 1 )*abs( incx ) ).
!           before entry, the incremented array x must contain the m
!           element vector x.
!           unchanged on exit.

!  incx   - integer.
!           on entry, incx specifies the increment for the elements of
!           x. incx must not be zero.
!           unchanged on exit.

!  y      - complex(kind=16)       array of dimension at least
!           ( 1 + ( n - 1 )*abs( incy ) ).
!           before entry, the incremented array y must contain the n
!           element vector y.
!           unchanged on exit.

!  incy   - integer.
!           on entry, incy specifies the increment for the elements of
!           y. incy must not be zero.
!           unchanged on exit.

!  a      - complex(kind=16)       array of dimension ( lda, n ).
!           before entry, the leading m by n part of the array a must
!           contain the matrix of coefficients. on exit, a is
!           overwritten by the updated matrix.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least
!           max( 1, m ).
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

complex(kind=16), parameter :: zero= (0.0d+0,0.0d+0)
!     ..
!     .. local scalars ..
complex(kind=16) :: temp
integer :: i,info,ix,j,jy,kx
!     ..
!     .. external subroutines ..
external xerbla
!     ..
!     .. intrinsic functions ..
intrinsic conjg,max
!     ..

!     test the input parameters.

info = 0
if (m < 0) then
  info = 1
else if (n < 0) then
  info = 2
else if (incx == 0) then
  info = 5
else if (incy == 0) then
  info = 7
else if (lda < max(1,m)) then
  info = 9
end if
if (info /= 0) then
  call xerbla('ogerc ',info)
  return
end if

!     quick return if possible.

if ((m == 0) .or. (n == 0) .or. (alpha == zero)) return

!     start the operations. in this version the elements of a are
!     accessed sequentially with one pass through a.

if (incy > 0) then
  jy = 1
else
  jy = 1 - (n-1)*incy
end if
if (incx == 1) then
  do  j = 1,n
    if (y(jy) /= zero) then
      temp = alpha*conjg(y(jy))
      do  i = 1,m
        a(i,j) = a(i,j) + x(i)*temp
      end do
    end if
    jy = jy + incy
  end do
else
  if (incx > 0) then
    kx = 1
  else
    kx = 1 - (m-1)*incx
  end if
  do  j = 1,n
    if (y(jy) /= zero) then
      temp = alpha*conjg(y(jy))
      ix = kx
      do  i = 1,m
        a(i,j) = a(i,j) + x(ix)*temp
        ix = ix + incx
      end do
    end if
    jy = jy + incy
  end do
end if

return

!     end of ogerc .

end subroutine ogerc
