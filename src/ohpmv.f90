subroutine ohpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:36

!     .. scalar arguments ..

character(len=*),intent(in)        :: uplo
integer, intent(in)                      :: n
complex(kind=16), intent(in)               :: alpha
complex(kind=16), intent(in)               :: ap(*)
complex(kind=16), intent(in)               :: x(*)
integer, intent(in)                      :: incx
complex(kind=16), intent(in)               :: beta
complex(kind=16), intent(out)              :: y(*)
integer, intent(in)                      :: incy



!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  ohpmv  performs the matrix-vector operation

!     y := alpha*a*x + beta*y,

!  where alpha and beta are scalars, x and y are n element vectors and
!  a is an n by n hermitian matrix, supplied in packed form.

!  arguments
!  ==========

!  uplo   - character*1.
!           on entry, uplo specifies whether the upper or lower
!           triangular part of the matrix a is supplied in the packed
!           array ap as follows:

!              uplo = 'u' or 'u'   the upper triangular part of a is
!                                  supplied in ap.

!              uplo = 'l' or 'l'   the lower triangular part of a is
!                                  supplied in ap.

!           unchanged on exit.

!  n      - integer.
!           on entry, n specifies the order of the matrix a.
!           n must be at least zero.
!           unchanged on exit.

!  alpha  - complex(kind=16)      .
!           on entry, alpha specifies the scalar alpha.
!           unchanged on exit.

!  ap     - complex(kind=16)       array of dimension at least
!           ( ( n*( n + 1 ) )/2 ).
!           before entry with uplo = 'u' or 'u', the array ap must
!           contain the upper triangular part of the hermitian matrix
!           packed sequentially, column by column, so that ap( 1 )
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
!           and a( 2, 2 ) respectively, and so on.
!           before entry with uplo = 'l' or 'l', the array ap must
!           contain the lower triangular part of the hermitian matrix
!           packed sequentially, column by column, so that ap( 1 )
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
!           and a( 3, 1 ) respectively, and so on.
!           note that the imaginary parts of the diagonal elements need
!           not be set and are assumed to be zero.
!           unchanged on exit.

!  x      - complex(kind=16)       array of dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ).
!           before entry, the incremented array x must contain the n
!           element vector x.
!           unchanged on exit.

!  incx   - integer.
!           on entry, incx specifies the increment for the elements of
!           x. incx must not be zero.
!           unchanged on exit.

!  beta   - complex(kind=16)      .
!           on entry, beta specifies the scalar beta. when beta is
!           supplied as zero then y need not be set on input.
!           unchanged on exit.

!  y      - complex(kind=16)       array of dimension at least
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

complex(kind=16), parameter :: one= (1.0d+0,0.0d+0)

complex(kind=16), parameter :: zero= (0.0d+0,0.0d+0)
!     ..
!     .. local scalars ..
complex(kind=16) :: temp1,temp2
integer :: i,info,ix,iy,j,jx,jy,k,kk,kx,ky
!     ..
!     .. external functions ..
logical :: lsame
external lsame
!     ..
!     .. external subroutines ..
external xerbla
!     ..
!     .. intrinsic functions ..
intrinsic dble,conjg
!     ..

!     test the input parameters.

info = 0
if (.not.lsame(uplo,'u') .and. .not.lsame(uplo,'l')) then
  info = 1
else if (n < 0) then
  info = 2
else if (incx == 0) then
  info = 6
else if (incy == 0) then
  info = 9
end if
if (info /= 0) then
  call xerbla('ohpmv ',info)
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

!     start the operations. in this version the elements of the array ap
!     are accessed sequentially with one pass through ap.

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
kk = 1
if (lsame(uplo,'u')) then
  
!        form  y  when ap contains the upper triangle.
  
  if ((incx == 1) .and. (incy == 1)) then
    do  j = 1,n
      temp1 = alpha*x(j)
      temp2 = zero
      k = kk
      do  i = 1,j - 1
        y(i) = y(i) + temp1*ap(k)
        temp2 = temp2 + conjg(ap(k))*x(i)
        k = k + 1
      end do
      y(j) = y(j) + temp1*real(ap(kk+j-1),kind=16) + alpha*temp2
      kk = kk + j
    end do
  else
    jx = kx
    jy = ky
    do  j = 1,n
      temp1 = alpha*x(jx)
      temp2 = zero
      ix = kx
      iy = ky
      do  k = kk,kk + j - 2
        y(iy) = y(iy) + temp1*ap(k)
        temp2 = temp2 + conjg(ap(k))*x(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      y(jy) = y(jy) + temp1*real(ap(kk+j-1),kind=16) + alpha*temp2
      jx = jx + incx
      jy = jy + incy
      kk = kk + j
    end do
  end if
else
  
!        form  y  when ap contains the lower triangle.
  
  if ((incx == 1) .and. (incy == 1)) then
    do  j = 1,n
      temp1 = alpha*x(j)
      temp2 = zero
      y(j) = y(j) + temp1*real(ap(kk),kind=16)
      k = kk + 1
      do  i = j + 1,n
        y(i) = y(i) + temp1*ap(k)
        temp2 = temp2 + conjg(ap(k))*x(i)
        k = k + 1
      end do
      y(j) = y(j) + alpha*temp2
      kk = kk + (n-j+1)
    end do
  else
    jx = kx
    jy = ky
    do  j = 1,n
      temp1 = alpha*x(jx)
      temp2 = zero
      y(jy) = y(jy) + temp1*real(ap(kk),kind=16)
      ix = jx
      iy = jy
      do  k = kk + 1,kk + n - j
        ix = ix + incx
        iy = iy + incy
        y(iy) = y(iy) + temp1*ap(k)
        temp2 = temp2 + conjg(ap(k))*x(ix)
      end do
      y(jy) = y(jy) + alpha*temp2
      jx = jx + incx
      jy = jy + incy
      kk = kk + (n-j+1)
    end do
  end if
end if

return

!     end of ohpmv .

end subroutine ohpmv
