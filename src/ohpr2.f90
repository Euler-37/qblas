subroutine ohpr2(uplo,n,alpha,x,incx,y,incy,ap)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:36

!     .. scalar arguments ..

character(len=1),intent(in)        :: uplo
integer, intent(in)                      :: n
complex(kind=16), intent(in)               :: alpha
complex(kind=16), intent(in)               :: x(*)
integer, intent(in)                      :: incx
complex(kind=16), intent(in)               :: y(*)
integer, intent(in)                      :: incy
complex(kind=16), intent(out)              :: ap(*)



!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  ohpr2  performs the hermitian rank 2 operation

!     a := alpha*x*y**h + conjg( alpha )*y*x**h + a,

!  where alpha is a scalar, x and y are n element vectors and a is an
!  n by n hermitian matrix, supplied in packed form.

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

!  x      - complex(kind=16)       array of dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ).
!           before entry, the incremented array x must contain the n
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

!  ap     - complex(kind=16)       array of dimension at least
!           ( ( n*( n + 1 ) )/2 ).
!           before entry with  uplo = 'u' or 'u', the array ap must
!           contain the upper triangular part of the hermitian matrix
!           packed sequentially, column by column, so that ap( 1 )
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
!           and a( 2, 2 ) respectively, and so on. on exit, the array
!           ap is overwritten by the upper triangular part of the
!           updated matrix.
!           before entry with uplo = 'l' or 'l', the array ap must
!           contain the lower triangular part of the hermitian matrix
!           packed sequentially, column by column, so that ap( 1 )
!           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
!           and a( 3, 1 ) respectively, and so on. on exit, the array
!           ap is overwritten by the lower triangular part of the
!           updated matrix.
!           note that the imaginary parts of the diagonal elements need
!           not be set, they are assumed to be zero, and on exit they
!           are set to zero.

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
  info = 5
else if (incy == 0) then
  info = 7
end if
if (info /= 0) then
  call xerbla('ohpr2 ',info)
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

!     start the operations. in this version the elements of the array ap
!     are accessed sequentially with one pass through ap.

kk = 1
if (lsame(uplo,'u')) then
  
!        form  a  when upper triangle is stored in ap.
  
  if ((incx == 1) .and. (incy == 1)) then
    do  j = 1,n
      if ((x(j) /= zero) .or. (y(j) /= zero)) then
        temp1 = alpha*conjg(y(j))
        temp2 = conjg(alpha*x(j))
        k = kk
        do  i = 1,j - 1
          ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
          k = k + 1
        end do
        ap(kk+j-1) = real(ap(kk+j-1),kind=16) + real(x(j)*temp1+y(j)*temp2,kind=16)
      else
        ap(kk+j-1) = real(ap(kk+j-1),kind=16)
      end if
      kk = kk + j
    end do
  else
    do  j = 1,n
      if ((x(jx) /= zero) .or. (y(jy) /= zero)) then
        temp1 = alpha*conjg(y(jy))
        temp2 = conjg(alpha*x(jx))
        ix = kx
        iy = ky
        do  k = kk,kk + j - 2
          ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
          ix = ix + incx
          iy = iy + incy
        end do
        ap(kk+j-1) = real(ap(kk+j-1),kind=16) + real(x(jx)*temp1+y(jy)*temp2,kind=16)
      else
        ap(kk+j-1) = real(ap(kk+j-1),kind=16)
      end if
      jx = jx + incx
      jy = jy + incy
      kk = kk + j
    end do
  end if
else
  
!        form  a  when lower triangle is stored in ap.
  
  if ((incx == 1) .and. (incy == 1)) then
    do  j = 1,n
      if ((x(j) /= zero) .or. (y(j) /= zero)) then
        temp1 = alpha*conjg(y(j))
        temp2 = conjg(alpha*x(j))
        ap(kk) = real(ap(kk),kind=16) + real(x(j)*temp1+y(j)*temp2,kind=16)
        k = kk + 1
        do  i = j + 1,n
          ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
          k = k + 1
        end do
      else
        ap(kk) = real(ap(kk),kind=16)
      end if
      kk = kk + n - j + 1
    end do
  else
    do  j = 1,n
      if ((x(jx) /= zero) .or. (y(jy) /= zero)) then
        temp1 = alpha*conjg(y(jy))
        temp2 = conjg(alpha*x(jx))
        ap(kk) = real(ap(kk),kind=16) + real(x(jx)*temp1+y(jy)*temp2,kind=16)
        ix = jx
        iy = jy
        do  k = kk + 1,kk + n - j
          ix = ix + incx
          iy = iy + incy
          ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
        end do
      else
        ap(kk) = real(ap(kk),kind=16)
      end if
      jx = jx + incx
      jy = jy + incy
      kk = kk + n - j + 1
    end do
  end if
end if

return

!     end of ohpr2 .

end subroutine ohpr2
