subroutine ohpr(uplo,n,alpha,x,incx,ap)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:36

!     .. scalar arguments ..

character(len=1),intent(in)        :: uplo
integer, intent(in)                      :: n
real(kind=16), intent(in)             :: alpha
complex(kind=16), intent(in)               :: x(*)
integer, intent(in)                      :: incx
complex(kind=16), intent(out)              :: ap(*)



!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  ohpr    performs the hermitian rank 1 operation

!     a := alpha*x*x**h + a,

!  where alpha is a real scalar, x is an n element vector and a is an
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

!  alpha  - real(kind=16).
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
complex(kind=16) :: temp
integer :: i,info,ix,j,jx,k,kk,kx
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
end if
if (info /= 0) then
  call xerbla('ohpr  ',info)
  return
end if

!     quick return if possible.

if ((n == 0) .or. (alpha == real(zero,kind=16))) return

!     set the start point in x if the increment is not unity.

if (incx <= 0) then
  kx = 1 - (n-1)*incx
else if (incx /= 1) then
  kx = 1
end if

!     start the operations. in this version the elements of the array ap
!     are accessed sequentially with one pass through ap.

kk = 1
if (lsame(uplo,'u')) then
  
!        form  a  when upper triangle is stored in ap.
  
  if (incx == 1) then
    do  j = 1,n
      if (x(j) /= zero) then
        temp = alpha*conjg(x(j))
        k = kk
        do  i = 1,j - 1
          ap(k) = ap(k) + x(i)*temp
          k = k + 1
        end do
        ap(kk+j-1) = real(ap(kk+j-1),kind=16) + real(x(j)*temp,kind=16)
      else
        ap(kk+j-1) = real(ap(kk+j-1),kind=16)
      end if
      kk = kk + j
    end do
  else
    jx = kx
    do  j = 1,n
      if (x(jx) /= zero) then
        temp = alpha*conjg(x(jx))
        ix = kx
        do  k = kk,kk + j - 2
          ap(k) = ap(k) + x(ix)*temp
          ix = ix + incx
        end do
        ap(kk+j-1) = real(ap(kk+j-1),kind=16) + real(x(jx)*temp,kind=16)
      else
        ap(kk+j-1) = real(ap(kk+j-1),kind=16)
      end if
      jx = jx + incx
      kk = kk + j
    end do
  end if
else
  
!        form  a  when lower triangle is stored in ap.
  
  if (incx == 1) then
    do  j = 1,n
      if (x(j) /= zero) then
        temp = alpha*conjg(x(j))
        ap(kk) = real(ap(kk),kind=16) + real(temp*x(j),kind=16)
        k = kk + 1
        do  i = j + 1,n
          ap(k) = ap(k) + x(i)*temp
          k = k + 1
        end do
      else
        ap(kk) = real(ap(kk),kind=16)
      end if
      kk = kk + n - j + 1
    end do
  else
    jx = kx
    do  j = 1,n
      if (x(jx) /= zero) then
        temp = alpha*conjg(x(jx))
        ap(kk) = real(ap(kk),kind=16) + real(temp*x(jx),kind=16)
        ix = jx
        do  k = kk + 1,kk + n - j
          ix = ix + incx
          ap(k) = ap(k) + x(ix)*temp
        end do
      else
        ap(kk) = real(ap(kk),kind=16)
      end if
      jx = jx + incx
      kk = kk + n - j + 1
    end do
  end if
end if

return

!     end of ohpr  .

end subroutine ohpr
