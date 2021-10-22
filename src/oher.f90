subroutine oher(uplo,n,alpha,x,incx,a,lda)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:36

!     .. scalar arguments ..

character(len=*),intent(in)        :: uplo
integer, intent(in)                      :: n
real(kind=16), intent(in)             :: alpha
complex(kind=16), intent(in)               :: x(*)
integer, intent(in)                      :: incx
complex(kind=16), intent(out)              :: a(lda,*)
integer, intent(in out)                  :: lda



!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  oher   performs the hermitian rank 1 operation

!     a := alpha*x*x**h + a,

!  where alpha is a real scalar, x is an n element vector and a is an
!  n by n hermitian matrix.

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

!  x      - complex(kind=16)       array of dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ).
!           before entry, the incremented array x must contain the n
!           element vector x.
!           unchanged on exit.

!  incx   - integer.
!           on entry, incx specifies the increment for the elements of
!           x. incx must not be zero.
!           unchanged on exit.

!  a      - complex(kind=16)       array of dimension ( lda, n ).
!           before entry with  uplo = 'u' or 'u', the leading n by n
!           upper triangular part of the array a must contain the upper
!           triangular part of the hermitian matrix and the strictly
!           lower triangular part of a is not referenced. on exit, the
!           upper triangular part of the array a is overwritten by the
!           upper triangular part of the updated matrix.
!           before entry with uplo = 'l' or 'l', the leading n by n
!           lower triangular part of the array a must contain the lower
!           triangular part of the hermitian matrix and the strictly
!           upper triangular part of a is not referenced. on exit, the
!           lower triangular part of the array a is overwritten by the
!           lower triangular part of the updated matrix.
!           note that the imaginary parts of the diagonal elements need
!           not be set, they are assumed to be zero, and on exit they
!           are set to zero.

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

complex(kind=16), parameter :: zero= (0.0d+0,0.0d+0)
!     ..
!     .. local scalars ..
complex(kind=16) :: temp
integer :: i,info,ix,j,jx,kx
!     ..
!     .. external functions ..
logical :: lsame
external lsame
!     ..
!     .. external subroutines ..
external xerbla
!     ..
!     .. intrinsic functions ..
intrinsic dble,conjg,max
!     ..

!     test the input parameters.

info = 0
if (.not.lsame(uplo,'u') .and. .not.lsame(uplo,'l')) then
  info = 1
else if (n < 0) then
  info = 2
else if (incx == 0) then
  info = 5
else if (lda < max(1,n)) then
  info = 7
end if
if (info /= 0) then
  call xerbla('oher  ',info)
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

!     start the operations. in this version the elements of a are
!     accessed sequentially with one pass through the triangular part
!     of a.

if (lsame(uplo,'u')) then
  
!        form  a  when a is stored in upper triangle.
  
  if (incx == 1) then
    do  j = 1,n
      if (x(j) /= zero) then
        temp = alpha*conjg(x(j))
        do  i = 1,j - 1
          a(i,j) = a(i,j) + x(i)*temp
        end do
        a(j,j) = real(a(j,j),kind=16) + real(x(j)*temp,kind=16)
      else
        a(j,j) = real(a(j,j),kind=16)
      end if
    end do
  else
    jx = kx
    do  j = 1,n
      if (x(jx) /= zero) then
        temp = alpha*conjg(x(jx))
        ix = kx
        do  i = 1,j - 1
          a(i,j) = a(i,j) + x(ix)*temp
          ix = ix + incx
        end do
        a(j,j) = real(a(j,j),kind=16) + real(x(jx)*temp,kind=16)
      else
        a(j,j) = real(a(j,j),kind=16)
      end if
      jx = jx + incx
    end do
  end if
else
  
!        form  a  when a is stored in lower triangle.
  
  if (incx == 1) then
    do  j = 1,n
      if (x(j) /= zero) then
        temp = alpha*conjg(x(j))
        a(j,j) = real(a(j,j),kind=16) + real(temp*x(j),kind=16)
        do  i = j + 1,n
          a(i,j) = a(i,j) + x(i)*temp
        end do
      else
        a(j,j) = real(a(j,j),kind=16)
      end if
    end do
  else
    jx = kx
    do  j = 1,n
      if (x(jx) /= zero) then
        temp = alpha*conjg(x(jx))
        a(j,j) = real(a(j,j),kind=16) + real(temp*x(jx),kind=16)
        ix = jx
        do  i = j + 1,n
          ix = ix + incx
          a(i,j) = a(i,j) + x(ix)*temp
        end do
      else
        a(j,j) = real(a(j,j),kind=16)
      end if
      jx = jx + incx
    end do
  end if
end if

return

!     end of oher  .

end subroutine oher
