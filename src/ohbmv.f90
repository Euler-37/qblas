subroutine ohbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:36

!     .. scalar arguments ..

character(len=*),intent(in)        :: uplo
integer, intent(in)                      :: n
integer, intent(in)                      :: k
complex(kind=16), intent(in)               :: alpha
complex(kind=16), intent(in)               :: a(lda,*)
integer, intent(in out)                  :: lda
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

!  ohbmv  performs the matrix-vector  operation

!     y := alpha*a*x + beta*y,

!  where alpha and beta are scalars, x and y are n element vectors and
!  a is an n by n hermitian band matrix, with k super-diagonals.

!  arguments
!  ==========

!  uplo   - character*1.
!           on entry, uplo specifies whether the upper or lower
!           triangular part of the band matrix a is being supplied as
!           follows:

!              uplo = 'u' or 'u'   the upper triangular part of a is
!                                  being supplied.

!              uplo = 'l' or 'l'   the lower triangular part of a is
!                                  being supplied.

!           unchanged on exit.

!  n      - integer.
!           on entry, n specifies the order of the matrix a.
!           n must be at least zero.
!           unchanged on exit.

!  k      - integer.
!           on entry, k specifies the number of super-diagonals of the
!           matrix a. k must satisfy  0 .le. k.
!           unchanged on exit.

!  alpha  - complex(kind=16)      .
!           on entry, alpha specifies the scalar alpha.
!           unchanged on exit.

!  a      - complex(kind=16)       array of dimension ( lda, n ).
!           before entry with uplo = 'u' or 'u', the leading ( k + 1 )
!           by n part of the array a must contain the upper triangular
!           band part of the hermitian matrix, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. the top left k by k triangle
!           of the array a is not referenced.
!           the following program segment will transfer the upper
!           triangular part of a hermitian band matrix from conventional
!           full matrix storage to band storage:

!                 do 20, j = 1, n
!                    m = k + 1 - j
!                    do 10, i = max( 1, j - k ), j
!                       a( m + i, j ) = matrix( i, j )
!              10    continue
!              20 continue

!           before entry with uplo = 'l' or 'l', the leading ( k + 1 )
!           by n part of the array a must contain the lower triangular
!           band part of the hermitian matrix, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. the bottom right k by k triangle of the
!           array a is not referenced.
!           the following program segment will transfer the lower
!           triangular part of a hermitian band matrix from conventional
!           full matrix storage to band storage:

!                 do 20, j = 1, n
!                    m = 1 - j
!                    do 10, i = j, min( n, j + k )
!                       a( m + i, j ) = matrix( i, j )
!              10    continue
!              20 continue

!           note that the imaginary parts of the diagonal elements need
!           not be set and are assumed to be zero.
!           unchanged on exit.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least
!           ( k + 1 ).
!           unchanged on exit.

!  x      - complex(kind=16)       array of dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ).
!           before entry, the incremented array x must contain the
!           vector x.
!           unchanged on exit.

!  incx   - integer.
!           on entry, incx specifies the increment for the elements of
!           x. incx must not be zero.
!           unchanged on exit.

!  beta   - complex(kind=16)      .
!           on entry, beta specifies the scalar beta.
!           unchanged on exit.

!  y      - complex(kind=16)       array of dimension at least
!           ( 1 + ( n - 1 )*abs( incy ) ).
!           before entry, the incremented array y must contain the
!           vector y. on exit, y is overwritten by the updated vector y.

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
integer :: i,info,ix,iy,j,jx,jy,kplus1,kx,ky,l
!     ..
!     .. external functions ..
logical :: lsame
external lsame
!     ..
!     .. external subroutines ..
external xerbla
!     ..
!     .. intrinsic functions ..
intrinsic conjg,max,min
!     ..

!     test the input parameters.

info = 0
if (.not.lsame(uplo,'u') .and. .not.lsame(uplo,'l')) then
  info = 1
else if (n < 0) then
  info = 2
else if (k < 0) then
  info = 3
else if (lda < (k+1)) then
  info = 6
else if (incx == 0) then
  info = 8
else if (incy == 0) then
  info = 11
end if
if (info /= 0) then
  call xerbla('ohbmv ',info)
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

!     start the operations. in this version the elements of the array a
!     are accessed sequentially with one pass through a.

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
  
!        form  y  when upper triangle of a is stored.
  
  kplus1 = k + 1
  if ((incx == 1) .and. (incy == 1)) then
    do  j = 1,n
      temp1 = alpha*x(j)
      temp2 = zero
      l = kplus1 - j
      do  i = max(1,j-k),j - 1
        y(i) = y(i) + temp1*a(l+i,j)
        temp2 = temp2 + conjg(a(l+i,j))*x(i)
      end do
      y(j) = y(j) + temp1*real(a(kplus1,j),kind=16) + alpha*temp2
    end do
  else
    jx = kx
    jy = ky
    do  j = 1,n
      temp1 = alpha*x(jx)
      temp2 = zero
      ix = kx
      iy = ky
      l = kplus1 - j
      do  i = max(1,j-k),j - 1
        y(iy) = y(iy) + temp1*a(l+i,j)
        temp2 = temp2 + conjg(a(l+i,j))*x(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      y(jy) = y(jy) + temp1*real(a(kplus1,j),kind=16) + alpha*temp2
      jx = jx + incx
      jy = jy + incy
      if (j > k) then
        kx = kx + incx
        ky = ky + incy
      end if
    end do
  end if
else
  
!        form  y  when lower triangle of a is stored.
  
  if ((incx == 1) .and. (incy == 1)) then
    do  j = 1,n
      temp1 = alpha*x(j)
      temp2 = zero
      y(j) = y(j) + temp1*real(a(1,j),kind=16)
      l = 1 - j
      do  i = j + 1,min(n,j+k)
        y(i) = y(i) + temp1*a(l+i,j)
        temp2 = temp2 + conjg(a(l+i,j))*x(i)
      end do
      y(j) = y(j) + alpha*temp2
    end do
  else
    jx = kx
    jy = ky
    do  j = 1,n
      temp1 = alpha*x(jx)
      temp2 = zero
      y(jy) = y(jy) + temp1*real(a(1,j),kind=16)
      l = 1 - j
      ix = jx
      iy = jy
      do  i = j + 1,min(n,j+k)
        ix = ix + incx
        iy = iy + incy
        y(iy) = y(iy) + temp1*a(l+i,j)
        temp2 = temp2 + conjg(a(l+i,j))*x(ix)
      end do
      y(jy) = y(jy) + alpha*temp2
      jx = jx + incx
      jy = jy + incy
    end do
  end if
end if

return

!     end of ohbmv .

end subroutine ohbmv
