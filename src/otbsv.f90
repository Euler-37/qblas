subroutine otbsv(uplo,trans,diag,n,k,a,lda,x,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:37

!     .. scalar arguments ..

character(len=1),intent(in)        :: uplo
character(len=1),intent(in)        :: trans
character(len=1),intent(in)        :: diag
integer, intent(in)                      :: n
integer, intent(in)                      :: k
complex(kind=16), intent(in)               :: a(lda,*)
integer, intent(in out)                  :: lda
complex(kind=16), intent(out)              :: x(*)
integer, intent(in)                      :: incx


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  otbsv  solves one of the systems of equations

!     a*x = b,   or   a**t*x = b,   or   a**h*x = b,

!  where b and x are n element vectors and a is an n by n unit, or
!  non-unit, upper or lower triangular band matrix, with ( k + 1 )
!  diagonals.

!  no test for singularity or near-singularity is included in this
!  routine. such tests must be performed before calling this routine.

!  arguments
!  ==========

!  uplo   - character*1.
!           on entry, uplo specifies whether the matrix is an upper or
!           lower triangular matrix as follows:

!              uplo = 'u' or 'u'   a is an upper triangular matrix.

!              uplo = 'l' or 'l'   a is a lower triangular matrix.

!           unchanged on exit.

!  trans  - character*1.
!           on entry, trans specifies the equations to be solved as
!           follows:

!              trans = 'n' or 'n'   a*x = b.

!              trans = 't' or 't'   a**t*x = b.

!              trans = 'c' or 'c'   a**h*x = b.

!           unchanged on exit.

!  diag   - character*1.
!           on entry, diag specifies whether or not a is unit
!           triangular as follows:

!              diag = 'u' or 'u'   a is assumed to be unit triangular.

!              diag = 'n' or 'n'   a is not assumed to be unit
!                                  triangular.

!           unchanged on exit.

!  n      - integer.
!           on entry, n specifies the order of the matrix a.
!           n must be at least zero.
!           unchanged on exit.

!  k      - integer.
!           on entry with uplo = 'u' or 'u', k specifies the number of
!           super-diagonals of the matrix a.
!           on entry with uplo = 'l' or 'l', k specifies the number of
!           sub-diagonals of the matrix a.
!           k must satisfy  0 .le. k.
!           unchanged on exit.

!  a      - complex(kind=16)       array of dimension ( lda, n ).
!           before entry with uplo = 'u' or 'u', the leading ( k + 1 )
!           by n part of the array a must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. the top left k by k triangle
!           of the array a is not referenced.
!           the following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:

!                 do 20, j = 1, n
!                    m = k + 1 - j
!                    do 10, i = max( 1, j - k ), j
!                       a( m + i, j ) = matrix( i, j )
!              10    continue
!              20 continue

!           before entry with uplo = 'l' or 'l', the leading ( k + 1 )
!           by n part of the array a must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. the bottom right k by k triangle of the
!           array a is not referenced.
!           the following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:

!                 do 20, j = 1, n
!                    m = 1 - j
!                    do 10, i = j, min( n, j + k )
!                       a( m + i, j ) = matrix( i, j )
!              10    continue
!              20 continue

!           note that when diag = 'u' or 'u' the elements of the array a
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           unchanged on exit.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least
!           ( k + 1 ).
!           unchanged on exit.

!  x      - complex(kind=16)       array of dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ).
!           before entry, the incremented array x must contain the n
!           element right-hand side vector b. on exit, x is overwritten
!           with the solution vector x.

!  incx   - integer.
!           on entry, incx specifies the increment for the elements of
!           x. incx must not be zero.
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
integer :: i,info,ix,j,jx,kplus1,kx,l
logical :: noconj,nounit
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
else if (.not.lsame(trans,'n') .and. .not.lsame(trans,'t') .and.  &
      .not.lsame(trans,'c')) then
  info = 2
else if (.not.lsame(diag,'u') .and. .not.lsame(diag,'n')) then
  info = 3
else if (n < 0) then
  info = 4
else if (k < 0) then
  info = 5
else if (lda < (k+1)) then
  info = 7
else if (incx == 0) then
  info = 9
end if
if (info /= 0) then
  call xerbla('otbsv ',info)
  return
end if

!     quick return if possible.

if (n == 0) return

noconj = lsame(trans,'t')
nounit = lsame(diag,'n')

!     set up the start point in x if the increment is not unity. this
!     will be  ( n - 1 )*incx  too small for descending loops.

if (incx <= 0) then
  kx = 1 - (n-1)*incx
else if (incx /= 1) then
  kx = 1
end if

!     start the operations. in this version the elements of a are
!     accessed by sequentially with one pass through a.

if (lsame(trans,'n')) then
  
!        form  x := inv( a )*x.
  
  if (lsame(uplo,'u')) then
    kplus1 = k + 1
    if (incx == 1) then
      do  j = n,1,-1
        if (x(j) /= zero) then
          l = kplus1 - j
          if (nounit) x(j) = x(j)/a(kplus1,j)
          temp = x(j)
          do  i = j - 1,max(1,j-k),-1
            x(i) = x(i) - temp*a(l+i,j)
          end do
        end if
      end do
    else
      kx = kx + (n-1)*incx
      jx = kx
      do  j = n,1,-1
        kx = kx - incx
        if (x(jx) /= zero) then
          ix = kx
          l = kplus1 - j
          if (nounit) x(jx) = x(jx)/a(kplus1,j)
          temp = x(jx)
          do  i = j - 1,max(1,j-k),-1
            x(ix) = x(ix) - temp*a(l+i,j)
            ix = ix - incx
          end do
        end if
        jx = jx - incx
      end do
    end if
  else
    if (incx == 1) then
      do  j = 1,n
        if (x(j) /= zero) then
          l = 1 - j
          if (nounit) x(j) = x(j)/a(1,j)
          temp = x(j)
          do  i = j + 1,min(n,j+k)
            x(i) = x(i) - temp*a(l+i,j)
          end do
        end if
      end do
    else
      jx = kx
      do  j = 1,n
        kx = kx + incx
        if (x(jx) /= zero) then
          ix = kx
          l = 1 - j
          if (nounit) x(jx) = x(jx)/a(1,j)
          temp = x(jx)
          do  i = j + 1,min(n,j+k)
            x(ix) = x(ix) - temp*a(l+i,j)
            ix = ix + incx
          end do
        end if
        jx = jx + incx
      end do
    end if
  end if
else
  
!        form  x := inv( a**t )*x  or  x := inv( a**h )*x.
  
  if (lsame(uplo,'u')) then
    kplus1 = k + 1
    if (incx == 1) then
      do  j = 1,n
        temp = x(j)
        l = kplus1 - j
        if (noconj) then
          do  i = max(1,j-k),j - 1
            temp = temp - a(l+i,j)*x(i)
          end do
          if (nounit) temp = temp/a(kplus1,j)
        else
          do  i = max(1,j-k),j - 1
            temp = temp - conjg(a(l+i,j))*x(i)
          end do
          if (nounit) temp = temp/conjg(a(kplus1,j))
        end if
        x(j) = temp
      end do
    else
      jx = kx
      do  j = 1,n
        temp = x(jx)
        ix = kx
        l = kplus1 - j
        if (noconj) then
          do  i = max(1,j-k),j - 1
            temp = temp - a(l+i,j)*x(ix)
            ix = ix + incx
          end do
          if (nounit) temp = temp/a(kplus1,j)
        else
          do  i = max(1,j-k),j - 1
            temp = temp - conjg(a(l+i,j))*x(ix)
            ix = ix + incx
          end do
          if (nounit) temp = temp/conjg(a(kplus1,j))
        end if
        x(jx) = temp
        jx = jx + incx
        if (j > k) kx = kx + incx
      end do
    end if
  else
    if (incx == 1) then
      do  j = n,1,-1
        temp = x(j)
        l = 1 - j
        if (noconj) then
          do  i = min(n,j+k),j + 1,-1
            temp = temp - a(l+i,j)*x(i)
          end do
          if (nounit) temp = temp/a(1,j)
        else
          do  i = min(n,j+k),j + 1,-1
            temp = temp - conjg(a(l+i,j))*x(i)
          end do
          if (nounit) temp = temp/conjg(a(1,j))
        end if
        x(j) = temp
      end do
    else
      kx = kx + (n-1)*incx
      jx = kx
      do  j = n,1,-1
        temp = x(jx)
        ix = kx
        l = 1 - j
        if (noconj) then
          do  i = min(n,j+k),j + 1,-1
            temp = temp - a(l+i,j)*x(ix)
            ix = ix - incx
          end do
          if (nounit) temp = temp/a(1,j)
        else
          do  i = min(n,j+k),j + 1,-1
            temp = temp - conjg(a(l+i,j))*x(ix)
            ix = ix - incx
          end do
          if (nounit) temp = temp/conjg(a(1,j))
        end if
        x(jx) = temp
        jx = jx - incx
        if ((n-j) >= k) kx = kx - incx
      end do
    end if
  end if
end if

return

!     end of otbsv .

end subroutine otbsv
