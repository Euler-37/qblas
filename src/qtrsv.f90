subroutine qtrsv(uplo,trans,diag,n,a,lda,x,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:32

!     .. scalar arguments ..

character(len=1),intent(in)        :: uplo
character(len=1),intent(in)        :: trans
character(len=1),intent(in)        :: diag
integer, intent(in)                      :: n
real(kind=16), intent(in)             :: a(lda,*)
integer, intent(in out)                  :: lda
real(kind=16), intent(out)            :: x(*)
integer, intent(in)                      :: incx


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  qtrsv  solves one of the systems of equations

!     a*x = b,   or   a**t*x = b,

!  where b and x are n element vectors and a is an n by n unit, or
!  non-unit, upper or lower triangular matrix.

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

!              trans = 'c' or 'c'   a**t*x = b.

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

!  a      - real(kind=16) array of dimension ( lda, n ).
!           before entry with  uplo = 'u' or 'u', the leading n by n
!           upper triangular part of the array a must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           a is not referenced.
!           before entry with uplo = 'l' or 'l', the leading n by n
!           lower triangular part of the array a must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           a is not referenced.
!           note that when  diag = 'u' or 'u', the diagonal elements of
!           a are not referenced either, but are assumed to be unity.
!           unchanged on exit.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. lda must be at least
!           max( 1, n ).
!           unchanged on exit.

!  x      - real(kind=16) array of dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ).
!           before entry, the incremented array x must contain the n
!           element right-hand side vector b. on exit, x is overwritten
!           with the solution vector x.

!  incx   - integer.
!           on entry, incx specifies the increment for the elements of
!           x. incx must not be zero.
!           unchanged on exit.


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
real(kind=16) :: temp
integer :: i,info,ix,j,jx,kx
logical :: nounit
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
else if (.not.lsame(trans,'n') .and. .not.lsame(trans,'t') .and.  &
      .not.lsame(trans,'c')) then
  info = 2
else if (.not.lsame(diag,'u') .and. .not.lsame(diag,'n')) then
  info = 3
else if (n < 0) then
  info = 4
else if (lda < max(1,n)) then
  info = 6
else if (incx == 0) then
  info = 8
end if
if (info /= 0) then
  call xerbla('qtrsv ',info)
  return
end if

!     quick return if possible.

if (n == 0) return

nounit = lsame(diag,'n')

!     set up the start point in x if the increment is not unity. this
!     will be  ( n - 1 )*incx  too small for descending loops.

if (incx <= 0) then
  kx = 1 - (n-1)*incx
else if (incx /= 1) then
  kx = 1
end if

!     start the operations. in this version the elements of a are
!     accessed sequentially with one pass through a.

if (lsame(trans,'n')) then
  
!        form  x := inv( a )*x.
  
  if (lsame(uplo,'u')) then
    if (incx == 1) then
      do  j = n,1,-1
        if (x(j) /= zero) then
          if (nounit) x(j) = x(j)/a(j,j)
          temp = x(j)
          do  i = j - 1,1,-1
            x(i) = x(i) - temp*a(i,j)
          end do
        end if
      end do
    else
      jx = kx + (n-1)*incx
      do  j = n,1,-1
        if (x(jx) /= zero) then
          if (nounit) x(jx) = x(jx)/a(j,j)
          temp = x(jx)
          ix = jx
          do  i = j - 1,1,-1
            ix = ix - incx
            x(ix) = x(ix) - temp*a(i,j)
          end do
        end if
        jx = jx - incx
      end do
    end if
  else
    if (incx == 1) then
      do  j = 1,n
        if (x(j) /= zero) then
          if (nounit) x(j) = x(j)/a(j,j)
          temp = x(j)
          do  i = j + 1,n
            x(i) = x(i) - temp*a(i,j)
          end do
        end if
      end do
    else
      jx = kx
      do  j = 1,n
        if (x(jx) /= zero) then
          if (nounit) x(jx) = x(jx)/a(j,j)
          temp = x(jx)
          ix = jx
          do  i = j + 1,n
            ix = ix + incx
            x(ix) = x(ix) - temp*a(i,j)
          end do
        end if
        jx = jx + incx
      end do
    end if
  end if
else
  
!        form  x := inv( a**t )*x.
  
  if (lsame(uplo,'u')) then
    if (incx == 1) then
      do  j = 1,n
        temp = x(j)
        do  i = 1,j - 1
          temp = temp - a(i,j)*x(i)
        end do
        if (nounit) temp = temp/a(j,j)
        x(j) = temp
      end do
    else
      jx = kx
      do  j = 1,n
        temp = x(jx)
        ix = kx
        do  i = 1,j - 1
          temp = temp - a(i,j)*x(ix)
          ix = ix + incx
        end do
        if (nounit) temp = temp/a(j,j)
        x(jx) = temp
        jx = jx + incx
      end do
    end if
  else
    if (incx == 1) then
      do  j = n,1,-1
        temp = x(j)
        do  i = n,j + 1,-1
          temp = temp - a(i,j)*x(i)
        end do
        if (nounit) temp = temp/a(j,j)
        x(j) = temp
      end do
    else
      kx = kx + (n-1)*incx
      jx = kx
      do  j = n,1,-1
        temp = x(jx)
        ix = kx
        do  i = n,j + 1,-1
          temp = temp - a(i,j)*x(ix)
          ix = ix - incx
        end do
        if (nounit) temp = temp/a(j,j)
        x(jx) = temp
        jx = jx - incx
      end do
    end if
  end if
end if

return

!     end of qtrsv .

end subroutine qtrsv
