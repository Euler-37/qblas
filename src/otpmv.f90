subroutine otpmv(uplo,trans,diag,n,ap,x,incx)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:37

!     .. scalar arguments ..

character(len=1),intent(in)        :: uplo
character(len=1),intent(in)        :: trans
character(len=1),intent(in)        :: diag
integer, intent(in)                      :: n
complex(kind=16), intent(in)               :: ap(*)
complex(kind=16), intent(in out)           :: x(*)
integer, intent(in)                      :: incx


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  otpmv  performs one of the matrix-vector operations

!     x := a*x,   or   x := a**t*x,   or   x := a**h*x,

!  where x is an n element vector and  a is an n by n unit, or non-unit,
!  upper or lower triangular matrix, supplied in packed form.

!  arguments
!  ==========

!  uplo   - character*1.
!           on entry, uplo specifies whether the matrix is an upper or
!           lower triangular matrix as follows:

!              uplo = 'u' or 'u'   a is an upper triangular matrix.

!              uplo = 'l' or 'l'   a is a lower triangular matrix.

!           unchanged on exit.

!  trans  - character*1.
!           on entry, trans specifies the operation to be performed as
!           follows:

!              trans = 'n' or 'n'   x := a*x.

!              trans = 't' or 't'   x := a**t*x.

!              trans = 'c' or 'c'   x := a**h*x.

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

!  ap     - complex(kind=16)       array of dimension at least
!           ( ( n*( n + 1 ) )/2 ).
!           before entry with  uplo = 'u' or 'u', the array ap must
!           contain the upper triangular matrix packed sequentially,
!           column by column, so that ap( 1 ) contains a( 1, 1 ),
!           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )
!           respectively, and so on.
!           before entry with uplo = 'l' or 'l', the array ap must
!           contain the lower triangular matrix packed sequentially,
!           column by column, so that ap( 1 ) contains a( 1, 1 ),
!           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )
!           respectively, and so on.
!           note that when  diag = 'u' or 'u', the diagonal elements of
!           a are not referenced, but are assumed to be unity.
!           unchanged on exit.

!  x      - complex(kind=16)       array of dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ).
!           before entry, the incremented array x must contain the n
!           element vector x. on exit, x is overwritten with the
!           tranformed vector x.

!  incx   - integer.
!           on entry, incx specifies the increment for the elements of
!           x. incx must not be zero.
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

complex(kind=16), parameter :: zero= (0.0d+0,0.0d+0)
!     ..
!     .. local scalars ..
complex(kind=16) :: temp
integer :: i,info,ix,j,jx,k,kk,kx
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
intrinsic conjg
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
else if (incx == 0) then
  info = 7
end if
if (info /= 0) then
  call xerbla('otpmv ',info)
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

!     start the operations. in this version the elements of ap are
!     accessed sequentially with one pass through ap.

if (lsame(trans,'n')) then
  
!        form  x:= a*x.
  
  if (lsame(uplo,'u')) then
    kk = 1
    if (incx == 1) then
      do  j = 1,n
        if (x(j) /= zero) then
          temp = x(j)
          k = kk
          do  i = 1,j - 1
            x(i) = x(i) + temp*ap(k)
            k = k + 1
          end do
          if (nounit) x(j) = x(j)*ap(kk+j-1)
        end if
        kk = kk + j
      end do
    else
      jx = kx
      do  j = 1,n
        if (x(jx) /= zero) then
          temp = x(jx)
          ix = kx
          do  k = kk,kk + j - 2
            x(ix) = x(ix) + temp*ap(k)
            ix = ix + incx
          end do
          if (nounit) x(jx) = x(jx)*ap(kk+j-1)
        end if
        jx = jx + incx
        kk = kk + j
      end do
    end if
  else
    kk = (n* (n+1))/2
    if (incx == 1) then
      do  j = n,1,-1
        if (x(j) /= zero) then
          temp = x(j)
          k = kk
          do  i = n,j + 1,-1
            x(i) = x(i) + temp*ap(k)
            k = k - 1
          end do
          if (nounit) x(j) = x(j)*ap(kk-n+j)
        end if
        kk = kk - (n-j+1)
      end do
    else
      kx = kx + (n-1)*incx
      jx = kx
      do  j = n,1,-1
        if (x(jx) /= zero) then
          temp = x(jx)
          ix = kx
          do  k = kk,kk - (n- (j+1)),-1
            x(ix) = x(ix) + temp*ap(k)
            ix = ix - incx
          end do
          if (nounit) x(jx) = x(jx)*ap(kk-n+j)
        end if
        jx = jx - incx
        kk = kk - (n-j+1)
      end do
    end if
  end if
else
  
!        form  x := a**t*x  or  x := a**h*x.
  
  if (lsame(uplo,'u')) then
    kk = (n* (n+1))/2
    if (incx == 1) then
      do  j = n,1,-1
        temp = x(j)
        k = kk - 1
        if (noconj) then
          if (nounit) temp = temp*ap(kk)
          do  i = j - 1,1,-1
            temp = temp + ap(k)*x(i)
            k = k - 1
          end do
        else
          if (nounit) temp = temp*conjg(ap(kk))
          do  i = j - 1,1,-1
            temp = temp + conjg(ap(k))*x(i)
            k = k - 1
          end do
        end if
        x(j) = temp
        kk = kk - j
      end do
    else
      jx = kx + (n-1)*incx
      do  j = n,1,-1
        temp = x(jx)
        ix = jx
        if (noconj) then
          if (nounit) temp = temp*ap(kk)
          do  k = kk - 1,kk - j + 1,-1
            ix = ix - incx
            temp = temp + ap(k)*x(ix)
          end do
        else
          if (nounit) temp = temp*conjg(ap(kk))
          do  k = kk - 1,kk - j + 1,-1
            ix = ix - incx
            temp = temp + conjg(ap(k))*x(ix)
          end do
        end if
        x(jx) = temp
        jx = jx - incx
        kk = kk - j
      end do
    end if
  else
    kk = 1
    if (incx == 1) then
      do  j = 1,n
        temp = x(j)
        k = kk + 1
        if (noconj) then
          if (nounit) temp = temp*ap(kk)
          do  i = j + 1,n
            temp = temp + ap(k)*x(i)
            k = k + 1
          end do
        else
          if (nounit) temp = temp*conjg(ap(kk))
          do  i = j + 1,n
            temp = temp + conjg(ap(k))*x(i)
            k = k + 1
          end do
        end if
        x(j) = temp
        kk = kk + (n-j+1)
      end do
    else
      jx = kx
      do  j = 1,n
        temp = x(jx)
        ix = jx
        if (noconj) then
          if (nounit) temp = temp*ap(kk)
          do  k = kk + 1,kk + n - j
            ix = ix + incx
            temp = temp + ap(k)*x(ix)
          end do
        else
          if (nounit) temp = temp*conjg(ap(kk))
          do  k = kk + 1,kk + n - j
            ix = ix + incx
            temp = temp + conjg(ap(k))*x(ix)
          end do
        end if
        x(jx) = temp
        jx = jx + incx
        kk = kk + (n-j+1)
      end do
    end if
  end if
end if

return

!     end of otpmv .

end subroutine otpmv
