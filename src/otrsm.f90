subroutine otrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:37

!     .. scalar arguments ..

character(len=*),intent(in)        :: side
character(len=*),intent(in)        :: uplo
character(len=*),intent(in)        :: transa
character(len=*),intent(in)        :: diag
integer, intent(in)                      :: m
integer, intent(in)                      :: n
complex(kind=16), intent(in)               :: alpha
complex(kind=16), intent(in)               :: a(lda,*)
integer, intent(in out)                  :: lda
complex(kind=16), intent(out)              :: b(ldb,*)
integer, intent(in out)                  :: ldb



!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  otrsm  solves one of the matrix equations

!     op( a )*x = alpha*b,   or   x*op( a ) = alpha*b,

!  where alpha is a scalar, x and b are m by n matrices, a is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( a )  is one  of

!     op( a ) = a   or   op( a ) = a**t   or   op( a ) = a**h.

!  the matrix x is overwritten on b.

!  arguments
!  ==========

!  side   - character*1.
!           on entry, side specifies whether op( a ) appears on the left
!           or right of x as follows:

!              side = 'l' or 'l'   op( a )*x = alpha*b.

!              side = 'r' or 'r'   x*op( a ) = alpha*b.

!           unchanged on exit.

!  uplo   - character*1.
!           on entry, uplo specifies whether the matrix a is an upper or
!           lower triangular matrix as follows:

!              uplo = 'u' or 'u'   a is an upper triangular matrix.

!              uplo = 'l' or 'l'   a is a lower triangular matrix.

!           unchanged on exit.

!  transa - character*1.
!           on entry, transa specifies the form of op( a ) to be used in
!           the matrix multiplication as follows:

!              transa = 'n' or 'n'   op( a ) = a.

!              transa = 't' or 't'   op( a ) = a**t.

!              transa = 'c' or 'c'   op( a ) = a**h.

!           unchanged on exit.

!  diag   - character*1.
!           on entry, diag specifies whether or not a is unit triangular
!           as follows:

!              diag = 'u' or 'u'   a is assumed to be unit triangular.

!              diag = 'n' or 'n'   a is not assumed to be unit
!                                  triangular.

!           unchanged on exit.

!  m      - integer.
!           on entry, m specifies the number of rows of b. m must be at
!           least zero.
!           unchanged on exit.

!  n      - integer.
!           on entry, n specifies the number of columns of b.  n must be
!           at least zero.
!           unchanged on exit.

!  alpha  - complex(kind=16)      .
!           on entry,  alpha specifies the scalar  alpha. when  alpha is
!           zero then  a is not referenced and  b need not be set before
!           entry.
!           unchanged on exit.

!  a      - complex(kind=16)       array of dimension ( lda, k ), where k is m
!           when  side = 'l' or 'l'  and is  n  when  side = 'r' or 'r'.
!           before entry  with  uplo = 'u' or 'u',  the  leading  k by k
!           upper triangular part of the array  a must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           a is not referenced.
!           before entry  with  uplo = 'l' or 'l',  the  leading  k by k
!           lower triangular part of the array  a must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           a is not referenced.
!           note that when  diag = 'u' or 'u',  the diagonal elements of
!           a  are not referenced either,  but are assumed to be  unity.
!           unchanged on exit.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program.  when  side = 'l' or 'l'  then
!           lda  must be at least  max( 1, m ),  when  side = 'r' or 'r'
!           then lda must be at least max( 1, n ).
!           unchanged on exit.

!  b      - complex(kind=16)       array of dimension ( ldb, n ).
!           before entry,  the leading  m by n part of the array  b must
!           contain  the  right-hand  side  matrix  b,  and  on exit  is
!           overwritten by the solution matrix  x.

!  ldb    - integer.
!           on entry, ldb specifies the first dimension of b as declared
!           in  the  calling  (sub)  program.   ldb  must  be  at  least
!           max( 1, m ).
!           unchanged on exit.

!  further details
!  ===============

!  level 3 blas routine.

!  -- written on 8-february-1989.
!     jack dongarra, argonne national laboratory.
!     iain duff, aere harwell.
!     jeremy du croz, numerical algorithms group ltd.
!     sven hammarling, numerical algorithms group ltd.

!  =====================================================================

!     .. external functions ..
logical :: lsame
external lsame
!     ..
!     .. external subroutines ..
external xerbla
!     ..
!     .. intrinsic functions ..
intrinsic conjg,max
!     ..
!     .. local scalars ..
complex(kind=16) :: temp
integer :: i,info,j,k,nrowa
logical :: lside,noconj,nounit,upper
!     ..
!     .. parameters ..

complex(kind=16), parameter :: one= (1.0d+0,0.0d+0)

complex(kind=16), parameter :: zero= (0.0d+0,0.0d+0)
!     ..

!     test the input parameters.

lside = lsame(side,'l')
if (lside) then
  nrowa = m
else
  nrowa = n
end if
noconj = lsame(transa,'t')
nounit = lsame(diag,'n')
upper = lsame(uplo,'u')

info = 0
if ((.not.lside) .and. (.not.lsame(side,'r'))) then
  info = 1
else if ((.not.upper) .and. (.not.lsame(uplo,'l'))) then
  info = 2
else if ((.not.lsame(transa,'n')) .and.  &
      (.not.lsame(transa,'t')) .and. (.not.lsame(transa,'c'))) then
  info = 3
else if ((.not.lsame(diag,'u')) .and. (.not.lsame(diag,'n'))) then
  info = 4
else if (m < 0) then
  info = 5
else if (n < 0) then
  info = 6
else if (lda < max(1,nrowa)) then
  info = 9
else if (ldb < max(1,m)) then
  info = 11
end if
if (info /= 0) then
  call xerbla('otrsm ',info)
  return
end if

!     quick return if possible.

if (m == 0 .or. n == 0) return

!     and when  alpha.eq.zero.

if (alpha == zero) then
  do  j = 1,n
    do  i = 1,m
      b(i,j) = zero
    end do
  end do
  return
end if

!     start the operations.

if (lside) then
  if (lsame(transa,'n')) then
    
!           form  b := alpha*inv( a )*b.
    
    if (upper) then
      do  j = 1,n
        if (alpha /= one) then
          do  i = 1,m
            b(i,j) = alpha*b(i,j)
          end do
        end if
        do  k = m,1,-1
          if (b(k,j) /= zero) then
            if (nounit) b(k,j) = b(k,j)/a(k,k)
            do  i = 1,k - 1
              b(i,j) = b(i,j) - b(k,j)*a(i,k)
            end do
          end if
        end do
      end do
    else
      do  j = 1,n
        if (alpha /= one) then
          do  i = 1,m
            b(i,j) = alpha*b(i,j)
          end do
        end if
        do  k = 1,m
          if (b(k,j) /= zero) then
            if (nounit) b(k,j) = b(k,j)/a(k,k)
            do  i = k + 1,m
              b(i,j) = b(i,j) - b(k,j)*a(i,k)
            end do
          end if
        end do
      end do
    end if
  else
    
!           form  b := alpha*inv( a**t )*b
!           or    b := alpha*inv( a**h )*b.
    
    if (upper) then
      do  j = 1,n
        do  i = 1,m
          temp = alpha*b(i,j)
          if (noconj) then
            do  k = 1,i - 1
              temp = temp - a(k,i)*b(k,j)
            end do
            if (nounit) temp = temp/a(i,i)
          else
            do  k = 1,i - 1
              temp = temp - conjg(a(k,i))*b(k,j)
            end do
            if (nounit) temp = temp/conjg(a(i,i))
          end if
          b(i,j) = temp
        end do
      end do
    else
      do  j = 1,n
        do  i = m,1,-1
          temp = alpha*b(i,j)
          if (noconj) then
            do  k = i + 1,m
              temp = temp - a(k,i)*b(k,j)
            end do
            if (nounit) temp = temp/a(i,i)
          else
            do  k = i + 1,m
              temp = temp - conjg(a(k,i))*b(k,j)
            end do
            if (nounit) temp = temp/conjg(a(i,i))
          end if
          b(i,j) = temp
        end do
      end do
    end if
  end if
else
  if (lsame(transa,'n')) then
    
!           form  b := alpha*b*inv( a ).
    
    if (upper) then
      do  j = 1,n
        if (alpha /= one) then
          do  i = 1,m
            b(i,j) = alpha*b(i,j)
          end do
        end if
        do  k = 1,j - 1
          if (a(k,j) /= zero) then
            do  i = 1,m
              b(i,j) = b(i,j) - a(k,j)*b(i,k)
            end do
          end if
        end do
        if (nounit) then
          temp = one/a(j,j)
          do  i = 1,m
            b(i,j) = temp*b(i,j)
          end do
        end if
      end do
    else
      do  j = n,1,-1
        if (alpha /= one) then
          do  i = 1,m
            b(i,j) = alpha*b(i,j)
          end do
        end if
        do  k = j + 1,n
          if (a(k,j) /= zero) then
            do  i = 1,m
              b(i,j) = b(i,j) - a(k,j)*b(i,k)
            end do
          end if
        end do
        if (nounit) then
          temp = one/a(j,j)
          do  i = 1,m
            b(i,j) = temp*b(i,j)
          end do
        end if
      end do
    end if
  else
    
!           form  b := alpha*b*inv( a**t )
!           or    b := alpha*b*inv( a**h ).
    
    if (upper) then
      do  k = n,1,-1
        if (nounit) then
          if (noconj) then
            temp = one/a(k,k)
          else
            temp = one/conjg(a(k,k))
          end if
          do  i = 1,m
            b(i,k) = temp*b(i,k)
          end do
        end if
        do  j = 1,k - 1
          if (a(j,k) /= zero) then
            if (noconj) then
              temp = a(j,k)
            else
              temp = conjg(a(j,k))
            end if
            do  i = 1,m
              b(i,j) = b(i,j) - temp*b(i,k)
            end do
          end if
        end do
        if (alpha /= one) then
          do  i = 1,m
            b(i,k) = alpha*b(i,k)
          end do
        end if
      end do
    else
      do  k = 1,n
        if (nounit) then
          if (noconj) then
            temp = one/a(k,k)
          else
            temp = one/conjg(a(k,k))
          end if
          do  i = 1,m
            b(i,k) = temp*b(i,k)
          end do
        end if
        do  j = k + 1,n
          if (a(j,k) /= zero) then
            if (noconj) then
              temp = a(j,k)
            else
              temp = conjg(a(j,k))
            end if
            do  i = 1,m
              b(i,j) = b(i,j) - temp*b(i,k)
            end do
          end if
        end do
        if (alpha /= one) then
          do  i = 1,m
            b(i,k) = alpha*b(i,k)
          end do
        end if
      end do
    end if
  end if
end if

return

!     end of otrsm .

end subroutine otrsm
