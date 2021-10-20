subroutine osymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:37

!     .. scalar arguments ..

character(len=1),intent(in)        :: side
character(len=1),intent(in)        :: uplo
integer, intent(in)                      :: m
integer, intent(in)                      :: n
complex(kind=16), intent(in)               :: alpha
complex(kind=16), intent(in)               :: a(lda,*)
integer, intent(in out)                  :: lda
complex(kind=16), intent(in)               :: b(ldb,*)
integer, intent(in out)                  :: ldb
complex(kind=16), intent(in)               :: beta
complex(kind=16), intent(out)              :: c(ldc,*)
integer, intent(in out)                  :: ldc



!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  osymm  performs one of the matrix-matrix operations

!     c := alpha*a*b + beta*c,

!  or

!     c := alpha*b*a + beta*c,

!  where  alpha and beta are scalars, a is a symmetric matrix and  b and
!  c are m by n matrices.

!  arguments
!  ==========

!  side   - character*1.
!           on entry,  side  specifies whether  the  symmetric matrix  a
!           appears on the  left or right  in the  operation as follows:

!              side = 'l' or 'l'   c := alpha*a*b + beta*c,

!              side = 'r' or 'r'   c := alpha*b*a + beta*c,

!           unchanged on exit.

!  uplo   - character*1.
!           on  entry,   uplo  specifies  whether  the  upper  or  lower
!           triangular  part  of  the  symmetric  matrix   a  is  to  be
!           referenced as follows:

!              uplo = 'u' or 'u'   only the upper triangular part of the
!                                  symmetric matrix is to be referenced.

!              uplo = 'l' or 'l'   only the lower triangular part of the
!                                  symmetric matrix is to be referenced.

!           unchanged on exit.

!  m      - integer.
!           on entry,  m  specifies the number of rows of the matrix  c.
!           m  must be at least zero.
!           unchanged on exit.

!  n      - integer.
!           on entry, n specifies the number of columns of the matrix c.
!           n  must be at least zero.
!           unchanged on exit.

!  alpha  - complex(kind=16)      .
!           on entry, alpha specifies the scalar alpha.
!           unchanged on exit.

!  a      - complex(kind=16)       array of dimension ( lda, ka ), where ka is
!           m  when  side = 'l' or 'l'  and is n  otherwise.
!           before entry  with  side = 'l' or 'l',  the  m by m  part of
!           the array  a  must contain the  symmetric matrix,  such that
!           when  uplo = 'u' or 'u', the leading m by m upper triangular
!           part of the array  a  must contain the upper triangular part
!           of the  symmetric matrix and the  strictly  lower triangular
!           part of  a  is not referenced,  and when  uplo = 'l' or 'l',
!           the leading  m by m  lower triangular part  of the  array  a
!           must  contain  the  lower triangular part  of the  symmetric
!           matrix and the  strictly upper triangular part of  a  is not
!           referenced.
!           before entry  with  side = 'r' or 'r',  the  n by n  part of
!           the array  a  must contain the  symmetric matrix,  such that
!           when  uplo = 'u' or 'u', the leading n by n upper triangular
!           part of the array  a  must contain the upper triangular part
!           of the  symmetric matrix and the  strictly  lower triangular
!           part of  a  is not referenced,  and when  uplo = 'l' or 'l',
!           the leading  n by n  lower triangular part  of the  array  a
!           must  contain  the  lower triangular part  of the  symmetric
!           matrix and the  strictly upper triangular part of  a  is not
!           referenced.
!           unchanged on exit.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in the  calling (sub) program. when  side = 'l' or 'l'  then
!           lda must be at least  max( 1, m ), otherwise  lda must be at
!           least max( 1, n ).
!           unchanged on exit.

!  b      - complex(kind=16)       array of dimension ( ldb, n ).
!           before entry, the leading  m by n part of the array  b  must
!           contain the matrix b.
!           unchanged on exit.

!  ldb    - integer.
!           on entry, ldb specifies the first dimension of b as declared
!           in  the  calling  (sub)  program.   ldb  must  be  at  least
!           max( 1, m ).
!           unchanged on exit.

!  beta   - complex(kind=16)      .
!           on entry,  beta  specifies the scalar  beta.  when  beta  is
!           supplied as zero then c need not be set on input.
!           unchanged on exit.

!  c      - complex(kind=16)       array of dimension ( ldc, n ).
!           before entry, the leading  m by n  part of the array  c must
!           contain the matrix  c,  except when  beta  is zero, in which
!           case c need not be set on entry.
!           on exit, the array  c  is overwritten by the  m by n updated
!           matrix.

!  ldc    - integer.
!           on entry, ldc specifies the first dimension of c as declared
!           in  the  calling  (sub)  program.   ldc  must  be  at  least
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
intrinsic max
!     ..
!     .. local scalars ..
complex(kind=16) :: temp1,temp2
integer :: i,info,j,k,nrowa
logical :: upper
!     ..
!     .. parameters ..

complex(kind=16), parameter :: one= (1.0d+0,0.0d+0)

complex(kind=16), parameter :: zero= (0.0d+0,0.0d+0)
!     ..

!     set nrowa as the number of rows of a.

if (lsame(side,'l')) then
  nrowa = m
else
  nrowa = n
end if
upper = lsame(uplo,'u')

!     test the input parameters.

info = 0
if ((.not.lsame(side,'l')) .and. (.not.lsame(side,'r'))) then
  info = 1
else if ((.not.upper) .and. (.not.lsame(uplo,'l'))) then
  info = 2
else if (m < 0) then
  info = 3
else if (n < 0) then
  info = 4
else if (lda < max(1,nrowa)) then
  info = 7
else if (ldb < max(1,m)) then
  info = 9
else if (ldc < max(1,m)) then
  info = 12
end if
if (info /= 0) then
  call xerbla('osymm ',info)
  return
end if

!     quick return if possible.

if ((m == 0) .or. (n == 0) .or. ((alpha == zero).and. (beta == one))) return

!     and when  alpha.eq.zero.

if (alpha == zero) then
  if (beta == zero) then
    do  j = 1,n
      do  i = 1,m
        c(i,j) = zero
      end do
    end do
  else
    do  j = 1,n
      do  i = 1,m
        c(i,j) = beta*c(i,j)
      end do
    end do
  end if
  return
end if

!     start the operations.

if (lsame(side,'l')) then
  
!        form  c := alpha*a*b + beta*c.
  
  if (upper) then
    do  j = 1,n
      do  i = 1,m
        temp1 = alpha*b(i,j)
        temp2 = zero
        do  k = 1,i - 1
          c(k,j) = c(k,j) + temp1*a(k,i)
          temp2 = temp2 + b(k,j)*a(k,i)
        end do
        if (beta == zero) then
          c(i,j) = temp1*a(i,i) + alpha*temp2
        else
          c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
        end if
      end do
    end do
  else
    do  j = 1,n
      do  i = m,1,-1
        temp1 = alpha*b(i,j)
        temp2 = zero
        do  k = i + 1,m
          c(k,j) = c(k,j) + temp1*a(k,i)
          temp2 = temp2 + b(k,j)*a(k,i)
        end do
        if (beta == zero) then
          c(i,j) = temp1*a(i,i) + alpha*temp2
        else
          c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
        end if
      end do
    end do
  end if
else
  
!        form  c := alpha*b*a + beta*c.
  
  do  j = 1,n
    temp1 = alpha*a(j,j)
    if (beta == zero) then
      do  i = 1,m
        c(i,j) = temp1*b(i,j)
      end do
    else
      do  i = 1,m
        c(i,j) = beta*c(i,j) + temp1*b(i,j)
      end do
    end if
    do  k = 1,j - 1
      if (upper) then
        temp1 = alpha*a(k,j)
      else
        temp1 = alpha*a(j,k)
      end if
      do  i = 1,m
        c(i,j) = c(i,j) + temp1*b(i,k)
      end do
    end do
    do  k = j + 1,n
      if (upper) then
        temp1 = alpha*a(j,k)
      else
        temp1 = alpha*a(k,j)
      end if
      do  i = 1,m
        c(i,j) = c(i,j) + temp1*b(i,k)
      end do
    end do
  end do
end if

return

!     end of osymm .

end subroutine osymm
