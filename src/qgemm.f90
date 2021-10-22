subroutine qgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:30

!     .. scalar arguments ..

character(len=*),intent(in)        :: transa
character(len=*),intent(in)        :: transb
integer, intent(in)                      :: m
integer, intent(in)                      :: n
integer, intent(in)                      :: k
real(kind=16), intent(in)             :: alpha
real(kind=16), intent(in)             :: a(lda,*)
integer, intent(in out)                  :: lda
real(kind=16), intent(in)             :: b(ldb,*)
integer, intent(in out)                  :: ldb
real(kind=16), intent(in)             :: beta
real(kind=16), intent(out)            :: c(ldc,*)
integer, intent(in out)                  :: ldc



!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  qgemm  performs one of the matrix-matrix operations

!     c := alpha*op( a )*op( b ) + beta*c,

!  where  op( x ) is one of

!     op( x ) = x   or   op( x ) = x**t,

!  alpha and beta are scalars, and a, b and c are matrices, with op( a )
!  an m by k matrix,  op( b )  a  k by n matrix and  c an m by n matrix.

!  arguments
!  ==========

!  transa - character*1.
!           on entry, transa specifies the form of op( a ) to be used in
!           the matrix multiplication as follows:

!              transa = 'n' or 'n',  op( a ) = a.

!              transa = 't' or 't',  op( a ) = a**t.

!              transa = 'c' or 'c',  op( a ) = a**t.

!           unchanged on exit.

!  transb - character*1.
!           on entry, transb specifies the form of op( b ) to be used in
!           the matrix multiplication as follows:

!              transb = 'n' or 'n',  op( b ) = b.

!              transb = 't' or 't',  op( b ) = b**t.

!              transb = 'c' or 'c',  op( b ) = b**t.

!           unchanged on exit.

!  m      - integer.
!           on entry,  m  specifies  the number  of rows  of the  matrix
!           op( a )  and of the  matrix  c.  m  must  be at least  zero.
!           unchanged on exit.

!  n      - integer.
!           on entry,  n  specifies the number  of columns of the matrix
!           op( b ) and the number of columns of the matrix c. n must be
!           at least zero.
!           unchanged on exit.

!  k      - integer.
!           on entry,  k  specifies  the number of columns of the matrix
!           op( a ) and the number of rows of the matrix op( b ). k must
!           be at least  zero.
!           unchanged on exit.

!  alpha  - real(kind=16).
!           on entry, alpha specifies the scalar alpha.
!           unchanged on exit.

!  a      - real(kind=16) array of dimension ( lda, ka ), where ka is
!           k  when  transa = 'n' or 'n',  and is  m  otherwise.
!           before entry with  transa = 'n' or 'n',  the leading  m by k
!           part of the array  a  must contain the matrix  a,  otherwise
!           the leading  k by m  part of the array  a  must contain  the
!           matrix a.
!           unchanged on exit.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in the calling (sub) program. when  transa = 'n' or 'n' then
!           lda must be at least  max( 1, m ), otherwise  lda must be at
!           least  max( 1, k ).
!           unchanged on exit.

!  b      - real(kind=16) array of dimension ( ldb, kb ), where kb is
!           n  when  transb = 'n' or 'n',  and is  k  otherwise.
!           before entry with  transb = 'n' or 'n',  the leading  k by n
!           part of the array  b  must contain the matrix  b,  otherwise
!           the leading  n by k  part of the array  b  must contain  the
!           matrix b.
!           unchanged on exit.

!  ldb    - integer.
!           on entry, ldb specifies the first dimension of b as declared
!           in the calling (sub) program. when  transb = 'n' or 'n' then
!           ldb must be at least  max( 1, k ), otherwise  ldb must be at
!           least  max( 1, n ).
!           unchanged on exit.

!  beta   - real(kind=16).
!           on entry,  beta  specifies the scalar  beta.  when  beta  is
!           supplied as zero then c need not be set on input.
!           unchanged on exit.

!  c      - real(kind=16) array of dimension ( ldc, n ).
!           before entry, the leading  m by n  part of the array  c must
!           contain the matrix  c,  except when  beta  is zero, in which
!           case c need not be set on entry.
!           on exit, the array  c  is overwritten by the  m by n  matrix
!           ( alpha*op( a )*op( b ) + beta*c ).

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
real(kind=16) :: temp
integer :: i,info,j,l,ncola,nrowa,nrowb
logical :: nota,notb
!     ..
!     .. parameters ..

real(kind=16), parameter :: one=1.0d+0
real(kind=16), parameter :: zero=0.0d+0
!     ..

!     set  nota  and  notb  as  true if  a  and  b  respectively are not
!     transposed and set  nrowa, ncola and  nrowb  as the number of rows
!     and  columns of  a  and the  number of  rows  of  b  respectively.

nota = lsame(transa,'n')
notb = lsame(transb,'n')
if (nota) then
  nrowa = m
  ncola = k
else
  nrowa = k
  ncola = m
end if
if (notb) then
  nrowb = k
else
  nrowb = n
end if

!     test the input parameters.

info = 0
if ((.not.nota) .and. (.not.lsame(transa,'c')) .and.  &
      (.not.lsame(transa,'t'))) then
  info = 1
else if ((.not.notb) .and. (.not.lsame(transb,'c')) .and.  &
      (.not.lsame(transb,'t'))) then
  info = 2
else if (m < 0) then
  info = 3
else if (n < 0) then
  info = 4
else if (k < 0) then
  info = 5
else if (lda < max(1,nrowa)) then
  info = 8
else if (ldb < max(1,nrowb)) then
  info = 10
else if (ldc < max(1,m)) then
  info = 13
end if
if (info /= 0) then
  call xerbla('qgemm ',info)
  return
end if

!     quick return if possible.

if ((m == 0) .or. (n == 0) .or.  &
    (((alpha == zero).or. (k == 0)).and. (beta == one))) return

!     and if  alpha.eq.zero.

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

if (notb) then
  if (nota) then
    
!           form  c := alpha*a*b + beta*c.
    
    do  j = 1,n
      if (beta == zero) then
        do  i = 1,m
          c(i,j) = zero
        end do
      else if (beta /= one) then
        do  i = 1,m
          c(i,j) = beta*c(i,j)
        end do
      end if
      do  l = 1,k
        if (b(l,j) /= zero) then
          temp = alpha*b(l,j)
          do  i = 1,m
            c(i,j) = c(i,j) + temp*a(i,l)
          end do
        end if
      end do
    end do
  else
    
!           form  c := alpha*a**t*b + beta*c
    
    do  j = 1,n
      do  i = 1,m
        temp = zero
        do  l = 1,k
          temp = temp + a(l,i)*b(l,j)
        end do
        if (beta == zero) then
          c(i,j) = alpha*temp
        else
          c(i,j) = alpha*temp + beta*c(i,j)
        end if
      end do
    end do
  end if
else
  if (nota) then
    
!           form  c := alpha*a*b**t + beta*c
    
    do  j = 1,n
      if (beta == zero) then
        do  i = 1,m
          c(i,j) = zero
        end do
      else if (beta /= one) then
        do  i = 1,m
          c(i,j) = beta*c(i,j)
        end do
      end if
      do  l = 1,k
        if (b(j,l) /= zero) then
          temp = alpha*b(j,l)
          do  i = 1,m
            c(i,j) = c(i,j) + temp*a(i,l)
          end do
        end if
      end do
    end do
  else
    
!           form  c := alpha*a**t*b**t + beta*c
    
    do  j = 1,n
      do  i = 1,m
        temp = zero
        do  l = 1,k
          temp = temp + a(l,i)*b(j,l)
        end do
        if (beta == zero) then
          c(i,j) = alpha*temp
        else
          c(i,j) = alpha*temp + beta*c(i,j)
        end if
      end do
    end do
  end if
end if

return

!     end of qgemm .

end subroutine qgemm
