subroutine oherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:36

!     .. scalar arguments ..

character(len=*),intent(in)        :: uplo
character(len=*),intent(in)        :: trans
integer, intent(in)                      :: n
integer, intent(in)                      :: k
real(kind=16), intent(in)             :: alpha
complex(kind=16), intent(in)               :: a(lda,*)
integer, intent(in out)                  :: lda
real(kind=16), intent(in)             :: beta
complex(kind=16), intent(out)              :: c(ldc,*)
integer, intent(in out)                  :: ldc



!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  oherk  performs one of the hermitian rank k operations

!     c := alpha*a*a**h + beta*c,

!  or

!     c := alpha*a**h*a + beta*c,

!  where  alpha and beta  are  real scalars,  c is an  n by n  hermitian
!  matrix and  a  is an  n by k  matrix in the  first case and a  k by n
!  matrix in the second case.

!  arguments
!  ==========

!  uplo   - character*1.
!           on  entry,   uplo  specifies  whether  the  upper  or  lower
!           triangular  part  of the  array  c  is to be  referenced  as
!           follows:

!              uplo = 'u' or 'u'   only the  upper triangular part of  c
!                                  is to be referenced.

!              uplo = 'l' or 'l'   only the  lower triangular part of  c
!                                  is to be referenced.

!           unchanged on exit.

!  trans  - character*1.
!           on entry,  trans  specifies the operation to be performed as
!           follows:

!              trans = 'n' or 'n'   c := alpha*a*a**h + beta*c.

!              trans = 'c' or 'c'   c := alpha*a**h*a + beta*c.

!           unchanged on exit.

!  n      - integer.
!           on entry,  n specifies the order of the matrix c.  n must be
!           at least zero.
!           unchanged on exit.

!  k      - integer.
!           on entry with  trans = 'n' or 'n',  k  specifies  the number
!           of  columns   of  the   matrix   a,   and  on   entry   with
!           trans = 'c' or 'c',  k  specifies  the number of rows of the
!           matrix a.  k must be at least zero.
!           unchanged on exit.

!  alpha  - real(kind=16)            .
!           on entry, alpha specifies the scalar alpha.
!           unchanged on exit.

!  a      - complex(kind=16)       array of dimension ( lda, ka ), where ka is
!           k  when  trans = 'n' or 'n',  and is  n  otherwise.
!           before entry with  trans = 'n' or 'n',  the  leading  n by k
!           part of the array  a  must contain the matrix  a,  otherwise
!           the leading  k by n  part of the array  a  must contain  the
!           matrix a.
!           unchanged on exit.

!  lda    - integer.
!           on entry, lda specifies the first dimension of a as declared
!           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
!           then  lda must be at least  max( 1, n ), otherwise  lda must
!           be at least  max( 1, k ).
!           unchanged on exit.

!  beta   - real(kind=16).
!           on entry, beta specifies the scalar beta.
!           unchanged on exit.

!  c      - complex(kind=16)          array of dimension ( ldc, n ).
!           before entry  with  uplo = 'u' or 'u',  the leading  n by n
!           upper triangular part of the array c must contain the upper
!           triangular part  of the  hermitian matrix  and the strictly
!           lower triangular part of c is not referenced.  on exit, the
!           upper triangular part of the array  c is overwritten by the
!           upper triangular part of the updated matrix.
!           before entry  with  uplo = 'l' or 'l',  the leading  n by n
!           lower triangular part of the array c must contain the lower
!           triangular part  of the  hermitian matrix  and the strictly
!           upper triangular part of c is not referenced.  on exit, the
!           lower triangular part of the array  c is overwritten by the
!           lower triangular part of the updated matrix.
!           note that the imaginary parts of the diagonal elements need
!           not be set,  they are assumed to be zero,  and on exit they
!           are set to zero.

!  ldc    - integer.
!           on entry, ldc specifies the first dimension of c as declared
!           in  the  calling  (sub)  program.   ldc  must  be  at  least
!           max( 1, n ).
!           unchanged on exit.

!  further details
!  ===============

!  level 3 blas routine.

!  -- written on 8-february-1989.
!     jack dongarra, argonne national laboratory.
!     iain duff, aere harwell.
!     jeremy du croz, numerical algorithms group ltd.
!     sven hammarling, numerical algorithms group ltd.

!  -- modified 8-nov-93 to set c(j,j) to real(kind=16, c(j,j) ) when beta = 1.
!     ed anderson, cray research inc.

!  =====================================================================

!     .. external functions ..
logical :: lsame
external lsame
!     ..
!     .. external subroutines ..
external xerbla
!     ..
!     .. intrinsic functions ..
intrinsic dble,dcmplx,conjg,max
!     ..
!     .. local scalars ..
complex(kind=16) :: temp
real(kind=16) :: rtemp
integer :: i,info,j,l,nrowa
logical :: upper
!     ..
!     .. parameters ..

real(kind=16), parameter :: one=1.0d+0
real(kind=16), parameter :: zero=0.0d+0
!     ..

!     test the input parameters.

if (lsame(trans,'n')) then
  nrowa = n
else
  nrowa = k
end if
upper = lsame(uplo,'u')

info = 0
if ((.not.upper) .and. (.not.lsame(uplo,'l'))) then
  info = 1
else if ((.not.lsame(trans,'n')) .and.  &
      (.not.lsame(trans,'c'))) then
  info = 2
else if (n < 0) then
  info = 3
else if (k < 0) then
  info = 4
else if (lda < max(1,nrowa)) then
  info = 7
else if (ldc < max(1,n)) then
  info = 10
end if
if (info /= 0) then
  call xerbla('oherk ',info)
  return
end if

!     quick return if possible.

if ((n == 0) .or. (((alpha == zero).or. (k == 0)).and. (beta == one))) return

!     and when  alpha.eq.zero.

if (alpha == zero) then
  if (upper) then
    if (beta == zero) then
      do  j = 1,n
        do  i = 1,j
          c(i,j) = zero
        end do
      end do
    else
      do  j = 1,n
        do  i = 1,j - 1
          c(i,j) = beta*c(i,j)
        end do
        c(j,j) = beta*real(c(j,j),kind=16)
      end do
    end if
  else
    if (beta == zero) then
      do  j = 1,n
        do  i = j,n
          c(i,j) = zero
        end do
      end do
    else
      do  j = 1,n
        c(j,j) = beta*real(c(j,j),kind=16)
        do  i = j + 1,n
          c(i,j) = beta*c(i,j)
        end do
      end do
    end if
  end if
  return
end if

!     start the operations.

if (lsame(trans,'n')) then
  
!        form  c := alpha*a*a**h + beta*c.
  
  if (upper) then
    do  j = 1,n
      if (beta == zero) then
        do  i = 1,j
          c(i,j) = zero
        end do
      else if (beta /= one) then
        do  i = 1,j - 1
          c(i,j) = beta*c(i,j)
        end do
        c(j,j) = beta*real(c(j,j),kind=16)
      else
        c(j,j) = real(c(j,j),kind=16)
      end if
      do  l = 1,k
        if (a(j,l) /= dcmplx(zero)) then
          temp = alpha*conjg(a(j,l))
          do  i = 1,j - 1
            c(i,j) = c(i,j) + temp*a(i,l)
          end do
          c(j,j) = real(c(j,j),kind=16) + real(temp*a(i,l),kind=16)
        end if
      end do
    end do
  else
    do  j = 1,n
      if (beta == zero) then
        do  i = j,n
          c(i,j) = zero
        end do
      else if (beta /= one) then
        c(j,j) = beta*real(c(j,j),kind=16)
        do  i = j + 1,n
          c(i,j) = beta*c(i,j)
        end do
      else
        c(j,j) = real(c(j,j),kind=16)
      end if
      do  l = 1,k
        if (a(j,l) /= dcmplx(zero)) then
          temp = alpha*conjg(a(j,l))
          c(j,j) = real(c(j,j),kind=16) + real(temp*a(j,l),kind=16)
          do  i = j + 1,n
            c(i,j) = c(i,j) + temp*a(i,l)
          end do
        end if
      end do
    end do
  end if
else
  
!        form  c := alpha*a**h*a + beta*c.
  
  if (upper) then
    do  j = 1,n
      do  i = 1,j - 1
        temp = zero
        do  l = 1,k
          temp = temp + conjg(a(l,i))*a(l,j)
        end do
        if (beta == zero) then
          c(i,j) = alpha*temp
        else
          c(i,j) = alpha*temp + beta*c(i,j)
        end if
      end do
      rtemp = zero
      do  l = 1,k
        rtemp = rtemp + conjg(a(l,j))*a(l,j)
      end do
      if (beta == zero) then
        c(j,j) = alpha*rtemp
      else
        c(j,j) = alpha*rtemp + beta*real(c(j,j),kind=16)
      end if
    end do
  else
    do  j = 1,n
      rtemp = zero
      do  l = 1,k
        rtemp = rtemp + conjg(a(l,j))*a(l,j)
      end do
      if (beta == zero) then
        c(j,j) = alpha*rtemp
      else
        c(j,j) = alpha*rtemp + beta*real(c(j,j),kind=16)
      end if
      do  i = j + 1,n
        temp = zero
        do  l = 1,k
          temp = temp + conjg(a(l,i))*a(l,j)
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

!     end of oherk .

end subroutine oherk
