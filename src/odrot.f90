subroutine odrot( n, cx, incx, cy, incy, c, s )

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:35

!     .. scalar arguments ..

integer, intent(in)                      :: n
complex(kind=16), intent(in out)               :: cx( * )
integer, intent(in)                      :: incx
complex(kind=16), intent(in out)               :: cy( * )
integer, intent(in)                      :: incy
real(kind=16), intent(in)             :: c
real(kind=16), intent(in)             :: s


!     ..
!     .. array arguments ..

!     ..

!  purpose
!  =======

!  applies a plane rotation, where the cos and sin (c and s) are real
!  and the vectors cx and cy are complex.
!  jack dongarra, linpack, 3/11/78.

!  arguments
!  ==========

!  n        (input) integer
!           on entry, n specifies the order of the vectors cx and cy.
!           n must be at least zero.
!           unchanged on exit.

!  cx       (input) complex(kind=16) array, dimension at least
!           ( 1 + ( n - 1 )*abs( incx ) ).
!           before entry, the incremented array cx must contain the n
!           element vector cx. on exit, cx is overwritten by the updated
!           vector cx.

!  incx     (input) integer
!           on entry, incx specifies the increment for the elements of
!           cx. incx must not be zero.
!           unchanged on exit.

!  cy       (input) complex(kind=16) array, dimension at least
!           ( 1 + ( n - 1 )*abs( incy ) ).
!           before entry, the incremented array cy must contain the n
!           element vector cy. on exit, cy is overwritten by the updated
!           vector cy.

!  incy     (input) integer
!           on entry, incy specifies the increment for the elements of
!           cy. incy must not be zero.
!           unchanged on exit.

!  c        (input) real(kind=16)
!           on entry, c specifies the cosine, cos.
!           unchanged on exit.

!  s        (input) real(kind=16)
!           on entry, s specifies the sine, sin.
!           unchanged on exit.

! =====================================================================

!     .. local scalars ..
integer :: i, ix, iy
complex(kind=16)         ctemp
!     ..
!     .. executable statements ..

if( n <= 0 ) return
if( incx == 1 .and. incy == 1 ) then
  
!        code for both increments equal to 1
  
  do i = 1, n
    ctemp = c*cx( i ) + s*cy( i )
    cy( i ) = c*cy( i ) - s*cx( i )
    cx( i ) = ctemp
  end do
else
  
!        code for unequal increments or equal increments not equal
!          to 1
  
  ix = 1
  iy = 1
  if( incx < 0 ) ix = ( -n+1 )*incx + 1
  if( incy < 0 ) iy = ( -n+1 )*incy + 1
  do i = 1, n
    ctemp = c*cx( ix ) + s*cy( iy )
    cy( iy ) = c*cy( iy ) - s*cx( ix )
    cx( ix ) = ctemp
    ix = ix + incx
    iy = iy + incy
  end do
end if
return
end subroutine odrot
