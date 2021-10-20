real(kind=16) function qsdot(n,sx,incx,sy,incy)

! code converted using to_f90 by alan miller
! date: 2021-09-18  time: 18:57:31

!     .. scalar arguments ..

integer, intent(in)                      :: n
real, intent(in out)                     :: sx(*)
integer, intent(in)                      :: incx
real, intent(in out)                     :: sy(*)
integer, intent(in)                      :: incy

!     ..
!     .. array arguments ..

!     ..

!  authors
!  =======
!  lawson, c. l., (jpl), hanson, r. j., (snla),
!  kincaid, d. r., (u. of texas), krogh, f. t., (jpl)

!  purpose
!  =======
!  compute the inner product of two vectors with extended
!  precision accumulation and result.

!  returns d.p. dot product accumulated in d.p., for s.p. sx and sy
!  qsdot = sum for i = 0 to n-1 of  sx(lx+i*incx) * sy(ly+i*incy),
!  where lx = 1 if incx .ge. 0, else lx = 1+(1-n)*incx, and ly is
!  defined in a similar way using incy.

!  arguments
!  =========

!  n      (input) integer
!         number of elements in input vector(s)

!  sx     (input) real array, dimension(n)
!         single precision vector with n elements

!  incx   (input) integer
!          storage spacing between elements of sx

!  sy     (input) real array, dimension(n)
!         single precision vector with n elements

!  incy   (input) integer
!         storage spacing between elements of sy

!  qsdot  (output) real(kind=16)
!         qsdot  double precision dot product (zero if n.le.0)

!  further details
!  ===============

!  references

!  c. l. lawson, r. j. hanson, d. r. kincaid and f. t.
!  krogh, basic linear algebra subprograms for fortran
!  usage, algorithm no. 539, transactions on mathematical
!  software 5, 3 (september 1979), pp. 308-323.

!  revision history  (yymmdd)

!  791001  date written
!  890831  modified array declarations.  (wrb)
!  890831  revision date from version 3.2
!  891214  prologue converted to version 4.0 format.  (bab)
!  920310  corrected definition of lx in description.  (wrb)
!  920501  reformatted the references section.  (wrb)
!  070118  reformat to lapack style (jl)

!  =====================================================================

!     .. local scalars ..
integer :: i,kx,ky,ns
!     ..
!     .. intrinsic functions ..
intrinsic dble
!     ..
qsdot = 0.0_16
if (n <= 0) return
if (incx == incy .and. incx > 0) then
  
!     code for equal, positive, non-unit increments.
  
  ns = n*incx
  do i = 1,ns,incx
    qsdot = qsdot + real(sx(i),kind=16)*real(sy(i),kind=16)
  end do
else
  
!     code for unequal or nonpositive increments.
  
  kx = 1
  ky = 1
  if (incx < 0) kx = 1 + (1-n)*incx
  if (incy < 0) ky = 1 + (1-n)*incy
  do i = 1,n
    qsdot = qsdot + real(sx(kx),kind=16)*real(sy(ky),kind=16)
    kx = kx + incx
    ky = ky + incy
  end do
end if
return
end function qsdot
