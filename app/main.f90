program main
    use qblas
    implicit none
    real(kind=16)::a(3)=1._16
    real(kind=16),external::qdot
    write(*,*)qdot(3, a, 1, a, 1)
end program main
