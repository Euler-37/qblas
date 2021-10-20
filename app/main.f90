program main
    use qblas
    implicit none
    block
        real(kind=16)::a(3)
        complex(kind=16)::b(3)
        real(kind=16),external::qdot
        complex(kind=16),external::odotc,odotu
        call random_number(a)
        call random_number(b%im)
        call random_number(b%re)
        write(*,*)"qblas:qdot ",qdot(3, a, 1, a, 1)
        write(*,*)"exact      ",dot_product(a, a)
        write(*,*)"qblas:odotc",odotc(3, b, 1, b, 1)
        write(*,*)"exact      ",dot_product(b, b)
        write(*,*)"qblas:odotu",odotu(3, b, 1, b, 1)
        write(*,*)"exact      ",dot_product(conjg(b), b)
    end block
end program main
