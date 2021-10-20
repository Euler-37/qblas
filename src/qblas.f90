module qblas
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, qblas!"
  end subroutine say_hello
end module qblas
