!> Reads the length (size) of the fields from the
!> Coulomb calculation
subroutine get_length(n)

    implicit none
    integer :: n
!f2py intent(out) :: n
    logical :: a, b, c

    open(2, file='pe_fock_1.bin', status='old', form='unformatted')
    read(2) a, b, c
    read(2) n
    close(2)

    return

end subroutine get_length

!> does the actual readin of the field to pass it onto the python program
subroutine get_field(n, field)

    implicit none
    integer :: n
    real*8, dimension(3*n) :: field
!f2py intent(in) n
!f2py intent(out) field

    integer :: ntmp
    logical :: a, b, c

    field = 0.0d0

    open(2, file='pe_fock_1.bin', status='old', form='unformatted')
    read(2) a, b, c
    read(2) ntmp

    if (n .ne. ntmp) then
        write(*,'(A)') "ERROR: dimensions do not match"
        goto 99
    end if
    read(2) field

 99 close(2)

    return

end subroutine get_field
