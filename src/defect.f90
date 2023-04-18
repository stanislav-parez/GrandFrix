module mod_defect
    use mycommons
    use mod_generate, only : rand1
    implicit none

    private
    
    public :: defect

    contains
!=======================================================================
Subroutine defect(percentdefects, sizevar)

    !     introducing the defects into initial hexagonal packing by changing sizes of grains
    integer indexdefect, Ninter, i, numberd
    
    real(8) sizevar
    real(8) percentdefects
    
    !  Nbound(2) - index of last wall grain
    ! N - total number of grains
    ! Ninter - number of interior grains 

    Ninter = N-Nbound(2)

    numberd = Ninter*nint(percentdefects)

    do i = 1,numberd

        indexdefect = nint(rand1(randum)*Ninter)+Nbound(2)

        if (gtype(indexdefect) == 1) then
            radius(indexdefect) = radius(indexdefect)*sizevar
            gtype(indexdefect) = 2
        else
            numberd = numberd+1
        endif

    enddo

end subroutine


end module
