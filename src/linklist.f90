module mod_linklist
    use mycommons
    use mod_inputReader, only: mx, my, ib
    implicit none

    private
    
    public :: links, maps, getcell

    contains
!    *******************************************************************
!    ** CONSTRUCTION OF CELL LINKED-LISTS AND USE IN FORCE ROUTINE.   **
!    **                                                               **
!    ** REFERENCES:                                                   **
!    **                                                               **
!    ** QUENTREC AND BROT, J. COMPUT. PHYS. 13, 430, 1975.            **
!    ** HOCKNEY AND EASTWOOD, COMPUTER SIMULATION USING PARTICLES,    **
!    **    MCGRAW HILL, 1981.                                         **
!    **                                                               **
!    ** ROUTINES SUPPLIED:                                            **
!    **                                                               **
!    ** SUBROUTINE MAPS                                               **
!    **    SETS UP MAP OF CELL STRUCTURE FOR USE IN FORCE             **
!    ** SUBROUTINE LINKS ( RCUT )                                     **
!    **    SETS UP HEAD OF CHAIN ARRAY AND LINKED LIST                **
!    ** SUBROUTINE FORCE ( SIGMA, RCUT, V, W )                        **
!    **    CALCULATES FORCES USING A LINKED LIST                      **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** SUBROUTINE MAPS IS CALLED ONCE AT THE START OF A SIMULATION   **
!    ** TO ESTABLISH CELL NEIGHBOUR IDENTITIES.  AT EACH TIMESTEP,    **
!    ** SUBROUTINE LINKS IS CALLED TO SET UP THE LINKED LIST AND THIS **
!    ** IS IMMEDIATELY USED BY SUBROUTINE FORCE.                      **
!    *******************************************************************


!=======================================================================
subroutine MAPS()

!    *******************************************************************
!    ** ROUTINE TO SET UP A LIST OF NEIGHBOURING CELLS                **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
!    ** INTEGER MAPSIZ             SIZE OF CELL-CELL MAP              **
!    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THIS SUBROUTINE SETS UP A LIST OF THE THIRTEEN NEIGHBOURING   **
!    ** CELLS OF EACH OF THE SMALL CELLS IN THE CENTRAL BOX. THE      **
!    ** EFFECTS OF THE PERIODIC BOUNDARY CONDITIONS ARE INCLUDED.     **
!    ** THE SUBROUTINE IS CALLED ONCE AT THE BEGINNING OF THE         **
!    ** SIMULATION AND THE MAP IS USED IN THE FORCE SUBROUTINE        **
!    *******************************************************************

    integer IX, IY, IMAP

    !    ** FIND HALF THE NEAREST NEIGHBOURS OF EACH CELL **

    DO IY = 1, my-1
        if (mx == 1) then
            IMAP = ( ICELL ( 1, IY ) - 1 ) * 4
            MAP( IMAP + 1 ) = 0
            MAP( IMAP + 2 ) = 0
            MAP( IMAP + 3 ) = ICELL( 1    , IY + 1)
            MAP( IMAP + 4 ) = 0
            else if (mx == 2) then
            IMAP = ( ICELL ( 1, IY ) - 1 ) * 4
            MAP( IMAP + 1 ) = ICELL( 2, IY    )
            MAP( IMAP + 2 ) = ICELL( 2, IY + 1)
            MAP( IMAP + 3 ) = ICELL( 1    , IY + 1)
            MAP( IMAP + 4 ) = 0

            IMAP = ( ICELL ( 2, IY ) - 1 ) * 4
            MAP( IMAP + 1 ) = 0
            MAP( IMAP + 2 ) = ICELL(1     , IY + 1)
            MAP( IMAP + 3 ) = ICELL( 2    , IY + 1)
            MAP( IMAP + 4 ) = 0
        else

            do IX = 1, mx-1
                IMAP = ( ICELL ( IX, IY ) - 1 ) * 4
                MAP( IMAP + 1 ) = ICELL( IX + 1, IY    )
                MAP( IMAP + 2 ) = ICELL( IX + 1, IY + 1)
                MAP( IMAP + 3 ) = ICELL( IX    , IY + 1)
                MAP( IMAP + 4 ) = ICELL( IX - 1, IY + 1)
            enddo
            
            IMAP = ( ICELL ( mx, IY ) - 1 ) * 4
            MAP( IMAP + 1 ) = ICELL( 1, IY    )
            MAP( IMAP + 2 ) = ICELL( 1, IY + 1)
            MAP( IMAP + 3 ) = ICELL( mx    , IY + 1)
            MAP( IMAP + 4 ) = ICELL( mx - 1, IY + 1)

        endif

    enddo


    !  top boundary
    iy=my
    if (mx == 1) then
        IMAP = ( ICELL ( 1, IY ) - 1 ) * 4
        MAP( IMAP + 1 ) = 0
        MAP( IMAP + 2 ) = 0
        MAP( IMAP + 3 ) = 0
        MAP( IMAP + 4 ) = 0
    else if (mx == 2) then
        IMAP = ( ICELL ( 1, IY ) - 1 ) * 4
        MAP( IMAP + 1 ) = ICELL( 2, IY    )
        MAP( IMAP + 2 ) = 0
        MAP( IMAP + 3 ) = 0
        MAP( IMAP + 4 ) = 0

        IMAP = ( ICELL ( 2, IY ) - 1 ) * 4
        MAP( IMAP + 1 ) = 0
        MAP( IMAP + 2 ) = 0
        MAP( IMAP + 3 ) = 0
        MAP( IMAP + 4 ) = 0
    else
        DO IX = 1, mx-1
            IMAP = ( ICELL ( IX, IY ) - 1 ) * 4
            MAP( IMAP + 1 ) = ICELL( IX + 1, IY    )
            MAP( IMAP + 2 ) = 0
            MAP( IMAP + 3 ) = 0
            MAP( IMAP + 4 ) = 0
        enddo
        
        IMAP = ( ICELL ( mx, IY ) - 1 ) * 4
        MAP( IMAP + 1 ) = ICELL( 1, IY    )
        MAP( IMAP + 2 ) = 0
        MAP( IMAP + 3 ) = 0
        MAP( IMAP + 4 ) = 0
    endif

end subroutine
!=======================================================================
subroutine LINKS 

!    *******************************************************************
!    ** ROUTINE TO SET UP LINKED LIST AND THE HEAD OF CHAIN ARRAYS    **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                  NUMBER OF ATOMS                    **
!    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
!    ** INTEGER NCELL              TOTAL NUMBER OF CELLS (M**3)       **
!    ** INTEGER LIST(NGRAINS+NBOUND)LINKED LIST OF ATOMS               **
!    ** INTEGER HEAD(NCELL)        HEAD OF CHAIN FOR EACH CELL        **
!    ** real(8)    RX(N),RY(N),RZ(N)  POSITIONS                          **
!    ** real(8)    RCUT               THE CUTOFF DISTANCE FOR THE FORCE  **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** EACH GRAIN IS SORTED INTO ONE OF THE MX*MY SMALL CELLS.        **
!    ** THE FIRST GRAIN IN EACH CELL IS PLACED IN THE HEAD ARRAY.      **
!    ** SUBSEQUENT GRAINS ARE PLACED IN THE LINKED LIST ARRAY.         **
!    ** THEN EACH OF THE BOUNDARY ATOMS IS PLACED IN THE LISTS        **
!    ** THE ROUTINE IS CALLED EVERY TIMESTEP BEFORE THE FORCE ROUTINE.**
!    *******************************************************************

    integer i,icell
    
    real(8) dxi,dyi,bx,by

    bx=(xright_out-xleft_out)
    by=(ytop_out-ybot_out)

    dxi = dble(mx)/bx
    dyi = dble(my)/by

    !    ** ZERO HEAD OF CHAIN ARRAY **
    do ICELL = 1, mx*my
        HEAD(ICELL) = 0
    enddo

    !    ** SORT ALL GRAINS **  
    ICELL = getcell(rx(1),ry(1),bdgrain(1),dxi,dyi)
    LIST(1)     = HEAD(ICELL)
    HEAD(ICELL) = 1

    do i = 2, n
        if (bdgrain(i) < 10) then
            ICELL = getcell(rx(i),ry(i),bdgrain(i),dxi,dyi)
            LIST(i)     = HEAD(ICELL)
            HEAD(ICELL) = i
        endif
    enddo


end subroutine
!=======================================================================
function getcell(x,y,ibd,dxi,dyi)

    integer ibd
    integer getcell  
      
    real(8) dxi,dyi,x,y

    if (ib(3) >= 0) then
        getcell = 1 + min(mx-1,INT((x-xleft_out)*dxi)) &
        + min(my-1,INT((y-ybot_out)*dyi)) * mx
    else
        if (ibd == 2) then
            getcell = 1 + INT ( ( x - xleft_out) * dxi ) &
                            + (my-1) * mx
            if (x >= xright_out) getcell = mx*my

        else if (ibd == 4) then
            getcell = mx *(1+ INT ( ( y - ybot_out) * dyi ) )
        else
            getcell = 1 + INT ( ( x - xleft_out) * dxi ) &
                           + INT ( ( y - ybot_out ) * dyi ) * mx
        endif
    endif

end function

!=======================================================================
    !    ** STATEMENT FUNCTION TO GIVE CELL INDEX **
integer function icell(ix , iy)
    integer,intent(in) :: ix, iy
    
    ICELL = 1 + MOD( IX - 1 + mx, mx ) + MOD( IY - 1 + my, my)* mx
end function

end module
