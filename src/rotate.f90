module mod_rotate
    use mycommons
    use mod_inputReader 
    implicit none

    private
    
    public :: rotatea, rotateb

    contains
!=======================================================================
subroutine rotatea()

    integer i
    
    real(8) dtsq2,dt2,twopi


    DT2 = 0.5*DT  
    DTSQ2 = DT * DT2
    twopi = 2.*pi

    ! do atom rotations
    do i=1,n
        if (gtype(i) /= 0) then 
            if (bdgrain(i) == 0) then 
                rnt(i)=rnt(i)+dt*w(i)+dtsq2*tq(i)*radius(i)/inertiamom(i)
                rnt(i)=mod(rnt(i),twopi)
                w(i) = w(i) + dt2*tq(i)*radius(i)/inertiamom(i)
            else
                rnt(i)=0.
                w(i)=0.
            endif
        endif
    enddo

end subroutine
!=======================================================================
subroutine rotateb()

    integer i
    
    real(8) dt2

    DT2 = 0.5*DT  

    do i=1,n
        if (gtype(i) /= 0) then 
            if (bdgrain(i) == 0) then 
                w(i) = w(i) + dt2*tq(i)*radius(i)/inertiamom(i)
            else
                w(i)=0.
            endif
        endif
    enddo

end subroutine
end module
