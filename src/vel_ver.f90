module mod_vel_ver
    use mycommons 
    use mod_inputReader
    implicit none

    private
    
    public :: MOVEA, MOVEB

    contains
!**************************************************************
!*  VELOCITY VERSION OF VERLET ALGORITHM                           **
!**************************************************************

!    *******************************************************************
!    ** TWO ROUTINES THAT TOGETHER IMPLEMENT VELOCITY VERLET METHOD.  **
!    **                                                               **
!    ** REFERENCE:                                                    **
!    **                                                               **
!    ** SWOPE ET AL., J. CHEM. PHYS. 76, 637, 1982.                   **
!    **                                                               **
!    ** ROUTINES SUPPLIED:                                            **
!    **                                                               **
!    ** SUBROUTINE MOVEA ( DT, N, NTOT )                                    **
!    **    MOVES POSITIONS AND PARTIALLY UPDATES VELOCITIES.          **
!    ** SUBROUTINE MOVEB ( DT, N, NTOT)                                **
!    **    COMPLETES VELOCITY MOVE AND CALCULATES KINETIC ENERGY.     **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                   NUMBER OF MOLECULES               **
!    ** real(8)    DT                  TIMESTEP                          **
!    ** real(8)    M                   ATOMIC MASS                       **
!    ** real(8)    RX(N),RY(N)   POSITIONS                         **
!    ** real(8)    VX(N),VY(N)   VELOCITIES                        **
!    ** real(8)    FX(N),FY(N)   FORCES                            **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** AT THE START OF A TIMESTEP, MOVEA IS CALLED TO ADVANCE THE    **
!    ** POSITIONS AND 'HALF-ADVANCE' THE VELOCITIES.  THEN THE FORCE  **
!    ** ROUTINE IS CALLED, AND THIS IS FOLLOWED BY MOVEB WHICH        **
!    ** COMPLETES THE ADVANCEMENT OF VELOCITIES.                      **
!    *******************************************************************


!=======================================================================
SUBROUTINE MOVEA ( FYTOP, FYBOT, fxbot,topmass)

!    *******************************************************************
!    ** FIRST PART OF VELOCITY VERLET ALGORITHM                       **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THE FIRST PART OF THE ALGORITHM IS A TAYLOR SERIES WHICH      **
!    ** ADVANCES POSITIONS FROM T TO T + DT AND VELOCITIES FROM       **
!    ** T TO T + DT/2.  AFTER THIS, THE FORCE ROUTINE IS CALLED.      **
!    *******************************************************************

    integer I
    
    real(8) DT2, DTSQ2
    real(8) topmass,tmass1
    real(8) fytop,fybot,fxbot,rms,xl,yb

    dt2   = 0.5*dt  
    dtsq2 = dt * dt2

    tmass1=1./topmass

    !  bottom boundary
    if (ib(1) == 0) then
        do concurrent (i=1:nbound(1))
            vx(i)=0.
            vy(i)=0.
        enddo
    else if (ib(1) == 1) then
        !          applied vertical force, free horizontal
        do concurrent (i=1:nbound(1))
            rx(i) = rx(i) + dt * vx(i) + dtsq2*fx(i)
            ry(i) = ry(i) + dt * vy(i) + dtsq2*fybot
            vx(i) = vx(i) + dt2* fx(i)
            vy(i) = vy(i) + dt2* fybot
        enddo
    else if (ib(1) == 2) then
        !          applied vertical force, applied horizontal velocity
        do concurrent (i=1:nbound(1))
            rx(i) = rx(i) + dt * vx(i)
            ry(i) = ry(i) + dt * vy(i) + dtsq2*fybot
            vx(i) = fb(1)
            vy(i) = vy(i) + dt2* fybot
        enddo
    else if (ib(1) == 3) then
        !          applied vertical velocity, free horizontal 
        do concurrent (i=1:nbound(1))
            rx(i) = rx(i) + dt * vx(i) + dtsq2*fx(i)
            ry(i) = ry(i) + dt * vy(i)
            vx(i) = vx(i) + dt2* fx(i)
            vy(i) = fb(1)
        enddo
    else if (ib(1) == 4) then
        !          zero vertical velocity, applied horizontal velocity
        do concurrent (i=1:nbound(1))
            rx(i) = rx(i) + dt * vx(i)
            vx(I) = fb(1)
            vy(I) = 0.
        enddo
    else if (ib(1) == 5) then
        !          applied vertical velocity, zero horizontal velocity
        do concurrent (i=1:nbound(1))
            ry(i) = ry(i) + dt * vy(i)
            vy(I) = fb(1)
            vx(I) = 0.
        enddo
    else if ((ib(1) == 6).or.(ib(1) == 7)) then
        !          zero vertical velocity, applied horizontal force
        do concurrent (i=1:nbound(1))
            rx(i) = rx(i) + dt * vx(i) + dtsq2*fxbot
            vx(i) = vx(i) + dt2* fxbot
            vy(i) = 0.
        enddo

    else if (ib(1) == 8) then
        !          prescribed oscilations
        do concurrent (i=1:nbound(1))
            rx(i) = rx(i)+dt * vsin
            vx(I) = vsin
            vy(I) = 0.
        enddo
    endif

    !  top boundary
    if (ib(2) == 0.) then
        do concurrent (i=nbound(1)+1:nbound(2))
            vx(i)=0.
            vy(i)=0.
        enddo
    else if (ib(2) == 1 .or. ib(2) == 5) then
        if (ib(2) == 5) then
          dwall=dwall+dt*vx(nbound(1)+1)+dtsq2*fx(nbound(1)+1)*tmass1
          dspring = dspring + dt*uspring
          !print*, 'XX v 1: dwall,dspring = ',dwall,dspring
        endif
        do concurrent (i = nbound(1)+1:nbound(2))
            rx(i) = rx(i) + dt * vx(i) + dtsq2*fx(i)*tmass1
            ry(i) = ry(i) + dt * vy(i) + dtsq2*fytop*tmass1
            vx(i) = vx(i) + dt2* fx(i)*tmass1
            vy(i) = vy(i) + dt2* fytop*tmass1
        enddo
        !print*, 'XX v 2: rx,ry,vx,vy = ',rx(nbound(1)+2),ry(nbound(1)+2),vx(nbound(1)+2),vy(nbound(1)+2)
    else if (ib(2) == 2) then
        do concurrent (i = nbound(1)+1:nbound(2))
            rx(i) = rx(i) + dt * vx(i)
            ry(i) = ry(i) + dt * vy(i) + dtsq2*fytop*tmass1
            vx(i) = fb(2)
            vy(i) = vy(i) + dt2* fytop*tmass1
        enddo
        !print*, 'XX v 2: rx,ry,vx,vy = ',rx(nbound(1)+2),ry(nbound(1)+2),vx(nbound(1)+2),vy(nbound(1)+2)
    else if (ib(2) == 3) then
        do concurrent (i = nbound(1)+1:nbound(2))
            rx(i) = rx(i) + dt * vx(i) + dtsq2*fx(i)*tmass1
            ry(i) = ry(i) + dt * vy(i)
            vx(i) = vx(i) + dt2* fx(i)*tmass1
            vy(i) = fb(2)
        enddo
    else if (ib(2) == 4) then
        do concurrent (i = nbound(1)+1:nbound(2))
            rx(i) = rx(i) + dt * vx(i)
            vx(i) = fb(2)
        enddo
    endif
          !print*, 'XX vel on top particle=',vy(nbound(1)+1)
          !print*, 'XX vel on bottom particle=',vy(nbound(1))
    ! interior grains
    do concurrent (i = nbound(4)+1:n)
        rms=radinv2(i)
        rx(i) = rx(i) + dt * vx(i) + dtsq2 * fx(i) *rms
        ry(i) = ry(i) + dt * vy(i) + dtsq2 * fy(i) *rms
        vx(i) = vx(i) + dt2 * fx(i) *rms
        vy(i) = vy(i) + dt2 * fy(i) *rms
    enddo

    ! UPDATE DOMAIN BOUNDARIES
    !ybot=ry(nbound(1))
    !ytop=ry(nbound(2))
    !print*, 'XX v 3-: ybot,ytop = ',ybot,ytop
     ybot=ry(2)
     ytop=ry(nbound(1)+2)
    !print*, 'XX v 3+: ybot,ytop = ',ybot,ytop
    if (ib(3) >= 0) then
        ! assuming static bottom wall
        !xleft=min(rx(1),rx(nbound(3)))
        xleft=rx(nbound(3))
        if (xleft < xleft_out) then
            print*, "ERROR: left wall flies away!"
            stop
        endif
        xright=rx(nbound(4))
        if (xright > xright_out) then
            print*, "ERROR: right wall flies away!"
            stop
        endif
        if (ybot < ybot_out) then
            print*, "ERROR: bottom wall flies away!"
            stop
        endif
        if (ytop > ytop_out) then
            print*, "ERROR: top wall flies away!"
            stop
        endif
    else
        ybot_out=ybot
        ytop_out=ytop
    endif

    !  ugly bit to take care of grain wraparound
    !  in case of periodic boundaries
    xl=xright_out-xleft_out
    do concurrent (i = 1:n)
        if (rx(i) >= xright_out) then
            rx(i)=rx(i)-xl
        else if (rx(i) < xleft_out) then
            rx(i)=rx(i)+xl
        endif
    enddo
    if (ib(3) >= 0) then
        yb=ytop_out-ybot_out
        do concurrent (i = 1:n)
            if (ry(i) >= ytop_out) then
                ry(i)=ry(i)-yb
            else if (ry(i) < ybot_out) then
                ry(i)=ry(i)+yb
            endif
        enddo
    endif
    
end subroutine
!=======================================================================
SUBROUTINE MOVEB ( FYTOP, FYBOT, fxbot,topmass)

!    *******************************************************************
!    ** SECOND PART OF VELOCITY VERLET ALGORITHM                      **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** THE SECOND PART OF THE ALGORITHM ADVANCES VELOCITIES FROM     **
!    ** T + DT/2 TO T + DT. THIS ASSUMES THAT FORCES HAVE BEEN        **
!    ** COMPUTED IN THE FORCE ROUTINE AND STORED IN FX, FY, FZ.       **
!    *******************************************************************

    integer I
    
    real(8) DT2
    real(8) fytop,fybot,rms,fxbot
    real(8) topmass,tmass1

    dt2 = 0.5*dt  

    tmass1=1./topmass

    !  bottom boundary
    if (ib(1) == 1) then
        do concurrent (i=1:nbound(1))
            vx(i) = vx(i) + dt2* fx(i)
            vy(i) = vy(i) + dt2* fybot
        enddo
    else if (ib(1) == 2) then
        do concurrent (i=1:nbound(1))
            vx(i) = fb(1)
            vy(i) = vy(i) + dt2* fybot
        enddo
    else if (ib(1) == 3) then
        do concurrent (i=1:nbound(1))
            vx(i) = vx(i) + dt2* fx(i)
            vy(i) = fb(1)
        enddo
    else if (ib(1) == 4) then
        do concurrent (i=1:nbound(1))
            vx(I) = fb(1)
        enddo
    else if (ib(1) == 5) then
        !          applied vertical velocity, zero horizontal velocity
        do concurrent (i=1:nbound(1))
            vy(I) = fb(1)
        enddo
    else if ((ib(1) == 6).or.(ib(1) == 7)) then
        !          zero vertical velocity, applied horizontal force
        do concurrent (i=1:nbound(1))
            vx(i) = vx(i) + dt2* fxbot
            vy(i) = 0.
        enddo
    else if (ib(1) == 8) then
        do concurrent (i=1:nbound(1))
            vx(I) = vsin
            vy(i)=0d0
        enddo
    endif

    !  top boundary
    if (ib(2) == 1 .or. ib(2) == 5) then
        do concurrent (i=nbound(1)+1:nbound(2))
            vx(i) = vx(i) + dt2* fx(i)*tmass1
            vy(i) = vy(i) + dt2* fytop*tmass1
        enddo
    else if (ib(2) == 2) then
        do concurrent (i=nbound(1)+1:nbound(2))
            vx(i) = fb(2)
            vy(i) = vy(i) + dt2* fytop*tmass1
        enddo
    else if (ib(2) == 3) then
        do concurrent (i=nbound(1)+1:nbound(2))
            vx(i) = vx(i) + dt2* fx(i)*tmass1
            vy(i) = fb(2)
        enddo
    else if (ib(2) == 4) then
        do concurrent (i=nbound(1)+1:nbound(2))
            vx(I) = fb(2)
        enddo
    endif

    ! interior grains
    do concurrent (i=nbound(4)+1:n)
        rms=radinv2(i)
        vx(i) = vx(i) + dt2 * fx(i) *rms
        vy(i) = vy(i) + dt2 * fy(i) *rms
    enddo

end subroutine
end module
