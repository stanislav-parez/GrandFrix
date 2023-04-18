module mod_init
    use mycommons
    use mod_generate, only: addgrain2, gauss, rand1
    use mod_inputReader
    implicit none
    
    private
    
    public :: INITBOUND

    contains
! ******************************************************************
! ** ROUTINES TO INITIATE VELOCITY AND POSITION VECTORS, **
! ******************************************************************
!=======================================================================
subroutine INITBOUND(grainarea)

    integer   i
    integer nt
    integer dflag
    
    real(8) dist,deviate,MyValue,logg
    real(8) xtogo
    real(8) grainarea
    real(8) rbb
    real(8) rxg,ryg,rad0
    
    xright_out  = boxx*.5 + w_offset
    xleft_out   = -xright_out
    ytop_out    = boxy*.5 + w_offset
    ybot_out    = -ytop_out

    ! BOUNDARY GRAINS ARE RAND0M-SIZE SINGLE ATOMS            

    ! for simplicity, make them doubly-periodic
    !  fill in from left to right and top to bottom
    print*, 'initialize boundaries'

    if (ib(1) >= 0) then
        !     ** FORM THE BOTTOM BOUNDARY ***
        rxg = xleft_out
        ryg = ybot_out + w_offset
        ! find grain size from set of gaussian distributions
        dist=rand1(randum)
        if (dist <= frac1) then 
            do
                deviate=gauss(randum)
                deviate=std1*deviate
                MyValue=mean1 + deviate 
                dflag=1
                if (.not.(MyValue < logr .or. MyValue > higr)) exit
            enddo
        else 
            do
                deviate=gauss(randum)
                deviate=std2*deviate
                MyValue=mean2 + deviate 
                dflag=2
                if (.not.(MyValue < logr2 .or. MyValue > higr2)) exit
            enddo
        endif
        rbb = sigb + sigb
        rad0=rbb*0.5*MyValue
        call addgrain2(rxg,ryg,rad0,1)
        gtype(n)=dflag
        xtogo=xright_out-2*radius(1)

        do i = 2, maxn
            dist=rand1(randum)
            if (dist <= frac1) then 
                dflag=1
                logg=logr
                do
                    deviate=gauss(randum)
                    deviate=std1*deviate
                    MyValue=mean1 + deviate 
                    if (.not.(MyValue < logr .or. MyValue > higr)) exit
                enddo
            else 
                dflag=2
                logg=logr2
                do
                    deviate=gauss(randum)
                    deviate=std2*deviate
                    MyValue=mean2 + deviate 
                    if (.not.(MyValue < logr2 .or. MyValue > higr2)) exit
                enddo
            endif
            rbb = sigb + mod(i,2)*sigb
            rad0=rbb*0.5*MyValue
            rxg = rx(i-1)+radius(i-1)+rad0
            if (2*rad0 > xtogo) then
                if (xtogo < logg) then
                    radius(i-1)=radius(i-1)+.5*xtogo
                    rx(i-1) = rx(i-1)+.5*xtogo
                    nbound(1)=i-1
                else
                    rad0=xtogo/2.
                    rxg = rx(i-1)+radius(i-1)+rad0
                    nbound(1)=i
                    call addgrain2(rxg,ryg,rad0,1)
                    gtype(n)=dflag
                endif
                exit
            else
                call addgrain2(rxg,ryg,rad0,1)
                gtype(n)=dflag
                xtogo = xright_out-radius(1)-rx(i)-radius(i)
            endif
        
        enddo
    endif

    nt=nbound(1)+1
    
    rxg = xleft_out
    ryg = ytop_out - w_offset
    
    ! find grain size from set of gaussian distributions
    dist=rand1(randum)
    
    if (dist <= frac1) then 
        dflag=1
        logg=logr
        do
            deviate=gauss(randum)
            deviate=std1*deviate
            MyValue=mean1 + deviate 
            if (.not.(MyValue < logr .or. MyValue > higr)) exit
        enddo
    else 
        dflag=2
        logg=logr2
        do
            deviate=gauss(randum)
            deviate=std2*deviate
            MyValue=mean2 + deviate 
            if (.not.(MyValue < logr2 .or. MyValue > higr2)) exit
        enddo
    endif
    rbb = sigb + mod(nt,2)*sigb
    rad0=rbb*0.5*MyValue
    call addgrain2(rxg,ryg,rad0,2)
    gtype(n)=dflag
    xtogo=xright_out-2*radius(nt)

    DO i = nt+1, maxn
        dist=rand1(randum)
        if (dist <= frac1) then 
            do
                deviate=gauss(randum)
                deviate=std1*deviate
                MyValue=mean1 + deviate 
                dflag=1
                if (.not.(MyValue < logr .or. MyValue > higr)) exit
            enddo
        else
            do
                deviate=gauss(randum)
                deviate=std2*deviate
                MyValue=mean2 + deviate 
                dflag=2
                if (.not.(MyValue < logr2 .or. MyValue > higr2)) exit
            enddo
        endif
        rbb = sigb + mod(i,2)*sigb
        rad0=rbb*0.5*MyValue
        rxg = rx(i-1)+radius(i-1)+rad0
        if (2*rad0 > xtogo) then
            if (xtogo < logg) then
                radius(i-1)=radius(i-1)+.5*xtogo
                rx(i-1) = rx(i-1)+.5*xtogo
                nbound(2)=i-1
            else
                rad0=xtogo/2.
                rxg = rx(i-1)+radius(i-1)+rad0
                nbound(2)=i
                call addgrain2(rxg,ryg,rad0,2)
                gtype(n)=dflag
            endif
            exit
        else
            call addgrain2(rxg,ryg,rad0,2)
            gtype(n)=dflag
            xtogo = xright_out-radius(nt)-rx(i)-radius(i)
        endif
    
    enddo

    ! ADD VERTICAL WALLS
    if (ib(3) >= 0) then

        ! LEFT WALL
        nt=nbound(2)+1
        rxg = xleft_out + w_offset
        ryg = ybot_out

        ! find grain size from set of gaussian distributions
        dist=rand1(randum)
        if (dist <= frac1) then
            dflag=1
            logg=logr
            do
                deviate=gauss(randum)
                deviate=std1*deviate
                MyValue=mean1 + deviate
                if (.not.(MyValue < logr .or. MyValue > higr)) exit
            enddo
        else
            dflag=2
            logg=logr2
            do
                deviate=gauss(randum)
                deviate=std2*deviate
                MyValue=mean2 + deviate
                if (.not.(MyValue < logr2 .or. MyValue > higr2)) exit
            enddo
        endif
        rbb = sigb + mod(nt,2)*sigb
        rad0=rbb*0.5*MyValue
        call addgrain2(rxg,ryg,rad0,3)
        gtype(n)=dflag
        xtogo=ytop_out-2*radius(nt)

        do i = nt+1, maxn
            dist=rand1(randum)
            if (dist <= frac1) then
                do
                    deviate=gauss(randum)
                    deviate=std1*deviate
                    MyValue=mean1 + deviate
                    dflag=1
                    if (.not.(MyValue < logr .or. MyValue > higr)) exit
                enddo
            else
                do
                    deviate=gauss(randum)
                    deviate=std2*deviate
                    MyValue=mean2 + deviate
                    dflag=2
                    if (.not.(MyValue < logr2 .or. MyValue > higr2)) exit
                enddo
            endif
            rbb = sigb + mod(i,2)*sigb
            rad0=rbb*0.5*MyValue
            ryg = ry(i-1)+radius(i-1)+rad0
            if (2*rad0 > xtogo) then
                if (xtogo < logg) then
                    radius(i-1)=radius(i-1)+.5*xtogo
                    ry(i-1) = ry(i-1)+.5*xtogo
                    nbound(3)=i-1
                else
                    rad0=xtogo/2.
                    ryg = ry(i-1)+radius(i-1)+rad0
                    nbound(3)=i
                    call addgrain2(rxg,ryg,rad0,3)
                    gtype(n)=dflag
                endif
                exit
            else
                call addgrain2(rxg,ryg,rad0,3)
                gtype(n)=dflag
                xtogo = ytop_out-radius(nt)-ry(i)-radius(i)
            endif
        enddo

        ! RIGHT WALL
        nt=nbound(3)+1
        rxg = xright_out - w_offset
        ryg = ybot_out

        ! find grain size from set of gaussian distributions
        dist=rand1(randum)
        if (dist <= frac1) then
            dflag=1
            logg=logr
            do
                deviate=gauss(randum)
                deviate=std1*deviate
                MyValue=mean1 + deviate
                if (.not.(MyValue < logr .or. MyValue > higr)) exit
            enddo
        else
            dflag=2
            logg=logr2
            do
                deviate=gauss(randum)
                deviate=std2*deviate
                MyValue=mean2 + deviate
                if (.not.(MyValue < logr2 .or. MyValue > higr2)) exit
            enddo
        endif
        rbb = sigb + mod(nt,2)*sigb
        rad0=rbb*0.5*MyValue
        call addgrain2(rxg,ryg,rad0,4)
        gtype(n)=dflag
        xtogo=ytop_out-2*radius(nt)

        do i = nt+1, maxn
            dist=rand1(randum)
            if (dist <= frac1) then
                do
                    deviate=gauss(randum)
                    deviate=std1*deviate
                    MyValue=mean1 + deviate
                    dflag=1
                    if (.not.(MyValue < logr .or. MyValue > higr)) exit
                enddo
            else
                do
                    deviate=gauss(randum)
                    deviate=std2*deviate
                    MyValue=mean2 + deviate
                    dflag=2
                    if (.not.(MyValue < logr2 .or. MyValue > higr2)) exit
                enddo
            endif
            rbb = sigb + mod(i,2)*sigb
            rad0=rbb*0.5*MyValue
            ryg = ry(i-1)+radius(i-1)+rad0
            if (2*rad0 > xtogo) then
                if (xtogo < logg) then
                    radius(i-1)=radius(i-1)+.5*xtogo
                    ry(i-1) = ry(i-1)+.5*xtogo
                    nbound(4)=i-1
                else
                    rad0=xtogo/2.
                    ryg = ry(i-1)+radius(i-1)+rad0
                    nbound(4)=i
                    call addgrain2(rxg,ryg,rad0,4)
                    gtype(n)=dflag
                endif
                exit
            else
                call addgrain2(rxg,ryg,rad0,4)
                gtype(n)=dflag
                xtogo = ytop_out-radius(nt)-ry(i)-radius(i)
            endif
        enddo

    else    ! (ib(3) < 0)
        nbound(3)=nbound(2)
        nbound(4)=nbound(3)
    endif

    do concurrent (i=1:nbound(4))
        grainarea=grainarea+.5*pi*radius(i)*radius(i)
    enddo

end subroutine

end module
