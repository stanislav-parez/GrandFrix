module mod_findforces
    use mycommons 
    use mod_inputReader
    implicit none

    private
    
    public :: getneighbors

    contains
!=======================================================================
    subroutine getneighbors(distint)

    integer ICELL, JCELL0, JCELL, L, NABOR
    integer ncell
    integer ixcell,iycell,ixcell2,ipf,iycell2
    integer i,k
    integer kindex(maxk),kk,kindex2,kntold,kmin,kmax,kspan
    
    real(8) temp(maxk),distint
    real(8) ryij,rxij,xl,rmaxsq

    ncell = mx*my

    ! save the shear forces on the contacts
    do k = 1,contactknt
        kindex(k)=(contacti(k)-1)*n+contactj(k)
        temp(k)=contft(k)
        !         print*, contacti(k),contactj(k)
    enddo
    
    kntold = contactknt

    !  zero out the neighbor lists
    do k = 1,contactknt
        contacti(k)=0
        contactj(k)=0
        contft(k)=0.    
    enddo

    ! **LOOP OVER ALL GRAINS IN THE CELL **
    icell = 0
    contactknt = 0
    
    do iycell = 1,my
        do ixcell = 1,mx
        
            icell = icell+1 
            K = HEAD(ICELL)
            do while (K > 0)

                !  ** LOOP OVER ALL GRAINS BELOW K IN THE CURRENT CELL **
                L = LIST(K)
                do while (L > 0)
                    !if ( (bdgrain(k)/=bdgrain(l)) .or. (bdgrain(k)*bdgrain(l)==0) ) then
                    if (bdgrain(k)*bdgrain(l)==0) then
                        if (gtype(k)*gtype(l) /= 0) then
                            call setneighbor(k,l,0,distint)
                        endif
                    endif
                    L = LIST(L)
                enddo

                !  ** LOOP OVER NEIGHBOURING CELLS **
                JCELL0 = 4 * (ICELL - 1)
                do NABOR = 1, 4
                    JCELL = MAP ( JCELL0 + NABOR )
                    ipf=0
                    if (mx == 2) then 
                        ipf=2
                    else
                        if (ixcell == 1) then 
                            iycell2=(jcell-1)/mx
                            ixcell2=jcell-(iycell2)*mx
                            if (ixcell2 == mx) ipf=1
                        else if (ixcell == mx) then
                            iycell2=(jcell-1)/mx
                            ixcell2=jcell-(iycell2)*mx
                            if (ixcell2 == 1) ipf=1
                        endif           
                    endif 

                    !     **  LOOP OVER ALL GRAINS IN NEIGHBOURING CELLS **
                    if (JCELL == 0) cycle
                    
                    L = HEAD(JCELL)
                    do while ( L /= 0 )

                        !if ( (bdgrain(k)/=bdgrain(l)) .or. &
                        if ((bdgrain(k)*bdgrain(l)==0) .and. &
                            (gtype(k)*gtype(l) /= 0)) then
                                call setneighbor(k,l,ipf,distint)
                        endif
                        
                        L = LIST(L)
                
                    enddo

                enddo

                K = LIST(K)

            enddo
        enddo
    enddo

    if(intruder==1)then
        xl=xright_out-xleft_out
        rmaxsq=(rintruder+1)**2
        do i=1,n-1
            ryij  = ry(n) - ry(i)
            rxij  = rx(n) - rx(i)
            if(abs(rxij)>(0.5d0*xl))then
                rxij=rxij-sign(xl,rxij)
            endif
            if((rxij**2+ryij**2)<rmaxsq)then
                contactknt=contactknt+1
                contacti(contactknt)=i
                contactj(contactknt)=n
            endif                        
        enddo        
    endif

    !  now reset the shear force array
    !    do k=1,contactknt
    !     contft(k)=temp(contacti(k),contactj(k))
    !    enddo

    kspan=3
outer:do k=1,contactknt
        ! index hasn't changed
        kindex2=(contacti(k)-1)*n+contactj(k)
        if (kindex2 == kindex(k)) then
            contft(k)=temp(k)
            cycle outer
        endif
        ! index is nearby in the list
        kmin=max(1,k-kspan)
        kmax=min(kntold,k+kspan)
        do kk=kmin,kmax
            if (kindex2 == kindex(kk)) then
                contft(k)=temp(kk)
                cycle outer
            endif
        enddo
        ! contact has moved a lot
        do kk=kmin-1,1,-1
            if (kindex2 == kindex(kk)) then
                contft(k)=temp(kk)
                cycle outer
            endif
        enddo
        
        do kk=kmax+1,kntold
            if (kindex2 == kindex(kk)) then
                contft(k)=temp(kk)
                cycle outer
            endif
        enddo
    
    enddo outer

end subroutine

!=======================================================================
subroutine setneighbor(i,j,ipf,distint)
    !**************************************************************
    !   check the neighbor list, add contact if particles are
    !   within distint of each other

    integer i,imx,imn,j
    integer ipf
    
    real(8) rijsq,distint,rxij2
    real(8) xl,ryij,rxij,rmaxsq,rmax

    if((intruder==1).and.((i==n).or.(j==n))) return

    xl=xright_out-xleft_out

    imn=min(i,j)
    imx=max(i,j)

    rmax = radius(imn)+radius(imx)+distint
    rmaxsq=rmax*rmax

    ryij  = ry(imx) - ry(imn)
    rxij  = rx(imx) - rx(imn)
    
    if (ipf == 1) then
        rxij = rxij-sign(xl,rxij)
    else if (ipf == 2) then
        rxij2 = rxij-sign(xl,rxij)
        if (abs(rxij2) < abs(rxij)) rxij=rxij2
    endif
    
    rijsq = (rxij*rxij + ryij*ryij)

    ! particles are interacting
    if (rijsq < rmaxsq) then
        contactknt=contactknt+1
        contacti(contactknt)=imn
        contactj(contactknt)=imx
    endif              

end subroutine


end module
