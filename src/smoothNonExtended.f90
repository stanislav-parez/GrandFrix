module mod_smoothNonExtended
    use mycommons 
    use mod_inputReader
    use mod_outputWriter
    implicit none

    private
    
    real(8), public, protected :: vc(MAXD+1,MAXD+1)
    real(8), public, protected :: siterd2(MAXD+1, MAXD+1)

    public :: smooth

    contains
!cccccccccccccccc 15/2/09 cccccccccccccccccccccccccccccccc


!      Here we assign values to grid points representing some continous properties of the fluid+grain system.
!       the values are assign using the halo function which is a 2D linear interpolation fuction which allows
!       calculating the relative weight of each particle in the system to the grid point arounf it.
!
!       We use staggared grid  so that 
!       P, k and phi are given at grid point (i,j)
!       ux is between x grid points (i+0.5,j)
!       uy is between y grid points (i,j+0.5)
!       
! 
!       vc -  is a 2D array of the area that each grain contributes to the grid point. On grid
!       siterd2 -  is a 2D array of the average squre of radiuses. On grid. 
!                  This is not a conservative quantity because of division by numpar 
!       sitevx and sitevy - are 2D arrays of grains velocity in each grid point. Staggared
!       divvel - is a 2D array of the divergence of the velocity. On grid
!      solfra - is a 2D array of the bulk density of system grains area/total area. On grid
!                Not a conserved quantity because of boundary effect.
!      perm - is a 2D array of the permeability. On grid
!      phimat - is a 2D array of the porosity. On grid
!       numpar - is a 2D array of the number of particles in each grid point, used for averaging radius^2. On grid 
!       numvx - is a 2D array of the number of particles in each vx grid point, used for averaging
!       numvy - is a 2D array of the number of particles in each vy grid point, used for averaging

!     Note that there are two more lines of grid point for sitevy and numvy at the top and at the bottom
!     and one more line of grid points for phimat,perm
!
!     This file should go together with updatePressure4
!     It was contracted from smoothExtended by removing the code of 
!     the top and bottom extra grains layers.
!=======================================================================
subroutine smooth()

    real(8),parameter :: rsmin=0.25d0
    real(8),parameter :: phiint=0.01d0
    real(8),parameter :: maxsiterd2 =4d0
    !,permint
    !minimum solfra is rsmin because below that the Carman-Kozeny permeability law fails.

    integer ix,iy,k
    integer ixvxh,iyvyh
    
    real(8) gs
    real(8) dxr,dyr
!    real(8) xb1,yb1
    real(8) xl
!    real(8) xbvx1,ybvy1
    real(8) dxrvx,dyrvy
    real(8) numvx(MAXD+1, MAXD+1)
    real(8) numvy(MAXD+1, MAXD+1)
    

    dxc=(xright-xleft)/mx
    dyc=(ytop-ybot)/my

    ! zeroing out the arrays
    ! my - # of cells in the y direction
    ! my+1 # of grid points lines inside the boc
    ! my+2 # number of grid point lines for y velocity

    do ix=1,mx+1
        do iy=1,my+2
            vc(ix,iy)=0.
            siterd2(ix,iy)=0.
            sitevx(ix,iy)=0.
            sitevy(ix,iy)=0.
            divvel(ix,iy)=0.
            solfra(ix,iy)=0.
            !perm(ix,iy) = 1.
            phimat(ix,iy) = 0.
            numpar(ix,iy)=0.
            numvx(ix,iy)=0.
            numvy(ix,iy)=0.
        enddo
    enddo
    !print*, 'XX s 1: vy = ',sitevy(1:mx,my+1)

    ! ix,iy for regular grid
    ! ixvxh for staggared vx
    ! iyvyh for staggared vy
    do k = 1,n
        if (rx(k) < xleft .or. rx(k) > xright &
          .or. ry(k) < (ybot-ftiny) .or. ry(k) > (ytop+ftiny)) then
            print*, '^^^ k',k
            print*,'*** rx, ry',rx(k),ry(k)
            cycle
        endif
        gs=4d0/3d0*pi*radius(k)*radius(k)*radius(k)
        if (bdgrain(k) > 0) gs=0.5d0*gs
        ix = int((rx(k) - xleft)/dxc) + 1
        if (ix > mx) ix = mx    !XX added
        if (rx(k) < (xleft + 0.5d0*dxc)) then
            ixvxh = mx
        else
            ixvxh = int((rx(k)-(xleft + 0.5d0*dxc))/dxc)+1
        endif
        iy = int((ry(k) - ybot)/dyc) + 1
        if (iy > my) iy = my    !XX added
        if (ry(k) < (ybot + 0.5d0*dyc)) then
            iyvyh = 1
        else
            iyvyh=int((ry(k)-(ybot + 0.5d0*dyc))/dyc)+2
        endif
        if (ixvxh <1 .or. ixvxh >mx .or. iyvyh <1 .or. iyvyh >my+1 ) then
            print*, '^^^',k
            print*,'***',rx(k),ix,ixvxh,dxc
            print*,'***',ry(k),iy,iyvyh,dyc
            print*,'^^^',ybot,ytop
        endif
        if (ix <1 .or. ix > mx .or. &
        iy < 1 .or. iy > my+1) then
            print*,'***',k
            print*,'***',rx(k),ix,dxc
            print*,'***',ry(k),iy,dyc
            print*,'---',ybot,ytop
        endif

        ! xb1, yb1 - coordinates of grid
        ! xbvx1,ybvy1 - coordinate of the relevant staggared 
!        xb1=xleft+(ix-1)*dxc
!        yb1=ybot+(iy-1)*dyc
!        if (ixvxh == mx) then
!            xbvx1 = xright - 0.5d0*dxc
!        else
!            xbvx1 = xleft+(ixvxh-0.5d0)*dxc
!        endif
!        ybvy1 = ybot + (iyvyh - 1.5d0)*dyc
!        if(iyvyh == 1) ybvy1 = ybot - 0.5*dyc
!        if(iyvyh == my+2) ybvy1 = ytop + 0.5*dyc

        ! scaled distance of grain from grid point
        
!        dxr=(rx(k)-xb1)/dxc
        dxr=(rx(k)-xleft)/dxc - (ix-1)
        if (dxr < 0.) print*, 'dxr: ',dxr
        if (dxr > 1.) print*, 'dxr: ',dxr   ! may happen when a grain is at x=xright because ix is still mx
!        dyr=(ry(k)-yb1)/dyc
        dyr=(ry(k)-ybot)/dyc - (iy-1)
        if (dyr > 1.d0+ftiny) print*, 'dyr: ',dyr, k
        if (dyr < 0.d0) print*, 'dyr: ',dyr, k
!        dxrvx = (rx(k) - xbvx1)/dxc
!        if (rx(k) < xleft + 0.5*dxc) then
!            dxrvx = (rx(k) - xleft + 0.5*dxc)/dxc
!        endif
        if (rx(k) < xleft + 0.5d0*dxc) then
           dxrvx = (rx(k) - xleft)/dxc + 0.5d0
        else
           dxrvx = (rx(k) - xleft)/dxc - (ixvxh - 0.5d0)
        endif
!        dyrvy = (ry(k) - ybvy1)/dyc
        dyrvy = (ry(k) - ybot)/dyc - (iyvyh - 1.5d0)

!     debugging
         if (dxrvx > 1. .or. dxrvx < 0.) then
            print*, 'dxrvx is : ',dxrvx, k
            print*, rx(k),xleft,dxc
         endif
        if (dyrvy > 1. .or. dyrvy < 0. ) then
            print*, 'dyrvy is : ',dyrvy, k
            print*, ry(k),ybot,dyc
        endif
        
        !          starting to calculate relevant quantities
        !          botom left
        vc(ix,iy)=vc(ix,iy)+((1.-dyr)*(1.-dxr)*gs)
        siterd2(ix,iy)=siterd2(ix,iy)+ &
        (1.-dyr)*(1.-dxr)*(radius(k)**2)
        
        if (gtype(k) /= 0) then
            numpar(ix,iy)=numpar(ix,iy)+(1.-dyr)*(1.-dxr)
        endif
        
        sitevx(ixvxh,iy)=sitevx(ixvxh,iy)+ &
        vx(k)*(1.-dyr)*(1.-dxrvx)
        numvx(ixvxh,iy)=numvx(ixvxh,iy)+ &
        (1.-dyr)*(1.-dxrvx)
        sitevy(ix,iyvyh)=sitevy(ix,iyvyh)+ &
        vy(k)*(1.-dyrvy)*(1.-dxr)
        !if(ix==22 .and. iyvyh==my+1) print*, 'XX s 2: k,vy = ',k,vy(k),dyrvy,dxr
        numvy(ix,iyvyh)=numvy(ix,iyvyh)+ &
        (1.-dyrvy)*(1.-dxr)

        !         top left
        vc(ix,iy+1)=vc(ix,iy+1)+(dyr*(1.-dxr)*gs)
        siterd2(ix,iy+1) = siterd2(ix,iy+1) +  &
        dyr*(1.-dxr)*(radius(k)**2)
        
        if (gtype(k) /= 0) then
            numpar(ix,iy+1) =numpar(ix,iy+1)+dyr*(1.-dxr) 
        endif
        
        sitevx(ixvxh,iy+1) = sitevx(ixvxh,iy+1) +  &
        vx(k)*dyr*(1.-dxrvx)
        numvx(ixvxh,iy+1) = numvx(ixvxh,iy+1) +  &
        dyr*(1.-dxrvx)

        sitevy(ix,iyvyh+1) = sitevy(ix,iyvyh+1) +  &
        vy(k)*dyrvy*(1.-dxr)
        !if(ix==22 .and. iyvyh==my) print*, 'XX s 2: k,vy = ',k,vy(k),dyrvy,dxr
        numvy(ix,iyvyh+1) = numvy(ix,iyvyh+1) +  &
        dyrvy*(1.-dxr)

        !          bottom right           
        if (ix < mx) then
            vc(ix+1,iy)=vc(ix+1,iy)+((1.-dyr)*dxr*gs)
            siterd2(ix+1,iy) = siterd2(ix+1,iy) + &
            (1.-dyr)*dxr*(radius(k)**2)
            
            if (gtype(k) /= 0) then
                numpar(ix+1,iy) =numpar(ix+1,iy) + &
                (1.-dyr)*dxr
            endif
            
            sitevy(ix+1,iyvyh) = sitevy(ix+1,iyvyh) + &
            vy(k)*(1.-dyrvy)*dxr
            !if(ix==21 .and. iyvyh==my+1) print*, 'XX s 2: k,vy = ',k,vy(k),dyrvy,dxr
            numvy(ix+1,iyvyh) = numvy(ix+1,iyvyh) + &
            (1.-dyrvy)*dxr
            !          top right
            vc(ix+1,iy+1)=vc(ix+1,iy+1)+(dyr*dxr*gs) 
            siterd2(ix+1,iy+1) = siterd2(ix+1,iy+1) + &
            dyr*dxr*(radius(k)**2)
            
            if (gtype(k) /= 0) then
                numpar(ix+1,iy+1)= &
                numpar(ix+1,iy+1)+dyr*dxr
            endif

            sitevy(ix+1,iyvyh+1) = sitevy(ix+1,iyvyh+1)+  &
            vy(k)*dyrvy*dxr
            !if(ix==21 .and. iyvyh==my) print*, 'XX s 2: k,vy = ',k,vy(k),dyrvy,dxr
            numvy(ix+1,iyvyh+1) = numvy(ix+1,iyvyh+1)+  &
            dyrvy*dxr              
            !C if ix = mx
        else
            !          bottom right
            vc(1,iy)=vc(1,iy)+((1.-dyr)*dxr*gs)
            siterd2(1,iy) = siterd2(1,iy) + &
            (1.-dyr)*dxr*(radius(k)**2)
            
            if (gtype(k) /= 0) then
                numpar(1,iy) =numpar(1,iy) + (1.-dyr)*dxr
            endif
            
            sitevy(1,iyvyh) = sitevy(1,iyvyh) + &
            vy(k)*(1.-dyrvy)*dxr
            numvy(1,iyvyh)=numvy(1,iyvyh)+(1.-dyrvy)*dxr
            !           top right
            vc(1,iy+1)=vc(1,iy+1)+(dyr*dxr*gs) 
            siterd2(1,iy+1) = siterd2(1,iy+1) +  &
            dyr*dxr*(radius(k)**2)
            
            if (gtype(k) /= 0) then
                numpar(1,iy+1)=numpar(1,iy+1)+dyr*dxr 
            endif
            
            sitevy(1,iyvyh+1) = sitevy(1,iyvyh+1) +  &
            vy(k)*dyrvy*dxr
            numvy(1,iyvyh+1)=numvy(1,iyvyh+1)+dyrvy*dxr
        endif

        if (ixvxh < mx) then
            sitevx(ixvxh+1,iy) = sitevx(ixvxh+1,iy) + &
            vx(k)*(1.-dyr)*dxrvx
            numvx(ixvxh+1,iy) = numvx(ixvxh+1,iy) + &
            (1.-dyr)*dxrvx
            sitevx(ixvxh+1,iy+1) = sitevx(ixvxh+1,iy+1)+  &
            vx(k)*dyr*dxrvx
            numvx(ixvxh+1,iy+1) = numvx(ixvxh+1,iy+1) +  &
            dyr*dxrvx
        else
            sitevx(1,iy)=sitevx(1,iy)+ &
            vx(k)*(1.-dyr)*dxrvx
            numvx(1,iy) = numvx(1,iy) + (1.-dyr)*dxrvx
            sitevx(1,iy+1)=sitevx(1,iy+1)+vx(k)*dyr*dxrvx
            numvx(1,iy+1) = numvx(1,iy+1) + dyr*dxrvx
        endif

    enddo
    !print*, 'XX s 2: vy = ',sitevy(1:mx,my+1)

!    do iy=1,my+1
!        vc(mx+1,iy)=vc(1,iy)
!    enddo


    do iy=1,my+1
        do ix =1,mx
            solfra(ix,iy) = vc(ix,iy)/(dxc*dyc)
            if(intruder/=1)then
                if (solfra(ix,iy) > 1.) &
                print*, 'solfra(',ix,iy,')=', &
                solfra(ix,iy), &
                vc(ix,iy),dxc,dyc,my
            endif
            if (solfra(ix,iy) < rsmin) then
                solfra(ix,iy) = rsmin
            endif
        enddo
        solfra(mx+1,iy) = solfra(1,iy) !This is because solfra(mx+1,iy) is needed in force.f90
    enddo

    do iy=1,my+1
        do ix=1,mx
            phimat(ix,iy)=1.- solfra(ix,iy)
            if (phimat(ix,iy) < 0.)then
                phimat(ix,iy) = 0.

            endif
        enddo
    enddo

!    do iy=1,my+2
!        phimat(mx+1,iy)=phimat(1,iy)
!        siterd2(mx+1,iy)=siterd2(1,iy)
!        numpar(mx+1,iy)=numpar(1,iy)
!        sitevy(mx+1,iy) = sitevy(1,iy)
!        numvy(mx+1,iy) = numvy(1,iy)
!    enddo


    ! the permeability is calculated in this way because rho_s{3D} = (2/3)rho_s{2D}
    ! the coefficient which enters the scaling of fluid equation is alpha*d^2
    ! where alpha = 1/45/12

    do ix=1,mx
        do iy=1,my+1
            if (numpar(ix,iy) /= 0.)then
                siterd2(ix,iy)=siterd2(ix,iy)/ &
                numpar(ix,iy)
            else
                siterd2(ix,iy) = 0.
            endif
        enddo
    enddo

    if (const_perm < 0)then
        do ix=1,mx
            do iy=1,my+1

                if (siterd2(ix,iy) > 0.) then
                    perm(ix,iy) = permfac*siterd2(ix,iy)*phimat(ix,iy)**3/ &
                    solfra(ix,iy)**2
                else
                    perm(ix,iy) = permfac*maxsiterd2*phimat(ix,iy)**3/ &
                    solfra(ix,iy)**2
                endif

            enddo
        enddo
    endif


    do ix=1,mx
        do iy=1,my+1
            if (numvx(ix,iy) < 0.) then
                print*,'numvx is: ',numvx(ix,iy), ix, iy
            endif
            
            if (numvx(ix,iy) > 0.)then
                sitevx(ix,iy) = sitevx(ix,iy)/numvx(ix,iy)
            else
                sitevx(ix,iy) = 0.
            endif
        enddo
    enddo
    
    do ix=1,mx
        do iy=1,my+2
            if (numvy(ix,iy) < 0.) then
                print*,'numvy is: ',numvy(ix,iy), ix, iy
            endif
            if (numvy(ix,iy) > 0.)then
                sitevy(ix,iy) = sitevy(ix,iy)/numvy(ix,iy)
            else
                sitevy(ix,iy) = 0.
            endif
        enddo
    enddo

    if(intruder==1)then       
        xl=xright-xleft 
        do ix=1,mx
            do iy=1,my+1
!                xb1=xleft+(ix-1)*dxc
!                yb1=ybot+(iy-1)*dyc
!                dxr=rx(n)-xb1
                dxr=rx(n)-xleft-(ix-1)*dxc
!                dyr=ry(n)-yb1
                dyr=ry(n)-ybot-(iy-1)*dyc
                if (ib(3) < 0) then
                    if(abs(dxr)>(0.5d0*xl)) dxr=dxr-sign(xl,dxr)
                endif
                
                if((dxr**2+dyr**2)<rintruder**2)then
                    phimat(ix,iy)=phiint
                    solfra(ix,iy)=1-phiint
                    perm(ix,iy) = permfac*maxsiterd2*phimat(ix,iy)**3/ &
                    solfra(ix,iy)**2
                endif    
                        
!                ybvy1=ybot+(iy-1.5)*dyc
!                dyrvy=ry(n)-ybvy1
                dyrvy=ry(n)-ybot-(iy-1.5)*dyc
                
                if((dxr**2+dyrvy**2)<(rintruder+1d0)**2)then
                    sitevy(ix,iy)=vy(n)
                endif
                
!                xbvx1=xleft+(ix-0.5)*dxc
!                dxrvx=rx(n)-xbvx1
                dxrvx=rx(n)-xleft-(ix-0.5)*dxc
                if (ib(3) < 0) then
                    if(abs(dxrvx)>(0.5d0*xl)) dxrvx=dxrvx-sign(xl,dxrvx)
                endif
                
                if((dxrvx**2+dyr**2)<(rintruder+1d0)**2)then
                    sitevx(ix,iy)=vx(n)
                endif
                
            enddo
        enddo
!        do iy=1,my+2
!            phimat(mx+1,iy)=phimat(1,iy)
!            perm(mx+1,iy)=perm(1,iy)
!            sitevy(mx+1,iy) = sitevy(1,iy)
!        enddo
    endif

    !       * CALCULATING THE DIVERGENCE OF THE VELOCITY *
    !       INNER GRID

    do ix=2,mx
        do iy=1,my+1
            divvel(ix,iy) = (sitevx(ix,iy) - &
            sitevx(ix-1,iy))/(dxc) +  &
            (sitevy(ix,iy+1) - sitevy(ix,iy))/(dyc)
        enddo
    enddo
    !       RIGHT AND LEFT

    do iy=1,my+1
        divvel(1,iy) = (sitevx(1,iy) - &
        sitevx(mx,iy))/(dxc) +  &
        (sitevy(1,iy+1) - sitevy(1,iy))/(dyc)
!        divvel(mx+1,iy) = divvel(1,iy)
    enddo



!    do ix = 2,mx
!        divvel(ix,1) = (sitevx(ix,1) - &
!        sitevx(ix-1,1))/(dxc) +  &
!        (sitevy(ix,2) - sitevy(ix,1)) /(dyc)
!        divvel(ix,my+1)=(sitevx(ix,my+1) - &
!        sitevx(ix-1,my+1))/(dxc) +  &
!        (sitevy(ix,my+2) - sitevy(ix,my+1))/(dyc)
!    enddo

    !       CORNERS
!    divvel(1,1) = (sitevx(1,1) - sitevx(mx,1))/(dxc) + &
!    (sitevy(1,2) - sitevy(1,1))/(dyc)
!    divvel(mx+1,1)=divvel(1,1)
!    divvel(1,my+1)=(sitevx(1,my+1) &
!    - sitevx(mx,my+1))/(dxc) +  &
!    (sitevy(1,my+2) - sitevy(1,my+1))/(dyc)
!    divvel(mx+1,my+1) = divvel(1,my+1)

    if(intruder==1)then       
        xl=xright-xleft 
        do ix=1,mx
            do iy=1,my+1
                !xb1=xleft+(ix-1)*dxc
                !yb1=ybot+(iy-1)*dyc
                !dxr=rx(n)-xb1
                dxr=rx(n)-xleft-(ix-1)*dxc
                !dyr=ry(n)-yb1
                dyr=ry(n)-ybot-(iy-1)*dyc
                if (ib(3) < 0) then
                    if(abs(dxr)>(0.5d0*xl)) dxr=dxr-sign(xl,dxr)
                endif
                
                if((dxr**2+dyr**2)<rintruder**2)then
                    if(abs(divvel(ix,iy))>ftiny)then
                        print*,'WARNING: divvel inside the intruder not zero'
                    endif
                endif            
            enddo
        enddo
    endif

    do iy = 1,my+1
        do ix = 1,mx
            permave(ix,iy)=permave(ix,iy)+perm(ix,iy)
            solfraave(ix,iy)=solfraave(ix,iy)+solfra(ix,iy)
            numparave(ix,iy)=numparave(ix,iy)+numpar(ix,iy)
            sitevxave(ix,iy)=sitevxave(ix,iy)+sitevx(ix,iy)
            sitevyave(ix,iy)=sitevyave(ix,iy)+sitevy(ix,iy)
            divvelave(ix,iy)=divvelave(ix,iy)+divvel(ix,iy)
        enddo
    enddo

end subroutine
end module
