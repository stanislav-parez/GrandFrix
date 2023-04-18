PROGRAM GranLayer
    use mycommons
    use mod_breaks
    use mod_force
    use mod_defect
    use mod_findforces
    use mod_generate, only:shearstep, partgen, strength, rand1
    use mod_init
    use mod_linklist
    use mod_rotate
    use mod_rw_conf
    use mod_vel_ver
    use mod_updatepressure
    use mod_smoothNonExtended
    use mod_principalstress, only : PRINCIPALSTRESS
    use mod_postScriptTool
    use mod_inputReader
    use mod_outputWriter
    IMPLICIT NONE

    !    *******************************************************************
    !    ** FORTRAN GRANULAR DYANMICS PROGRAM*         
    !    **  
    !    **   BUILDS AND/OR SHEARS A LAYER OF CIRCULAR GRAINS
    !    **     (NOTE NOT SET UP FOR MULTI-PARTICLE GRAINS)
    !    **     PERIODIC IN THE HORIZONTAL DIRECTION
    !    **     VARIOUS CONDITIONS ON THE TOP AND BOTTOM WALLS
    !    **     C    *******************************************************************
    !    ** PRINCIPAL VARIABLES:                   
    !    **                                             
    !    **  N            NO. OF GRAINS
    !    **  RX,RY      GRAIN POSITIONS    
    !    **  VX,VY      GRAIN VELOCITIES
    !    **  FX,FY      GRAIN FORCES
    !    **  RADIUS GRAIN RADII


    !    ** BOUNDARIES ARE LISTED IN THE ORDER: BOTTOM,TOP,LEFT,RIGHT
    !    **  NBOUND(4)      NO. OF ATOMS IN EACH BOUNDARY
    !    **  IB(4)      TYPE OF EACH BOUNDARY

    !    **  HEAD      HEAD OF CHAIN FOR EACH CELL
    !    **  LIST      LINKED LIST OF GRAINS
    !    **  MAP            LIST OF NEIGHBORING CELLS (ONLY UPPER-RIGHT HALF)
    !    **  MX             TOTAL NUMBER OF CELLS IN X-DIRECTION
    !    **  MY             TOTAL NUMBER OF CELLS IN Y-DIRECTION

    !    **  SIG            AVERAGE ATOM DIAMETER (SET TO 1.)

    !    ** **INPUT PARAMETERS**
    !    **  BOXX, BOXY  SIZE OF BOX IN X AND Y (RELATIVE TO AVERAGE
    !    **          ATOM DIAMETER)
    !    **  SIGB      AVERAGE BOUNDARY ATOM DIAMETER
    !    **  SKSHEAR     SPRING CONSTANT FOR LEAF SPRING
    !    **  GAMMA       DAMPING PARAMETER FOR NORMAL VELOCITY
    !    **  GAMSHEAR    DAMPING PARAMETER FOR TANGENTIAL VELOCITY
    !    **  FRICTION   COEFFICIENT OF FRICTION


    !  DWS 11/9 
    !  ^^^^non-dimensionalization^^^^
    ! choosing values for :
    !  k = stiffness between particles 
    !  m = mass of standard particle
    !  length scale , x0, is one standard atom diameter

    ! the following scales apply:
    !  velocity scale,v0 = x0 * sqrt(k/m)
    !  time scale, t0 = x0/v0
    !  resulting force scale, f0 = k*x0
    !  resulting energy scale, e0 = k*x0*x0

    !  non-dim equations become
    !    f = abs(overlap) - gamma*(velocity difference) 
    !                     - rho*(unit vertical vector)
    !
    !    x(t+dt) = x(t) + v(t)*dt  + .5*dt*dt*f(t)
    !
    !    v(t+dt) = v(t) + .5*dt*(f(t) + f(t+dt))
    !
    !   where
    !   damping parameter gamma = gamma0/sqrt(k*m)
    !                      = gamma0/grainradius/sqrt(density*pi*k)
    !   buoyancy parameter rho  = g*m/(k*x0)
    !                           = g *(grainradius)**2/(acoustic velocity)**2/x0
    !
    !   g=acceleration of gravity
    !   gamma0 = some inherent damping(kg-m/s)
    !
    !      note: rho (input as gx and gy, will be small (order 10-3 or less),
    !            gamma will be of order 1 for well-damped systems
    !       

    !  critical time step = grain radius/acoustic velocity (from Mora and Place)
    !                     ~ r/( sqrt(9/8 * k/m) *r)
    !
    !  critical time step scaled by t0 = sqrt(8/9)
    !  Mora and Place recommend a time step of 0.2*critical, or about .15
    !

    !  DWS 11/9 

    INTEGER(8) :: STEP    
    integer i,k
    integer j
    integer icell
    integer ipfc
    integer ixcell,iycell
    integer nbcheck
    integer(8) :: step0
    integer icrush,nevent,ibreak,inext
    integer ibtotal
    integer indint,nintn,nshort,nintshort
    
    real(8) FYTOP, FYBOT
    real(8) mvx, mvy
    real(8) grainarea,ev,ek,ekr
    real(8) grainarea0
    real(8) fxbot,fxtop
    real(8) xl,boxvol
    real(8) angmom,yl
    real(8) lavevx(maxd),lavevy(maxd),layervx(maxd),layervy(maxd)
    real(8) px(maxd),py(maxd),lavepx(maxd),lavepy(maxd)
    real(8) lavew(maxd),layerw(maxd)
    real(8) heightave,height
    real(8) ldens(maxd),lavedens(maxd)
    real(8) avesxxt(maxd),avesxyt(maxd),avesyxt(maxd),avesyyt(maxd)
    real(8) aveskxxt(maxd),aveskxyt(maxd),aveskyxt(maxd),aveskyyt(maxd)
    real(8) lavetempx(maxd),lavetempy(maxd),layertempx(maxd),layertempy(maxd)
    real(8) layerpor(maxd),lavepor(maxd),laverho(maxd),lrho(maxd)
    real(8) presyfinal
    real(8) vm,vmax
    real(8) dtsave
    real(8) presy0
    real(8) rsmall,randnm
    real(8) distint,tresume
    real(8) wmax,topg
    real(8) averradius2
    real(8) curpresy
    real(8) fgrain
    real(8) histdy
    real(8) tmp
!    real(8) :: phigoal0 = 0.
    real(8) :: start_time = 0.d0
    real(8) :: masstot = 0d0
    real(8) :: alap1=0d0
    
    character(len=4) fnn

    call cpu_time(start_time)
    
    call random_seed()

    call read_input_file()
    
    call postScriptTool_init(jpeg_resolution)
  
    call zero_variables()

    presy0 = presy
    !  for building a new set
    !   start at confining pressure of 2.0e-5, increase pressure gradually
    presyfinal = 0.
    if (phigoal > 0.) then
        presyfinal=presy
        presy=1.0d-5
    endif

    !  adjust number of cells to account for final box size and grain sizes

    ! input damping parameters are scaled by gamma
    !        gammamol=gamma*gammamol
    !        gamshear=gamma*gamshear
  
    if (rand_seed == 0) then
        call random_number(tmp)
        randum = int(10000000*tmp)
    else
        randum = rand_seed
    endif

    ! &&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT
    !  INITIALIZATION SECTION: BUILD OR READ IN GRAINS

    ! irstart=0, build all the grains
    if (irstart == 0) then
        !if (ifill == 0) ib(2)=-1
        tau=0.
        grainarea=0d0
        call initbound(grainarea)
        if (ifill < 0) then
            ! if gfac > 1. --> sedimentation test with sparse distribution of grains.
            ! for sedimentation ratio of boxy to boxfy is about 2:3
            if (gfac > 1.) then
                boxvol=boxx*boxy
                call partgen(phigoal,boxvol,grainarea,alap1,gfac)
                phigoal = 0.0
            else
                boxvol=boxx*boxfy
                call partgen(phigoal,boxvol,grainarea,alap1,1.d0)
            endif
            print*, 'finished generating ',n,' grains'
            call dumpgrains(afrac1)
        endif

        ! initializing the fluid preesure 2D array (liran)

        if (ifluid == 0) then

            print*, 'Initializing pressure to:" ', presy*fstress  
            pout =  0.
            ! 0.379 is the ratio between fluid density and grains density 
            print*, 'Pressure outside the domain is: ',pout
            do i=1,MAXD+1
                do j=1,MAXD+1
                    press(i,j) = presy*fstress
                    gradpx(i,j) = 0.
                    gradpy(i,j) = 0.
                    phimat(i,j) = 0.
                    perm(i,j) = 0.
                    sitevx(i,j)= 0.
                    sitevy(i,j)=0.
                    divvel(i,j)=0.
                    solfra(i,j)=1.
                    numpar(i,j)=1.
                enddo
            enddo

            call read_fbcinput(fbc,irstart)

        endif
    else
        ! irstart>1, read a restart file
        if (ibw > 0) then 
            print*, 'IBW>0, WILL ONLY MAKE ONE PLOT OF '
            print*, title
            call readrestart(title,ifluid)
        else
            !print*, 'XX m 8-: ytop = ',ytop
            call readrestart(title,ifluid)
            !print*, 'XX m 8+: ytop = ',ytop
!            print*, 'XX m 8+: perm = ',perm(1:mx,my)
!            print*, 'XX m 8+: coef2 = ',phimat(1:mx,my)
!            print*, 'XX m 8+: vy = ',sitevy(1:mx,my+1)
!            print*, 'XX m 8+: rx = ',rx(1:n)
!            print*, 'XX m 8+: ry = ',ry(1:n)
!            print*, 'XX m 8+: radius = ',radius(1:n)
!            print*, 'XX m 8+: vx = ',vx(1:n)
!            print*, 'XX m 8+: vy = ',vy(1:n)
        endif

        if (ifluid == 0) call read_fbcinput(fbc,irstart)
        
        if (idefect) call defect(percentdefects, sizevar)

        call dumpgrains(afrac1)

        !  adjust number of cells to account for final box size and grain sizes
!        bg=0.
!        do i=1,n
!            bg=max(bg,radius(i))
!        enddo
!        bg=bg*2.
        !boxxg=(xright-xleft)/bg
        !     ******************************************************************
        !boxx=boxx*(1-presy)

        ! maximum number of cells allowed
        !mxmax=int((xright-xleft)/bg)-1
        !mymax=max(1,int((ytop-ybot)/bg))

        print*, 'read restart file', title

        write(*,*) 'tau',tau

        !xl=xright-xleft

        ! zero out tau at beginning of run
        if (iztau == 0) tau=0.

    endif

    do  ! LOOP FOR DOUBLE RUN WHEN IRSTART=0; ENDS BEFORE "WRITE(*,'(/1X,''**** END OF DYNAMICS **** ''//)')"

    ! AFTER GENERATING OR READING INITIAL CONFIGURATION UPDATE DOMAIN BOUNDARIES
    ybot=ry(nbound(1))
    ytop=ry(nbound(2))
    if (ib(3) >= 0) then
        xleft=rx(nbound(3))
        xright=rx(nbound(4))
    else
        xleft=xleft_out
        xright=xright_out
        ybot_out=ybot
        ytop_out=ytop
    endif


    if(ifluid==0) then
        if(wl<0) then 
            wl=ytop
        else
            wl=wl+ybot
        endif
        print*, 'water level', wl
        dyc=(ytop-ybot)/my
        mywl=nint((wl-ybot)/dyc)
        if(mywl>my) mywl=my
        print*, 'mywl', mywl
        
        do iycell=mywl+1,my+2
            do ixcell=1,mx+1
                press(ixcell,iycell)=0d0
            enddo
        enddo
        
        if(mywl<my) then
            if(fbc==2 .or. fbc == 4 .or. fbc == 6) then
                print*, 'Neumann bc is not compatible with partially saturated systems'
                stop 
            endif
        endif
        
        !  re-initialize the fluid pressure array (liran)
        if (irstart == 2) then
            print*, 'Initializing pressure to:" ', presy*fstress   
            pout =  0.
            ! 0.379 is the ratio between fluid density and grains density 
            print*, 'Pressure outside the domain is: ',pout
            
            press = presy*fstress
            gradpx = 0.
            gradpy = 0.
            phimat = 0.
            perm = 0.
            sitevx= 0.
            sitevy=0.
            divvel=0.
            solfra=1.
            numpar=1.
            
            call read_fbcinput(fbc,irstart)
        endif

    endif
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! AFTER GENERATING SYSTEM OR READING RSTART FILE:
    ! CHECK EXISTANCE AND READ FILE WITH INFORMATION ABOUT GRAINS TO BE REMOVED

    nshort = n

    !  set up factor for getting 1/mass of a particle (particle of
    !    radius=0.5 has mass=1
    do i=1,n
        if (gtype(i) /= 0) then
            if((intruder==1).and.(i==n))then
                radinv2(i)=1d0/(6d0*radius(i)*radius(i)*densratio)
                inertiamom(i)=0.5d0/radinv2(i)*radius(i)*radius(i)
            else
                radinv2(i)=.125/(radius(i)*radius(i)*radius(i))
                inertiamom(i)=0.4d0/radinv2(i)*radius(i)*radius(i)
            endif
            mass2D(i)=4d0*radius(i)*radius(i)
            masstot=masstot+1d0/radinv2(i)
        else
            radinv2(i) = 0
        endif
    enddo

    call calculate_area_and_porosity()
    print*, 'initial porosity = ', 1. - grainarea/xl/yl

    ! find smallest grain for calculating time step
    rsmall=radius(1)
    do i=2,n
        if (gtype(i) /= 0) then
            if (radius(i) < rsmall) then
            rsmall= radius(i)
            endif
        endif
    enddo

    ! Calculate dt
    print*, "smallest grain's diameter =",2d0*rsmall
    dt=dt*((2*rsmall)**1.5)    ! perhaps due to T~sqrt(m/k)

    WRITE(*,*) '*********************************'
    WRITE(*, *) 'NUMBER OF ATOMS = ',N
    WRITE(*, *)' NUMBER OF BOUNDARY ATOMS:' 
    WRITE(*, *)'  BOTTOM: ',NBOUND(1)
    WRITE(*, *)'     TOP: ',NBOUND(2)-NBOUND(1)
    if (ib(3) >= 0) then
        WRITE(*, *)'    LEFT: ',NBOUND(3)-NBOUND(2)
        WRITE(*, *)'   RIGHT: ',NBOUND(4)-NBOUND(3)
    endif
    WRITE(*, *)' BOTTOM TOP LEFT RIGHT = ',ybot,ytop,xleft,xright
    if (ib(3) >= 0) then
        WRITE(*, *)' BOX BTLR = ',ybot_out,ytop_out,xleft_out,xright_out
    endif
    WRITE(*, *)' NO. OF CELLS IN LINKLIST =', MX,    ' x ', MY
    WRITE(*, *)' TIME STEP        = ', dt
    WRITE(*, *)'OUTPUT FREQUENCY = ',  IPRINT, iprintf
    WRITE(*, *)' '
    WRITE(*, *)'CONTACT MODEL = ', contact_model
    WRITE(*, *)'NORMAL STIFFNESS = ', skmol
    WRITE(*, *)'TANGENTIAL STIFFNESS = ', skshear
    WRITE(*, *)'NORMAL DAMPING COEF= ', gamma
    WRITE(*, *)'FRICTION = ', friction
    WRITE(*, *)' EXTERNAL LOADING = ', presy
    if (ib(2) == 5) then
        WRITE(*, *)' SPRING: kspring, dspring, dwall =', kspring,dspring,dwall
    endif
    WRITE(*, *)' BODY FORCES = ', gx,gy
    WRITE(*, *)'FSTRESS =', fstress
    WRITE(*, *)'TOPMASS =', topmass
    WRITE(*, *)'FLUID MODE =', ifluid
    if (ifluid == 0) then
        WRITE(*, *)'FBC: ', fbc
        if (const_perm > 0d0) then
            WRITE(*, *)'perm = ',const_perm
        else
            WRITE(*, *)'PERMFAC = ',permfac
        endif
        WRITE(*, *)'viscosity = ',visco
        WRITE(*, *)'acoef = ',acoef
        WRITE(*, *)'bcoef = ',bcoef
        WRITE(*, *)'Pressure advection = ',yn_advect
        WRITE(*, *)'Explicit scheme = ',yn_explicit
    endif
    if(intruder==1)then
    WRITE(*, *)'INTRUDER, RINTRUDER = ', intruder,rintruder
    WRITE(*, *)'DENSRATIO = ',densratio
    endif

    !   index of first interior atom
    indint=nbound(4)+1
    !   number of interior atoms
    nintn=n-nbound(4)        
    nintshort = nintn

    ! END OF BUILDING AND READING
    ! &&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT&&&&&INIT

    !    ** CALCULATE INITIAL VALUES, NEED TO FORM INITIAL  LIST **
    !    ** AND GET FORCES AT T=0 **


    !      ** INITIALIZE THE MAP OF NEIGHBOR CELLS
    !print*, 'XX m 5-: ytop = ',ytop
    call maps()
    call links()
    !print*, 'XX m 5+: ytop = ',ytop
    
    ! &&&&&BOUNDARY&&&&BOUNDARY&&&&&BOUNDARY&&&&&BOUNDARY&&&&&BOUNDARY&&&&&
    ! &&  SET UP THE BOUNDARY AND INITIAL CONDITIONS    

    !  zero out the initial momentum, for a dead restart
    if (izmom == 0) then
        do i= 1,n
            vx(i)=0.
            vy(i)=0.
            w(i)=0.
        enddo

    else if (izmom == 1) then
        do i= nbound(1),n
            vx(i)=fb(2)
            vy(i)=0.
            w(i)=0.
        enddo
    endif

    !  zero out the slip info, to avoid crazy springs during
    !   no friction relaxation
    if (izslip == 0) then
        do k=1,contactknt
            contft(k)=0.
        enddo
    endif
    ! if freezing(or resetting) boundaries, also need to freeze velocities,
    !  otherwise, there will be slight movement during first 
    !  time step at the previous velocity
    if (ib(1)==3) then
        do concurrent (i=1:nbound(1))
            vy(i)=fb(1)
        enddo
    else if (ib(1) == 4) then
        do concurrent (i=1:nbound(1))
            vx(i)=fb(1)
        enddo
    endif
    
    if (ib(2)==3) then
        do concurrent (i=nbound(1)+1:nbound(2))
            vy(i)=fb(2)
        enddo
    else if (ib(2) == 4) then
        do concurrent (i=nbound(1)+1:nbound(2))
            vx(i)=fb(2)
        enddo
    endif

    !  top wall pulled by shear spring
    !   initialize the distance moved by the wall, dwall
    !   and the pulling spring, dspring 
    !dtsq2 = 0.5*dt * dt
    if (ib(2) == 5) then
!        if ((irstart == 0).or.(iztau == 0)) then
!            dwall=0.
!            dspring=0.
!        endif
        uspring=fb(2)
    endif

    !  treat the top row as a rigid block of some weight:
    !   applied boundary force (presy)/gravity (gy) gives the
    !   height of the overlying mass in particles. Give the
    !   boundary particles effectively the mass of the overlying
    !  column
!    if (itmass) then
!        topmass=0.00463/kspring
!    else
!        topmass=1.
!    endif

    ! &&&&&BOUNDARY&&&&BOUNDARY&&&&&BOUNDARY&&&&&BOUNDARY&&&&&BOUNDARY&&&&&

    if (const_perm > 0) then
        do ixcell=1,mx+1
            do iycell=1,my+2
                perm(ixcell,iycell) = const_perm
            enddo
        enddo
    endif

    ! 11/30/01 - forces are now read in from restart file
    ! GET INITIAL FORCES AND PLOT
    print*, 'calling first force'
    call links
    distint=0.1
    print*, 'number of contacts =',contactknt
    call getneighbors(distint)
    
    if (ifluid == 0) then
        !print*, 'XX m 2-: gradpy = ',gradpy(:,my+1)
        call smooth()
        !print*, 'XX m 2+: gradpy = ',gradpy(:,my+1)
!        print*, 'XX m 2+: perm = ',perm(1:mx,my)
!        print*, 'XX m 2+: coef2 = ',phimat(1:mx,my)
!        print*, 'XX m 2+: vy = ',sitevy(1:mx,my+1)
        do i = 1,TRATIO
            !print*, 'XX m 3-: gradpy = ',gradpy(:,my+1)
            call updatepressure(fbc,i)
            !print*, 'XX m 3+: gradpy = ',gradpy(:,my+1)
        enddo
    endif

    curpresy = presy

    !print*, 'XX m 1-: ytop = ',ytop
    call force( gx, gy, ev, presx, curpresy, presy0, fytop, fybot, &
    fxtop,fxbot,1,crush1,ifluid,fgrain,fbc)
    !print*, 'XX m 1+: ytop = ',ytop
    
    call write_to('monitor.out', tau,maxvn,sumvn/max(ncontacts,1),maxcomp, &
        sumcomp/max(ncontacts,1),refi,refj,ncontacts,max_flux_x,max_flux_y)


    ! &&&SPECIAL RUN MODES &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ! RUN PROGRAM IN SINGLE PLOTTING MODE (ibw>0, just do first plot and quit)
    if (ibw == 1) then
        fnn=title(8:11)
        print*, fnn
        call pscriptplot(next,igif,fnn,0)
        !     call bwplot(fnn)       !alternate plotting routine
        print*, 'MADE COLOR PLOT'
        do k=1,contactknt
            i=contacti(k)
            j=contactj(k)
            if (contfn(k) < 0.0)  call write_to('breakevents.out', i,j,contfn(k), &
                                                contft(k),contang(k))
        enddo

        stop
    else if (ibw == 2) then

        if (fbc /= 3) then
            if(ib(1) == 8) then
                !rsin = amp*(1-cos(omega*(tau+dt)))
                !vsin = amp*omega*cos(omega*(tau+0.5d0*dt))
                vsin = amp*omega*sin(omega*(tau+0.5d0*dt))
            endif
            CALL MOVEA ( fytop, fybot,fxbot,topmass)

            CALL ROTATEA
        endif


        if (ifluid == 0) then
            call smooth()
            do i = 1,TRATIO
                call updatepressure(fbc,i)
            enddo
        endif

        CALL FORCE ( gx, gy, ev,presx,presy,presy0,fytop,fybot,   &
        fxtop,fxbot,ipfc,1,ifluid,fgrain,fbc)

        call principalstress()
        call pscriptplot2(next,igif,fnn)
        print*, 'MADE COLOR PLOT2'
        stop

    else 
        !print*, 'XX m 10-: ytop = ',ytop
        call pscriptplot(next,igif,fnn,0)
        !print*, 'XX m 10+: ytop = ',ytop

    endif

    !  DO A SHEAR STEP
    !    shift the x position of all particles by some amount
    if (gstep > 0.) then 
        call shearstep(gstep)
        call pscriptplot(next,igif,fnn,0)
    endif
    

    ! &&&SPECIAL RUN MODES &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&        
    nbcheck=20

    !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
    ! FOR BREAKING GRAINS RUNS
    ! If using grains without strength assigned,
    ! or you want to change the strength of the grains 
    ! change the switch restrg to 1 
    
    if (restrg) then
        randnm=rand1(randum)
        do i=1,n
            call strength(i)
        enddo
    endif

    ! Initialize the counters used for breakage events 
    icrush=crush1
    nevent=0

    !   ZERO OUT FLUID RELATED FIELDS THAT ARE PRINTED SO THAT THEY DO NOT CONTAIN
    !   GARBAGGE FROM PREVIOUS CALLS OF UPDATEPRESSURE OR SMOOTH
    call zero_fluid_outputs()
    !   THE SAME TO REMOVE GARBAGGE IN SXXT,... FROM PREVIOUS CALLS OF FORCE
    call zero_variables()

    !    *******MAIN*******MAIN*******MAIN*******MAIN*******MAIN*******MAIN
    !    ** MAIN LOOP BEGINS                                            **
    !    *******************************************************************

    WRITE(*,'(//1X,''**** START OF DYNAMICS ****'')')

        step0=0
        DO STEP = 1, NSTEP

             if (yn_explicit) then
                 if (ifluid == 0) then
                    call smooth
                    do i = 1,TRATIO
                       call updatepressure(fbc,i)
                    enddo
                 endif
             endif

            !  GET NEW POSITION AND VELOCITY AND ROTATION
            if (fbc /=3) then
                if(ib(1) == 8) then
                    !vsin = amp*omega*cos(omega*(tau+0.5d0*dt))
                    !rsin = amp*(1-cos(omega*(tau+dt)))
                    vsin = amp*omega*sin(omega*(tau+0.5d0*dt))
                    !write(98,*) tau, vsin, amp
                endif
                !print*, 'XX m 11-: ytop = ',ytop
                CALL MOVEA ( fytop, fybot,fxbot,topmass)
                !print*, 'XX m 11+: ytop = ',ytop
                !     COMMENT NEXT LINE TO STOP ROTATION
                if (iroll == 0) CALL ROTATEA
            endif

!            if (ib(2) == 5) then
!                dwall = dwall+dt*vx(nbound(1)+1)+dtsq2*fx(nbound(1)+1)
!            endif
            ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
            ! TIME-DEPENDENT BOUNDARY FORCES
            ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

            ! spring-driven run
!            if (ib(2) == 5) then
!                dspring=dspring + dt*uspring
!                fspring=kspring*(dspring-dwall)
!            else
!                fspring=0.
!            endif

            ! Horizonal periodic acceleration
            if (ib(1) == 7) then
                fsin = amp*sin(omega*tau)
            else
                fsin = 0.
            endif

            ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

            !  IF NECESSARY BUILD NEIGHBOR LIST
            !    INCLUDES ALL PARTICLES THAT ARE distint FROM OVERLAP OR LESS.
            !   the maximum number of time steps between rebuilding the list (nbcheck) is
            !   calculated based on the time step and particle velocities. 
            !   however (based on tests with similar sized grains),
            !   once nbcheck gets over 100 it has minimal effect on 
            !   running time, and at 500 the solution begins to change.

            if ( step == step0+nbcheck) then
                step0=step
                !print*, 'XX m 10-: ytop = ',ytop
                call links()
                !print*, 'XX m 10+: ytop = ',ytop
                distint=0.1
                !print*, 'XX m 9-: ytop = ',ytop
                call getneighbors(distint)
                !print*, 'XX m 9+: ytop = ',ytop
                vmax=0.
                do i=1,n
                    vm=vx(i)*vx(i)+vy(i)*vy(i)
                    vmax=max(vmax,vm)
                enddo
                vmax=sqrt(vmax)
                nbcheck=int(distint/vmax/dt/5)
                nbcheck=min(100,nbcheck)
            endif

            !  CALCULATE NEW FORCES
            if ( mod( step, iprint ) == 0 ) then
                ipfc=1
            else
                ipfc=0
            endif

            if (.not.(yn_explicit)) then
                if (ifluid == 0) then
                    call smooth()
                    do i = 1,TRATIO
                        call updatepressure(fbc,i)
                    enddo
                endif
            endif
            
            curpresy = presy
            if (IPULSE) then  
                curpresy = presy + amp*cos(omega*tau)
                if (ipfc==1) then 
                    call write_to('coswave.out', tau, curpresy)
                endif
            endif
            
            !print*, 'XX m 8-: ytop = ',ytop
            CALL FORCE ( gx, gy, ev,presx,curpresy,presy0,fytop,fybot,   &
            fxtop,fxbot,ipfc,crush1,ifluid,fgrain,fbc)
            !print*, 'XX m 8+: ytop = ',ytop
            
            if (ipfc == 1) then
                call write_to('monitor.out', tau,maxvn,sumvn/max(ncontacts,1),maxcomp, &
                    sumcomp/max(ncontacts,1),refi,refj,ncontacts,max_flux_x,max_flux_y)
                if (maxcomp > 1d-2) then
                    print*, 'WARNING: max overlap > 1e-2'
                endif
            endif

            ! add a viscous fluid into the pore space
            !         visc=gammamol
!            do i=1,n
!                fx(i)=fx(i)-gammamol*vx(i)
!                fy(i)=fy(i)-gammamol*vy(i)
!                tq(i)=tq(i)-gammamol*w(i)
!            enddo

            !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
            !     BREAKING GRAINS                             *
            !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
            if (crush1 > 0) then
                if (icrush == 0) then
                    if (tau > tresume) then
                        icrush=crush1
                    endif
                else 
                    dtsave=dt
                    call breaking(rsmall,ibreak,tresume,grainarea,icrush,inext)
                    if (ibreak > 0) then
                        icrush=0
                        ibtotal=ibtotal+1
                        call write_to('breakevents.out', tau,ibtotal)
                        call links
                        distint=0.1
                        call getneighbors(distint)
                        call force(gx, gy, ev, presx, presy, presy0, fytop, fybot, &
                        fxtop,fxbot,ipfc,icrush,ifluid,fgrain,fbc)
                        print*, 'grain break at time ', tau, ' step ',step
                        call pscriptplot(next,.true.,fnn,inext)
                        if (dt < dtsave) then
                            iprint=nint(iprint*dtsave/dt)
                            iprintf=iprint/20
                            print*, 'iprint, iprintf=',iprint,iprintf
                        endif
                    endif
                endif
            endif

            !BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

            !  UPDATE VELOCITIES
            if (fbc /=3) then
                if(ib(1) == 8) then
                    !vsin = amp*omega*cos(omega*(tau+dt))
                    vsin = amp*omega*sin(omega*(tau+0.5d0*dt))
                endif
                !print*, 'XX m 12-: ytop = ',ytop
                CALL MOVEB ( fytop, fybot,fxbot, topmass)
                !print*, 'XX m 12+: ytop = ',ytop
                !     COMMENT NEXT LINE TO STOP ROTATION
                if (iroll == 0) CALL ROTATEB
            endif

            tau = tau + dt

            !     IF CRUSHING, MODIFY WALL STRESS
            if (icompact) then
                presy=presy+dt*prate
            endif

            ! also sum up the average velocities into lavevx and lavevy
            do iycell=1,my
                do ixcell=1,mx
                    icell=ixcell+mx*(iycell-1)
                    k=head(icell)
                    do while (k > 0)
                        !if (bdgrain(k) == 0) then
                        layervx(iycell) =layervx(iycell) + vx(k)
                        layervy(iycell) =layervy(iycell) + vy(k)
                        layertempx(iycell)=layertempx(iycell)+vx(k)*vx(k)
                        layertempy(iycell)=layertempy(iycell)+vy(k)*vy(k)
                        layerw(iycell) =layerw(iycell) + w(k)
                        ldens(iycell)=ldens(iycell)+1d0
!                            layertemp(iycell)=layertemp(iycell)+ &
!                                (vx(k)*vx(k)+vy(k)*vy(k))
                        layerpor(iycell)=layerpor(iycell)+mass2D(k)
                        lrho(iycell)=lrho(iycell)+1d0/radinv2(k)
                        skxxt(iycell) = skxxt(iycell) + vx(k)*vx(k)/radinv2(k)
                        skxyt(iycell) = skxyt(iycell) + vx(k)*vy(k)/radinv2(k)
                        skyxt(iycell) = skyxt(iycell) + vy(k)*vx(k)/radinv2(k)
                        skyyt(iycell) = skyyt(iycell) + vy(k)*vy(k)/radinv2(k)
                        px(iycell) = px(iycell) + vx(k)/radinv2(k)
                        py(iycell) = py(iycell) + vy(k)/radinv2(k)
                        !endif
                        k=list(k)
                    enddo
                enddo
            enddo
            
            if(ib(2)==0)then
                topg=ry(1)+radius(1)
                do i=nbound(4)+1,n
                    topg=max(topg,ry(i)+radius(i))
                enddo
            else
                topg=ytop
            endif
            
            height=height+topg-ybot

            !  ##### OUTPUT ##### OUTPUT ##### OUTPUT ##### OUTPUT ##### OUTPUT 
            ! ****************************************************************
            ! &&&  FREQUENT OUTPUT &&&&&&&  FREQUENT OUTPUT &&&&
            IF ( MOD( STEP, iprintf ) == 0 ) call iprintf_sub()
            
            IF ( MOD( STEP, IPRINT ) == 0 ) call iprint_sub()

            if ( mod(step,iprint_restart) == 0) call write_restart()

            if( any(ib(1)==[7,8]) .and. mod(tau, period) < dt) then
                call write_restart("period")
                write(*,*) "write period"
            endif
            ! &&&  INFREQUENT OUTPUT &&&&&&&  INFREQUENT OUTPUT &&&&
            !  ##### OUTPUT ##### OUTPUT ##### OUTPUT ##### OUTPUT ##### OUTPUT 
            
!            if (phigoal0 > 0.) then
!                call calculate_area_and_porosity()
!                if (1. - grainarea/xl/yl <= phigoal0) exit
!            endif

            if(sim_duration()) exit

        enddo

        !    *******MAIN*******MAIN*******MAIN*******MAIN*******MAIN*******MAIN
        !    ** MAIN LOOP ENDS                                                **
        !    *******************************************************************

        if(sim_duration()) exit

        ! find the top of the grains
        if (phigoal > 0.) then
            call set_phi()
            call calculate_area_and_porosity()
            write(*,*) "phi set to ", 1. - grainarea/xl/yl
    
            cycle
        endif
        exit
    enddo

    WRITE(*,'(/1X,''**** END OF DYNAMICS **** ''//)')

    call pscriptplot(next+1,igif,fnn,0)
    call write_restart("final")

    write(*,*) 'masstot = ',masstot

    ! finalized output

    
!    call postScriptTool_clean()
 
contains   
!=======================================================================
logical function sim_duration() result(res)
    real(8) :: actual_time
    
    res = .false.
    call cpu_time(actual_time)
    if(actual_time - start_time > max_time) res = .true.

end function
!=======================================================================
subroutine iprintf_sub()

    real(8) sxt,syt,sxb,syb
    real(8) area, area2, phi, phi0,strain,rmass,coordave
    
    !  * coordination number
    coordave=0.
    do i=indint,n
        coordave=coordave+coordn(i)
    enddo
            
    coordave=coordave/dble(nintshort)
    call write_to('coordination.out', tau, coordave)

    if (icompact) then
        strain=-((ytop-ybot)-h0)/h0
        call write_to('strain', tau,strain, presy)
    endif

    ! distance between the two walls
    call write_to('dilatation.out', tau,topg-ybot)

    !  * porosity

    area2 = (xright-xleft)*(topg-ybot)
    area = (xright-xleft)*(ytop-ybot)

    phi = 1.-grainarea/area
    phi0 = 1.-grainarea0/area2
    call write_to('porosity.out', tau,phi,phi0)
    
    ! these are the average forces on each grain. 
    !  convert to stresses
    sxt=fxtop*(nbound(2)-nbound(1))/(xright-xleft)
    syt=fytop*(nbound(2)-nbound(1))/(xright-xleft)
    sxb=fxbot*nbound(1)/(xright-xleft)
    syb=fybot*nbound(1)/(xright-xleft)
    
    !if (ifluid == 0) then
    !    call write_to('stresses.out', tau,sxt,syt,sxb,syb,fgrain)
    !else
        call write_to('stresses.out', tau,sxt,syt,sxb,syb)
    !endif
    
    !  waves at walls
    if (ipulse) then
        if (ib(2) == 5) then
            call write_to('waves.out', tau, fspring,fxbot,fybot)
        endif
    endif

    !  spring-driven wall
    !  if with fluid top friction also contains the ratio shear stress/effective normal
    !  when fgrain contains the effective normal stress
    if (ib(2) == 5) then
!        if (ifluid == 0) then
!            call write_to('topfriction.out', tau, fspring/presy,  &
!                fspring/fgrain, vx(nbound(1)+1)/fb(2))
!        else
!            call write_to('topfriction.out', tau, fspring/presy,  &
!                vx(nbound(1)+1)/fb(2))
!        endif
        
        !call write_to('spring.out', tau, fx(nbound(1)+1),dwall,dspring)
        call write_to('spring.out', tau, fspring,dwall,dspring)
    endif

!        if (ifluid == 0) then
!            call write_to('topfriction.out', tau, -sxt/presy, -sxt/fgrain)
!        else
!            call write_to('topfriction.out', tau, -sxt/presy)
!        endif
    call write_to('topfriction.out', tau, fxtop,fytop,fgrain, &
        vx(nbound(2)),vx(1),fb(2))

    if ((ib(1) == 7).or.(ib(1) == 8)) then
        i=int(nbound(1)/2d0)
        call write_to('periodicf.out', tau, fsin, vx(i), rx(i))
    endif

    !  * kinetic energy and momentum
    ek = 0.0
    ekr = 0.0
    mvx=0.
    mvy=0.
    angmom=0.
    do i=indint,n
        rmass=1d0/radinv2(i)
        ek = ek +  &
        rmass*(vx(i)*vx(i) + vy(i)*vy(i)) 
        ekr = ekr + inertiamom(i)*w(i)*w(i)
        mvx = mvx + rmass*vx(i) 
        mvy = mvy + rmass*vy(i) 
        angmom = angmom + inertiamom(i)*w(i)
    enddo
    ek = 0.5  * ek
    ekr = 0.5  * ekr

    ev   = ev  / real ( nintshort )
    ek   = ek  / real ( nintshort )
    ekr  = ekr  / real ( nintshort )
    mvx  = mvx / real(nintshort)
    mvy  = mvy / real(nintshort)
    angmom  = angmom / real(nintshort) 

    call write_to('momentum.out', tau, mvx,mvy,angmom)
    
    call write_to('energy.out', tau, ek, ev, ekr)

    do iycell=1,my
      avesxxt(iycell)=sxxt(iycell)/iprintf
      avesxyt(iycell)=sxyt(iycell)/iprintf
      avesyxt(iycell)=syxt(iycell)/iprintf
      avesyyt(iycell)=syyt(iycell)/iprintf
      aveskxxt(iycell)=skxxt(iycell)/iprintf
      aveskxyt(iycell)=skxyt(iycell)/iprintf
      aveskyxt(iycell)=skyxt(iycell)/iprintf
      aveskyyt(iycell)=skyyt(iycell)/iprintf
      lavevx(iycell)=layervx(iycell)/iprintf
      lavevy(iycell)=layervy(iycell)/iprintf
      lavepx(iycell)=px(iycell)/iprintf
      lavepy(iycell)=py(iycell)/iprintf
      lavetempx(iycell)=layertempx(iycell)/iprintf
      lavetempy(iycell)=layertempy(iycell)/iprintf
      lavew(iycell)=layerw(iycell)/iprintf
      lavedens(iycell)=ldens(iycell)/iprintf
      lavepor(iycell)=layerpor(iycell)/iprintf
      laverho(iycell)=lrho(iycell)/iprintf
      sxxt(iycell)=0d0
      sxyt(iycell)=0d0
      syxt(iycell)=0d0
      syyt(iycell)=0d0
      skxxt(iycell)=0d0
      skxyt(iycell)=0d0
      skyxt(iycell)=0d0
      skyyt(iycell)=0d0
      layervx(iycell)=0d0
      layervy(iycell)=0d0
      px(iycell)=0d0
      py(iycell)=0d0
      layertempx(iycell)=0d0
      layertempy(iycell)=0d0
      layerw(iycell)=0d0
      ldens(iycell)=0d0
      layerpor(iycell)=0d0
      lrho(iycell)=0d0
    enddo
    heightave=height/iprintf
    height=0d0

    do iycell=1,my+1
        do ixcell=1,mx  
            permavetp(ixcell,iycell)=permave(ixcell,iycell)/iprintf
            sitevxavetp(ixcell,iycell)=sitevxave(ixcell,iycell)/iprintf
            sitevyavetp(ixcell,iycell)=sitevyave(ixcell,iycell)/iprintf
            dpdttp(ixcell,iycell)=pressave(ixcell,iycell)/iprintf - &
                pressavetp(ixcell,iycell)
            pressavetp(ixcell,iycell)=pressave(ixcell,iycell)/iprintf
            divvelavetp(ixcell,iycell)=divvelave(ixcell,iycell)/iprintf
            solfraavetp(ixcell,iycell)=solfraave(ixcell,iycell)/iprintf
            numparavetp(ixcell,iycell)=numparave(ixcell,iycell)/iprintf
            permave(ixcell,iycell)=0d0
            sitevxave(ixcell,iycell)=0d0
            sitevyave(ixcell,iycell)=0d0
            pressave(ixcell,iycell)=0d0
            divvelave(ixcell,iycell)=0d0
            solfraave(ixcell,iycell)=0d0
            numparave(ixcell,iycell)=0d0
        enddo
    enddo       
end subroutine
!=======================================================================
subroutine iprint_sub()

    !heightave=heightave+height

    print*,'writing postscript page ',next+1,' at step ',step
    next=next+1
    inext=0
    call pscriptplot(next,igif,fnn,0)

    if(ib(2)==0)then
        histdy=(ytop-ybot)/my
    else
        histdy=heightave/my
    endif

    call write_to('shearvel.out', "# "//fnn)
    do iycell=1,my
        call write_to('shearvel.out', histdy*(iycell-0.5d0), &
            lavevx(iycell)/max(lavedens(iycell),1d0), &
            lavevy(iycell)/max(lavedens(iycell),1d0), &
            (lavetempx(iycell)-lavevx(iycell)**2/max(lavedens(iycell),1d0))/max(lavedens(iycell),1d0), &
            (lavetempy(iycell)-lavevy(iycell)**2/max(lavedens(iycell),1d0))/max(lavedens(iycell),1d0), &
            lavew(iycell)/max(lavedens(iycell),1d0))
    enddo
    call write_to('shearvel.out', " ")

    call write_to('dens.out', "# "//fnn)
    do iycell=1,my
        call write_to('dens.out', histdy*(iycell-0.5d0), &
            lavedens(iycell)/(histdy*(xright-xleft)), &
            pi/4d0*lavepor(iycell)/(histdy*(xright-xleft)), &
            laverho(iycell)/(histdy*(xright-xleft)))
    enddo
    call write_to('dens.out', " ")

    call write_to('stress.out', "# "//fnn)
    do iycell=1,my
        ek=histdy*(xright-xleft)
        call write_to('stress.out', histdy*(iycell-0.5d0), &
            avesxxt(iycell)/ek, avesxyt(iycell)/ek, &
            avesyxt(iycell)/ek, avesyyt(iycell)/ek)
    enddo
    call write_to('stress.out', " ")

    call write_to('kin_stress.out', "# "//fnn)
    do iycell=1,my
        ek=histdy*(xright-xleft)
        call write_to('kin_stress.out', histdy*(iycell-0.5d0), &
         (aveskxxt(iycell) - 2d0*lavepx(iycell)*lavevx(iycell)/max(lavedens(iycell),1d0) + &
          laverho(iycell)*lavevx(iycell)*lavevx(iycell)/max(lavedens(iycell)*lavedens(iycell),1d0))/ek, &
         (aveskxyt(iycell) - lavepx(iycell)*lavevy(iycell)/max(lavedens(iycell),1d0) - &
          lavepy(iycell)*lavevx(iycell)/max(lavedens(iycell),1d0) + &
          laverho(iycell)*lavevx(iycell)*lavevy(iycell)/max(lavedens(iycell)*lavedens(iycell),1d0))/ek, &
         (aveskyxt(iycell) - lavepx(iycell)*lavevy(iycell)/max(lavedens(iycell),1d0) - &
          lavepy(iycell)*lavevx(iycell)/max(lavedens(iycell),1d0) + &
          laverho(iycell)*lavevx(iycell)*lavevy(iycell)/max(lavedens(iycell)*lavedens(iycell),1d0))/ek, &
         (aveskyyt(iycell) - 2d0*lavepy(iycell)*lavevy(iycell)/max(lavedens(iycell),1d0) + &
          laverho(iycell)*lavevy(iycell)*lavevy(iycell)/max(lavedens(iycell)*lavedens(iycell),1d0))/ek)
    enddo
    call write_to('kin_stress.out', " ")

    !height=0d0

    if (ifluid==0) then
        call write_to("pgradfile.out", "# "//fnn)
        call write_to("press.out", "# "//fnn)
        do i=1,mx
            do j = 1,my+1
                call write_to("pgradfile.out", xleft+(i-1)*dxc, (j-1)*dyc, &
                    gradpx(i,j), gradpy(i,j), vol_flux_x(i,j), vol_flux_y(i,j), &
                    diff_vol_flux_x(i,j), diff_vol_flux_y(i,j), &
                    sitevx(i,j), sitevy(i,j))
                call write_to("press.out", xleft+(i-1)*dxc, (j-1)*dyc, &
                    pressavetp(i,j), dpdttp(i,j)/(iprintf*dt))
            enddo
        enddo
        call write_to("pgradfile.out", " ")
        call write_to("press.out", " ")

        ! files:
        ! phipermfile - porosity | permeability
        ! areasfile - total grain area | weighted radiuses^2 | volume fraction (area/cellsize) | number grains in a cell
        ! velocityfile - vx | vy | div(v)
        ! numvelofile - numvx | numvy
        call write_to("phipermfile.out", "# "//fnn)
        call write_to("areasfile.out", "# "//fnn)
        call write_to("velocityfile.out", "# "//fnn)
!        call write_to("numvelofile.out", "# "//fnn)
        do i=1,mx
            do j = 1,my+1
                call write_to("phipermfile.out", xleft+(i-1)*dxc, (j-1)*dyc, &
                    1d0-solfraavetp(i,j), permavetp(i,j))
                call write_to("areasfile.out", xleft+(i-1)*dxc,(j-1)*dyc, vc(i,j), &
                    siterd2(i,j), solfraavetp(i,j), numparavetp(i,j))
                call write_to("velocityfile.out", xleft+(i-1)*dxc, (j-1)*dyc, &
                    sitevxavetp(i,j), sitevyavetp(i,j), divvelavetp(i,j))
!                call write_to("numvelofile.out", xleft+(i-1)*dxc, (j-1)*dyc, &
!                    numvx(i,j), numvy(i,j))
            enddo
        enddo
        call write_to("phipermfile.out", " ")
        call write_to("areasfile.out", " ")
        call write_to("velocityfile.out", " ")
!        call write_to("numvelofile.out", " ")
    endif   ! ifluid eq. 0
    
end subroutine
!=======================================================================
subroutine set_phi()
    topg=ry(1)+radius(1)
    do i=nbound(4)+1,n
        topg=max(topg,ry(i)+radius(i))
    enddo
    
    wmax=0.
    do i=nbound(1)+1,nbound(2)
        wmax=max(wmax,radius(i))
    enddo
    
!    write(*,*) "ib(2) = ",ib(2)
    
!    tmp = grainarea / (1. - phigoal) / (xright-xleft) + ybot
!    if(tmp+wmax > topg+wmax) then
!         topg = tmp - wmax
!         ib(2) = 0
!         presy=0.
!     else
!         ib(2) = 2
!         phigoal0 = phigoal
!         presy=presyfinal
!     endif
    ib(2) = 1
!    phigoal0 = phigoal
    presy=presyfinal
    gy=0.
    phigoal = 0.

    write(*,*) "ib(2) = ",ib(2)
    
    !write(*,*) "set_phi()", topg, grainarea, phigoal
                
    do i=nbound(1)+1,nbound(2)
        ry(i)=topg+wmax
    enddo

end subroutine
!=======================================================================
subroutine zero_variables()
! DO NOT INCLUDE VARIABLES THAT ARE READ FROM RESTART FILE, SO THAT THEY
! ARE NOT OVERWRITTEN BEFORE THE MAIN LOOP
    do iycell=1,my
        ldens(iycell)=0d0
        lavedens(iycell)=0d0
        lrho(iycell)=0d0
        laverho(iycell)=0d0
        layerpor(iycell)=0d0
        lavepor(iycell)=0d0
        layervx(iycell)=0.
        layervy(iycell)=0.
        lavevx(iycell)=0.
        lavevy(iycell)=0.
        px(iycell)=0.
        py(iycell)=0.
        lavepx(iycell)=0.
        lavepy(iycell)=0.
        lavetempx(iycell)=0d0
        lavetempy(iycell)=0d0
        layertempx(iycell)=0d0
        layertempy(iycell)=0d0
        layerw(iycell)=0d0
        lavew(iycell)=0d0
        sxxt(iycell)=0d0
        sxyt(iycell)=0d0
        syxt(iycell)=0d0
        syyt(iycell)=0d0
        avesxxt(iycell)=0d0
        avesxyt(iycell)=0d0
        avesyxt(iycell)=0d0
        avesyyt(iycell)=0d0
        skxxt(iycell)=0d0
        skxyt(iycell)=0d0
        skyxt(iycell)=0d0
        skyyt(iycell)=0d0
        aveskxxt(iycell)=0d0
        aveskxyt(iycell)=0d0
        aveskyxt(iycell)=0d0
        aveskyyt(iycell)=0d0
    enddo
    
    height=0d0
    heightave=0d0

end subroutine
!=======================================================================
subroutine zero_fluid_outputs()

    do iycell=1,my+1
        do ixcell=1,mx  
            permave(ixcell,iycell)=0d0
            sitevxave(ixcell,iycell)=0d0
            sitevyave(ixcell,iycell)=0d0
            pressave(ixcell,iycell)=0d0
            divvelave(ixcell,iycell)=0d0
            solfraave(ixcell,iycell)=0d0
            numparave(ixcell,iycell)=0d0
            pressavetp(ixcell,iycell)=0d0
        enddo
    enddo

end subroutine
!=======================================================================
subroutine calculate_area_and_porosity()
    grainarea=0.
    averradius2 = 0.
    
    do i=1,n
        averradius2 = averradius2 + radius(i)*radius(i);
        if (bdgrain(i) == 0) then
            grainarea = grainarea+radius(i)*radius(i)*pi
        else
            if ((rx(i) >= xleft).and.(rx(i) <= xright)) then
                grainarea = grainarea + .5*radius(i)*radius(i)*pi
            endif
        endif
    enddo
    
    grainarea0=grainarea
    averradius2 =  averradius2/n
    !write(6,*) 'initial grain area in main=',grainarea
    xl=xright-xleft
    yl=ytop-ybot
end subroutine

end program
!=======================================================================
subroutine dumpgrains(afrac1)
    use mycommons
    use mod_inputReader
    use mod_outputWriter
    implicit none
    
    integer n1sum,n2sum
    integer i
    
    real(8) a1sum,a2sum,r1sum,r2sum,s1sum,s2sum,afrac1
    real(8) amean1,amean2,astd1,astd2,alogr2,ahigr2,alogr,ahigr
    real(8) a1,a2,tmean,tstd
    
    emod = 1.
    do i=1,n
        call write_to('grainsizes.out', gtype(i),2*radius(i),bdgrain(i))
    enddo
    call close_file('grainsizes.out')

    !  get stats about grains
    a1sum=0.
    a2sum=0.
    n1sum=0
    n2sum=0
    r1sum=0.
    r2sum=0.
    s1sum=0.
    s2sum=0.
    alogr=higr
    ahigr=logr
    alogr2=higr2
    ahigr2=logr2
    
    do i=1,n
        if (gtype(i) == 1) then
            r1sum=r1sum+radius(i)
            a1=pi*radius(i)*radius(i)
            if (bdgrain(i) == 0) then
                a1sum=a1sum+a1
            else
                a1sum=a1sum+0.5*a1
            endif
            n1sum=n1sum+1
            alogr=min(alogr , 2.d0*radius(i))
            ahigr=max(ahigr , 2.d0*radius(i))
        else
            r2sum=r2sum+radius(i)
            a2=pi*radius(i)*radius(i)
            if (bdgrain(i) == 0) then
                a2sum=a2sum+a2
            else
                a2sum=a2sum+0.5*a2
            endif
            n2sum=n2sum+1
            alogr2=min(alogr2,2d0*radius(i))
            ahigr2=max(ahigr2,2d0*radius(i))
        endif
    enddo
    
    amean1=2.*r1sum/n1sum
    if (n2sum > 0) amean2=2.*r2sum/n2sum

    do i=1,n
        if (gtype(i) == 1) then
            s1sum=s1sum+(2.*radius(i)-amean1)**2
        else
            s2sum=s2sum+(2.*radius(i)-amean2)**2
        endif
    enddo
    
    astd1=s1sum/n1sum
    astd2=s2sum/max(1,n2sum)

    call write_to('graininfo.out', 'DISTRIBUTION 1: ')
    call write_to('graininfo.out', '   total grains:', n1sum)
    call write_to('graininfo.out', '   total area:', a1sum)
    call write_to('graininfo.out', '            target values     actual values')
    call write_to('graininfo.out', '   mean     ', mean1,amean1)
    call write_to('graininfo.out', '   std. dev.', std1,astd1)
    call write_to('graininfo.out', '   smallest ', logr,alogr)
    call write_to('graininfo.out', '   largest  ', higr,ahigr)
    call write_to('graininfo.out', " ")
    call write_to('graininfo.out', 'DISTRIBUTION 2: ')
    call write_to('graininfo.out', '   total grains:', n2sum)
    call write_to('graininfo.out', '   total area:', a2sum)
    call write_to('graininfo.out', '              target values     actual values')
    call write_to('graininfo.out', '   mean     ', mean2,amean2)
    call write_to('graininfo.out', '   std. dev.', std2,astd2)
    call write_to('graininfo.out', '   smallest ', logr2,alogr2)
    call write_to('graininfo.out', '   largest  ', higr2,ahigr2)
    call write_to('graininfo.out', " ")

    afrac1=a1sum/(a1sum+a2sum)
    call write_to('graininfo.out', 'fractional area of Dist. 1:', afrac1)

    call write_to('graininfo.out', " ")
    tmean = (amean1*n1sum+amean2*n2sum)/dble(n)
    call write_to('graininfo.out', 'mean diameter of entire assemblage:', tmean)
    do i=1,n
        tstd=tstd+(2.*radius(i)-tmean)**2
    enddo
    tstd=tstd/dble(n)
    call write_to('graininfo.out', 'std. dev. of entire assemblage:', tstd)

    call close_file("graininfo.out")
    
end subroutine
