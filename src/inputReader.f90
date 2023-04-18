module mod_inputReader
    use mycommons
!***********************************************************************
!
! MODULE FOR INPUT FILE READING
!
! Martin Svoboda
! svobod.martin@gmail.com .or. svobodam@icpf.cas.cz 
! 2019
!
!***********************************************************************

    implicit none
   
    interface read_val
        ! overload interface for statements reading function
        procedure :: read_char, read_int, read_intarr, read_real, &
                    read_realarr, read_logic, read_int_realarr, &
                    read_int8
    end interface
       
    ! DURATION AND OUTPUT FREQUENCY ------------------------------------
    character(len=50), protected :: title = "unnamed"     ! title
    
    integer(8),protected :: NSTEP = 0              ! number of time steps to do
    integer(8),public    :: iprint = huge(0)       ! postscript and profiles page plotted every iprint time steps
    integer(8),public    :: iprintf = huge(0)      ! frequent output dumped every iprintf time steps
    integer(8),public    :: iprint_restart = huge(0)   ! restart printed  every iprint_restart time steps
    real(8),public    :: max_time = huge(0.d0)   ! maximum seconds of computing
    
    
    ! PARAMETERS FOR BUILDING A NEW SYSTE2M ----------------------------
    integer,protected :: irstart = 1 ! =0 build a new set of grains,
                                     ! =1 read file 'restart'
                                     ! =2 initialize fluid pressure with some fraction of presy
    integer,protected :: ifill = 0  ! =0 don't make any new particles (not useful)
                                     ! =-1 pack in a box full of grains at the start
    integer,protected :: rand_seed = 0         ! random generator seed if 0 then generate random seed at start
    real(8),public    :: boxx = 0., boxy = 0.  ! horizontal and (initial) vertical size of box
    real(8),protected :: boxfy = tiny(0.)      !  >0, final size of box in y direction, where dispalcment is held const. 
    real(8),public    :: phigoal = 0.          ! the target porosity to try to compress
    
    real(8),protected :: frac1 = 1. ! the number fraction of the grains in population 1
    
    real(8),protected :: mean1 = 1.   ! mean diameter of population 1
    real(8),protected :: std1 = 1.    ! standard deviation of population 1
    real(8),protected :: logr = 0.8   ! smallest diameter of population 1
    real(8),protected :: higr = 1.2   ! largest diameter of population 1
    real(8),protected :: e1 = 1.      ! ??? maybe Young's Mod  for gaussian 1
                                            
    real(8),protected :: mean2 = 1.    ! mean diameter of population 2
    real(8),protected :: std2 = 1.     ! standard deviation of population 2
    real(8),protected :: logr2 = 0.8   ! smallest diameter of population 2
    real(8),protected :: higr2 = 1.2   ! largest diameter of population 2
    real(8),protected :: e2 = 1        ! ??? maybe Young's Mod  for gaussian 2
    
    real(8),protected :: sigb = 1. ! mean diameter of boundary grains
    real(8),protected :: w_offset = 0. ! length of walls that is not part of the box with interior grains. Used when ib(3)>=0
    
    ! BOUNDARY AND INITIAL CONDITIONS ----------------------------------  
    integer,protected :: iztau = 0           ! = 0, zero time counter tau
                                             ! = 0, use tau from restart file
    integer,protected :: izmom = 2           ! = 0: zero out all velocities
                                             ! = 1: give all particles fb(2) horizontal velocity
                                             ! other: don't touch velocities, use restart
    integer,protected :: izslip = 2          ! = 0: zero out slip array (shear forces between particles)
                                             ! other: use restart
    integer,public    :: ib(4) = 0           ! boundary type switch (array order: bottom,top,left,right)
                                             ! <0 no wall (for side walls this makes it periodic)
                                             ! =0 stationary wall  (both velocities =0.)
                                             ! =1 no fixed velocities
                                             ! =2 fixed shear velocity
                                             ! =3 fixed normal velocity
                                             ! =4 fixed shear velocity, zero normal velocity

    ! APPLIED FORCES AND VELOCITIES ------------------------------------            
    real(8),protected :: fb(4) = 0      ! boundary velocity to impose on the 4 walls    
    real(8),public    :: presx = 0.,&   ! size of external force applied to left and right, top and bottom
                         presy = 0.        ! expressed as the fraction of an atom diameter that would be 
                                        ! overlapping 
    real(8),public    :: gx = 0.,&      ! body force acting in x and y directions (g/sqrt(mass/k))
                         gy = 0.
    
    ! PHYSICAL PARAMETERS OF GRAINS ------------------------------------
    real(8),protected :: skmol = 1.    ! spring constant for interactions within a molecule, scaled
                                       ! by intermolecule interactions 
    real(8),protected :: skshear = 0.5  ! spring constant for shear interactions within a molecule scaled to skmol
    real(8),protected :: friction = 0.5 ! ratio of shear to normal stress at which a grain begins to slide
    real(8),protected :: gamma = 0.5   ! damping parameter for atom/atom interactions (generally between
                                       ! 0 and 1, 0.8 is a good number)
    real(8),protected :: gammamol = 0. ! damping by pore fluid (viscous damping, 0.01 is a good number?) 
    
    ! FLOW CONTROL AND EFFICIENCY PARAMETERS ---------------------------
    real(8),public    :: dt = 0.1           !time step factor:  
                                  !    maximum stable time step will be calculated, 
                                  !    then multiplied by this factor: should be <= 0.15  
    integer,protected :: mx = 0,& ! number of grid cells in x and y directions
                         my = 0   ! (usually mx < boxx/(2*natom), my < boxy/(2*natom))
    
    ! PLOTTING SWITCHES ------------------------------------------------
    integer,protected :: ibw = 0            ! =1 just read a restart file and make a plot, then quit
                                            ! = 2 ???
    logical,protected :: igif = .false.     ! = false just make postscript output
                                            ! = true convert postscript to jpeg
    integer,public    :: next = 0              ! file extension for first plot from this run
    
    ! SPECIAL RUN SWITCHES ---------------------------------------------
    logical,protected :: ipulse = .false.  ! if true send a seismic pulse through the system
    real(8),protected :: gstep = 0.        ! >0 do a uniform horizontal displacement of all
                                           ! particles from 0 at bottom to gstep at top
                                           ! =0, ignore
    real(8),protected :: kspring = 0.      ! for spring tied to top wall driven runs, used with ib(2)=5
    real(8),protected :: topmass = 1d0     ! mass of each grain in the top wall
    real(8),public :: dwall = 0d0
    real(8),public :: dspring = 0d0
    
    ! BREAKING GRAINS PARAMETERS ---------------------------------------
    integer,protected :: crush1 = 0     ! switch for breakage (0=don't break,1=allow breaking)
    integer,protected :: pieces = 7     ! number of pieces into which a particle breaks (3, 7, or 13)
    
    integer,protected :: typest = 2   ! switch for strength.(0=all particles have same strength equal to purest
                                      ! 1=strength inversly depends on grain size with random  noise added
                                      ! 2=strength inversly depends on grain size. No random noise added
                                      ! 3=strength is equal to purest plus random noise)

    real(8),protected :: purest = 0.0004  ! max value of strength
    integer,protected :: esp = 1          ! exponent used in grain size dependence relation with strength. Could be any
                                          ! integer value, but we can limit our values to 1 or 2.
    logical,protected :: restrg = .false. ! switch for re-assigning strengts to grains    
                                          ! false=strengths are the same as the ones assigned at the creation of grains
                                          ! true=re-assignes the strengths to the grains 
    
    logical,protected :: icompact = .false.   ! ???
    real(8),protected :: prate = 0.           ! ???
    real(8),protected :: h0 = 0.              ! ???
    
    logical,protected :: idefect = .false.    ! ???
    real(8),protected :: percentdefects = 0.  ! ???
    real(8),protected :: sizevar = 0.         ! ???
    
    real(8),protected :: permfac = 1.d-5      ! ???
    
    ! gfac is used only in sedimentation test to determine the distance between 
    ! grains in the initial configuration. gfac = 1. means no free space. 
    real(8),protected :: gfac = 1.
    
    ! fstress is a factor determining the initial fluid pressure in the system
    ! the normal stress (pressy) will be multiplied by this factor.     
    real(8),protected :: fstress = 0.
    
    ! fbc are boundary conditions for the fluid pressure at the top and bottom of the box. possibilities:
    ! 1 = Dirichlet, 2 = Neumann with derivative = 0, 
    ! 3 = Dirichlet with permeability test i.e. not moving the grains in responce to fluid pressure.
    ! 4 = Mixed boundary conditions. Constant pressure gradient at the top and constant pressure at bottom
    ! 5 = Mixed boundary conditions. Constant pressure gradient at the bottom and constant pressure at top
    ! 6 = Uniform pressure
    integer,protected :: fbc = 1
    
    integer,protected :: contact_model = 0   ! ???
    integer,protected :: ifluid = 1          ! ifluid: 1 - no fluid; 0 - by default with fluid
    integer,protected :: intruder = 0        ! ???
    integer,protected :: iroll = 0           ! iroll: 1 - no roll; 0 - by default with rolling

    real(8),protected :: densratio = 1.      ! ???
    real(8),protected :: rintruder = 0.      ! ???
    real(8),public    :: wl = -1.            ! ???
    
    real(8),protected :: acoef = 1d-1 !3.33d-2    ! for d=1cm,rho=2640kg/m3,E=1e10Pa^M
    real(8),protected :: bcoef = 3.72d6
    real(8),protected :: const_perm = -1d0  ! switch and value (if positive) for uniform permeability
    
    ! for periodic acceleration of the the bottom wall 
    ! SEISMIC WAVE PROPAGATION      
    integer(4),protected :: pall = 0
    real(8),protected :: period = 0.d0 ,&
                         omega = 0.d0 ,&
                         amp = 0.d0

    integer,protected :: jpeg_resolution = 450

    integer,protected :: ext_force_to = 0
    real(8),protected :: ext_force(2) = 0.d0
    logical,protected :: use_ext_force = .false.
    logical,protected :: yn_advect = .true.  ! switch for including advection term into pressure equation
    logical,protected :: yn_explicit = .false.  ! switch for employing the explicit scheme for solving the pressure equation
    
    private :: calculate_omited, read_intarr, read_realarr, read_char, &
               read_int, read_logic, read_real, read_val, read_int_realarr, &
               read_int8

    contains
!=======================================================================
subroutine read_input_file()

    integer istat
    
    character(len=100) one_line
    character(len=20) cmd

    print*, ' **  PROGRAM UniGranDam **'
    WRITE(*,*) ' MOLECULAR DYNAMICS OF UNIFORM GRANULAR MEDIA'
    WRITE(*,*) 'VELOCITY-VERLET ALGORITHM, TRINGULAR GRAINS'

    istat=0
    
    do
        read(*,'(a100)',iostat=istat) one_line
        
        if(istat/=0) exit
        
        if(len_trim(one_line)==0) cycle
        
        if(one_line(1:1) == "#") cycle ! note

        read(one_line,*) cmd
        select case(trim(cmd))
            case("title")
                call read_val(title, one_line)
            case("nstep")
                call read_val(nstep, one_line)
            case("iprint")
                call read_val(iprint, one_line)
            case("iprintf")
                call read_val(iprintf, one_line)
            case("iprint_restart")
                call read_val(iprint_restart, one_line)
            case("irstart")
                call read_val(irstart, one_line)
            case("ifill")
                call read_val(ifill, one_line)
            case("boxx")
                call read_val(boxx, one_line)
            case("boxy")
                call read_val(boxy, one_line)
            case("boxfy")
                call read_val(boxfy, one_line)
            case("phigoal")
                call read_val(phigoal, one_line)
            case("frac1")
                call read_val(frac1, one_line)
            case("mean1")
                call read_val(mean1, one_line)
            case("std1")
                call read_val(std1, one_line)
            case("logr")
                call read_val(logr, one_line)
            case("higr")
                call read_val(higr, one_line)
            case("e1")
                call read_val(e1, one_line)
            case("mean2")
                call read_val(mean2, one_line)
            case("std2")
                call read_val(std2, one_line)
            case("logr2")
                call read_val(logr2, one_line)
            case("higr2")
                call read_val(higr2, one_line)
            case("e2")
                call read_val(e2, one_line)
            case("sigb")
                call read_val(sigb, one_line)
            case("w_offset")
                call read_val(w_offset, one_line)
            case("seed")
                call read_val(rand_seed, one_line)
            case("iztau")
                call read_val(iztau, one_line)
            case("izmom")
                call read_val(izmom, one_line)
            case("izslip")
                call read_val(izslip, one_line)
            case("ib")
                call read_val(ib, one_line)
            case("fb")
                call read_val(fb, one_line)
            case("presy")
                call read_val(presy, one_line)
            case("presx")
                call read_val(presx, one_line)
            case("gx")
                call read_val(gx, one_line)
            case("gy")
                call read_val(gy, one_line)
            case("skmol")
                call read_val(skmol, one_line)
            case("skshear")
                call read_val(skshear, one_line)
            case("friction")
                call read_val(friction, one_line)
            case("gamma")
                call read_val(gamma, one_line)
            case("gammamol")
                call read_val(gammamol, one_line)
            case("dt")
                call read_val(dt, one_line)
            case("mx")
                call read_val(mx, one_line)
            case("my")
                call read_val(my, one_line)
            case("ibw")
                call read_val(ibw, one_line)
            case("ijpeg","igif")
                call read_val(igif, one_line)
            case("next")
                call read_val(next, one_line)
            case("ipulse")
                call read_val(ipulse, one_line)
            case("gstep")
                call read_val(gstep, one_line)
            case("kspring")
                call read_val(kspring, one_line)
            case("topmass")
                call read_val(topmass, one_line)
            case("crush1")
                call read_val(crush1, one_line)
            case("pieces")
                call read_val(pieces, one_line)
            case("typest")
                call read_val(typest, one_line)
            case("purest")
                call read_val(purest, one_line)
            case("esp")
                call read_val(esp, one_line)
            case("restrg")
                call read_val(restrg, one_line)
            case("idefect")
                call read_val(idefect, one_line)
            case("percentdefects")
                call read_val(percentdefects, one_line)
            case("sizevar")
                call read_val(sizevar, one_line)
            case("permfac")
                call read_val(permfac, one_line)
            case("contact_model")
                call read_val(contact_model, one_line)
            case("iroll")
                call read_val(iroll, one_line)
            case("ifluid")
                call read_val(ifluid, one_line)
            case("fstress")
                call read_val(fstress, one_line)
            case("gfac")
                call read_val(gfac, one_line)
            case("intruder")
                call read_val(intruder, one_line)
            case("rintruder")
                call read_val(rintruder, one_line)
            case("densratio")
                call read_val(densratio, one_line)
            case("fbc")
                call read_val(fbc, one_line)
            case("wl")
                call read_val(wl, one_line)
            case("icompact")
                call read_val(icompact, one_line)
            case("prate")
                call read_val(prate, one_line)
            case("h0")
                call read_val(h0, one_line)
            case("acoef")
                call read_val(acoef, one_line)
            case("bcoef")
                call read_val(bcoef, one_line)
            case("max_time")
                call read_val(max_time, one_line)
            case("jpeg_resolution","resolution")
                call read_val(jpeg_resolution, one_line)
            case("ext_force_to")
                call read_val(use_ext_force, one_line)
                call read_val(ext_force_to, ext_force, one_line)
            case("yn_advect")
                call read_val(yn_advect, one_line)
            case("yn_explicit")
                call read_val(yn_explicit, one_line)
            case("const_perm")
                call read_val(const_perm, one_line)
            case default
                write(*,*) "Unknown input command: ",cmd
                stop
        end select
    enddo
    
    ! for periodic acceleration of the the bottom wall 
    ! SEISMIC WAVE PROPAGATION      
    if (ib(1) == 7 .or. ib(1) == 8 .or. ipulse) call read_waveinput()

    call calculate_omited()

end subroutine
!=======================================================================
subroutine read_fbcinput(fbc,irstart)

    integer(4),intent(in) :: fbc,irstart
    integer(4) :: nu, i, j

    open(newunit=nu,file='fbcinput.dat')
    if (any(fbc == [1,3,6])) then
        read(nu,*) pbot,ptop,pall
        print*,  "pbot, ptop, pall =", pbot,ptop,pall
        press_bot(:)=pbot
        
        if(pall==0)then
            read(nu,*) (press_bot(i),i=1,mx+1)
            if(irstart /= 1)then
                do j=2,my
                    read(nu,*) (press(i,j),i=1,mx+1)
                enddo
            endif
        endif
        if(irstart /= 1)then
            if (pall == 1) then
                do i=1,MAXD+1
                    do j=1,MAXD+1
                        press(i,j) = pbot
                    enddo
                enddo
            elseif(pall == 2) then
                do i=1,MAXD+1
                    do j=1,MAXD+1
                        press(i,j) = ptop
                    enddo
                enddo
            endif
        endif
    elseif (any(fbc == [4,5])) then
        read(nu,*) fluxbot,ptop
        print*,  'fluxbot, ptop =', fluxbot,ptop
    endif
    close(nu)
    
end subroutine
!=======================================================================
subroutine read_waveinput()
    open(unit=7,file='waveinput') 
        read(7,*) period,amp
        print*,  "period, amp = ", period,amp
    close(7)
    omega=2.*pi/period
end subroutine
!=======================================================================
subroutine calculate_omited()

  if ( iprint <= 0 ) iprint = huge(0)
  if ( iprintf <= 0 ) iprintf = huge(0)
  
  if ( iprint_restart <= 0 ) iprint_restart = huge(0)
  if ( iprint_restart == huge(0) ) iprint_restart = iprint

  if(mx == 0) then
      mx = floor(boxx*0.5)
      write(*,*) "mx set to", mx
  endif
  
  if(my == 0) then
      my = floor(boxy*0.5)
      write(*,*) "my set to", my
  endif

  visco = acoef/bcoef

end subroutine
!=======================================================================
subroutine read_char(var, one_line)
    character(len=100),intent(in) :: one_line
    character(len=20) :: cmd
    character(len=*),intent(out) :: var
    
    cmd = ""
    
    read(one_line,*) cmd, var
    write(*,*) cmd, var
  
end subroutine
!=======================================================================
subroutine read_int(var, one_line)
    character(len=*),intent(in) :: one_line
    character(len=20) :: cmd
    integer,intent(out) :: var

    cmd = ""
    
    read(one_line,*) cmd, var
    write(*,*) cmd, var
  
end subroutine
!=======================================================================
subroutine read_int8(var, one_line)
    character(len=*),intent(in) :: one_line
    character(len=20) :: cmd
    integer(8),intent(out) :: var

    cmd = ""
    
    read(one_line,*) cmd, var
    write(*,*) cmd, var
  
end subroutine
!=======================================================================
subroutine read_intarr(var, one_line)
    character(len=*),intent(in) :: one_line
    character(len=20) :: cmd
    integer,dimension(:),intent(out) :: var

    cmd = ""
    
    read(one_line,*) cmd, var
    write(*,*) cmd, var
  
end subroutine
!=======================================================================
subroutine read_realarr(var, one_line)
    character(len=*),intent(in) :: one_line
    character(len=20) :: cmd
    real(8),dimension(:),intent(out) :: var

    cmd = ""
    
    read(one_line,*) cmd, var
    write(*,*) cmd, var
  
end subroutine
!=======================================================================
subroutine read_real(var, one_line)
    character(len=*),intent(in) :: one_line
    character(len=20) :: cmd
    real(8),intent(out) :: var

    cmd = ""
    
    read(one_line,*) cmd, var
    write(*,*) cmd, var
  
end subroutine
!=======================================================================
subroutine read_logic(var, one_line)
    character(len=*),intent(in) :: one_line
    character(len=20) :: cmd
    logical,intent(out) :: var
    
    cmd = ""
    !var = .true.
    
    read(one_line,*) cmd, var
    write(*,*) cmd, var
  
end subroutine
!=======================================================================
subroutine read_int_realarr(var_int, var_real, one_line)
    character(len=*),intent(in) :: one_line
    character(len=20) :: cmd
    integer,intent(out) :: var_int
    real(8),dimension(:),intent(out) :: var_real
    
    cmd = ""
    
    read(one_line,*) cmd, var_int, var_real
    write(*,*) cmd, var_int, var_real
  
end subroutine

end module
