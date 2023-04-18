module mycommons 
    implicit none

    public

    integer, parameter :: maxn    = 40000
    integer, parameter :: maxcell = int(maxn/3d0)
    integer, parameter :: maxk    = 4*maxn
    integer, parameter :: maxd    = 140
    integer, parameter :: tratio  = 1

    real(8), parameter :: pi = 4.*atan(1.)
    real(8), parameter :: ftiny = 1d-13

    ! CURRENT TALLY OF ARRAYS AND SIZES
    !  SIZE MAXN: 29 real(8) + 5 integer 
    !  SIZE MAXK=MAXN*25: 4 real(8) + 3 integer 
    !  SIZE MAXCELL: 2
    !  TOTAL: 1352*MAXN bytes

    !  ATOM ARRAYS    
    integer :: n = 0
    integer gtype(maxn),bdgrain(maxn)
    integer color(maxn)
    
    integer randum
    
    real(8) fspring,uspring
 
    real(8) rx(MAXN), ry(MAXN)
    real(8) vx(MAXN), vy(MAXN)
    real(8) fx(MAXN), fy(MAXN)
    real(8) rnt(MAXN),w(MAXN),tq(MAXN)
    real(8) radius(MAXN),radinv2(MAXN),inertiamom(MAXN)
    real(8) mass2D(MAXN),coordn(MAXN)
    real(8) emod(maxn)
    
    !    **  N        NO. OF ATOMS
    !    **  INDINT    INDEX OF FIRST NON-BOUNDARY ATOM
    !    **  NINTN    NO. OF ATOMS (NOT IN BOUNDARY)
    !    **  RX,RY    ATOMIC POSITIONS
    !    **  VX,VY    ATOMIC VELOCITIES
    !    **  FX,FY    ATOMIC FORCES
    !    **  RNT        ORIENTATION OF THE ATOM (DEGREES CCW FROM EAST)
    !    **  W        ANGULAR VELOCITY OF ATOM
    !    **  TQ        TORQUE ON ATOM
    !    **  RADIUS     ATOMIC RADIUS
    !    **  RADINV2    INVERSE SQUARE OF RADIUS
    !       **  ONBD    TRUE IF ATOM IS TOUCHING A BOUNDARY ATOM
    !       **  SX,SY,SXY   NORMAL AND SHEAR STRESSES IN X-Y SYSTEM
    !       **  S1,S2       PRINCIPLE STRESSES
    !       **  THETA       ORIENTATION OF PRINCIPLE COMPRESSIVE STRESS
    !       **  COORDN      COORDINATION NUMBER OF GRAIN
    !       **  DEDGE    DISTANCE TO PERIODIC BOUNDARY
    !    **  COLOR    COLOR OF ATOM IN LEFT-HAND PICTURE
    !       **  GTYPE    INDEX OF THE DISTRIBUTION (GRAIN TYPE)
    !       **  EMOD        YOUNG'S MODULUS OF THE GRAIN

    !  BOUNDARY CONDITION ARRAYS
    integer nbound(4)
    
    real(8) vsin

    !   LINK-LIST ARRAYS
    integer list(2*MAXN),head(MAXCELL),map(4*MAXCELL)
    integer mywl
    
    !  FRICTION ARRAYS
    !      integer numcont(maxn),contact(20,2,maxn)
    !      real(8) slip(20,maxn)
    !      COMMON/contacts/ numcont, contact, slip

    ! CONTACT ARRAYS
    
    integer ncontacts,refi,refj
    integer contacti(maxk),contactj(maxk),contactknt
    
    real(8) contfn(maxk),contft(maxk),contang(maxk)
    real(8) sumcomp,sumvn,maxcomp,maxvn 
    real(8) sxxt(maxd),sxyt(maxd),syxt(maxd),syyt(maxd), skxxt(maxd),skxyt(maxd),skyxt(maxd),skyyt(maxd)
    
    real(8) ftij

    !  SCALAR PARAMETERS  
    integer natoms
    real(8) dmol,sig
    real(8) xleft,xright,ytop,ybot
    real(8) xleft_out,xright_out,ytop_out,ybot_out
    real(8) tau,tauadd 
    real(8) cutlo,cuthi 
    

    !    **  SIG        AVERAGE ATOM DIAMETER (DEFINED AS 1.)
    !       ** **INPUT PARAMETERS**
    !    **  BOXX, BOXY  SIZE OF BOX IN X AND Y (RELATIVE TO AVERAGE
    !    **    ATOM DIAMETER
    !    **  SIGB    AVERAGE BOUNDARY ATOM DIAMETER 
    !    **  DMOL    SETS DISTANCE FROM MOLECULE CENTER TO CENTER OF ATOMS
    !    **        1.=ATOMS ARE TANGENT, < 1. ATOMS OVERLAP
    !    **  NATOMS    NUMBER OF ATOMS PER MOLECULE
    !    **  SKMOL    SPRING CONSTANT FOR INTRAMOLECULAR FORCES
    !    **  SKSHEAR     SPRING CONSTANT FOR LEAF SPRING
    !       **  GAMMA     DAMPING PARAMETER FOR NORMAL VELOCITY
    !       **  GAMMAMOL    DAMPING PARAMETER FOR INTERACTIONS BETWEEN ATOMS
    !       **              IN A MOLECULE 
    !       **  GAMSHEAR    DAMPING PARAMETER FOR TANGENTIAL VELOCITY
    !    **  TAU        CURRENT TIME
    !    **  TAUADD    TIME AT WHICH NEXT PARTICLE ADDED (FOR GROWING
    !            SANDPILE PROBLEM)
    !
    !    ** **CALCULATED PARAMETERS**
    !    ** XLEFT, XRIGHT, YTOP,YBOT
    !     **         CURRENT COORDINATES FOR BOX EDGES


    ! BREAKING GRAINS
    
    integer newgrn(13)
    
    real(8) teindx(MAXN)
    real(8) stot(MAXN),sx(MAXN),sy(MAXN),sxy(MAXN)
    real(8) s1(MAXN),s2(MAXN),theta(MAXN)
    real(8) u1(MAXN),u2(MAXN),u3(MAXN)
    real(8) sigmax(MAXN),sigmin(MAXN),strgth(MAXN)
    real(8) l1(MAXN),l2(MAXN),cron(MAXN)
    real(8) bvx(MAXN),bvy(MAXN)

    ! WITH THE ADDITION OF FLUID

    real(8) dxc,dyc
    real(8) phimat(MAXD+1,MAXD+1)
    real(8) perm(MAXD+1, MAXD+1)
    real(8) sitevx(MAXD+1, MAXD+1)
    real(8) sitevy(MAXD+1,MAXD+1)
    real(8) divvel(MAXD+1,MAXD+1)        
    real(8) press(MAXD+1,MAXD+1)    
    real(8) gradpx(MAXD+1,MAXD+1)
    real(8) gradpy(MAXD+1,MAXD+1)    
    real(8) solfra(MAXD+1,MAXD+1)    
    real(8) numpar(MAXD+1,MAXD+1)          
    real(8) pbot,ptop,fluxbot
    real(8) famp, fomega,pout
    real(8) permave(MAXD+1, MAXD+1),sitevxave(MAXD+1, MAXD+1)
    real(8) sitevyave(MAXD+1, MAXD+1),pressave(MAXD+1, MAXD+1)
    real(8) divvelave(MAXD+1,MAXD+1),solfraave(MAXD+1,MAXD+1)
    real(8) numparave(MAXD+1,MAXD+1)
    real(8) permavetp(MAXD+1,MAXD+1),sitevxavetp(MAXD+1,MAXD+1)
    real(8) sitevyavetp(MAXD+1,MAXD+1),pressavetp(MAXD+1,MAXD+1)
    real(8) divvelavetp(MAXD+1,MAXD+1),solfraavetp(MAXD+1,MAXD+1)
    real(8) numparavetp(MAXD+1,MAXD+1)
    real(8) press_bot(MAXD+1)
    real(8) visco
    real(8) vol_flux_x(MAXD+1,MAXD+1),vol_flux_y(MAXD+1,MAXD+1)
    real(8) diff_vol_flux_x(MAXD+1,MAXD+1)
    real(8) diff_vol_flux_y(MAXD+1,MAXD+1)
    real(8) max_flux_x,max_flux_y
    real(8) dpdttp(MAXD+1,MAXD+1)

!    integer,parameter :: MAXWORK = 9999999
!
!    real work(MAXWORK)
!
!    integer iprm(17),mgopt(4)
!    integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny
!    integer iguess,maxcy,method,nwork,lwrkqd,itero
!    integer ierror
!    integer imaxcy
!
!    real fprm(6)
!    real xa,xb,yc,yd,tolmax,relmax
!    real alfxa,gbdxa,alfxb,gbdxb,alfyc,gbdyc,alfyd,gbdyd
!    real bcxa,bcxb,bcyc,bcyd
!    real rtolmax

    contains

end module
