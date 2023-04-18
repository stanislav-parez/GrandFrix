# GranFrix

Fortran 90 program for granular + fluid simulations, using DEM + fluid solver.

## How to install

1. go to `src`
2. open Makefile
	1. change compiler's name (`FC=gfortran`), if necessery
	2. change execution file's name (`GranFrixrm`), if necessery
3. `make`

Exe file is in src directory.
	
**or**

1. `sh cmake.sh`

Exe file is now in bin directory.

## How to run it

`GranFrixrm < {inputFile}`

## Input file
Input file is a list of settings in the form `variable value`. Default values of variables can be found in _./src/inputReader.f90_. and _./src/mycommons.f90_.
Lines starting with `#` are interpreted as comments. For example input files see _./examples/_.
Units used in the inputs as well as output files are non-dimensional, such that diameter and mass of a mean-size grain are both 1, and Young's modulus of a grain is 1.

### Statements
####  DURATION AND OUTPUT FREQUENCY
variable | type  | default value | meaning
-| - | - | -
 title | character(len=12) | "unnamed     "  | title of a restart file
 nstep | integer |  0              | number of time steps
 iprint | integer |  huge(0)       | postscript page plotted every iprint time steps
 iprintf | integer |  huge(0)      | frequent output dumped every iprintf time steps
 iprint_restart | integer |  huge(0)      | restart printed  every iprint_restart time steps
 max_time | real(8) |  huge(0)      | maximum duration of the simulation
    
    
####  PARAMETERS FOR BUILDING A NEW SYSTE2M
variable 	| type  | default value | meaning
-| - | - | -
 irstart 	| integer |  0 			| = 0  build a new set of grains,
 |||                     			| =1 read file 'restart'
 |||                     			| =2 initialize fluid pressure with some fraction of presy
 ifill 		| integer |	-1			| =0 don't  make any new particles
 |||                  				| =-1 create a packing of grains by pouring them into the box
 seed 		| integer |  0       	| random generator seed; if 0 then generate random seed at start
 boxx		| real |  0.   			| horizontal  size of box (active for ifill=-1)
 boxy 		| real |  0.  			| (initial) vertical size of box (active for ifill=-1)
 boxfy		| real |  tiny(0.)      |  >0, final vertical size of box (active for ifill=-1) 
 phigoal 	| real |  0.2         	| target porosity (active for ifill=-1)
 frac1 		| real |  1. 			| the number fraction of the grains in population 1
 mean1 		| real |  1.  	 		| mean diameter of population 1
 std1 		| real |   1.    		| standard deviation of diameters in population 1
 logr 		| real |   0.8   		| smallest diameter of population 1
 higr 		| real |   1.2   		| largest diameter of population 1
 e1  		| real |  1.      		| Young's modulus for population 1
 mean2 		| real |  1.    		| mean diameter of population 2
 std2 		| real |  1.     		| standard deviation of diameters in population 2
 logr2 		| real |  0.8   		| smallest diameter of population 2
 higr2 		| real |  1.2   		| largest diameter of population 2
 e2         | real |	1.			| Young's modulus for population 2
 sigb 		| real |  1. 			| mean diameter of boundary grains
    
####  BOUNDARY AND INITIAL CONDITIONS 
variable 	| type  | default value | meaning
-| - | - | -
 iztau  | integer |    0   | = 0, zero out the time counter
 |||                       | other, read initial time from restart file
 izmom  | integer |    0   | = 0: zero out all velocities
 |||                       | = 1: give all particles fb(2) horizontal velocity
 |||                       | other: read velocities from the restart file
 izslip | integer |  0     | = 0: zero out slip array (shear forces between particles)
 |||                        | other: use the restart file
 ib 	| integer |   0    | boundary conition switch (array order: bottom,top,left,right wall)
 |||                       | <0 no wall (for side walls this makes it periodic)
 |||                       | =0 stationary wall  (both velocities =0.)
 |||                       | =1 no fixed velocities
 |||                       | =2 fixed shear velocity
 |||                       | =3 fixed normal velocity
 |||                       | =4 fixed shear velocity, zero normal velocity

####  APPLIED FORCES AND VELOCITIES
variable 	| type  | default value | meaning  
-| - | - | -       
 fb 	| real 	| 0     | boundary velocity on each of the 4 walls    
 presx 	| real 	|  0.   | size of external shear force
 presy 	| real	| 0.  	| size of external normal force
 gx		| real	| 0.	| body force acting in x directions
 gy 	| real	| 0.	| body force acting in y directions
    
####  PHYSICAL PARAMETERS OF GRAINS
variable 	| type  | default value | meaning  
-| - | - | -       
 contact_model | integer |  1   | 0 - linear, 1 - HERTZ-MINDLIN
 skmol | real |  1.    | spring constant for elastic repulsion between grains
 skshear | real |  0.  | spring constant for elastic repulsion between grains in the shear direction
 friction | real |  0. | ratio of shear to normal stress at which a grain begins to slide 
 gamma | real |  0.8   | damping parameter for grain collisions, characterizing inelasticity (generally between
 |||                                      | 0 and 1, 0.8 is a good number)
 gammamol | real |  0. | damping by pore fluid (not necessary even for fluid runs)
 kspring | real |  0.      | spring constant for a spring tied to top wall, used with ib(2)=5
 itmass | logical | .false.  | switch for a heavy block of mass topmass, placed on the top wall (true or false)
 topmass | real |  1.   | mass of the heavy block
    
####  FLOW CONTROL AND EFFICIENCY PARAMETERS
variable 	| type  | default value | meaning  
-| - | - | -       
 dt | real | 0.1           				 | time step factor:  
|||                                  |    maximum stable time step will be calculated, 
|||                                  |    then multiplied by this factor: should be <= 0.15  
 mx | integer |  0 		| number of grid cells in x direction for solving the pore pressure equation
 my | integer |   0   	| number of grid cells in y direction for solving the pore pressure equation
         
#### FLUID       
variable 	| type  | default value | meaning  
-| - | - | -                         
 ifluid | integer |  1          | ifluid: 1 - no fluid; 0 - with fluid             
 fstress | real |  0.   | fstress is a factor in proportion to presy determining the initial fluid pressure in the system
 fbc | integer |  2     | fluid boundary conditions at the top and bottom walls:
|||    | 1 = Dirichlet, 
|||    | 2 = Neumann with derivative = 0, 
|||    | 3 = Dirichlet with permeability test i.e. not moving the grains in response to fluid pressure.
 |||   | 4 = Dirichlet with moving grains
 |||   | 5 = Mixed boundary conditions. Constant pressure gradient at the bottom and constant pressure at top
 |||   | 6 = Sin/Cos wave - time dependent pressure on same b.c.
 |||   | 7 = Undrained BC that are implemented
 permfac | real |  1.d-5    | a factor used in Karman-Kozeny relation to scale permeability  
 acoef | real |  1.d-1      | coefficient in the pore pressure equation, 1/(fluid compressibility * grain Young's modulus)
 bcoef | real |  3.72d4     | coefficient in the pore pressure equation, 1/(fluid compressibility * viscosity)
 