# hydro3D
Repository for hydro3D LES code

Hydro3D – version 6.0				February/2018

The version 6.0 of Hydro3D gathers the newest versions of Immersed Boundary Method, Lagrangian Particle Tracking, Level-Set Method, Synthetic Eddy Method and Scalar transport equation. 

For new users or beginners, note some of the key considerations when starting to use Hydro3D:

- COMPILATION: type >> make -j 16, in order to compile the code both in your workstation or Supercomputing Wales (hereafter SCW). To recompile type >> make clean, and then >> make -j 16.  To add/remove compilation flags, see “makefile” file.

- Supercomputing Wales: there are two modules that need to be loaded in order to compile Hydro3D:
 	>> module load compiler/intel/15.0/3.187
	>> module load mpi/intel/5.1/3.210 

- TECPLOT: output files of Hydro3D are in tec360 format. In order to get the license, once you install tec360, select “Network license” and type: 10.74.7.131 and 27100 in the IP and Port numbers respectively. This has to be done at School of Engineering connected to internet. If you are going to use tec360 away from the University, once you activate the license, select in Help>Roaming the latest day available. This option will give you access to tec360 without internet connection.

- REPOSITORY: Version 6.0 of Hydro3D is available through a Git-hub repository which will be modified every time a bug is corrected. https://github.com/OuroPablo/Hydro3D. Need to ask for permission to access the GitHub repository to Pablo Ouro: ourobarbap@cardiff.ac.uk 

- REFERENCES: There is a list of papers in which some chunks of Hydro3D are described and validated. In this document, those referring to IBM, LSM or LPT are included. At the end there is also a list of other papers from people who has used Hydro3D and thus some of them must be acknowledge so we all can foster the impact of Hydro3D.

- ACKNOWLEDGEMENTS: it is vital to add a sentence on the support we receive when running in SCW. Hence, please you add in your journal/conference papers the following sentence:
“We acknowledge the support of the Supercomputing Wales project, which is part-funded by the European Regional Development Fund (ERDF) via Welsh Government”. Include the SCW logo in your presentations.

- MAIN FILES: the files fdstag.for (i.e. Finite Differences Stages) and flosol.for (Flow Solver) indicate the sequence followed by Hydro3D to allocate variables, advance the Navier-Stokes equations in time, etc.

In the following, the input files “.cin” are briefly described, and these are:
1.	control.cin: fluid solver variables and choice of additional tools.
2.	mdmap.cin: sub-domain descriptions.
3.	infodom.cin: assignation of sub-domains to CPUs.
4.	geom.cin: immersed boundary method variables.
5.	LPT.cin: lagrangian particle tracking variables.
6.	rough_info.cin: rough bed variables.


Log of changes made to this report are indicated at the end of this report.

 
1. CONTROL.CIN 
keyword = type of flow. Usually “channel” flow.
Ubulk = Inlet velocity.
dx,dy,dz = Resolution of the mesh. If using Local-Mesh Refinement, this corresponds to the coarser level.
1/vis = inverse of viscosity. If normalising with Reynolds number, set this as the corresponding Re.
Prandtl = Parameter used in the Energy equation.
Turbulent Schmidt number = Parameter used in the Scalar equation
Beta = Parameter used in the Energy equation.
Convection scheme = 
Diffusion scheme =
Differencing = 2nd and 4th order Central Differences, or 5th order WENO.
Solver (Poisson equation) = Pressure equation solver: 1- SIP solver, 2- Multi-grid technique.
Multigrid steps = 
Multigrid iteration scheme = 
MG maximum it per time step =
Restriction iter = 
Prolongation iter =
Fixed time step = “T” sets fixed time step, “F” for variable. 
dt = value of time step if deemed fixed. 
sweeps = number of sweeps in SIP solver
safety_factor = CFL condition, always <1.0
itime_end = number of total iterations. 
restart = “T” if the simulation restarts from a previous one. Needs tecbin files and final_ctime.dat.
Reinitmean = “T” restart the computation of time-averaged quantities. Sets iterations counter to 0.
n_out = indicates number of iterations frequency to write the output files. 
Number iter = Number of iterations to be used in the pressure solver. 
Residual (eps) = Maximum residual (error) to consider the pressure correction converged.
nswp(1-4) = parameters used in SIP solver.
Flow boundary conditions = explained below.
Friction coefficient for rough BC =
Save inflow data = write planes to be used as inlet condition in other simulations. Location of this x-plane is set in write_inflow in timesig.for as “xmapp”.
Number of inlets = Number of planes to write.
Synthetic Eddy Method = artificial turbulence imposed at the inlet 
time averaging = set it to “T” if time-averaged quantities need to be calculated.
t_start_averaging1 = time since first-order mean quantities start to be quantified. 
t_start_averaging2 = same for second-order quantities. 
Noise = introduce artificial noise in the initial velocity field. 
SGS models = sub-grid scale model. 
Local Mesh Refinement (LMR) = set which approach to be used. 
Number of Ghost cells (pl_ex) = Number of ghost cells used.
LIBM = “T” when using the immersed boundary method.
LROUGH = “T” when using rough beds. 
LSCALAR = “T” when computing Scalar transport. 
LPT = “T” when computing Lagrangian Particle Tracking.
OpenMP threads = number of OpenMP threads. It can be used in the IBM or LPT. 
LENERGY = “T” when using Energy equation.
Th, Tc = hot and cold temperatures in the Energy equation.
Temperature/Scalar BC = Boundary conditions in the Scalar/Energy equations. 

FREE SURFACE PARAMETERS
Phi_depth = depth of the free-surface initially.
L_LSM = “T” when computing the free-surface level set equation. 
re-initialisation =  “T” to use the re-initialisation of the level set function. 
Number_iter_reinit = number of iterations used during the re-initialisation.
Residual_LSM = residual to consider re-initialisation converged. 
CFL_LSM = CFL condition used during the re-initialisation
LEND = fixed level set function at the streamwise start and end of the domain. 
L_LSMbase = “T” when running the simulation with the level set function fixed, e.g. precursor simulations.
Animation = output files to generate videos. 
DensLIQUID, ,DensGAS,NuLIQUIS,NuGAS = densities and viscosities of the liquid and gas phases. 
Gravity components = magnitude of gravity.
Slope of channel = 

OUTPUT OPTIONS  
Turbulent quantities = writes instantaneous and mean velocity and Reynolds Stresses in “tecturb_XXXX.dat”. Overwritten at every nout. 
Tecinstantaneous = write velocity values in “tecinst_XXXX_it.dat” at each iteration “it”.
Planes outlet = writes instantaneous velocities at a x- or y- or z-plane set in post.dat for videos in “tecplane_XXXX_it.dat”.
Pressure = output of inst and mean pressure and viscosity values, ksgs and eps from RANS in “tecout_p_XXXX.dat”. Overwritten at every nout. 
Binary files = binary files containing relevant fluid variables in “tecbins_XXXX.dat”. Overwritten at every nout and required to restart.
Phi files = write out level set variables inst and mean phi (LSM), density and mu distribution in “tecout_phi_XXXX. Dat”. Overwritten at every nout and required to restart level set method.
Temperature file = write out energy equation variables in “tecout_T_XXXX.dat”. Overwritten at every nout and required to restart Energy equation.
Scalar transport file = write out scalar transport equation variables in “tecout_S_XXXX.dat”. Overwritten at every nout and required to restart scalar transport.

TIME SERIES 
Number of time series points = 
Points = in each line indicate: block number and I,j,k indices corresponding to the location from which fluid variables will be written in a “unst_xx.dat” file. 

BOUNDARY CONDITIONS:
Inflow conditions:
1 –   Uniform distribution
12 – 1/7th Power Law distribution both in vertical and horizontal
13 – 1/7th Power Law distribution only in horizontal
14 – 1/7th Power Law distribution only in vertical 
15 – Logarithmic distribution in the vertical directions
7   – Read inlet files from a precursor simulations stored in “inflow” folder.
77 – Mapping inflow (x-plane distance “xmapp” set in timesig.for). “inflow” folder required.
8   – Synthetic Eddy Method. Valid with velocity distributions 1, 12, 13, 14 and 15. In SEM.for the turbulence length scale is defined in SIGMA(1-3,:,:) and turbulence is considered homogeneous in the 3 directions, i.e. same turbulence intensity. If isotropic turbulence, set SIGMA equal in the 3 components. 
Outflow:
2 –  Neumann BC.
21 – Convective BC
- In bounds.for you could also set a constant outlet velocity (i.e. Dirichlet condition).
South, north, bottom and top:
3   – Slip BC
4   – No-slip BC (wall)
61 – Smooth log-law
62 – Rough log-law
63 – 1/6th law
64 – 1/7th law
65 – 1/8th law
5 - Periodic BC: valid in every of the 3 directions but must be set at opposite faces:
	West and East → when running periodically infinite channel flow
	South and North 
	Top and Bottom

2. mdmap.cin
Assign the sub-domains to CPUs.

3. infodom.cin
Specify the sub-domain dimensions in the following manner: introduce by z-layers beginning with strips in x-direction. 

4. geom.cin
This file generates geometries using the immersed boundary method, which needs to be activated in control.cin “LIMB” to make it work. 

Output forces = set if forces at the IB bodies should be written. If YES, the files will be named “F_Body_XX.dat” where XX corresponds to the number of the IB. This forces are computed in “caldrag” subroutine in IMB.for
Delta function = There are 7 delta functions in Hydro3D. These are used in the interpolation processes and adopt a different number of neighbours to do so. 2, 4 and 7 use 125 neighbours, 3 and 6 use 28 neighbours and 1 and 5 use 27. Note that the larger the number of neighbours the larger is the computational time.
Multi-Direct Forcing steps = extra IB loops to better enforce the no-slip condition at the solid boundaries. 
NIBMs = Number of immersed boundary bodies used. If more than one, copy the lines below as many time as bodies to be used.
Shape = 1 – square, 2- cylinder, 3- cube, 4- sphere, 11-cone, 12-dune, 13-hemisphere, 5-read geometry from file.
Length = Objects such as cylinders or others geometries can be extruded along the axis direction over the whole domain length or just between zini and zend. 
Centre = Centre of the IB geometry.
Radius = radius of the geometry for cylinders or spheres while characteristic side length for square or cubes.
Cmax = number of additional internal layers of Lagrangian markers.
Axis of extrusion = to be used if object if finite along the preferential direction. For instance, if reading geometry from file, this should be “-1”. 
File defining IB geometry = this file should have a first line with the number of markers to be read, and hereafter each line will have three columns defining x, y, and z coordinates of each marker. When generating this external geometry, this can be in local coordinates, i.e. centred in 0,0,0, and then using Cx, Cy, Cz to move it within the computational domain. 
Movement = set if the geometry is static or dynamic.
Reduction factor = usually leave it as 1.0. This should be used when IB geometries have a resolution larger than the fluid one (Reduction factor < 1.0) or lower resolution (RF > 1.0). It only adds stability to the method but it should be avoided whenever possible. 

– Turbine application –
Is this a turbine = set yes or no. If these are blades, then yes; if it is a hub, then no.
Axis of the turbine = specify if it is a 1-VAT, 2-HAT or 3-Actuator Line Model. All these require to read geometry from file!
Gravity centre of hydrofoils = move the airfoil gravity centre if needed. Only to be used for VATs.
Pitch Angle of hydrofoils = additional pitch angle to airfoils. Only to be used in VATs
Number of hydrofoils = number of blades for VATs only.
Rotational velocity (rad/s) = rotational speed of the turbine. Note this is prescribed and thus the turbine will spin always at this velocity.

SUBROUTINES: imb.for computes the immersed boundary method subroutine, DeltaF_MLS.for for shape functions, and shapes.for generates the geometries and also computes the new coordinates when the solid moves.

REFERENCES: To read/cite papers about the IBM and its implementation in Hydro3D:
•	Ouro P., Stoesser T. 2017. An immersed boundary-based large-eddy simulation approach to predict the performance of vertical axis tidal turbines. Computers and Fluids. 152: 74—87.
•	Ouro P, Harrold M, Stoesser T, Bromley P. 2017. Hydrodynamic loadings on a horizontal axis tidal turbine prototype. Journal of Fluids and Structures. 71: 78—95.
•	Ouro P. 2017. Large-Eddy Simulation of Tidal Turbines. PhD thesis. Cardiff University.
•	Kara M, Stoesser T, McSherry R. 2015. Calculation of fluid-structure interaction: methods, refinements, applications. Proceedings of the ICE - Engineering and Computational Mechanics. 168(2): 59—78.
•	Uhlmann M. 2005. An immersed boundary method with direct forcing for the simulation of particulate flows. Journal of Computational Physics. 209: 448-476.
•	Wang Z, Fan J, Luo K. 2008. Combined multi-direct forcing and immersed boundary method for simulating flows with moving particles. International Journal of Multiphase Flow. 34: 283—302. 


5. LPT.cin

PSIcell (T)/PSIball (F): the user may choose between 
-	PSI-cell approach (Crowe et al., 1977): the coupling is done within the cell where the particle is.
-	PSI-ball (Hu & Celik, 2008): the coupling is done within a local group of neighbour cells.
Delta function: choice of the delta function type and order of interpolation (). 
Results output: number of time steps between Lagrangian outputs (“tecout_#_pt.dat”).
Release frec: number of time steps between particle releases.
Number of parts/release: number of new particles per release.
Diameter of particles: 1) TRUE: random distribution of particles following a normal distribution whose average is given in 2); FALSE: constant diameter for all particles; 2) particle diameter in meters.
Diffusor width in cm.
LPT (F)/DF (T): Lagrangian Particle Tracking algorithm (recommended) or Direct Forcing.
Random release: if TRUE, only the first set of coordinates is read and the remaining bubbles are created in a random distribution around it, depending on the diffuser width; if FALSE, the initial set of coordinates for every bubble must be specified.
xp, yp, zp, up, vp, wp = coordinates of position and velocity for every particle if Random release=FALSE; if Random release=TRUE, only the first one is read and the rest will be created in a random field around it.


SUBROUTINES: LPT.for, DeltaF_MLS.for

REFERENCES: To read/cite papers about LPT and its implementation in Hydro3D:

•	Fraga, Stoesser, Lai, Socolofsky. 2016. Ocean Modelling. 87:27—36. 
•	Fraga, Stoesser. 2016. Journal of Geophysical Research: Oceans. 121: 3887—3904.


6. ROUGH_INFO.cin

Generation of rough bed.

SUBROUTINES: roughness_function.for

REFERENCES: 
•	Stoesser, McSherry, Fraga. 2016. Secondary currents and turbulence over a non-uniformly roughened open-channel bed. Water, 7: 4896-4913.



7. Benchmark cases:

1. Lid-driven cavity flow:

•	Description: Three-dimensional lid-driven cavity flow reported in Ouro, Fraga, Lopez-Novoa, Stoesser. Computers and Fluids (under preparation). 
•	Comparison datasets: Ghia et al. 1982, Wang et al. 2013 and Ouro et al. XXX.
•	Aims: (1) setup the code with different mdmap and infodom configuration, (2) use Tecplot to plot 2D contours and XY-lines comparing the computed results with data available.





8. List of Hydro3D references:



2018
•	Ouro, Stoesser, Ramirez. 2018. ASME Journal of Fluids Engineering. 
•	Ouro, Fraga, Viti, Angeloudis, Stoesser, Gualtieri. 2018. Environmental Fluid Mechanics.

2017
•	Liu, Stoesser, Fan, Papanicolau, Tsakiris. 2017. Computers and Fluids. 
•	McSherry, Chua, Stoesser. 2017. Journal of Hydrodynamics. 29:1—12. 
•	Ouro, Stoesser. 2017. Computers and Fluids. 152:74—87.
•	Ouro, Harrold, Stoesser, Bromley. 2017. Journal of Fluids and Structures. 71:78—95.
•	Ouro. Large-Eddy Simulation of Tidal Turbines. 2017. Phd Thesis. Cardiff University.

2016
•	Cevheri, McSherry, Stoesser. 2016. International Journal of Numerical Methods in Fluids. 82:261—285. 
•	Fraga, Stoesser, Lai, Socolofsky. 2016. Ocean Modelling. 87:27—36. 
•	Fraga, Stoesser. 2016. Journal of Geophysical Research: Oceans. 121: 3887—3904.

2015
•	Kara, Stoesser, McSherry. 2015. Proceedings of the ICE- Engineering and Computational Mechanics. 168:59—78.
•	Kara, Kara, Stoesser, Sturm. 2015. Journal of Hydraulic Engineering. 141: 04015019.
•	Kara, Stoesser, Sturm, Mulahasan. 2015. Journal of Hydraulic Research. 53: 186—195.
•	Stoesser, McSherry, Fraga. 2015. Water.

2014
•	Stoesser. 2014. Journal of Hydraulic Research. 52: 441—452.

2013
•	Kim, Stoesser, Kim. 2013. Journal of Hydraulic Research. 51: 558-568.
•	Kim, Stoesser, Kim. 2013. Applied Mathematical Modelling. 37:8029—8050.

2012
•	Kara, Stoesser, Sturm. 2012. Journal of Hydraulic Research. 50: 482—493. 

2011
•	Bomminayuni, Stoesser. 2011. Journal of Hydraulic Engineering. 11:1347—1358. 

2010
•	Kim, Kim, Kim, Stoesser. 2010. Journal of Environmental Engineering. 136:22--31. 
•	Stoesser, Kim, Diplas. 2010. Journal of Hydraulic Engineering. 12:1003--1017
•	Stoesser. 2010. Journal of Hydraulic Engineering. 136:812—819.

2008
•	Stoesser, Braun, Garcia-Villalba, Rodi. 2008. Journal of Hydraulic Engineering. 134:42-55.
•	Stoesser, Nikora. 2008. Acta Geophysica. 56:876—893.



Document updates log:

Version of the report	Author	Date	Details
1.0	Pablo Ouro	19/10/2017	•	Description of input files
•	List of references
•	Tasks to do to run Hydro3D
1.1	Pablo Ouro	21/02/2018	•	Update references
•	Benchmark cases section included.

