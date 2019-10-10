# hydro3D
Repository for Hydro3D.
Last updated: 10/Oct/2019 by Dr Pablo Ouro

Hydro3D – version 7.0				October/2018

The version 7.0 of Hydro3D has the following updates:

- Improved modularity. Less input data in control.cin with new files for LSM: in_LSM.cin; Energy and scalar transport: in_energy.cin.

- All secondary input files have been renamed starting with "in_" to ease the finding of the files.

- LPT: all the files have been included in LPT.for for convenience. 

- New domain setup: both infodom.cin and mdmap.cin have been changed with two options to generate the domain: option 1 is the classic "by sub-domains" way, whilst option 2 is a homogeneous "whole domain" partition. When using LRM, only option 1 works. 

- LSM: there has been some updates from the work of Aristos and this is expected to change.

- actuator.for: the ALM for HAT and VAT have been updated.

- SEM: The input files are now generated at every time step which is more efficient. Hence, no need to write the number of planes to be used.

- FORTRAN legacy: We have removed two issues in the syntax that appeared when compiling the code with a different compiler than Intel, e.g. gnu or Cray. Now in module_mpi.for the tag_offset have been modified. This is volatile so there could be some issues in the communication with the simulation not starting to run. If you use LMR, need to stick to the old tags. 

- COMPILATION: when using SCW use intel compilers 2017.


*********************************************************************** 
*********************************************************************** 

Hydro3D – version 6.0				February/2018

The version 6.0 of Hydro3D gathers the newest versions of Immersed Boundary Method, Lagrangian Particle Tracking, Level-Set Method, Synthetic Eddy Method and Scalar transport equation. 

For new users or beginners, note some of the key considerations when starting to use Hydro3D:

- COMPILATION: type >> make -j 16, in order to compile the code both in your workstation or Supercomputing Wales (hereafter SCW). To recompile type >> make clean, and then >> make -j 16.  To add/remove compilation flags, see “makefile” file.

- Supercomputing Wales: there are two modules that need to be loaded in order to compile Hydro3D:
   	>> module load compiler/intel/15.0/3.187
	>> module load mpi/intel/5.1/3.210 

- TECPLOT: output files of Hydro3D are in tec360 format. In order to get the license, once you install tec360, select “Network license” and type: 10.74.7.131 and 27100 in the IP and Port numbers respectively. This has to be done at School of Engineering connected to internet. If you are going to use tec360 away from the University, once you activate the license, select in Help>Roaming the latest day available. This option will give you access to tec360 without internet connection.

- REPOSITORY: Version 6.0 of Hydro3D is available through a Git-hub repository which will be modified every time a bug is corrected. https://github.com/OuroPablo/3DFD_v60 . Need to ask for permission to access the GitHub repository to Pablo Ouro: ourobarbap@cardiff.ac.uk 

- REFERENCES: There is a list of papers in which some chunks of Hydro3D are described and validated. In this document, those referring to IBM, LSM or LPT are included. At the end there is also a list of other papers from people who has used Hydro3D and thus some of them must be acknowledge so we all can foster the impact of Hydro3D.

- ACKNOWLEDGEMENTS: it is vital to add a sentence on the support we receive when running in SCW. Hence, please you add in your journal/conference papers the following sentence:
“We acknowledge the support of the Supercomputing Wales project, which is part-funded by the European Regional Development Fund (ERDF) via Welsh Government”. 
And perhaps the SCW logo:
 

- MAIN FILES: the files fdstag.for (i.e. Finite Differences Stages) and flosol.for (Flow Solver) indicate the sequence followed by Hydro3D to allocate variables, advance the Navier-Stokes equations in time, etc.

In the following, the input files “.cin” are briefly described, and these are:
1.	control.cin: fluid solver variables and choice of additional tools.
2.	mdmap.cin: sub-domain descriptions.
3.	infodom.cin: assignation of sub-domains to CPUs.
4.	geom.cin: immersed boundary method variables.
5.	LPT.cin: lagrangian particle tracking variables.
6.	rough_info.cin: rough bed variables.

