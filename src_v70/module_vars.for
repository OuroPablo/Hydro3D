!##########################################################################
        module vars
!##########################################################################
	  SAVE
        double precision g_dx,g_dy,g_dz,dens,Re,eps,fac,Pr,beta,Th,Tc 
        double precision qzero,qstpn,qstpn_y,qstpn_z
        double precision Sc_t,forcn,forcn_y,forcn_z
        double precision flomas,rmax,alfapr,resor,fric
        double precision ctime,dt,dtavg,dtsum,safety_factor,noise,Mdef
        double precision t_start_averaging1,t_start_averaging2,ubulk
	  double precision xst,xen,yst,yen,zst,zen,wtime_cfl,dtinit
        integer alfabc,niter,nswp(4),nsweep,iter,nsweeps,ntime
        integer ngg,nge,ngc,n_unstpt,OMP_threads,iaddinlet,ireadinlet
        integer ngrid_input,ngrd_gl,maxcy,iproln,irestr
        integer itmax,sweeps,numfile,numfile1,numfile2,numfile4,numfile5
        integer conv_sch,diff_sch,sgs_model,solver,mg_itrsch
        integer bc_w,bc_e,bc_s,bc_n,bc_b,bc_t
        integer Tbc_w,Tbc_e,Tbc_s,Tbc_n,Tbc_b,Tbc_t
        integer itime,itime_start,itime_end,n_out,ITMAX_PI
        integer ipref,jpref,kpref,prefdom
	  integer pl,pl_ex,differencing,LMR,normal_inter,order
        character*80 keyword,L_n
	  logical :: LRESTART,LIMB,SGS,PERIODIC,LENERGY,LROUGH
	  logical :: pressureforce,pressureforce_y,pressureforce_z
	  logical :: time_averaging,reinitmean
	  logical :: LPT,save_inflow,read_inflow,L_dt,LSCALAR
	  logical :: LTECP,LTECBIN,LTURB,LINST,LPLAN,variTS,LPAR,LLSM
!====================   SEM variables read in control.cin
	  INTEGER:: UPROF
	  DOUBLE PRECISION :: TI_SEM
!============================== LSM VARIABLES =============================
        double precision reldif_LSM,length,densl,densg,nul,nug,mul,mug
        double precision cfl_lsm,uprev,grx,gry,grz,slope
        double precision densip12,densjp12,denskp12
	  double precision densim12,densjm12,denskm12
        double precision muip12,mujp12,mukp12,muim12,mujm12,mukm12
        integer ntime_reinit,ngrid,numfile3
	  logical :: L_LSM,REINIT,L_LSMbase,L_LSMinit,LENDS
	  logical :: L_anim_phi,L_anim_grd
!=========================================================================
        end module
!##########################################################################
      module vars_pt
!##########################################################################
	  SAVE
      integer np,tsteps_pt,tsnr,cnt_pt,np_loc,ptnr
	logical DF,PSIcell,random,random_dp	
	integer,allocatable,dimension(:)::	ptsinproc
!	integer,allocatable,dimension(:)::  ipux,jpvy,kpwz
	integer,allocatable,dimension(:)::  id
!	double precision, allocatable, dimension(:)::  dh_acumu,dh_acumv
!	double precision, allocatable, dimension(:)::  dh_acumw
      double precision, allocatable, dimension(:):: xp_pt,yp_pt,zp_pt
      double precision, allocatable, dimension(:):: xp_loc,yp_loc,zp_loc
      double precision, allocatable, dimension(:):: uop_pt,vop_pt,wop_pt
      double precision, allocatable, dimension(:):: uop_loc,vop_loc
      double precision, allocatable, dimension(:):: wop_loc,dp_loc
	double precision, allocatable, dimension(:):: Fu,Fv,Fw
	double precision, allocatable, dimension(:):: Fpu,Fpv,Fpw
!	double precision, allocatable, dimension(:):: Fsu,Fau,Fdu,Flu
!	double precision, allocatable, dimension(:):: Fsv,Fav,Fdv,Flv
!	double precision, allocatable, dimension(:):: Fgw,Fsw,Faw,Fdw,Flw
      double precision, allocatable, dimension(:):: xpold,ypold,zpold
      double precision, allocatable, dimension(:):: uopold,vopold,wopold
      double precision, allocatable, dimension(:):: dp_pt,dp_old
	double precision xp,yp,zp,uop,vop,wop,div,Dp,sigma

      end module
