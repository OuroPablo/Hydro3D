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

