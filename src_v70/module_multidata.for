!##########################################################################
        module multidata
!##########################################################################
	  SAVE
        integer :: nbp,nbpmax,num_domains
        integer :: rdivmax,idom,jdom,kdom
        integer,allocatable,dimension(:) :: rdiv
        integer,allocatable,dimension(:,:) :: rdv
        integer,allocatable,dimension(:) :: dom_id,dom_ad,dom_indid
        integer,allocatable,dimension(:) :: i_unst,j_unst,k_unst
        integer,allocatable,dimension(:) :: id_unst
	  integer,allocatable,dimension(:) :: imbinblk	!Pablo
        double precision, allocatable,dimension(:,:)::xcor,ycor,zcor !Pablo

        type multidom
           integer :: inext,iprev,jnext,jprev,knext,kprev
           integer :: corprev1,corprev2,corprev3,corprev4
           integer :: cornext1,cornext2,cornext3,cornext4
           integer :: edgprev1,edgprev2,edgprev3
           integer :: edgprev4,edgprev5,edgprev6
           integer :: edgnext1,edgnext2,edgnext3
           integer :: edgnext4,edgnext5,edgnext6
           integer :: per_in,per_ip,per_jn,per_jp,per_kn,per_kp
           integer :: ttc_i,ttc_j,ttc_k,ttc_ijk,ngrid
           integer :: isp,iep,jsp,jep,ksp,kep
	   integer :: isu,ieu,jsu,jeu,ksu,keu
	   integer :: isv,iev,jsv,jev,ksv,kev
	   integer :: isw,iew,jsw,jew,ksw,kew
           integer :: nwork,nvars
           integer :: mximb
           integer :: rq_m1,rq_p1,rq_m2,rq_p2,rq_m3,rq_p3
           integer :: rq_c1m,rq_c2m,rq_c3m,rq_c4m
           integer :: rq_c1p,rq_c2p,rq_c3p,rq_c4p
           integer :: rq_e1m,rq_e2m,rq_e3m,rq_e4m,rq_e5m,rq_e6m
           integer :: rq_e1p,rq_e2p,rq_e3p,rq_e4p,rq_e5p,rq_e6p
           double precision    :: xsl,ysl,zsl,xel,yel,zel,dx,dy,dz
	     double precision,pointer,dimension(:) :: tauw
           double precision, pointer, dimension(:,:,:) :: S,So,Sm,Stm
           double precision, pointer, dimension(:,:,:) :: SUtm,SVtm,SWtm
	     double precision, pointer, dimension(:,:,:) :: sfactor
           double precision, pointer, dimension(:) :: x,y,z,xc,yc,zc
           double precision, pointer, dimension(:,:,:) :: u,v,w,p,pp
           double precision, pointer, dimension(:,:,:) :: ksgs,ksgso
           double precision, pointer, dimension(:,:,:) :: eps,epso
           double precision, pointer, dimension(:,:,:) :: T,To,Tm,Ttm
           double precision, pointer, dimension(:,:,:) :: su,ap,sup
           double precision, pointer, dimension(:,:,:) :: ae,aw,as,an
           double precision, pointer, dimension(:,:,:) :: at,ab
           double precision, pointer, dimension(:,:,:) :: um,vm,wm
           double precision, pointer, dimension(:,:,:) :: pm,ppm
           double precision, pointer, dimension(:,:,:) :: vis,vism,epsm
           double precision, pointer, dimension(:,:,:) :: uum,vvm,wwm
           double precision, pointer, dimension(:,:,:) :: uvm,uwm,vwm
           double precision, pointer, dimension(:,:,:) :: ustar,vstar
           double precision, pointer, dimension(:,:,:) :: wstar
           double precision, pointer, dimension(:,:,:) :: uo,uoo,vo
           double precision, pointer, dimension(:,:,:) :: voo,wo,woo
           double precision, pointer, dimension(:,:,:) :: stfcinf
           double precision, pointer, dimension(:,:,:) :: facp1,facp2
           double precision, pointer, dimension(:,:,:) :: facm1,facm2
           double precision, pointer, dimension (:) :: dh1,dh2,dh3 
           integer, pointer, dimension (:,:,:,:) :: ndimb,imbbodynum
           integer, pointer, dimension (:) :: faz,cntp
           integer, pointer, dimension (:,:,:) :: ntav1,ntav2
           integer, dimension (26) :: tg
	     integer, pointer, dimension (:,:,:) :: ibfactor 
           double precision, pointer, dimension(:)     :: cof
           double precision, pointer, dimension(:,:)   :: tauwe,tauww
           double precision, pointer, dimension(:,:)   :: tauws,tauwn
           double precision, pointer, dimension(:,:)   :: tauwt,tauwb
           double precision, pointer, dimension(:,:)   :: tauwe2,tauww2
           double precision, pointer, dimension(:,:)   :: tauws2,tauwn2
           double precision, pointer, dimension(:,:)   :: tauwt2,tauwb2
           double precision, pointer, dimension(:) :: sendb_m1,sendb_p1
           double precision, pointer, dimension(:) :: recvb_m1,recvb_p1
           double precision, pointer, dimension(:) :: sendb_m2,sendb_p2
           double precision, pointer, dimension(:) :: recvb_m2,recvb_p2
           double precision, pointer, dimension(:) :: sendb_m3,sendb_p3
           double precision, pointer, dimension(:) :: recvb_m3,recvb_p3
           double precision, pointer, dimension(:)::sc1m,sc1p,rc1m,rc1p
           double precision, pointer, dimension(:)::sc2m,sc2p,rc2m,rc2p
           double precision, pointer, dimension(:)::sc3m,sc3p,rc3m,rc3p
           double precision, pointer, dimension(:)::sc4m,sc4p,rc4m,rc4p
           double precision, pointer, dimension(:)::se1m,se1p,re1m,re1p
           double precision, pointer, dimension(:)::se2m,se2p,re2m,re2p
           double precision, pointer, dimension(:)::se3m,se3p,re3m,re3p
           double precision, pointer, dimension(:)::se4m,se4p,re4m,re4p
           double precision, pointer, dimension(:)::se5m,se5p,re5m,re5p
           double precision, pointer, dimension(:)::se6m,se6p,re6m,re6p
           double precision, pointer, dimension(:,:,:) ::d1,dphi_dxplus
           double precision, pointer, dimension(:,:,:) ::dphi_dyplus,
     &dphi_dzplus,dphi_dxminus,dphi_dyminus,dphi_dzminus
!============================== LSM VARIABLES =============================
           double precision, pointer, dimension(:,:,:) :: phi_init,
     &phi_new,phi_reinit,phi,dphi_dx,dphi_dy,dphi_dz,s_phi0,h_phi,
     &dens,mu,phim,abs_dphi_check
           double precision, pointer, dimension(:)     :: dens_mg
           real, pointer, dimension(:) :: sendb_m,sendb_p
           real, pointer, dimension(:) :: recvb_m,recvb_p
           integer, pointer, dimension (:) :: ijkp_lsm
	     integer :: tot
           integer :: niul,njul,nkul,nivl,njvl,nkvl,niwl,njwl,nkwl
           integer :: nipl,njpl,nkpl,nigl,njgl,nkgl
           integer :: nipl2,njpl2,nkpl2
           double precision, pointer, dimension(:,:,:) :: resmax
           double precision, pointer, dimension(:,:,:) :: resfact
           double precision, pointer, dimension(:,:,:) :: resfact1
!==========================================================================
           integer bc_west,bc_east,bc_south,bc_north,bc_bottom,bc_top 
           integer Tbc_west,Tbc_east,Tbc_south,Tbc_north
           integer Tbc_bottom,Tbc_top 
	   logical :: coarse_ng,fine_ng
        end type multidom

        type (multidom), pointer, dimension(:) :: dom

        end module multidata
!##########################################################################
