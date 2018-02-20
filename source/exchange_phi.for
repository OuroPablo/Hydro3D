!##########################################################################
        subroutine  exchange_phi(op,nly)
!##########################################################################
        use multidata
        use mpi
        use vars
        implicit none
        integer :: i,j,k,ijk,ib,op,nly,chck_f,chck_c,no1,no2
        integer :: is,ie,js,je,ks,ke,ief,jef,kef,iec,jec,kec
        integer :: isf,jsf,ksf,isc,jsc,ksc
        integer :: st1,st2,st3,pl1,pl2,pl3,pll,nn
        integer :: iepr,jepr,kepr,ispr,jspr,kspr
        integer :: ieprf,jeprf,keprf,isprf,jsprf,ksprf
        integer :: ieprc,jeprc,keprc,isprc,jsprc,ksprc
        integer :: i1,i2,j1,j2,k1,k2,ii,jj,kk
        integer :: tag,ta,ly
        integer :: tsend,trecv
        integer :: ni,nj,nk,nic,nif,njc,njf,nkc,nkf
        integer :: nix,njx,nkx,nixc,njxc,nkxc,nixf,njxf,nkxf
        double precision, allocatable,dimension(:,:,:,:) :: fif
        double precision, allocatable,dimension(:,:,:,:) :: fic
        double precision, pointer, dimension(:,:,:) :: fi
        double precision, pointer, dimension(:) :: sbuf

        pll=(pl+1)*(pl+1)

        MPI_FLT   = MPI_DOUBLE_PRECISION

        chck_f=0; chck_c=0
        do ib=1,nbp
           if(dom(ib)%fine_ng) chck_f=1
           if(dom(ib)%coarse_ng) chck_c=1
        end do

        nix=0; njx=0; nkx=0
        nixc=0; njxc=0; nkxc=0
        nixf=0; njxf=0; nkxf=0
        do ib=1,nbp
           ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k 
           nif=2*(ni-2*pl)+2*pl
           njf=2*(nj-2*pl)+2*pl
           nkf=2*(nk-2*pl)+2*pl
           if(rdiv(dom_id(ib)).eq.1) then
              nic=ni; njc=nj; nkc=nk
           else
              nic=int((ni-2*pl)/2)+2*pl
              njc=int((nj-2*pl)/2)+2*pl 
              nkc=int((nk-2*pl)/2)+2*pl
           end if
           nix=max(ni,nix); njx=max(nj,njx); nkx=max(nk,nkx)
           nixc=max(nic,nixc); njxc=max(njc,njxc); nkxc=max(nkc,nkxc)
           nixf=max(nif,nixf); njxf=max(njf,njxf); nkxf=max(nkf,nkxf)
        end do

        if(chck_c.eq.1) allocate(fic(nbp,nixc,njxc,nkxc))
        if(chck_f.eq.1) allocate(fif(nbp,nixf,njxf,nkxf))    

        do ib=1,nbp
           ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k 
           select case (op)
              case (14) 
                 fi => dom(ib)%phi
              case (15)
                 fi => dom(ib)%phi_reinit
              case (16) 
                 fi => dom(ib)%phi_new
              case (17) 
                 fi => dom(ib)%phi_init
              case (18) 
                 fi => dom(ib)%dens
              case (19) 
                 fi => dom(ib)%mu
           end select

        if(dom(ib)%fine_ng) then

!        if(pl.eq.1) then
!           no1=pl; no2=-pl
!        else
!           no1=pl-1; no2=1-pl
!        end if
!
!        do k=no1,nk+no2; do j=no1,nj+no2; do i=no1,ni+no2
!           ii=2*i-pl; jj=2*j-pl; kk=2*k-pl
!           fif(ib,ii,jj,kk)=
!     & (27.0d0*fi(i,j,k)+9.0d0*fi(i,j+1,k)+9.0d0*fi(i,j,k+1)+
!     & 3.0d0*fi(i,j+1,k+1)+9.0d0*fi(i+1,j,k)+
!     & 3.0d0*fi(i+1,j+1,k)+3.0d0*fi(i+1,j,k+1)+
!     & 1.0d0*fi(i+1,j+1,k+1))/64.0d0
!           ii=2*i-pl+1; jj=2*j-pl; kk=2*k-pl
!           fif(ib,ii,jj,kk)=
!     & (9.0d0*fi(i,j,k)+3.0d0*fi(i,j+1,k)+3.0d0*fi(i,j,k+1)+
!     & 1.0d0*fi(i,j+1,k+1)+27.0d0*fi(i+1,j,k)+
!     & 9.0d0*fi(i+1,j+1,k)+9.0d0*fi(i+1,j,k+1)+
!     & 3.0d0*fi(i+1,j+1,k+1))/64.0d0
!           ii=2*i-pl; jj=2*j-pl+1; kk=2*k-pl
!           fif(ib,ii,jj,kk)=
!     & (9.0d0*fi(i,j,k)+27.0d0*fi(i,j+1,k)+3.0d0*fi(i,j,k+1)+
!     & 9.0d0*fi(i,j+1,k+1)+3.0d0*fi(i+1,j,k)+
!     & 9.0d0*fi(i+1,j+1,k)+1.0d0*fi(i+1,j,k+1)+
!     & 3.0d0*fi(i+1,j+1,k+1))/64.0d0
!           ii=2*i-pl+1; jj=2*j-pl+1; kk=2*k-pl
!           fif(ib,ii,jj,kk)=
!     & (3.0d0*fi(i,j,k)+9.0d0*fi(i,j+1,k)+1.0d0*fi(i,j,k+1)+
!     & 3.0d0*fi(i,j+1,k+1)+9.0d0*fi(i+1,j,k)+
!     & 27.0d0*fi(i+1,j+1,k)+3.0d0*fi(i+1,j,k+1)+
!     & 9.0d0*fi(i+1,j+1,k+1))/64.0d0
!           ii=2*i-pl; jj=2*j-pl; kk=2*k-pl+1
!           fif(ib,ii,jj,kk)=
!     & (9.0d0*fi(i,j,k)+3.0d0*fi(i,j+1,k)+27.0d0*fi(i,j,k+1)+
!     & 9.0d0*fi(i,j+1,k+1)+3.0d0*fi(i+1,j,k)+
!     & 1.0d0*fi(i+1,j+1,k)+9.0d0*fi(i+1,j,k+1)+
!     & 3.0d0*fi(i+1,j+1,k+1))/64.0d0
!           ii=2*i-pl+1; jj=2*j-pl; kk=2*k-pl+1
!           fif(ib,ii,jj,kk)=
!     & (3.0d0*fi(i,j,k)+1.0d0*fi(i,j+1,k)+9.0d0*fi(i,j,k+1)+
!     & 3.0d0*fi(i,j+1,k+1)+9.0d0*fi(i+1,j,k)+
!     & 3.0d0*fi(i+1,j+1,k)+27.0d0*fi(i+1,j,k+1)+
!     & 9.0d0*fi(i+1,j+1,k+1))/64.0d0
!           ii=2*i-pl; jj=2*j-pl+1; kk=2*k-pl+1
!           fif(ib,ii,jj,kk)=
!     & (3.0d0*fi(i,j,k)+9.0d0*fi(i,j+1,k)+9.0d0*fi(i,j,k+1)+
!     & 27.0d0*fi(i,j+1,k+1)+1.0d0*fi(i+1,j,k)+
!     & 3.0d0*fi(i+1,j+1,k)+3.0d0*fi(i+1,j,k+1)+
!     & 9.0d0*fi(i+1,j+1,k+1))/64.0d0
!           ii=2*i-pl+1; jj=2*j-pl+1; kk=2*k-pl+1
!           fif(ib,ii,jj,kk)=
!     & (1.0d0*fi(i,j,k)+3.0d0*fi(i,j+1,k)+3.0d0*fi(i,j,k+1)+
!     & 9.0d0*fi(i,j+1,k+1)+3.0d0*fi(i+1,j,k)+
!     & 9.0d0*fi(i+1,j+1,k)+9.0d0*fi(i+1,j,k+1)+
!     & 27.0d0*fi(i+1,j+1,k+1))/64.0d0
!
!        end do; end do; end do; end if

        if(pl.eq.1) then
           no1=pl+1; no2=-pl
        else
           no1=pl; no2=1-pl
        end if

        do k=no1,nk+no2; do j=no1,nj+no2; do i=no1,ni+no2

        kk=2*k-pl-1; ii=2*i-pl-1; jj=2*j-pl-1; !(i-1/4,j-1/4,k-1/4) 
        fif(ib,ii,jj,kk)=(125.0*fi(i-1,j-1,k-1)+750.0*fi(i,j-1,k-1)-
     &75.0*fi(i+1,j-1,k-1)+750.0*fi(i-1,j,k-1)+4500.0*fi(i,j,k-1)-
     &450.0*fi(i+1,j,k-1)-75.0*fi(i-1,j+1,k-1)-450.0*fi(i,j+1,k-1)+
     &45.0*fi(i+1,j+1,k-1)+750.0*fi(i-1,j-1,k)+4500.0*fi(i,j-1,k)-
     &450.0*fi(i+1,j-1,k)+4500.0*fi(i-1,j,k)+27000.0*fi(i,j,k)-
     &2700.0*fi(i+1,j,k)-450.0*fi(i-1,j+1,k)-2700.0*fi(i,j+1,k)+
     &270.0*fi(i+1,j+1,k)-75.0*fi(i-1,j-1,k+1)-450.0*fi(i,j-1,k+1)+
     &45.0*fi(i+1,j-1,k+1)-450.0*fi(i-1,j,k+1)-2700.0*fi(i,j,k+1)+
     &270.0*fi(i+1,j,k+1)+45.0*fi(i-1,j+1,k+1)+270.0*fi(i,j+1,k+1)-
     &27.0*fi(i+1,j+1,k+1))/32768.0

        kk=2*k-pl-1; ii=2*i-pl; jj=2*j-pl-1;   !(i+1/4,j-1/4,k-1/4)
        fif(ib,ii,jj,kk)=(-75.0*fi(i-1,j-1,k-1)+750.0*fi(i,j-1,k-1)+
     &125.0*fi(i+1,j-1,k-1)-450.0*fi(i-1,j,k-1)+4500.0*fi(i,j,k-1)+
     &750.0*fi(i+1,j,k-1)+45.0*fi(i-1,j+1,k-1)-450.0*fi(i,j+1,k-1)-
     &75.0*fi(i+1,j+1,k-1)-450.0*fi(i-1,j-1,k)+4500.0*fi(i,j-1,k)+
     &750.0*fi(i+1,j-1,k)-2700.0*fi(i-1,j,k)+27000.0*fi(i,j,k)+
     &4500.0*fi(i+1,j,k)+270.0*fi(i-1,j+1,k)-2700.0*fi(i,j+1,k)-
     &450.0*fi(i+1,j+1,k)+45.0*fi(i-1,j-1,k+1)-450.0*fi(i,j-1,k+1)-
     &75.0*fi(i+1,j-1,k+1)+270.0*fi(i-1,j,k+1)-2700.0*fi(i,j,k+1)-
     &450.0*fi(i+1,j,k+1)-27.0*fi(i-1,j+1,k+1)+270.0*fi(i,j+1,k+1)+
     &45.0*fi(i+1,j+1,k+1))/32768.0

        kk=2*k-pl-1; ii=2*i-pl-1; jj=2*j-pl;   !(i-1/4,j+1/4,k-1/4)
        fif(ib,ii,jj,kk)=(-75.0*fi(i-1,j-1,k-1)-450.0*fi(i,j-1,k-1)+
     &45.0*fi(i+1,j-1,k-1)+750.0*fi(i-1,j,k-1)+4500.0*fi(i,j,k-1)-
     &450.0*fi(i+1,j,k-1)+125.0*fi(i-1,j+1,k-1)+750.0*fi(i,j+1,k-1)-
     &75.0*fi(i+1,j+1,k-1)-450.0*fi(i-1,j-1,k)-2700.0*fi(i,j-1,k)+
     &270.0*fi(i+1,j-1,k)+4500.0*fi(i-1,j,k)+27000.0*fi(i,j,k)-
     &2700.0*fi(i+1,j,k)+750.0*fi(i-1,j+1,k)+4500.0*fi(i,j+1,k)-
     &450.0*fi(i+1,j+1,k)+45.0*fi(i-1,j-1,k+1)+270.0*fi(i,j-1,k+1)-
     &27.0*fi(i+1,j-1,k+1)-450.0*fi(i-1,j,k+1)-2700.0*fi(i,j,k+1)+
     &270.0*fi(i+1,j,k+1)-75.0*fi(i-1,j+1,k+1)-450.0*fi(i,j+1,k+1)+
     &45.0*fi(i+1,j+1,k+1))/32768.0

        kk=2*k-pl-1; ii=2*i-pl; jj=2*j-pl;     !(i+1/4,j+1/4,k-1/4)
        fif(ib,ii,jj,kk)=(45.0*fi(i-1,j-1,k-1)-450.0*fi(i,j-1,k-1)-
     &75.0*fi(i+1,j-1,k-1)-450.0*fi(i-1,j,k-1)+4500.0*fi(i,j,k-1)+
     &750.0*fi(i+1,j,k-1)-75.0*fi(i-1,j+1,k-1)+750.0*fi(i,j+1,k-1)+
     &125.0*fi(i+1,j+1,k-1)+270.0*fi(i-1,j-1,k)-2700.0*fi(i,j-1,k)-
     &450.0*fi(i+1,j-1,k)-2700.0*fi(i-1,j,k)+27000.0*fi(i,j,k)+
     &4500.0*fi(i+1,j,k)-450.0*fi(i-1,j+1,k)+4500.0*fi(i,j+1,k)+
     &750.0*fi(i+1,j+1,k)-27.0*fi(i-1,j-1,k+1)+270.0*fi(i,j-1,k+1)+
     &45.0*fi(i+1,j-1,k+1)+270.0*fi(i-1,j,k+1)-2700.0*fi(i,j,k+1)-
     &450.0*fi(i+1,j,k+1)+45.0*fi(i-1,j+1,k+1)-450.0*fi(i,j+1,k+1)-
     &75.0*fi(i+1,j+1,k+1))/32768.0

        kk=2*k-pl; ii=2*i-pl-1; jj=2*j-pl-1;   !(i-1/4,j-1/4,k+1/4)
        fif(ib,ii,jj,kk)=(-75.0*fi(i-1,j-1,k-1)-450.0*fi(i,j-1,k-1)+
     &45.0*fi(i+1,j-1,k-1)-450.0*fi(i-1,j,k-1)-2700.0*fi(i,j,k-1)+
     &270.0*fi(i+1,j,k-1)+45.0*fi(i-1,j+1,k-1)+270.0*fi(i,j+1,k-1)-
     &27.0*fi(i+1,j+1,k-1)+750.0*fi(i-1,j-1,k)+4500.0*fi(i,j-1,k)-
     &450.0*fi(i+1,j-1,k)+4500.0*fi(i-1,j,k)+27000.0*fi(i,j,k)-
     &2700.0*fi(i+1,j,k)-450.0*fi(i-1,j+1,k)-2700.0*fi(i,j+1,k)+
     &270.0*fi(i+1,j+1,k)+125.0*fi(i-1,j-1,k+1)+750.0*fi(i,j-1,k+1)-
     &75.0*fi(i+1,j-1,k+1)+750.0*fi(i-1,j,k+1)+4500.0*fi(i,j,k+1)-
     &450.0*fi(i+1,j,k+1)-75.0*fi(i-1,j+1,k+1)-450.0*fi(i,j+1,k+1)+
     &45.0*fi(i+1,j+1,k+1))/32768.0

        kk=2*k-pl; ii=2*i-pl; jj=2*j-pl-1;     !(i+1/4,j-1/4,k+1/4)
        fif(ib,ii,jj,kk)=(45.0*fi(i-1,j-1,k-1)-450.0*fi(i,j-1,k-1)-
     &75.0*fi(i+1,j-1,k-1)+270.0*fi(i-1,j,k-1)-2700.0*fi(i,j,k-1)-
     &450.0*fi(i+1,j,k-1)-27.0*fi(i-1,j+1,k-1)+270.0*fi(i,j+1,k-1)+
     &45.0*fi(i+1,j+1,k-1)-450.0*fi(i-1,j-1,k)+4500.0*fi(i,j-1,k)+
     &750.0*fi(i+1,j-1,k)-2700.0*fi(i-1,j,k)+27000.0*fi(i,j,k)+
     &4500.0*fi(i+1,j,k)+270.0*fi(i-1,j+1,k)-2700.0*fi(i,j+1,k)-
     &450.0*fi(i+1,j+1,k)-75.0*fi(i-1,j-1,k+1)+750.0*fi(i,j-1,k+1)+
     &125.0*fi(i+1,j-1,k+1)-450.0*fi(i-1,j,k+1)+4500.0*fi(i,j,k+1)+
     &750.0*fi(i+1,j,k+1)+45.0*fi(i-1,j+1,k+1)-450.0*fi(i,j+1,k+1)-
     &75.0*fi(i+1,j+1,k+1))/32768.0

        kk=2*k-pl; ii=2*i-pl-1; jj=2*j-pl;     !(i-1/4,j+1/4,k+1/4)
        fif(ib,ii,jj,kk)=(45.0*fi(i-1,j-1,k-1)+270.0*fi(i,j-1,k-1)-
     &27.0*fi(i+1,j-1,k-1)-450.0*fi(i-1,j,k-1)-2700.0*fi(i,j,k-1)+
     &270.0*fi(i+1,j,k-1)-75.0*fi(i-1,j+1,k-1)-450.0*fi(i,j+1,k-1)+
     &45.0*fi(i+1,j+1,k-1)-450.0*fi(i-1,j-1,k)-2700.0*fi(i,j-1,k)+
     &270.0*fi(i+1,j-1,k)+4500.0*fi(i-1,j,k)+27000.0*fi(i,j,k)-
     &2700.0*fi(i+1,j,k)+750.0*fi(i-1,j+1,k)+4500.0*fi(i,j+1,k)-
     &450.0*fi(i+1,j+1,k)-75.0*fi(i-1,j-1,k+1)-450.0*fi(i,j-1,k+1)+
     &45.0*fi(i+1,j-1,k+1)+750.0*fi(i-1,j,k+1)+4500.0*fi(i,j,k+1)-
     &450.0*fi(i+1,j,k+1)+125.0*fi(i-1,j+1,k+1)+750.0*fi(i,j+1,k+1)-
     &75.0*fi(i+1,j+1,k+1))/32768.0

        kk=2*k-pl; ii=2*i-pl; jj=2*j-pl;       !(i+1/4,j+1/4,k+1/4)
        fif(ib,ii,jj,kk)=(-27.0*fi(i-1,j-1,k-1)+270.0*fi(i,j-1,k-1)+
     &45.0*fi(i+1,j-1,k-1)+270.0*fi(i-1,j,k-1)-2700.0*fi(i,j,k-1)-
     &450.0*fi(i+1,j,k-1)+45.0*fi(i-1,j+1,k-1)-450.0*fi(i,j+1,k-1)-
     &75.0*fi(i+1,j+1,k-1)+270.0*fi(i-1,j-1,k)-2700.0*fi(i,j-1,k)-
     &450.0*fi(i+1,j-1,k)-2700.0*fi(i-1,j,k)+27000.0*fi(i,j,k)+
     &4500.0*fi(i+1,j,k)-450.0*fi(i-1,j+1,k)+4500.0*fi(i,j+1,k)+
     &750.0*fi(i+1,j+1,k)+45.0*fi(i-1,j-1,k+1)-450.0*fi(i,j-1,k+1)-
     &75.0*fi(i+1,j-1,k+1)-450.0*fi(i-1,j,k+1)+4500.0*fi(i,j,k+1)+
     &750.0*fi(i+1,j,k+1)-75.0*fi(i-1,j+1,k+1)+750.0*fi(i,j+1,k+1)+
     &125.0*fi(i+1,j+1,k+1))/32768.0

        end do; end do; end do; end if

        if(dom(ib)%coarse_ng) then
        nic=int((ni-2*pl)/2)+2*pl
        njc=int((nj-2*pl)/2)+2*pl 
        nkc=int((nk-2*pl)/2)+2*pl

        if(pl.eq.1) then
           no1=pl+1; no2=-pl
        else
           no1=pl; no2=1-pl
        end if

        do k=no1,nkc+no2; do j=no1,njc+no2; do i=no1,nic+no2
           i1=2*i-pl-1; i2=2*i-pl; j1=2*j-pl-1; j2=2*j-pl
           k1=2*k-pl-1; k2=2*k-pl
           fic(ib,i,j,k)=0.125*(fi(i1,j1,k1)+fi(i1,j1,k2)+
     & fi(i1,j2,k1)+fi(i1,j2,k2)+fi(i2,j1,k1)+fi(i2,j1,k2)+
     & fi(i2,j2,k1)+fi(i2,j2,k2))
        end do; end do; end do; end if

        end do
!==========================================================================
!==========================================================================
        do ly=0,nly

           do ib=1,nbp
              ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k 
              select case (op)
              case (14) 
                 fi => dom(ib)%phi
              case (15)
                 fi => dom(ib)%phi_reinit
              case (16) 
                 fi => dom(ib)%phi_new
              case (17) 
                 fi => dom(ib)%phi_init
              case (18) 
                 fi => dom(ib)%dens
              case (19) 
                 fi => dom(ib)%mu
              end select
              is=dom(ib)%isp; ie=dom(ib)%iep
              js=dom(ib)%jsp; je=dom(ib)%jep
              ks=dom(ib)%ksp; ke=dom(ib)%kep

              if(dom(ib)%fine_ng) then
                 nif=2*(ni-2*pl)+2*pl
                 njf=2*(nj-2*pl)+2*pl
                 nkf=2*(nk-2*pl)+2*pl

                 if(LMR.eq.2) then
                    isf=pl+2;     jsf=pl+2;     ksf=pl+2
                    ief=nif-pl-1; jef=njf-pl-1; kef=nkf-pl-1
                 else
                    isf=pl+1;     jsf=pl+1;     ksf=pl+1
                    ief=nif-pl;   jef=njf-pl;   kef=nkf-pl
                 end if

                 isprf=pl+1; ieprf=nif-pl
                 jsprf=pl+1; jeprf=njf-pl
                 ksprf=pl+1; keprf=nkf-pl
              end if

              if(dom(ib)%coarse_ng) then
                 nic=int((ni-2*pl)/2)+2*pl
                 njc=int((nj-2*pl)/2)+2*pl 
                 nkc=int((nk-2*pl)/2)+2*pl
                 isc=pl+1;   jsc=pl+1;   ksc=pl+1
                 iec=nic-pl; jec=njc-pl; kec=nkc-pl
                 isprc=pl+1; ieprc=nic-pl
                 jsprc=pl+1; jeprc=njc-pl
                 ksprc=pl+1; keprc=nkc-pl
              end if

!..........................................................................
!=== Previous Neighbor  ===> 
!..........................................................................
           if (dom(ib)%iprev.ge.0)  then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%iprev)) then
              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%iprev)) then
                 sbuf => dom(dom_indid(dom(ib)%iprev))%recvb_p1
              else
                 sbuf => dom(ib) % sendb_m1
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%iprev)) then
              tsend=njc*nkc; trecv=nj*nk
              do k=1,nkc; do j=1,njc; ijk=(k-1)*njc+j
              sbuf(ijk)=fic(ib,isc+ly,j,k)
              end do; end do

              else !(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%iprev))
              tsend=njf*nkf; trecv=nj*nk
              do k=1,nkf; do j=1,njf; ijk=(k-1)*njf+j
              sbuf(ijk)=fif(ib,isf+ly,j,k)
              end do; end do
              end if

              if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%iprev)) then
                 tag=1*10**5+dom(ib)%iprev
                 ta=2
                 call MPI_IRECV  (dom(ib)%recvb_m1(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%iprev),dom(ib)%tg(ta),
     &MPI_COMM_WORLD,dom(ib)%rq_m1,ierr)
                 call MPI_SEND (dom(ib)%sendb_m1(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%iprev),tag,MPI_COMM_WORLD,ierr)

              end if

           end if
           end if

           if (dom(ib)%jprev.ge.0)  then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%jprev)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%jprev)) then
                 sbuf => dom(dom_indid(dom(ib)%jprev))%recvb_p2
              else
                 sbuf => dom(ib) % sendb_m2
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%jprev)) then
              tsend=nic*nkc; trecv=ni*nk
              do k=1,nkc; do i=1,nic; ijk=(k-1)*nic+i
              sbuf(ijk)=fic(ib,i,jsc+ly,k)
              end do; end do

              else !(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%jprev))
              tsend=nif*nkf; trecv=ni*nk
              do k=1,nkf; do i=1,nif; ijk=(k-1)*nif+i
              sbuf(ijk)=fif(ib,i,jsf+ly,k)
              end do; end do
              end if

              if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%jprev)) then
                 tag=3*10**5+dom(ib)%jprev
                 ta=4
                 call MPI_IRECV  (dom(ib)%recvb_m2(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%jprev),dom(ib)%tg(ta),
     &MPI_COMM_WORLD,dom(ib)%rq_m2,ierr)
                 call MPI_SEND (dom(ib)%sendb_m2(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%jprev),tag,MPI_COMM_WORLD,ierr)

              end if

           end if
           end if

           if (dom(ib)%kprev.ge.0)  then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%kprev)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%kprev)) then
                 sbuf => dom(dom_indid(dom(ib)%kprev))%recvb_p3
              else
                 sbuf => dom(ib) % sendb_m3
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%kprev)) then
              tsend=nic*njc; trecv=ni*nj
              do i=1,nic; do j=1,njc; ijk=(i-1)*njc+j
              sbuf(ijk)=fic(ib,i,j,ksc+ly)
              end do; end do

              else !(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%kprev))
              tsend=nif*njf; trecv=ni*nj
              do i=1,nif; do j=1,njf; ijk=(i-1)*njf+j
              sbuf(ijk)=fif(ib,i,j,ksf+ly)
              end do; end do
              end if

              if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%kprev)) then
                 tag=5*10**5+dom(ib)%kprev
                 ta=6
                 call MPI_IRECV  (dom(ib)%recvb_m3(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%kprev),dom(ib)%tg(ta),
     &MPI_COMM_WORLD,dom(ib)%rq_m3,ierr)
                 call MPI_SEND (dom(ib)%sendb_m3(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%kprev),tag,MPI_COMM_WORLD,ierr)

              end if

           end if
           end if
!..........................................................................
!=== Next Neighbor ===> 
!..........................................................................
           if (dom(ib)%inext.ge.0)  then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%inext)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%inext)) then
                 sbuf => dom(dom_indid(dom(ib)%inext))%recvb_m1
              else
                 sbuf => dom(ib) % sendb_p1
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%inext)) then
              tsend=njc*nkc; trecv=nj*nk
              do k=1,nkc; do j=1,njc; ijk=(k-1)*njc+j
              sbuf(ijk)=fic(ib,iec-ly,j,k)
              end do; end do

              else !(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%inext))
              tsend=njf*nkf; trecv=nj*nk
              do k=1,nkf; do j=1,njf; ijk=(k-1)*njf+j
              sbuf(ijk)=fif(ib,ief-ly,j,k)
              end do; end do
              end if

              if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%inext)) then
                 tag=2*10**5+dom(ib)%inext
                 ta=1
                 call MPI_IRECV  (dom(ib)%recvb_p1(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%inext),dom(ib)%tg(ta),
     &MPI_COMM_WORLD,dom(ib)%rq_p1,ierr)
                 call MPI_SEND (dom(ib)%sendb_p1(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%inext),tag,MPI_COMM_WORLD,ierr)

              end if

           end if
           end if

           if (dom(ib)%jnext.ge.0)  then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%jnext)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%jnext)) then
                 sbuf => dom(dom_indid(dom(ib)%jnext))%recvb_m2
              else
                 sbuf => dom(ib) % sendb_p2
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%jnext)) then
              tsend=nic*nkc; trecv=ni*nk
              do k=1,nkc; do i=1,nic; ijk=(k-1)*nic+i;
              sbuf(ijk)=fic(ib,i,jec-ly,k)
              end do; end do

              else !(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%jnext))
              tsend=nif*nkf; trecv=ni*nk
              do k=1,nkf; do i=1,nif; ijk=(k-1)*nif+i;
              sbuf(ijk)=fif(ib,i,jef-ly,k)
              end do; end do
              end if

              if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%jnext)) then
                 tag=4*10**5+dom(ib)%jnext
                 ta=3
                 call MPI_IRECV  (dom(ib)%recvb_p2(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%jnext),dom(ib)%tg(ta),
     &MPI_COMM_WORLD,dom(ib)%rq_p2,ierr)
                 call MPI_SEND (dom(ib)%sendb_p2(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%jnext),tag,MPI_COMM_WORLD,ierr)

              end if

           end if
           end if

           if (dom(ib)%knext.ge.0)  then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%knext)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%knext)) then
                 sbuf => dom(dom_indid(dom(ib)%knext))%recvb_m3
              else
                 sbuf => dom(ib) % sendb_p3
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%knext)) then
              tsend=nic*njc; trecv=ni*nj
              do i=1,nic; do j=1,njc; ijk=(i-1)*njc+j;
              sbuf(ijk)=fic(ib,i,j,kec-ly)
              end do; end do

              else !(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%knext))
              tsend=nif*njf; trecv=ni*nj
              do i=1,nif; do j=1,njf; ijk=(i-1)*njf+j;
              sbuf(ijk)=fif(ib,i,j,kef-ly)
              end do; end do
              end if

              if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%knext)) then
                 tag=6*10**5+dom(ib)%knext
                 ta=5
                 call MPI_IRECV  (dom(ib)%recvb_p3(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%knext),dom(ib)%tg(ta),
     &MPI_COMM_WORLD,dom(ib)%rq_p3,ierr)
                 call MPI_SEND (dom(ib)%sendb_p3(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%knext),tag,MPI_COMM_WORLD,ierr)

              end if

           end if
           end if

!======================================================================
        if (ly.eq.0)  then
!======================================================================
!=====> previous cor #1
        if (dom(ib)%corprev1.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%corprev1)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%corprev1)) then
                 sbuf => dom(dom_indid(dom(ib)%corprev1))%rc1p
              else
                 sbuf => dom(ib) % sc1m
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%corprev1)) then
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=isprc-1+pl1; j=jsprc-1+pl2; k=ksprc-1+pl3
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=isprf-1+pl1; j=jsprf-1+pl2; k=ksprf-1+pl3
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev1)) then
                 tag=7*10**5+dom(ib)%corprev1
               call MPI_IRECV  (dom(ib)%rc1m(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%corprev1),dom(ib)%tg(8),
     &MPI_COMM_WORLD,dom(ib)%rq_c1m,ierr)
                 call MPI_SEND (dom(ib)%sc1m(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%corprev1),tag,MPI_COMM_WORLD,ierr)
  
              end if
           end if
        end if
!=====> previous cor #2
        if (dom(ib)%corprev2.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%corprev2)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%corprev2)) then
                 sbuf => dom(dom_indid(dom(ib)%corprev2))%rc2p
              else
                 sbuf => dom(ib) % sc2m
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%corprev2)) then
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=isprc-1+pl1; j=jeprc+1-pl2; k=ksprc-1+pl3 
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=isprf-1+pl1; j=jeprf+1-pl2; k=ksprf-1+pl3 
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev2)) then
                 tag=9*10**5+dom(ib)%corprev2
                 call MPI_IRECV  (dom(ib)%rc2m(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%corprev2),dom(ib)%tg(10),
     &MPI_COMM_WORLD,dom(ib)%rq_c2m,ierr)
                 call MPI_SEND (dom(ib)%sc2m(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%corprev2),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> previous cor #3
        if (dom(ib)%corprev3.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%corprev3)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%corprev3)) then
                 sbuf => dom(dom_indid(dom(ib)%corprev3))%rc3p
              else
                 sbuf => dom(ib) % sc3m
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%corprev3)) then
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ieprc+1-pl1; j=jeprc+1-pl2; k=ksprc-1+pl3;
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ieprf+1-pl1; j=jeprf+1-pl2; k=ksprf-1+pl3;
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev3)) then
                 tag=11*10**5+dom(ib)%corprev3
                 call MPI_IRECV  (dom(ib)%rc3m(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%corprev3),dom(ib)%tg(12),
     &MPI_COMM_WORLD,dom(ib)%rq_c3m,ierr)
                 call MPI_SEND (dom(ib)%sc3m(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%corprev3),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> previous cor #4
        if (dom(ib)%corprev4.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%corprev4)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%corprev4)) then
                 sbuf => dom(dom_indid(dom(ib)%corprev4))%rc4p
              else
                 sbuf => dom(ib) % sc4m
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%corprev4)) then
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ieprc+1-pl1; j=jsprc-1+pl2; k=ksprc-1+pl3;
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ieprf+1-pl1; j=jsprf-1+pl2; k=ksprf-1+pl3;
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev4)) then
                 tag=13*10**5+dom(ib)%corprev4
                 call MPI_IRECV  (dom(ib)%rc4m(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%corprev4),dom(ib)%tg(14),
     &MPI_COMM_WORLD,dom(ib)%rq_c4m,ierr)
                 call MPI_SEND (dom(ib)%sc4m(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%corprev4),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> next cor #1
        if (dom(ib)%cornext1.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%cornext1)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%cornext1)) then
                 sbuf => dom(dom_indid(dom(ib)%cornext1))%rc1m
              else
                 sbuf => dom(ib) % sc1p
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%cornext1)) then
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ieprc+1-pl1; j=jeprc+1-pl2; k=keprc+1-pl3
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ieprf+1-pl1; j=jeprf+1-pl2; k=keprf+1-pl3
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext1)) then
                 tag=8*10**5+dom(ib)%cornext1
                 call MPI_IRECV  (dom(ib)%rc1p(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%cornext1),dom(ib)%tg(7),
     &MPI_COMM_WORLD,dom(ib)%rq_c1p,ierr)
                 call MPI_SEND (dom(ib)%sc1p(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%cornext1),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> next cor #2
        if (dom(ib)%cornext2.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%cornext2)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%cornext2)) then
                 sbuf => dom(dom_indid(dom(ib)%cornext2))%rc2m
              else
                 sbuf => dom(ib) % sc2p
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%cornext2)) then
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ieprc+1-pl1; j=jsprc-1+pl2; k=keprc+1-pl3
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ieprf+1-pl1; j=jsprf-1+pl2; k=keprf+1-pl3
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext2)) then
                 tag=10*10**5+dom(ib)%cornext2
                 call MPI_IRECV  (dom(ib)%rc2p(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%cornext2),dom(ib)%tg(9),
     &MPI_COMM_WORLD,dom(ib)%rq_c2p,ierr)
                 call MPI_SEND (dom(ib)%sc2p(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%cornext2),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> next cor #3
        if (dom(ib)%cornext3.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%cornext3)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%cornext3)) then
                 sbuf => dom(dom_indid(dom(ib)%cornext3))%rc3m
              else
                 sbuf => dom(ib) % sc3p
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%cornext3)) then
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=isprc-1+pl1; j=jsprc-1+pl2; k=keprc+1-pl3
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=isprf-1+pl1; j=jsprf-1+pl2; k=keprf+1-pl3
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext3)) then
                 tag=12*10**5+dom(ib)%cornext3
                 call MPI_IRECV  (dom(ib)%rc3p(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%cornext3),dom(ib)%tg(11),
     &MPI_COMM_WORLD,dom(ib)%rq_c3p,ierr)
                 call MPI_SEND (dom(ib)%sc3p(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%cornext3),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> next cor #4
        if (dom(ib)%cornext4.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%cornext4)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%cornext4)) then
                 sbuf => dom(dom_indid(dom(ib)%cornext4))%rc4m
              else
                 sbuf => dom(ib) % sc4p
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%cornext4)) then
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=isprc-1+pl1; j=jeprc+1-pl2; k=keprc+1-pl3
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=isprf-1+pl1; j=jeprf+1-pl2; k=keprf+1-pl3
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext4)) then
                 tag=14*10**5+dom(ib)%cornext4
                 call MPI_IRECV  (dom(ib)%rc4p(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%cornext4),dom(ib)%tg(13),
     &MPI_COMM_WORLD,dom(ib)%rq_c4p,ierr)
                 call MPI_SEND (dom(ib)%sc4p(1),(pl+1)**3,MPI_FLT,
     &dom_ad(dom(ib)%cornext4),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> previous edge #1
        if (dom(ib)%edgprev1.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev1)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev1)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev1))%re1p
              else
                 sbuf => dom(ib) % se1m
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev1)) then
              tsend=njc*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,njc
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=isprc-1+pl1; j=nn; k=ksprc-1+pl2
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=njf*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,njf
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=isprf-1+pl1; j=nn; k=ksprf-1+pl2
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev1)) then
                 tag=15*10**5+dom(ib)%edgprev1
                 call MPI_IRECV  (dom(ib)%re1m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev1),dom(ib)%tg(16),
     &MPI_COMM_WORLD,dom(ib)%rq_e1m,ierr)
                 call MPI_SEND (dom(ib)%se1m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev1),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> previous edge #2
        if (dom(ib)%edgprev2.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev2)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev2)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev2))%re2p
              else
                 sbuf => dom(ib) % se2m
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev2)) then
              tsend=nic*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nic
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jeprc+1-pl1; k=ksprc-1+pl2
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=nif*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nif
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jeprf+1-pl1; k=ksprf-1+pl2
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev2)) then
                 tag=17*10**5+dom(ib)%edgprev2
                 call MPI_IRECV  (dom(ib)%re2m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev2),dom(ib)%tg(18),
     &MPI_COMM_WORLD,dom(ib)%rq_e2m,ierr)
                 call MPI_SEND (dom(ib)%se2m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev2),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> previous edge #3
        if (dom(ib)%edgprev3.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev3)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev3)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev3))%re3p
              else
                 sbuf => dom(ib) % se3m
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev3)) then
              tsend=njc*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,njc
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ieprc+1-pl1; j=nn; k=ksprc-1+pl2
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=njf*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,njf
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ieprf+1-pl1; j=nn; k=ksprf-1+pl2
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev3)) then
                 tag=19*10**5+dom(ib)%edgprev3
                 call MPI_IRECV  (dom(ib)%re3m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev3),dom(ib)%tg(20),
     &MPI_COMM_WORLD,dom(ib)%rq_e3m,ierr)
                 call MPI_SEND (dom(ib)%se3m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev3),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> previous edge #4
        if (dom(ib)%edgprev4.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev4)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev4)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev4))%re4p
              else
                 sbuf => dom(ib) % se4m
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev4)) then
              tsend=nic*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nic
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jsprc-1+pl1; k=ksprc-1+pl2
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=nif*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nif
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jsprf-1+pl1; k=ksprf-1+pl2
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev4)) then
                 tag=21*10**5+dom(ib)%edgprev4
                 call MPI_IRECV  (dom(ib)%re4m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev4),dom(ib)%tg(22),
     &MPI_COMM_WORLD,dom(ib)%rq_e4m,ierr)
                 call MPI_SEND (dom(ib)%se4m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev4),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> previous edge #5
        if (dom(ib)%edgprev5.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev5)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev5)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev5))%re5p
              else
                 sbuf => dom(ib) % se5m
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev5)) then
              tsend=nkc*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nkc
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=isprc-1+pl1; j=jsprc-1+pl2; k=nn
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=nkf*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nkf
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=isprf-1+pl1; j=jsprf-1+pl2; k=nn
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev5)) then
                 tag=23*10**5+dom(ib)%edgprev5
                 call MPI_IRECV  (dom(ib)%re5m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev5),dom(ib)%tg(24),
     &MPI_COMM_WORLD,dom(ib)%rq_e5m,ierr)
                 call MPI_SEND (dom(ib)%se5m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev5),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> previous edge #6
        if (dom(ib)%edgprev6.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev6)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev6)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev6))%re6p
              else
                 sbuf => dom(ib) % se6m
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev6)) then
              tsend=nkc*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nkc
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=isprc-1+pl1; j=jeprc+1-pl2; k=nn
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=nkf*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nkf
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=isprf-1+pl1; j=jeprf+1-pl2; k=nn
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev6)) then
                 tag=25*10**5+dom(ib)%edgprev6
                 call MPI_IRECV  (dom(ib)%re6m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev6),dom(ib)%tg(26),
     &MPI_COMM_WORLD,dom(ib)%rq_e6m,ierr)
                 call MPI_SEND (dom(ib)%se6m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev6),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> next edge #1
        if (dom(ib)%edgnext1.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext1)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext1)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext1))%re1m
              else
                 sbuf => dom(ib) % se1p
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext1)) then
              tsend=njc*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,njc
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ieprc+1-pl1; j=nn; k=keprc+1-pl2
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=njf*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,njf
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ieprf+1-pl1; j=nn; k=keprf+1-pl2
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext1)) then
                 tag=16*10**5+dom(ib)%edgnext1
                 call MPI_IRECV  (dom(ib)%re1p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext1),dom(ib)%tg(15),
     &MPI_COMM_WORLD,dom(ib)%rq_e1p,ierr)
                 call MPI_SEND (dom(ib)%se1p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext1),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> next edge #2
        if (dom(ib)%edgnext2.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext2)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext2)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext2))%re2m
              else
                 sbuf => dom(ib) % se2p
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext2)) then
              tsend=nic*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nic
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jsprc-1+pl1; k=keprc+1-pl2
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=nif*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nif
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jsprf-1+pl1; k=keprf+1-pl2
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext2)) then
                 tag=18*10**5+dom(ib)%edgnext2
                 call MPI_IRECV  (dom(ib)%re2p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext2),dom(ib)%tg(17),
     &MPI_COMM_WORLD,dom(ib)%rq_e2p,ierr)
                 call MPI_SEND (dom(ib)%se2p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext2),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> next edge #3
        if (dom(ib)%edgnext3.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext3)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext3)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext3))%re3m
              else
                 sbuf => dom(ib) % se3p
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext3)) then
              tsend=njc*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,njc
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=isprc-1+pl1; j=nn; k=keprc+1-pl2
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=njf*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,njf
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=isprf-1+pl1; j=nn; k=keprf+1-pl2
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext3)) then
                 tag=20*10**5+dom(ib)%edgnext3
                 call MPI_IRECV  (dom(ib)%re3p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext3),dom(ib)%tg(19),
     &MPI_COMM_WORLD,dom(ib)%rq_e3p,ierr)
                 call MPI_SEND (dom(ib)%se3p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext3),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> next edge #4
        if (dom(ib)%edgnext4.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext4)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext4)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext4))%re4m
              else
                 sbuf => dom(ib) % se4p
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext4)) then
              tsend=nic*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nic
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jeprc+1-pl1; k=keprc+1-pl2
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=nif*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nif
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jeprf+1-pl1; k=keprf+1-pl2
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext4)) then
                 tag=22*10**5+dom(ib)%edgnext4
                 call MPI_IRECV  (dom(ib)%re4p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext4),dom(ib)%tg(21),
     &MPI_COMM_WORLD,dom(ib)%rq_e4p,ierr)
                 call MPI_SEND (dom(ib)%se4p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext4),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> next edge #5
        if (dom(ib)%edgnext5.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext5)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext5)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext5))%re5m
              else
                 sbuf => dom(ib) % se5p
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext5)) then
              tsend=nkc*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nkc
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ieprc+1-pl1; j=jeprc+1-pl2; k=nn
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=nkf*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nkf
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ieprf+1-pl1; j=jeprf+1-pl2; k=nn
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext5)) then
                 tag=24*10**5+dom(ib)%edgnext5
                 call MPI_IRECV  (dom(ib)%re5p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext5),dom(ib)%tg(23),
     &MPI_COMM_WORLD,dom(ib)%rq_e5p,ierr)
                 call MPI_SEND (dom(ib)%se5p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext5),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!=====> next edge #6
        if (dom(ib)%edgnext6.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext6)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext6)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext6))%re6m
              else
                 sbuf => dom(ib) % se6p
              end if

              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext6)) then
              tsend=nkc*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nkc
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ieprc+1-pl1; j=jsprc-1+pl2; k=nn
                 sbuf(ijk)=fic(ib,i,j,k)
              end do; end do; end do
              else
              tsend=nkf*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nkf
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ieprf+1-pl1; j=jsprf-1+pl2; k=nn
                 sbuf(ijk)=fif(ib,i,j,k)
              end do; end do; end do
              end if

              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext6)) then
                 tag=26*10**5+dom(ib)%edgnext6
                 call MPI_IRECV  (dom(ib)%re6p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext6),dom(ib)%tg(25),
     &MPI_COMM_WORLD,dom(ib)%rq_e6p,ierr)
                 call MPI_SEND (dom(ib)%se6p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext6),tag,MPI_COMM_WORLD,ierr)

              end if
           end if
        end if
!======================================================================
        end if
!======================================================================
           end do
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
           do ib=1,nbp

              ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k

              ispr=pl+1; iepr=ni-pl
              jspr=pl+1; jepr=nj-pl
              kspr=pl+1; kepr=nk-pl

              select case (op)
              case (14) 
                 fi => dom(ib)%phi
              case (15)
                 fi => dom(ib)%phi_reinit
              case (16) 
                 fi => dom(ib)%phi_new
              case (17) 
                 fi => dom(ib)%phi_init
              case (18) 
                 fi => dom(ib)%dens
              case (19) 
                 fi => dom(ib)%mu
              end select
              is=dom(ib)%isp; ie=dom(ib)%iep
              js=dom(ib)%jsp; je=dom(ib)%jep
              ks=dom(ib)%ksp; ke=dom(ib)%kep

!======================================================================
        if (ly.eq.0)  then
!======================================================================
!=====> previous cor #1
        if (dom(ib)%corprev1.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%corprev1)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev1)) then
                 call MPI_WAIT(dom(ib)%rq_c1m,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1; st3=1
              do pl1=st1,pl; do pl2=st2,pl; do pl3=st3,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-pl1; j=jspr-pl2; k=kspr-pl3
                 if(i.lt.is .or. j.lt.js .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % rc1m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous cor #2
        if (dom(ib)%corprev2.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%corprev2)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev2)) then
                 call MPI_WAIT(dom(ib)%rq_c2m,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1; st3=1
              do pl1=st1,pl; do pl2=st2,pl; do pl3=st3,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-pl1; j=jepr+pl2; k=kspr-pl3 
                 if(i.lt.is .or. j.gt.je .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % rc2m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous cor #3
        if (dom(ib)%corprev3.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%corprev3)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev3)) then
                 call MPI_WAIT(dom(ib)%rq_c3m,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1; st3=1
              do pl1=st1,pl; do pl2=st2,pl; do pl3=st3,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+pl1; j=jepr+pl2; k=kspr-pl3
                 if(i.gt.ie .or. j.gt.je .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % rc3m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous cor #4
        if (dom(ib)%corprev4.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%corprev4)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev4)) then
                 call MPI_WAIT(dom(ib)%rq_c4m,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1; st3=1
              do pl1=st1,pl; do pl2=st2,pl; do pl3=st3,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+pl1; j=jspr-pl2; k=kspr-pl3 
                 if(i.gt.ie .or. j.lt.js .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % rc4m(ijk)
              end do; end do; end do
           end if
        end if
!=====> next cor #1
        if (dom(ib)%cornext1.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%cornext1)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext1)) then
                 call MPI_WAIT(dom(ib)%rq_c1p,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1; st3=1
              do pl1=st1,pl; do pl2=st2,pl; do pl3=st3,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+pl1; j=jepr+pl2; k=kepr+pl3
                 if(i.gt.ie .or. j.gt.je .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % rc1p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next cor #2
        if (dom(ib)%cornext2.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%cornext2)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext2)) then
                 call MPI_WAIT(dom(ib)%rq_c2p,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1; st3=1
              do pl1=st1,pl; do pl2=st2,pl; do pl3=st3,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+pl1; j=jspr-pl2; k=kepr+pl3
                 if(i.gt.ie .or. j.lt.js .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % rc2p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next cor #3
        if (dom(ib)%cornext3.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%cornext3)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext3)) then
                 call MPI_WAIT(dom(ib)%rq_c3p,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1; st3=1
              do pl1=st1,pl; do pl2=st2,pl; do pl3=st3,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-pl1; j=jspr-pl2; k=kepr+pl3
                 if(i.lt.is .or. j.lt.js .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % rc3p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next cor #4
        if (dom(ib)%cornext4.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%cornext4)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext4)) then
                 call MPI_WAIT(dom(ib)%rq_c4p,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1; st3=1
              do pl1=st1,pl; do pl2=st2,pl; do pl3=st3,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-pl1; j=jepr+pl2; k=kepr+pl3 
                 if(i.lt.is .or. j.gt.je .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % rc4p(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #1
        if (dom(ib)%edgprev1.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev1)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev1)) then
                 call MPI_WAIT(dom(ib)%rq_e1m,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=jspr,jepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-pl1; j=nn; k=kspr-pl2 
                 if(i.lt.is .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % re1m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #2
        if (dom(ib)%edgprev2.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev2)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev2)) then
                 call MPI_WAIT(dom(ib)%rq_e2m,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=ispr,iepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jepr+pl1; k=kspr-pl2
                 if(j.gt.je .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % re2m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #3
        if (dom(ib)%edgprev3.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev3)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev3)) then
                 call MPI_WAIT(dom(ib)%rq_e3m,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=jspr,jepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+pl1; j=nn; k=kspr-pl2
                 if(i.gt.ie .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % re3m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #4
        if (dom(ib)%edgprev4.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev4)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev4)) then
                 call MPI_WAIT(dom(ib)%rq_e4m,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=ispr,iepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jspr-pl1; k=kspr-pl2
                 if(j.lt.js .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % re4m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #5
        if (dom(ib)%edgprev5.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev5)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev5)) then
                 call MPI_WAIT(dom(ib)%rq_e5m,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=kspr,kepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-pl1 ; j=jspr-pl2; k=nn
                 if(i.lt.is .or. j.lt.js)
     & fi(i,j,k)=dom(ib) % re5m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #6
        if (dom(ib)%edgprev6.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgprev6)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev6)) then
                 call MPI_WAIT(dom(ib)%rq_e6m,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=kspr,kepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-pl1 ; j=jepr+pl2; k=nn
                 if(i.lt.is .or. j.gt.je)
     & fi(i,j,k)=dom(ib) % re6m(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #1
        if (dom(ib)%edgnext1.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext1)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext1)) then
                 call MPI_WAIT(dom(ib)%rq_e1p,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=jspr,jepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+pl1; j=nn; k=kepr+pl2
                 if(i.gt.ie .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % re1p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #2
        if (dom(ib)%edgnext2.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext2)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext2)) then
                 call MPI_WAIT(dom(ib)%rq_e2p,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=ispr,iepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jspr-pl1; k=kepr+pl2
                 if(j.lt.js .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % re2p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #3
        if (dom(ib)%edgnext3.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext3)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext3)) then
                 call MPI_WAIT(dom(ib)%rq_e3p,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=jspr,jepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-pl1; j=nn; k=kepr+pl2
                 if(i.lt.is .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % re3p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #4
        if (dom(ib)%edgnext4.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext4)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext4)) then
                 call MPI_WAIT(dom(ib)%rq_e4p,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=ispr,iepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jepr+pl1; k=kepr+pl2
                 if(j.gt.je .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % re4p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #5
        if (dom(ib)%edgnext5.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext5)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext5)) then
                 call MPI_WAIT(dom(ib)%rq_e5p,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=kspr,kepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+pl1 ; j=jepr+pl2; k=nn
                 if(i.gt.ie .or. j.gt.je)
     & fi(i,j,k)=dom(ib) % re5p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #6
        if (dom(ib)%edgnext6.ge.0) then
           if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%edgnext6)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext6)) then
                 call MPI_WAIT(dom(ib)%rq_e6p,MPI_STATUS_IGNORE,ierr)
              end if
              st1=1; st2=1
              do pl1=st1,pl; do pl2=st2,pl; do nn=kspr,kepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+pl1 ; j=jspr-pl2; k=nn
                 if(i.gt.ie .or. j.lt.js)
     & fi(i,j,k)=dom(ib) % re6p(ijk)
              end do; end do; end do
           end if
        end if
!======================================================================
        end if
!======================================================================
!..............................................................................
!=== Previous Neighbor  ===> 
!..............................................................................
              if (dom(ib)%iprev.ge.0)  then
              if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%iprev)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%iprev)) then
                    call MPI_WAIT(dom(ib)%rq_m1,MPI_STATUS_IGNORE,ierr)
                 end if
                 do k=kspr,kepr; do j=jspr,jepr; ijk=(k-1)*nj+j
                    fi(is-1-ly,j,k)=dom(ib) % recvb_m1(ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%jprev.ge.0)  then
              if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%jprev)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%jprev)) then
                    call MPI_WAIT(dom(ib)%rq_m2,MPI_STATUS_IGNORE,ierr)
                 end if
                 do k=kspr,kepr; do i=ispr,iepr; ijk=(k-1)*ni+i
                    fi(i,js-1-ly,k)=dom(ib) % recvb_m2(ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%kprev.ge.0)  then
              if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%kprev)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%kprev)) then
                    call MPI_WAIT(dom(ib)%rq_m3,MPI_STATUS_IGNORE,ierr)
                 end if
                 do i=ispr,iepr; do j=jspr,jepr; ijk=(i-1)*nj+j
                    fi(i,j,ks-1-ly)=dom(ib) % recvb_m3(ijk)
                 end do; end do
              end if
              end if

!..............................................................................
!=== Next Neighbor  ===> 
!..............................................................................
              if (dom(ib)%inext.ge.0)  then
              if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%inext)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%inext)) then
                    call MPI_WAIT(dom(ib)%rq_p1,MPI_STATUS_IGNORE,ierr)
                 end if
                 do k=kspr,kepr; do j=jspr,jepr; ijk=(k-1)*nj+j
                    fi(ie+1+ly,j,k)=dom(ib) % recvb_p1(ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%jnext.ge.0)  then
              if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%jnext)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%jnext)) then
                    call MPI_WAIT(dom(ib)%rq_p2,MPI_STATUS_IGNORE,ierr)
                 end if
                 do k=kspr,kepr; do i=ispr,iepr; ijk=(k-1)*ni+i
                    fi(i,je+1+ly,k)=dom(ib) % recvb_p2(ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%knext.ge.0)  then
              if(rdiv(dom_id(ib)).ne.rdiv(dom(ib)%knext)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%knext)) then
                    call MPI_WAIT(dom(ib)%rq_p3,MPI_STATUS_IGNORE,ierr)
                 end if
                 do i=ispr,iepr; do j=jspr,jepr; ijk=(i-1)*nj+j
                    fi(i,j,ke+1+ly)=dom(ib) % recvb_p3(ijk)
                 end do; end do
              end if
              end if

!--------------------------------------------------------------------------
           end do
!==========================================================================
        end do

        if(chck_c.eq.1) deallocate(fic)
        if(chck_f.eq.1) deallocate(fif)

        return
        end subroutine exchange_phi
!########################################################################## 
