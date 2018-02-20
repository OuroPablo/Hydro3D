!##########################################################################
        subroutine mom_buo
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib
        integer :: is,ie,js,je,ks,ke
 

        do ib=1,nbp

        is=dom(ib)%isu; ie=dom(ib)%ieu
        js=dom(ib)%jsu; je=dom(ib)%jeu
        ks=dom(ib)%ksu; ke=dom(ib)%keu

        do k=ks,ke
           do i=is,ie
              do j=js,je
                 dom(ib)%ustar(i,j,k)=
     & dom(ib)%ustar(i,j,k)+dt*alfapr*grx*
     & (1.d0-0.5*beta*(dom(ib)%T(i+1,j,k)+dom(ib)%T(i,j,k)))
              end do
           end do
        end do

        is=dom(ib)%isv; ie=dom(ib)%iev
        js=dom(ib)%jsv; je=dom(ib)%jev
        ks=dom(ib)%ksv; ke=dom(ib)%kev

        do k=ks,ke
           do i=is,ie
              do j=js,je
                 dom(ib)%vstar(i,j,k)=
     & dom(ib)%vstar(i,j,k)+dt*alfapr*gry*
     & (1.d0-0.5*beta*(dom(ib)%T(i,j+1,k)+dom(ib)%T(i,j,k)))
              end do
           end do
        end do

        is=dom(ib)%isw; ie=dom(ib)%iew
        js=dom(ib)%jsw; je=dom(ib)%jew
        ks=dom(ib)%ksw; ke=dom(ib)%kew

        do k=ks,ke
           do i=is,ie
              do j=js,je
                 dom(ib)%wstar(i,j,k)=
     & dom(ib)%wstar(i,j,k)+dt*alfapr*grz*
     & (1.d0-0.5*beta*(dom(ib)%T(i,j,k+1)+dom(ib)%T(i,j,k)))

              end do
           end do
        end do

        end do

        return
        end subroutine mom_buo
!##########################################################################
        subroutine energy
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,k,ib
        double precision :: dxx,dyy,dzz
        double precision :: conv,diff
        double precision :: duTdx,dvTdy,dwTdz
        double precision :: awT,aeT,asT,anT,abT,atT,apT
        double precision :: kp,km,ku,kc,kd,b_r

        do ib=1,nbp

        dxx=dom(ib)%dx*dom(ib)%dx
        dyy=dom(ib)%dy*dom(ib)%dy
        dzz=dom(ib)%dz*dom(ib)%dz

           do k=dom(ib)%ksp,dom(ib)%kep
              do i=dom(ib)%isp,dom(ib)%iep
                 do j=dom(ib)%jsp,dom(ib)%jep


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if(dom(ib)%u(i-1,j,k).gt.0.0) then
           ku=dom(ib)%To(i-2,j,k)
           kc=dom(ib)%To(i-1,j,k)
           kd=dom(ib)%To(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%u(i-1,j,k).lt.0.0) then
           ku=dom(ib)%To(i+1,j,k)
           kc=dom(ib)%To(i,j,k)
           kd=dom(ib)%To(i-1,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else
           km=0.5*(dom(ib)%To(i,j,k)+dom(ib)%To(i-1,j,k))
        end if
        if(dom(ib)%u(i,j,k).gt.0.0) then
           ku=dom(ib)%To(i-1,j,k)
           kc=dom(ib)%To(i,j,k)
           kd=dom(ib)%To(i+1,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%u(i,j,k).lt.0.0) then
           ku=dom(ib)%To(i+2,j,k)
           kc=dom(ib)%To(i+1,j,k)
           kd=dom(ib)%To(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else
           kp=0.5*(dom(ib)%To(i,j,k)+dom(ib)%To(i+1,j,k))
        end if
        duTdx=(dom(ib)%u(i,j,k)*kp-dom(ib)%u(i-1,j,k)*km)/dom(ib)%dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(dom(ib)%v(i,j-1,k).gt.0.0) then
           ku=dom(ib)%To(i,j-2,k)
           kc=dom(ib)%To(i,j-1,k)
           kd=dom(ib)%To(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%v(i,j-1,k).lt.0.0) then
           ku=dom(ib)%To(i,j+1,k)
           kc=dom(ib)%To(i,j,k)
           kd=dom(ib)%To(i,j-1,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else
           km=0.5*(dom(ib)%To(i,j,k)+dom(ib)%To(i,j-1,k))
        end if
        if(dom(ib)%v(i,j,k).gt.0.0) then
           ku=dom(ib)%To(i,j-1,k)
           kc=dom(ib)%To(i,j,k)
           kd=dom(ib)%To(i,j+1,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%v(i,j,k).lt.0.0) then
           ku=dom(ib)%To(i,j+2,k)
           kc=dom(ib)%To(i,j+1,k)
           kd=dom(ib)%To(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else
           kp=0.5*(dom(ib)%To(i,j,k)+dom(ib)%To(i,j+1,k))
        end if
        dvTdy=(dom(ib)%v(i,j,k)*kp-dom(ib)%v(i,j-1,k)*km)/dom(ib)%dy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(dom(ib)%w(i,j,k-1).gt.0.0) then
           ku=dom(ib)%To(i,j,k-2)
           kc=dom(ib)%To(i,j,k-1)
           kd=dom(ib)%To(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%w(i,j,k-1).lt.0.0) then
           ku=dom(ib)%To(i,j,k+1)
           kc=dom(ib)%To(i,j,k)
           kd=dom(ib)%To(i,j,k-1)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           km=kc+0.5*b_r*(kc-ku)
        else
           km=0.5*(dom(ib)%To(i,j,k)+dom(ib)%To(i,j,k-1))
        end if
        if(dom(ib)%w(i,j,k).gt.0.0) then
           ku=dom(ib)%To(i,j,k-1)
           kc=dom(ib)%To(i,j,k)
           kd=dom(ib)%To(i,j,k+1)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else if(dom(ib)%w(i,j,k).lt.0.0) then
           ku=dom(ib)%To(i,j,k+2)
           kc=dom(ib)%To(i,j,k+1)
           kd=dom(ib)%To(i,j,k)
           b_r=max(0.0,
     & min(2.0*((kd-kc)/(kc-ku)),0.75*((kd-kc)/(kc-ku))+0.25,4.0))
           kp=kc+0.5*b_r*(kc-ku)
        else
           kp=0.5*(dom(ib)%To(i,j,k)+dom(ib)%To(i,j,k+1))
        end if
        dwTdz=(dom(ib)%w(i,j,k)*kp-dom(ib)%w(i,j,k-1)*km)/dom(ib)%dz
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        conv=(duTdx+dvTdy+dwTdz)


                 awT=-(1.0/(Re*Pr)+(dom(ib)%vis(i,j,k)-1.0/Re)/Sc_t)/dxx
                 aeT=-(1.0/(Re*Pr)+(dom(ib)%vis(i,j,k)-1.0/Re)/Sc_t)/dxx
                 anT=-(1.0/(Re*Pr)+(dom(ib)%vis(i,j,k)-1.0/Re)/Sc_t)/dyy
                 asT=-(1.0/(Re*Pr)+(dom(ib)%vis(i,j,k)-1.0/Re)/Sc_t)/dyy
                 atT=-(1.0/(Re*Pr)+(dom(ib)%vis(i,j,k)-1.0/Re)/Sc_t)/dzz
                 abT=-(1.0/(Re*Pr)+(dom(ib)%vis(i,j,k)-1.0/Re)/Sc_t)/dzz

                 apT = -1.0*(awT+aeT+asT+anT+abT+atT)

                 diff=(apT*dom(ib)%To(i,j,k)+
     & anT*dom(ib)%To(i,j+1,k) + asT*dom(ib)%To(i,j-1,k)+
     & aeT*dom(ib)%To(i+1,j,k) + awT*dom(ib)%To(i-1,j,k)+
     & atT*dom(ib)%To(i,j,k+1) + abT*dom(ib)%To(i,j,k-1))
 
                 dom(ib)%T(i,j,k)=(dom(ib)%To(i,j,k)-
     & dt*(conv+diff)) 

                 end do
              end do
           end do
        end do

        call exchange(5)
        call boundT

        return
        end subroutine energy
!##########################################################################
        subroutine boundT
!##########################################################################
        use vars
        use multidata
        implicit none
        integer :: i,j,k,ib,ni,nj,nk,ly
        integer :: is,ie,js,je,ks,ke

        do ly=0,pl_ex

        do ib=1,nbp
           ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k
           is=dom(ib)%isp; ie=dom(ib)%iep
           js=dom(ib)%jsp; je=dom(ib)%jep
           ks=dom(ib)%ksp; ke=dom(ib)%kep
	
! Boundary Conditions for T
!..............................................................................
!=== West ===> 
!..............................................................................
        if (dom(ib)%iprev.lt.0) then
           if (dom(ib)%Tbc_west.eq.7) then
              do k=1,nk; do j=1,nj
                 dom(ib)%T(is-1-ly,j,k)= dom(ib)%T(is+ly,j,k)	
              end do; end do
           else if (dom(ib)%Tbc_west.eq.8) then
              do k=1,nk; do j=1,nj
                 dom(ib)%T(is-1-ly,j,k)= 2.d0*Tc-dom(ib)%T(is,j,k)
              end do; end do
           else if (dom(ib)%Tbc_west.eq.9) then
              do k=1,nk; do j=1,nj
                 dom(ib)%T(is-1-ly,j,k)= 2.d0*Th-dom(ib)%T(is,j,k)
              end do; end do
           end if
        end if
!...............................................................................
!=== East ===> 
!...............................................................................
        if (dom(ib)%inext.lt.0) then
           if (dom(ib)%Tbc_east.eq.7) then
              do k=1,nk; do j=1,nj
                 dom(ib)%T(ie+1+ly,j,k)= dom(ib)%T(ie-ly,j,k)	
              end do; end do
           else if (dom(ib)%Tbc_east.eq.8) then
              do k=1,nk; do j=1,nj
                 dom(ib)%T(ie+1+ly,j,k)= 2.d0*Tc-dom(ib)%T(ie,j,k)	
              end do; end do
           else if (dom(ib)%Tbc_east.eq.9) then
              do k=1,nk; do j=1,nj
                 dom(ib)%T(ie+1+ly,j,k)= 2.d0*Th-dom(ib)%T(ie,j,k)	
              end do; end do
           end if
        end if
!...............................................................................
!=== South ===>
!...............................................................................
        if (dom(ib)%jprev.lt.0) then
           if (dom(ib)%Tbc_south.eq.7) then 
              do k=1,nk; do i=1,ni
                 dom(ib)%T(i,js-1-ly,k)= dom(ib)%T(i,js+ly,k)	
              end do; end do
           else if (dom(ib)%Tbc_south.eq.8) then 
              do k=1,nk; do i=1,ni
                 dom(ib)%T(i,js-1-ly,k)= 2.d0*Tc-dom(ib)%T(i,js,k)	
              end do; end do
           else if (dom(ib)%Tbc_south.eq.9) then 
              do k=1,nk; do i=1,ni
                 dom(ib)%T(i,js-1-ly,k)= 2.d0*Th-dom(ib)%T(i,js,k)	
              end do; end do
           end if
        end if
!.............................................................................
!=== North ===>
!.............................................................................
        if (dom(ib)%jnext.lt.0) then
           if (dom(ib)%Tbc_north.eq.7) then
              do k=1,nk; do i=1,ni
                 dom(ib)%T(i,je+1+ly,k) = dom(ib)%T(i,je-ly,k) 	
              end do; end do
           else if (dom(ib)%Tbc_north.eq.8) then
              do k=1,nk; do i=1,ni
                 dom(ib)%T(i,je+1+ly,k) = 2.d0*Tc-dom(ib)%T(i,je,k) 	
              end do; end do
           else if (dom(ib)%Tbc_north.eq.9) then
              do k=1,nk; do i=1,ni
                 dom(ib)%T(i,je+1+ly,k) = 2.d0*Th-dom(ib)%T(i,je,k) 	
              end do; end do
           end if
        end if
!...............................................................................
!=== Bottom ===>
!...............................................................................
        if (dom(ib)%kprev.lt.0) then
           if (dom(ib)%Tbc_bottom.eq.7) then
              do j=1,nj; do i=1,ni
                 dom(ib)%T(i,j,ks-1-ly)= dom(ib)%T(i,j,ks+ly)	
              end do; end do
           else if (dom(ib)%Tbc_bottom.eq.8) then
              do j=1,nj; do i=1,ni
                 dom(ib)%T(i,j,ks-1-ly)= 2.d0*Tc-dom(ib)%T(i,j,ks)	
              end do; end do
           else if (dom(ib)%Tbc_bottom.eq.9) then
              do j=1,nj; do i=1,ni
                 dom(ib)%T(i,j,ks-1-ly)= 2.d0*Th-dom(ib)%T(i,j,ks)	
              end do; end do
           end if
        end if
!.............................................................................
!=== Top ===>
!.............................................................................
        if (dom(ib)%knext.lt.0) then
           if (dom(ib)%Tbc_top.eq.7) then
              do j=1,nj; do i=1,ni
                 dom(ib)%T(i,j,ke+1+ly) = dom(ib)%T(i,j,ke-ly)	
              end do; end do
           else if (dom(ib)%Tbc_top.eq.8) then
              do j=1,nj; do i=1,ni
                 dom(ib)%T(i,j,ke+1+ly) = 2.d0*Tc-dom(ib)%T(i,j,ke)	
              end do; end do
           else if (dom(ib)%Tbc_top.eq.9) then
              do j=1,nj; do i=1,ni
                 dom(ib)%T(i,j,ke+1+ly) = 2.d0*Th-dom(ib)%T(i,j,ke)	
              end do; end do
           end if
        end if

!==============================================================================
        end do
        end do

        end subroutine boundT
!#############################################################################
