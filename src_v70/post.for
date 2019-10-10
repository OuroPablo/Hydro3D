!##########################################################################
        subroutine tecplot_T(ki)
!##########################################################################
        use multidata
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecout_T_'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep
        toti=(ie+1)-(is-1)+1
        totj=(je+1)-(js-1)+1
        totk=(ke+1)-(ks-1)+1

        write (88,*) 'title = ',gf
        write (88,*)'variables=x,y,z,T,Tm,Ttm'
        write (88,*)'zone ', ' i=',toti,', ',
     &  ' j=',totj,', k= ',totk,' f=point'
        do k=ks-1,ke+1
           do j=js-1,je+1
              do i=is-1,ie+1
                 write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%zc(k)
     & ,dom(ib)%T(i,j,k),dom(ib)%Tm(i,j,k),dom(ib)%Ttm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do
88      format (15e25.8)

        end
!##########################################################################
        subroutine tecplot_S(ki)
!##########################################################################
        use multidata
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecout_S_'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep
        toti=(ie+1)-(is-1)+1
        totj=(je+1)-(js-1)+1
        totk=(ke+1)-(ks-1)+1

        write (88,*) 'title = ',gf
        write (88,*)'variables=x,y,z,S,Sm,Stm,SUM,SVM,SWM'
        write (88,*)'zone  i=',toti,',  j=',totj,', k= ',totk,' f=point'

        do k=ks-1,ke+1
           do j=js-1,je+1
              do i=is-1,ie+1
                 write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%zc(k)
     & ,dom(ib)%S(i,j,k),dom(ib)%Sm(i,j,k),dom(ib)%Stm(i,j,k)
     & ,dom(ib)%SUtm(i,j,k),dom(ib)%SVtm(i,j,k),dom(ib)%SWtm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do
88      format (15e25.8)
        end
!##########################################################################
        subroutine tecplot_p(ki)
!##########################################################################
        use multidata
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecout_p'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep
        toti=(ie+1)-(is-1)+1
        totj=(je+1)-(js-1)+1
        totk=(ke+1)-(ks-1)+1

        write (88,*) 'title = ',gf
        write (88,*)'variables=x,y,z,P,PM,ppm,vis,visM,ksgs,eps,epsm'
        write (88,*)'zone ', ' i=',toti,', ',
     &  ' j=',totj,', k= ',totk,' f=point'

        do k=ks-1,ke+1
           do j=js-1,je+1
              do i=is-1,ie+1
                 write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),
     & dom(ib)%zc(k),dom(ib)%p(i,j,k),dom(ib)%pm(i,j,k),
     & dom(ib)%ppm(i,j,k),dom(ib)%vis(i,j,k),dom(ib)%vism(i,j,k),
     & dom(ib)%ksgs(i,j,k),dom(ib)%eps(i,j,k),dom(ib)%epsm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do
88      format (15e25.8)

        end
!##########################################################################
        subroutine tec_turb(ki)
!##########################################################################
        use multidata
        use vars
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        double precision u_cn,v_cn,w_cn,um_cn,vm_cn,wm_cn
        double precision uum_cn,vvm_cn,wwm_cn,uvml,uwml,vwml
        integer :: itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecturb'//trim(adjustl(chb))//'.plt'
        open (unit=88, file=gf)

        is=pl+1; ie=dom(ib)%ttc_i-pl
        js=pl+1; je=dom(ib)%ttc_j-pl
        ks=pl+1; ke=dom(ib)%ttc_k-pl
        toti=ie-(is-1)+1
        totj=je-(js-1)+1
        totk=ke-(ks-1)+1

!Plot every X points in every direction: change itotX respectively
	  itoti=1 	; ntoti=toti/itoti+1	
	  itotj=1 	; ntotj=totj/itotj+1
	  itotk=1 	; ntotk=totk/itotk+1


	IF(L_LSM.EQ..FALSE.) THEN
        write (88,*) 'title = ',gf
        write (88,*)'variables=x,y,z,U,V,W,UM,VM,WM'
!     &, ,'uuM,vvM,wwM,uvM,uwM,vwM'
        write (88,*)'zone ', ' i=',ntoti,', ',
     &  ' j=',ntotj,', k= ',ntotk,' f=point'
	ELSE !WHEN LSM IS ALSO ACTIVATED
        write (88,*) 'title = ',gf
        write (88,*)'variables=x,y,z,U,V,W,phi,p'
        write (88,*)'zone ', ' i=',ntoti,', ',
     &  ' j=',ntotj,', k= ',ntotk,' f=point'
	ENDIF

        do k=ks-1,ke+1,itotk	
           do j=js-1,je+1,itotj
              do i=is-1,ie+1,itoti    

                 u_cn  =0.25*(dom(ib)%u(i,j,k)+
     &dom(ib)%u(i,j+1,k)+dom(ib)%u(i,j,k+1)+
     &dom(ib)%u(i,j+1,k+1))
                 um_cn  =0.25*(dom(ib)%um(i,j,k)+
     &dom(ib)%um(i,j+1,k)+dom(ib)%um(i,j,k+1)+
     &dom(ib)%um(i,j+1,k+1))
                 uum_cn  =0.25*(dom(ib)%uum(i,j,k)+
     &dom(ib)%uum(i,j+1,k)+dom(ib)%uum(i,j,k+1)+
     &dom(ib)%uum(i,j+1,k+1))

                 v_cn  =0.25*(dom(ib)%v(i,j,k)+
     &dom(ib)%v(i+1,j,k)+dom(ib)%v(i,j,k+1)+
     &dom(ib)%v(i+1,j,k+1))
                 vm_cn  =0.25*(dom(ib)%vm(i,j,k)+
     &dom(ib)%vm(i+1,j,k)+dom(ib)%vm(i,j,k+1)+
     &dom(ib)%vm(i+1,j,k+1))
                 vvm_cn  =0.25*(dom(ib)%vvm(i,j,k)+
     &dom(ib)%vvm(i+1,j,k)+dom(ib)%vvm(i,j,k+1)+
     &dom(ib)%vvm(i+1,j,k+1))

                 w_cn  =0.25*(dom(ib)%w(i,j,k)+
     &dom(ib)%w(i+1,j,k)+dom(ib)%w(i,j+1,k)+
     &dom(ib)%w(i+1,j+1,k)) 
                 wm_cn  =0.25*(dom(ib)%wm(i,j,k)+
     &dom(ib)%wm(i+1,j,k)+dom(ib)%wm(i,j+1,k)+
     &dom(ib)%wm(i+1,j+1,k)) 
                 wwm_cn  =0.25*(dom(ib)%wwm(i,j,k)+
     &dom(ib)%wwm(i+1,j,k)+dom(ib)%wwm(i,j+1,k)+
     &dom(ib)%wwm(i+1,j+1,k)) 

                 uvml  =0.125*(dom(ib)%uvm(i,j,k)+
     &dom(ib)%uvm(i+1,j,k)    +dom(ib)%uvm(i,j+1,k)+
     &dom(ib)%uvm(i+1,j+1,k)  +dom(ib)%uvm(i,j,k+1)+
     &dom(ib)%uvm(i+1,j,k+1)  +dom(ib)%uvm(i,j+1,k+1)+
     &dom(ib)%uvm(i+1,j+1,k+1))
                 uwml  =0.125*(dom(ib)%uwm(i,j,k)+
     &dom(ib)%uwm(i+1,j,k)    +dom(ib)%uwm(i,j+1,k)+
     &dom(ib)%uwm(i+1,j+1,k)  +dom(ib)%uwm(i,j,k+1)+
     &dom(ib)%uwm(i+1,j,k+1)  +dom(ib)%uwm(i,j+1,k+1)+
     &dom(ib)%uwm(i+1,j+1,k+1))
                 vwml  =0.125*(dom(ib)%vwm(i,j,k)+
     &dom(ib)%vwm(i+1,j,k)    +dom(ib)%vwm(i,j+1,k)+
     &dom(ib)%vwm(i+1,j+1,k)  +dom(ib)%vwm(i,j,k+1)+
     &dom(ib)%vwm(i+1,j,k+1)  +dom(ib)%vwm(i,j+1,k+1)+
     &dom(ib)%vwm(i+1,j+1,k+1))

	IF(L_LSM.EQ..FALSE.) THEN
          write (88,88) dom(ib)%x(i),dom(ib)%y(j),
     &  dom(ib)%z(k),u_cn,v_cn,w_cn,um_cn,vm_cn,wm_cn
!     &  ,uum_cn,vvm_cn,wwm_cn,uvml,uwml,vwml
	else
          write (88,88) dom(ib)%x(i),dom(ib)%y(j), dom(ib)%z(k)
     & ,u_cn,v_cn,w_cn,dom(ib)%phi(i,j,k),dom(ib)%p(i,j,k)
	ENDIF
              end do
           end do
        end do

        close (88)

        end do
88      format (20e25.8)
        end
!##########################################################################
        subroutine tecplot_phi(ki)
!##########################################################################
        use vars
        use multidata
        implicit none

        integer :: sn,sn1,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib
        character*8 :: chb,chb1
        character*27 :: gf

        do ib=1,nbp

	  if (L_anim_phi) then
          write(chb,'(i4)') dom_id(ib)
          write(chb1,'(i6)') ki
          sn=len(trim(adjustl(chb)))
          sn1=len(trim(adjustl(chb1)))
          chb=repeat('0',(4-sn))//trim(adjustl(chb))
          chb1=repeat('0',(6-sn1))//trim(adjustl(chb1))
          gf='tecout_phi_'//trim(adjustl(chb))//'_'//
     & trim(adjustl(chb1))//'.plt'
	  else
          write(chb,'(i4)') dom_id(ib)
          sn=len(trim(adjustl(chb)))
          chb=repeat('0',(4-sn))//trim(adjustl(chb))
          gf='tecout_phi_'//trim(adjustl(chb))//'.plt'
	  endif      

	  open (unit=88, file=gf)

        is=dom(ib)%isp; ie=dom(ib)%iep
        js=dom(ib)%jsp; je=dom(ib)%jep
        ks=dom(ib)%ksp; ke=dom(ib)%kep
        toti=ie-(is-1)+1
        totj=je-(js-1)+1
        totk=ke-(ks-1)+1

        write (88,*)'title = ',gf
        write (88,*)'variables=x,y,z,phi,phim,dens,mu'
	  if (L_anim_phi) then
          write (88,*)'zone ','STRANDID=', 1, 'SOLUTIONTIME=', ctime,
     &    ' i=',toti,', ',' j=',totj,', k= ',totk,
     &    'zonetype=', 'ordered',', DATAPACKING=point'
	  else
          write(88,*)'zone  i=',toti,', j=',totj,', k= ',totk,' f=point'
	  endif

        do k=ks-1,ke
           do j=js-1,je
              do i=is-1,ie
                write (88,88) dom(ib)%xc(i),dom(ib)%yc(j),dom(ib)%zc(k)
     & ,dom(ib)%phi(i,j,k),dom(ib)%phim(i,j,k),dom(ib)%dens(i,j,k),
     & dom(ib)%mu(i,j,k)
              end do
           end do
        end do
        close (88)

        end do

88      format (8e20.8)

        end
!##########################################################################
        subroutine tecbin(ki)
!##########################################################################
        use multidata
        use vars
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: ki,ib,inind,jnind,knind
        character*8 :: chb
        character*25 :: gf

        do ib=1,nbp

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        gf='tecbin'//trim(adjustl(chb))//'.bin'
        open (unit=88, file=gf, form='unformatted')

        toti=dom(ib)%ttc_i
        totj=dom(ib)%ttc_j
        totk=dom(ib)%ttc_k

        write (88) toti,totj,totk
        write (88) pl
!====================================================================
        inind=0; jnind=0; knind=0
        if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5) inind=-1
        if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5) jnind=-1
        if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5) knind=-1
        write (88) inind,jnind,knind
!====================================================================

        do k=1,totk
           do j=1,totj
              do i=1,toti

        write (88) dom(ib)%x(i),dom(ib)%y(j),dom(ib)%z(k),
     & dom(ib)%p(i,j,k),dom(ib)%pm(i,j,k),dom(ib)%ppm(i,j,k),
     & dom(ib)%vis(i,j,k),dom(ib)%vism(i,j,k),
     & dom(ib)%u(i,j,k),dom(ib)%um(i,j,k),dom(ib)%uum(i,j,k),
     & dom(ib)%v(i,j,k),dom(ib)%vm(i,j,k),dom(ib)%vvm(i,j,k),
     & dom(ib)%w(i,j,k),dom(ib)%wm(i,j,k),dom(ib)%wwm(i,j,k),
     & dom(ib)%uvm(i,j,k),dom(ib)%uwm(i,j,k),dom(ib)%vwm(i,j,k),
     & dom(ib)%ksgs(i,j,k),dom(ib)%eps(i,j,k),dom(ib)%epsm(i,j,k)
              end do
           end do
        end do
        close (88)

        end do
        end
!##########################################################################
        subroutine tec_instant(ki)
!##########################################################################
        use multidata
        use vars
	use mpi
        implicit none
        integer :: sn,i,j,k,toti,totj,totk
        integer :: is,ie,js,je,ks,ke
        integer :: ki,ib,itoti,itotj,itotk,ntoti,ntotj,ntotk
        character*8 :: chb,chb2
        character*35 :: gf

!	if (nbp.gt.1) RETURN

        do ib=1,nbp

!	  if(rdiv(dom_id(ib)).le.2) RETURN

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecinst'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
        open (unit=95, file=gf)        

        is=pl+1; ie=dom(ib)%ttc_i-pl
        js=pl+1; je=dom(ib)%ttc_j-pl
        ks=pl+1; ke=dom(ib)%ttc_k-pl
        toti=ie-(is-1)+1	!toti=nicell
        totj=je-(js-1)+1
        totk=ke-(ks-1)+1

!Plot every X points in every direction: change itotX respectively
	   itoti=1 	; ntoti=toti/itoti+1
	   itotj=1 	; ntotj=totj/itotj+1
	   itotk=1 	; ntotk=totk/itotk+1
	
        write(95,*)'variables=x,y,z,U,V,W'
        write(95,89)'zone i=',ntoti,', j=',ntotj,', k=',ntotk,' f=point'

        do k=ks-1,ke+1,itotk	
           do j=js-1,je+1,itotj
              do i=is-1,ie+1,itoti     
             write (95,88) dom(ib)%x(i),dom(ib)%y(j),dom(ib)%z(k),
     &  dom(ib)%u(i,j,k),dom(ib)%v(i,j,k),dom(ib)%w(i,j,k)
              end do
           end do
        end do

        close (95)
        end do

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end
!##########################################################################
        subroutine tec_inst_plane(ki)
!##########################################################################
        use multidata
        use vars
	  use mpi
        implicit none
        integer :: sn,i,j,k,toti,totj,totk,is,ie,js,je,ks,ke,ki,ib
        integer :: itoti,itotj,itotk,ntoti,ntotj,ntotk
	  double precision :: xval,yval,zval
	  integer :: iplane
        character*8 :: chb,chb2
        character*35 :: gf

        do ib=1,nbp

        is=pl+1; ie=dom(ib)%ttc_i-pl
        js=pl+1; je=dom(ib)%ttc_j-pl
        ks=pl+1; ke=dom(ib)%ttc_k-pl
        toti=ie-(is-1)+1 ; itoti=1 	; ntoti=toti/itoti+1
        totj=je-(js-1)+1 ; itotj=1 	; ntotj=totj/itotj+1
        totk=ke-(ks-1)+1 ; itotk=1 	; ntotk=totk/itotk+1
	  
!!!! Set the plane orientation and location at xval, yval or zval
	  iplane=2 			!1=x-plane, 2=y-plane, 3=z-plane
	  xval=(xen-xst)/2 ;  yval=(yen-yst)/2  ;  zval=(zen-zst)/2 

	 if(iplane.eq.1) then
        do i=is-1,ie
	    if(dom(ib)%x(i).le.xval .and. dom(ib)%x(i+1).gt.xval) then

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecplane'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
           open (unit=95, file=gf)     
        write(95,*)'variables=y,z,u,v,w'
        write(95,89)'zone i= 1, j=',ntotj,', k=',ntotk,' f=point'		
           do k=ks-1,ke+1,itotk	!k=pl,nicell+pl
             do j=js-1,je+1,itotj
             write (95,88) dom(ib)%y(j),dom(ib)%z(k),
     &  dom(ib)%u(i,j,k),dom(ib)%v(i,j,k),dom(ib)%w(i,j,k)
             end do
           end do
           close (95)
	    endif ! within the interval
	  enddo
	 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	 if(iplane.eq.2) then
        do j=js-1,je
	    if(dom(ib)%y(j).le.yval .and. dom(ib)%y(j+1).gt.yval) then

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecplane'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
           open (unit=95, file=gf)     
        write(95,*)'variables=x,z,u,v,w'
        write(95,89)'zone i=',ntoti,', j=1 , k=',ntotk,' f=point'		
           do k=ks-1,ke+1,itotk	
             do i=is-1,ie+1,itoti
             write (95,88) dom(ib)%x(i),dom(ib)%z(k),
     &  dom(ib)%u(i,j,k),dom(ib)%v(i,j,k),dom(ib)%w(i,j,k)
!     &  ,dom(ib)%phi(i,j,k)
             end do
           end do
           close (95)
	    endif ! within the interval
	  enddo
	 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	 if(iplane.eq.3) then
        do k=ks-1,ke
	    if(dom(ib)%z(k).le.zval .and. dom(ib)%z(k+1).gt.zval) then

        write(chb,'(i4)') dom_id(ib)
        sn=len(trim(adjustl(chb)))
        chb=repeat('0',(4-sn))//trim(adjustl(chb))
        write(chb2,'(i6)') ki
        sn=len(trim(adjustl(chb2)))
        chb2=repeat('0',(6-sn))//trim(adjustl(chb2))
        gf='tecplane'//trim(adjustl(chb2))//
     &		'_'//trim(adjustl(chb))//'.plt'
           open (unit=95, file=gf)     
        write(95,*)'variables=x,y,u,v,w'
        write(95,89)'zone i=',ntoti,', j=',ntotj,', k=1 f=point'		
           do j=js-1,je+1,itotj
            do i=is-1,ie+1,itoti
             write (95,88) dom(ib)%x(i),dom(ib)%y(j),
     &  dom(ib)%u(i,j,k),dom(ib)%v(i,j,k),dom(ib)%w(i,j,k)
             end do
           end do
           close (95)
	    endif ! within the interval
	  enddo
	 endif

	enddo !ib

88      format (10e15.5)
89      format (a7,i4,a4,i4,a4,i4,a8)

        end
