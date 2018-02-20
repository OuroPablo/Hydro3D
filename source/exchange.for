!##########################################################################
        subroutine  exchange(op)
!##########################################################################
        use vars
        use multidata
        implicit none
        integer :: op,ly,ib,ni,nj,nk,ijk,nly
        integer :: i,j,k,is,ie,js,je,ks,ke
        integer :: ispr,iepr,jspr,jepr,kspr,kepr
        double precision, pointer, dimension(:,:,:) :: fi

! stfcinf: 1=iprev, 3=jprev, 5=kprev, 2=inext, 4=jnext, 6=knext

        if (op.eq.4 .or. op.eq.44 .or. op.eq.7) then! .or. op.eq.14 .or.
!     & op.eq.15 .or. op.eq.16 .or. op.eq.17 .or. op.eq.18 .or. 
!     & op.eq.19) then
           nly=0
        else 
           nly=pl_ex
        end if

        if (PERIODIC) call exchange_bc(op,nly)
        call exchange_smlvl(op,nly)

        if(rdivmax.gt.1) then

        select case (op)
           case (1)  
              call exchangeu(op,nly)
           case (2)  
              call exchangev(op,nly)
           case (3)  
              call exchangew(op,nly)
           case (11)  
              call exchangeu(op,nly)
           case (22)  
              call exchangev(op,nly)
           case (33)  
              call exchangew(op,nly)
           case (4)  
              call exchangep(op)
           case (44)  
              call exchangep(op)
           case (5)  
              call exchangesca(op,nly)
           case (6)  
              call exchangesca(op,nly)
           case (7)  
              call exchangesca(op,nly)
           case (8)
              call exchangesca(op,nly)
           case (9)  
              call exchangesca(op,nly)
	     case (14)
		  call exchange_phi(op,nly) !phi
	     case (15)
		  call exchange_phi(op,nly) !phi_reinit
	     case (16)
		  call exchange_phi(op,nly) !phi_new
	     case (17)
		  call exchange_phi(op,nly) !phi_init
	     case (18)
		  call exchange_phi(op,nly)	!dens
	     case (19)
		  call exchange_phi(op,nly)	!mu
        end select

        if(LMR.eq.2) then 

        do ib=1,nbp

           ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k  

           select case (op)
              case (1)  
                 fi => dom(ib)%u
                 is=dom(ib)%isu; ie=dom(ib)%ieu
                 js=dom(ib)%jsu; je=dom(ib)%jeu
                 ks=dom(ib)%ksu; ke=dom(ib)%keu
              case (2)  
                 fi => dom(ib)%v
                 is=dom(ib)%isv; ie=dom(ib)%iev
                 js=dom(ib)%jsv; je=dom(ib)%jev
                 ks=dom(ib)%ksv; ke=dom(ib)%kev
              case (3)  
                 fi => dom(ib)%w
                 is=dom(ib)%isw; ie=dom(ib)%iew
                 js=dom(ib)%jsw; je=dom(ib)%jew
                 ks=dom(ib)%ksw; ke=dom(ib)%kew
              case (11)  
                 fi => dom(ib)%ustar
                 is=dom(ib)%isu; ie=dom(ib)%ieu
                 js=dom(ib)%jsu; je=dom(ib)%jeu
                 ks=dom(ib)%ksu; ke=dom(ib)%keu
              case (22)  
                 fi => dom(ib)%vstar
                 is=dom(ib)%isv; ie=dom(ib)%iev
                 js=dom(ib)%jsv; je=dom(ib)%jev
                 ks=dom(ib)%ksv; ke=dom(ib)%kev
              case (33)  
                 fi => dom(ib)%wstar
                 is=dom(ib)%isw; ie=dom(ib)%iew
                 js=dom(ib)%jsw; je=dom(ib)%jew
                 ks=dom(ib)%ksw; ke=dom(ib)%kew
              case (4) 
                 fi => dom(ib)%p
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              case (44) 
                 fi => dom(ib)%pp
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              case (5) 
                 fi => dom(ib)%T
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              case (6) 
                 fi => dom(ib)%ksgs
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              case (7) 
                 fi => dom(ib)%vis
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              case (8)
                 fi => dom(ib)%S
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep  
              case (9)
                 fi => dom(ib)%eps
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep	
	        case (14)
		     fi => dom(ib)%phi
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
	        case (15)
		     fi => dom(ib)%phi_reinit
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
	        case (16)
		     fi => dom(ib)%phi_new
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
	        case (17)
		     fi => dom(ib)%phi_init
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
	        case (18)
		     fi => dom(ib)%dens
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
	        case (19)
		     fi => dom(ib)%mu
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
           end select

           do ly=0,nly

!..............................................................................
!=== Previous Neighbor  ===> 
!..............................................................................
              if (dom(ib)%iprev.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%iprev)) then
!===>> iprev !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.2 .or. op.eq.22) then
        if(dom(ib)%edgprev5.ge.0 .and. 
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev5)) jspr=pl
                 else if (op.eq.3 .or. op.eq.33) then
        if(dom(ib)%edgprev1.ge.0 .and. 
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev1)) kspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do k=kspr,kepr; do j=jspr,jepr; ijk=(k-1)*nj+j
                    fi(is-1-ly,j,k)=dom(ib)%stfcinf(1,ly+1,ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%jprev.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%jprev)) then
!===>> jprev !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.1 .or. op.eq.11) then
        if(dom(ib)%edgprev5.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev5)) ispr=pl
                 else if (op.eq.3 .or. op.eq.33) then
        if(dom(ib)%edgprev4.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev4)) kspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do k=kspr,kepr; do i=ispr,iepr; ijk=(k-1)*ni+i
                    fi(i,js-1-ly,k)=dom(ib)%stfcinf(3,ly+1,ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%kprev.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%kprev)) then
!===>> kprev !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.1 .or. op.eq.11) then
        if(dom(ib)%edgprev1.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev1)) ispr=pl
                 else if (op.eq.2 .or. op.eq.22) then
        if(dom(ib)%edgprev4.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev4)) jspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do i=ispr,iepr; do j=jspr,jepr; ijk=(i-1)*nj+j
                    fi(i,j,ks-1-ly)=dom(ib)%stfcinf(5,ly+1,ijk)
                 end do; end do
              end if
              end if
!..............................................................................
!=== Next Neighbor  ===> 
!..............................................................................
              if (dom(ib)%inext.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%inext)) then
!===>> inext !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.2 .or. op.eq.22) then
        if(dom(ib)%edgnext6.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext6)) jspr=pl
                 else if (op.eq.3 .or. op.eq.33) then
        if(dom(ib)%edgprev3.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev3)) kspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do k=kspr,kepr; do j=jspr,jepr; ijk=(k-1)*nj+j
                    fi(ie+1+ly,j,k)=dom(ib)%stfcinf(2,ly+1,ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%jnext.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%jnext)) then
!===>> jnext !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.1 .or. op.eq.11) then
        if(dom(ib)%edgprev6.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev6)) ispr=pl
                 else if (op.eq.3 .or. op.eq.33) then
        if(dom(ib)%edgprev2.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgprev2)) kspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do k=kspr,kepr; do i=ispr,iepr; ijk=(k-1)*ni+i
                    fi(i,je+1+ly,k)=dom(ib)%stfcinf(4,ly+1,ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%knext.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%knext)) then
!===>> knext !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 ispr=pl+1;	jspr=pl+1;	kspr=pl+1
                 iepr=ni-pl;	jepr=nj-pl;	kepr=nk-pl
                 if (op.eq.1 .or. op.eq.11) then
        if(dom(ib)%edgnext3.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext3)) ispr=pl
                 else if (op.eq.2 .or. op.eq.22) then
        if(dom(ib)%edgnext2.ge.0 .and.
     &rdiv(dom_id(ib)).gt.rdiv(dom(ib)%edgnext2)) jspr=pl
                 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 do i=ispr,iepr; do j=jspr,jepr; ijk=(i-1)*nj+j
                    fi(i,j,ke+1+ly)=dom(ib)%stfcinf(6,ly+1,ijk)
                 end do; end do
              end if
              end if

           end do
        end do
        end if
        end if

        return
        end subroutine exchange
!##########################################################################
        subroutine  exchange_smlvl(op,nly)
!##########################################################################
        use multidata
        use mpi
        use vars
        implicit none
        integer :: i,j,k,ijk,ib,op,nly
        integer :: is,ie,js,je,ks,ke
        integer :: iepr,jepr,kepr,ispr,jspr,kspr
        integer :: tag,ta,ly
        integer :: tsend,trecv
        integer :: ni,nj,nk,pl1,pl2,pl3,pll,nn
        double precision, pointer, dimension(:,:,:) :: fi
        double precision, pointer, dimension(:) :: sbuf

        pll=(pl+1)*(pl+1)

        MPI_FLT   = MPI_DOUBLE_PRECISION

        do ly=0,nly
!==========================================================================
!==========================================================================
           do ib=1,nbp

              ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k  

              ispr=pl+1; iepr=ni-pl
              jspr=pl+1; jepr=nj-pl
              kspr=pl+1; kepr=nk-pl

              select case (op)
                 case (1)  
                    fi => dom(ib)%u
                    is=dom(ib)%isu; ie=dom(ib)%ieu
                    js=dom(ib)%jsu; je=dom(ib)%jeu
                    ks=dom(ib)%ksu; ke=dom(ib)%keu
                 case (2)  
                    fi => dom(ib)%v
                    is=dom(ib)%isv; ie=dom(ib)%iev
                    js=dom(ib)%jsv; je=dom(ib)%jev
                    ks=dom(ib)%ksv; ke=dom(ib)%kev
                 case (3)  
                    fi => dom(ib)%w
                    is=dom(ib)%isw; ie=dom(ib)%iew
                    js=dom(ib)%jsw; je=dom(ib)%jew
                    ks=dom(ib)%ksw; ke=dom(ib)%kew
                 case (11)  
                    fi => dom(ib)%ustar
                    is=dom(ib)%isu; ie=dom(ib)%ieu
                    js=dom(ib)%jsu; je=dom(ib)%jeu
                    ks=dom(ib)%ksu; ke=dom(ib)%keu
                 case (22)  
                    fi => dom(ib)%vstar
                    is=dom(ib)%isv; ie=dom(ib)%iev
                    js=dom(ib)%jsv; je=dom(ib)%jev
                    ks=dom(ib)%ksv; ke=dom(ib)%kev
                 case (33)  
                    fi => dom(ib)%wstar
                    is=dom(ib)%isw; ie=dom(ib)%iew
                    js=dom(ib)%jsw; je=dom(ib)%jew
                    ks=dom(ib)%ksw; ke=dom(ib)%kew
                 case (4) 
                    fi => dom(ib)%p
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (44) 
                    fi => dom(ib)%pp
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (5) 
                    fi => dom(ib)%T
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (6) 
                    fi => dom(ib)%ksgs
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (7) 
                    fi => dom(ib)%vis
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (8)
                    fi => dom(ib)%S
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (9)
                    fi => dom(ib)%eps
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (14)
			  fi => dom(ib)%phi
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (15)
		  	  fi => dom(ib)%phi_reinit
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (16)
		        fi => dom(ib)%phi_new
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (17)
		        fi => dom(ib)%phi_init
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (18)
		        fi => dom(ib)%dens
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (19)
		        fi => dom(ib)%mu
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
              end select
!..........................................................................
!=== Previous Neighbor  ===> 
!..........................................................................
           if (dom(ib)%iprev.ge.0)  then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%iprev)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%iprev)) then
                 sbuf => dom(dom_indid(dom(ib)%iprev))%recvb_p1
              else
                 sbuf => dom(ib) % sendb_m1
              end if

              tsend=nj*nk; trecv=nj*nk
              do k=1,nk; do j=1,nj; ijk=(k-1)*nj+j
              sbuf(ijk)=fi(is+ly,j,k)
              end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%jprev)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%jprev)) then
                 sbuf => dom(dom_indid(dom(ib)%jprev))%recvb_p2
              else
                 sbuf => dom(ib) % sendb_m2
              end if

              tsend=ni*nk; trecv=ni*nk
              do k=1,nk; do i=1,ni; ijk=(k-1)*ni+i
              sbuf(ijk)=fi(i,js+ly,k)
              end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%kprev)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%kprev)) then
                 sbuf => dom(dom_indid(dom(ib)%kprev))%recvb_p3
              else
                 sbuf => dom(ib) % sendb_m3
              end if

              tsend=ni*nj; trecv=ni*nj
              do i=1,ni; do j=1,nj; ijk=(i-1)*nj+j
              sbuf(ijk)=fi(i,j,ks+ly)
              end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%inext)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%inext)) then
                 sbuf => dom(dom_indid(dom(ib)%inext))%recvb_m1
              else
                 sbuf => dom(ib) % sendb_p1
              end if

              tsend=nj*nk; trecv=nj*nk
              do k=1,nk; do j=1,nj; ijk=(k-1)*nj+j
              sbuf(ijk)=fi(ie-ly,j,k)
              end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%jnext)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%jnext)) then
                 sbuf => dom(dom_indid(dom(ib)%jnext))%recvb_m2
              else
                 sbuf => dom(ib) % sendb_p2
              end if

              tsend=ni*nk; trecv=ni*nk
              do k=1,nk; do i=1,ni; ijk=(k-1)*ni+i;
              sbuf(ijk)=fi(i,je-ly,k)
              end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%knext)) then

              if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%knext)) then
                 sbuf => dom(dom_indid(dom(ib)%knext))%recvb_m3
              else
                 sbuf => dom(ib) % sendb_p3
              end if

              tsend=ni*nj; trecv=ni*nj
              do i=1,ni; do j=1,nj; ijk=(i-1)*nj+j;
              sbuf(ijk)=fi(i,j,ke-ly)
              end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%corprev1)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%corprev1)) then
                 sbuf => dom(dom_indid(dom(ib)%corprev1))%rc1p
              else
                 sbuf => dom(ib) % sc1m
              end if

              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-1+pl1; j=jspr-1+pl2; k=kspr-1+pl3
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%corprev2)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%corprev2)) then
                 sbuf => dom(dom_indid(dom(ib)%corprev2))%rc2p
              else
                 sbuf => dom(ib) % sc2m
              end if

              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-1+pl1; j=jepr+1-pl2; k=kspr-1+pl3 
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%corprev3)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%corprev3)) then
                 sbuf => dom(dom_indid(dom(ib)%corprev3))%rc3p
              else
                 sbuf => dom(ib) % sc3m
              end if

              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+1-pl1; j=jepr+1-pl2; k=kspr-1+pl3;
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%corprev4)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%corprev4)) then
                 sbuf => dom(dom_indid(dom(ib)%corprev4))%rc4p
              else
                 sbuf => dom(ib) % sc4m
              end if

              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+1-pl1; j=jspr-1+pl2; k=kspr-1+pl3;
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%cornext1)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%cornext1)) then
                 sbuf => dom(dom_indid(dom(ib)%cornext1))%rc1m
              else
                 sbuf => dom(ib) % sc1p
              end if

              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+1-pl1; j=jepr+1-pl2; k=kepr+1-pl3
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%cornext2)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%cornext2)) then
                 sbuf => dom(dom_indid(dom(ib)%cornext2))%rc2m
              else
                 sbuf => dom(ib) % sc2p
              end if

              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+1-pl1; j=jspr-1+pl2; k=kepr+1-pl3
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%cornext3)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%cornext3)) then
                 sbuf => dom(dom_indid(dom(ib)%cornext3))%rc3m
              else
                 sbuf => dom(ib) % sc3p
              end if

              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-1+pl1; j=jspr-1+pl2; k=kepr+1-pl3
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%cornext4)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%cornext4)) then
                 sbuf => dom(dom_indid(dom(ib)%cornext4))%rc4m
              else
                 sbuf => dom(ib) % sc4p
              end if

              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-1+pl1; j=jepr+1-pl2; k=kepr+1-pl3
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev1)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev1)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev1))%re1p
              else
                 sbuf => dom(ib) % se1m
              end if

              tsend=nj*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nj
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-1+pl1; j=nn; k=kspr-1+pl2
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev2)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev2)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev2))%re2p
              else
                 sbuf => dom(ib) % se2m
              end if

              tsend=ni*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,ni
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jepr+1-pl1; k=kspr-1+pl2
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev3)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev3)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev3))%re3p
              else
                 sbuf => dom(ib) % se3m
              end if

              tsend=nj*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nj
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+1-pl1; j=nn; k=kspr-1+pl2
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev4)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev4)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev4))%re4p
              else
                 sbuf => dom(ib) % se4m
              end if

              tsend=ni*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,ni
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jspr-1+pl1; k=kspr-1+pl2
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev5)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev5)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev5))%re5p
              else
                 sbuf => dom(ib) % se5m
              end if

              tsend=nk*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nk
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-1+pl1; j=jspr-1+pl2; k=nn
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev6)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgprev6)) then
                 sbuf => dom(dom_indid(dom(ib)%edgprev6))%re6p
              else
                 sbuf => dom(ib) % se6m
              end if

              tsend=nk*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nk
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-1+pl1; j=jepr+1-pl2; k=nn
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext1)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext1)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext1))%re1m
              else
                 sbuf => dom(ib) % se1p
              end if

              tsend=nj*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nj
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+1-pl1; j=nn; k=kepr+1-pl2
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext2)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext2)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext2))%re2m
              else
                 sbuf => dom(ib) % se2p
              end if

              tsend=ni*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,ni
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jspr-1+pl1; k=kepr+1-pl2
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext3)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext3)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext3))%re3m
              else
                 sbuf => dom(ib) % se3p
              end if

              tsend=nj*pll; trecv=nj*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nj
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-1+pl1; j=nn; k=kepr+1-pl2
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext4)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext4)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext4))%re4m
              else
                 sbuf => dom(ib) % se4p
              end if

              tsend=ni*pll; trecv=ni*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,ni
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jepr+1-pl1; k=kepr+1-pl2
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext5)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext5)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext5))%re5m
              else
                 sbuf => dom(ib) % se5p
              end if

              tsend=nk*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nk
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+1-pl1; j=jepr+1-pl2; k=nn
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext6)) then

              if (dom_ad(dom_id(ib)).eq.dom_ad(dom(ib)%edgnext6)) then
                 sbuf => dom(dom_indid(dom(ib)%edgnext6))%re6m
              else
                 sbuf => dom(ib) % se6p
              end if

              tsend=nk*pll; trecv=nk*pll
              do pl1=0,pl; do pl2=0,pl; do nn=1,nk
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+1-pl1; j=jspr-1+pl2; k=nn
                 sbuf(ijk)=fi(i,j,k)
              end do; end do; end do

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
                 case (1)  
                    fi => dom(ib)%u
                    is=dom(ib)%isu; ie=dom(ib)%ieu
                    js=dom(ib)%jsu; je=dom(ib)%jeu
                    ks=dom(ib)%ksu; ke=dom(ib)%keu
                 case (2)  
                    fi => dom(ib)%v
                    is=dom(ib)%isv; ie=dom(ib)%iev
                    js=dom(ib)%jsv; je=dom(ib)%jev
                    ks=dom(ib)%ksv; ke=dom(ib)%kev
                 case (3)  
                    fi => dom(ib)%w
                    is=dom(ib)%isw; ie=dom(ib)%iew
                    js=dom(ib)%jsw; je=dom(ib)%jew
                    ks=dom(ib)%ksw; ke=dom(ib)%kew
                 case (11)  
                    fi => dom(ib)%ustar
                    is=dom(ib)%isu; ie=dom(ib)%ieu
                    js=dom(ib)%jsu; je=dom(ib)%jeu
                    ks=dom(ib)%ksu; ke=dom(ib)%keu
                 case (22)  
                    fi => dom(ib)%vstar
                    is=dom(ib)%isv; ie=dom(ib)%iev
                    js=dom(ib)%jsv; je=dom(ib)%jev
                    ks=dom(ib)%ksv; ke=dom(ib)%kev
                 case (33)  
                    fi => dom(ib)%wstar
                    is=dom(ib)%isw; ie=dom(ib)%iew
                    js=dom(ib)%jsw; je=dom(ib)%jew
                    ks=dom(ib)%ksw; ke=dom(ib)%kew
                 case (4) 
                    fi => dom(ib)%p
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (44) 
                    fi => dom(ib)%pp
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (5) 
                    fi => dom(ib)%T
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (6) 
                    fi => dom(ib)%ksgs
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (7) 
                    fi => dom(ib)%vis
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (8)
                    fi => dom(ib)%S
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
                 case (9)
                    fi => dom(ib)%eps
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (14)
			  fi => dom(ib)%phi
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (15)
		  	  fi => dom(ib)%phi_reinit
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (16)
		        fi => dom(ib)%phi_new
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (17)
		        fi => dom(ib)%phi_init
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (18)
		        fi => dom(ib)%dens
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
	           case (19)
		        fi => dom(ib)%mu
                    is=dom(ib)%isp; ie=dom(ib)%iep
                    js=dom(ib)%jsp; je=dom(ib)%jep
                    ks=dom(ib)%ksp; ke=dom(ib)%kep
              end select

!======================================================================
        if (ly.eq.0)  then
!======================================================================
!=====> previous cor #1
        if (dom(ib)%corprev1.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%corprev1)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev1)) then
                 call MPI_WAIT(dom(ib)%rq_c1m,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-pl1; j=jspr-pl2; k=kspr-pl3
                 if(i.lt.is .or. j.lt.js .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % rc1m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous cor #2
        if (dom(ib)%corprev2.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%corprev2)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev2)) then
                 call MPI_WAIT(dom(ib)%rq_c2m,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-pl1; j=jepr+pl2; k=kspr-pl3 
                 if(i.lt.is .or. j.gt.je .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % rc2m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous cor #3
        if (dom(ib)%corprev3.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%corprev3)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev3)) then
                 call MPI_WAIT(dom(ib)%rq_c3m,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+pl1; j=jepr+pl2; k=kspr-pl3
                 if(i.gt.ie .or. j.gt.je .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % rc3m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous cor #4
        if (dom(ib)%corprev4.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%corprev4)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%corprev4)) then
                 call MPI_WAIT(dom(ib)%rq_c4m,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+pl1; j=jspr-pl2; k=kspr-pl3 
                 if(i.gt.ie .or. j.lt.js .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % rc4m(ijk)
              end do; end do; end do
           end if
        end if
!=====> next cor #1
        if (dom(ib)%cornext1.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%cornext1)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext1)) then
                 call MPI_WAIT(dom(ib)%rq_c1p,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+pl1; j=jepr+pl2; k=kepr+pl3
                 if(i.gt.ie .or. j.gt.je .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % rc1p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next cor #2
        if (dom(ib)%cornext2.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%cornext2)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext2)) then
                 call MPI_WAIT(dom(ib)%rq_c2p,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=iepr+pl1; j=jspr-pl2; k=kepr+pl3
                 if(i.gt.ie .or. j.lt.js .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % rc2p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next cor #3
        if (dom(ib)%cornext3.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%cornext3)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext3)) then
                 call MPI_WAIT(dom(ib)%rq_c3p,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-pl1; j=jspr-pl2; k=kepr+pl3
                 if(i.lt.is .or. j.lt.js .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % rc3p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next cor #4
        if (dom(ib)%cornext4.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%cornext4)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%cornext4)) then
                 call MPI_WAIT(dom(ib)%rq_c4p,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do pl3=0,pl
                 ijk=pl1*pll+pl2*(pl+1)+pl3+1
                 i=ispr-pl1; j=jepr+pl2; k=kepr+pl3 
                 if(i.lt.is .or. j.gt.je .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % rc4p(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #1
        if (dom(ib)%edgprev1.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev1)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev1)) then
                 call MPI_WAIT(dom(ib)%rq_e1m,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=jspr,jepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-pl1; j=nn; k=kspr-pl2 
                 if(i.lt.is .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % re1m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #2
        if (dom(ib)%edgprev2.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev2)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev2)) then
                 call MPI_WAIT(dom(ib)%rq_e2m,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=ispr,iepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jepr+pl1; k=kspr-pl2
                 if(j.gt.je .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % re2m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #3
        if (dom(ib)%edgprev3.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev3)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev3)) then
                 call MPI_WAIT(dom(ib)%rq_e3m,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=jspr,jepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+pl1; j=nn; k=kspr-pl2
                 if(i.gt.ie .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % re3m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #4
        if (dom(ib)%edgprev4.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev4)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev4)) then
                 call MPI_WAIT(dom(ib)%rq_e4m,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=ispr,iepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jspr-pl1; k=kspr-pl2
                 if(j.lt.js .or. k.lt.ks)
     & fi(i,j,k)=dom(ib) % re4m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #5
        if (dom(ib)%edgprev5.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev5)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev5)) then
                 call MPI_WAIT(dom(ib)%rq_e5m,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=kspr,kepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-pl1 ; j=jspr-pl2; k=nn
                 if(i.lt.is .or. j.lt.js)
     & fi(i,j,k)=dom(ib) % re5m(ijk)
              end do; end do; end do
           end if
        end if
!=====> previous edge #6
        if (dom(ib)%edgprev6.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgprev6)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgprev6)) then
                 call MPI_WAIT(dom(ib)%rq_e6m,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=kspr,kepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-pl1 ; j=jepr+pl2; k=nn
                 if(i.lt.is .or. j.gt.je)
     & fi(i,j,k)=dom(ib) % re6m(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #1
        if (dom(ib)%edgnext1.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext1)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext1)) then
                 call MPI_WAIT(dom(ib)%rq_e1p,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=jspr,jepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+pl1; j=nn; k=kepr+pl2
                 if(i.gt.ie .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % re1p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #2
        if (dom(ib)%edgnext2.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext2)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext2)) then
                 call MPI_WAIT(dom(ib)%rq_e2p,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=ispr,iepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jspr-pl1; k=kepr+pl2
                 if(j.lt.js .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % re2p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #3
        if (dom(ib)%edgnext3.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext3)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext3)) then
                 call MPI_WAIT(dom(ib)%rq_e3p,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=jspr,jepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=ispr-pl1; j=nn; k=kepr+pl2
                 if(i.lt.is .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % re3p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #4
        if (dom(ib)%edgnext4.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext4)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext4)) then
                 call MPI_WAIT(dom(ib)%rq_e4p,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=ispr,iepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=nn; j=jepr+pl1; k=kepr+pl2
                 if(j.gt.je .or. k.gt.ke)
     & fi(i,j,k)=dom(ib) % re4p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #5
        if (dom(ib)%edgnext5.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext5)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext5)) then
                 call MPI_WAIT(dom(ib)%rq_e5p,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=kspr,kepr
                 ijk=(nn-1)*pll+pl1*(pl+1)+pl2+1
                 i=iepr+pl1 ; j=jepr+pl2; k=nn
                 if(i.gt.ie .or. j.gt.je)
     & fi(i,j,k)=dom(ib) % re5p(ijk)
              end do; end do; end do
           end if
        end if
!=====> next edge #6
        if (dom(ib)%edgnext6.ge.0) then
           if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%edgnext6)) then
              if (dom_ad(dom_id(ib)).ne.dom_ad(dom(ib)%edgnext6)) then
                 call MPI_WAIT(dom(ib)%rq_e6p,MPI_STATUS_IGNORE,ierr)
              end if
              do pl1=0,pl; do pl2=0,pl; do nn=kspr,kepr
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
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%iprev)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%iprev)) then
                    call MPI_WAIT(dom(ib)%rq_m1,MPI_STATUS_IGNORE,ierr)
                 end if
                 do k=1,nk; do j=1,nj; ijk=(k-1)*nj+j
                    dom(ib)%stfcinf(1,ly+1,ijk)=dom(ib) % recvb_m1(ijk)
                 end do; end do
                 do k=kspr,kepr; do j=jspr,jepr; ijk=(k-1)*nj+j
                    fi(is-1-ly,j,k)=dom(ib) % recvb_m1(ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%jprev.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%jprev)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%jprev)) then
                    call MPI_WAIT(dom(ib)%rq_m2,MPI_STATUS_IGNORE,ierr)
                 end if
                 do k=1,nk; do i=1,ni; ijk=(k-1)*ni+i
                    dom(ib)%stfcinf(3,ly+1,ijk)=dom(ib) % recvb_m2(ijk)
                 end do; end do
                 do k=kspr,kepr; do i=ispr,iepr; ijk=(k-1)*ni+i
                    fi(i,js-1-ly,k)=dom(ib) % recvb_m2(ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%kprev.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%kprev)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%kprev)) then
                    call MPI_WAIT(dom(ib)%rq_m3,MPI_STATUS_IGNORE,ierr)
                 end if
                 do i=1,ni; do j=1,nj; ijk=(i-1)*nj+j
                    dom(ib)%stfcinf(5,ly+1,ijk)=dom(ib) % recvb_m3(ijk)
                 end do; end do
                 do i=ispr,iepr; do j=jspr,jepr; ijk=(i-1)*nj+j
                    fi(i,j,ks-1-ly)=dom(ib) % recvb_m3(ijk)
                 end do; end do
              end if
              end if

!..............................................................................
!=== Next Neighbor  ===> 
!..............................................................................
              if (dom(ib)%inext.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%inext)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%inext)) then
                    call MPI_WAIT(dom(ib)%rq_p1,MPI_STATUS_IGNORE,ierr)
                 end if
                 do k=1,nk; do j=1,nj; ijk=(k-1)*nj+j
                    dom(ib)%stfcinf(2,ly+1,ijk)=dom(ib) % recvb_p1(ijk)
                 end do; end do
                 do k=kspr,kepr; do j=jspr,jepr; ijk=(k-1)*nj+j
                    fi(ie+1+ly,j,k)=dom(ib) % recvb_p1(ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%jnext.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%jnext)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%jnext)) then
                    call MPI_WAIT(dom(ib)%rq_p2,MPI_STATUS_IGNORE,ierr)
                 end if
                 do k=1,nk; do i=1,ni; ijk=(k-1)*ni+i
                    dom(ib)%stfcinf(4,ly+1,ijk)=dom(ib) % recvb_p2(ijk)
                 end do; end do
                 do k=kspr,kepr; do i=ispr,iepr; ijk=(k-1)*ni+i
                    fi(i,je+1+ly,k)=dom(ib) % recvb_p2(ijk)
                 end do; end do
              end if
              end if

              if (dom(ib)%knext.ge.0)  then
              if(rdiv(dom_id(ib)).eq.rdiv(dom(ib)%knext)) then
                 if (dom_ad(dom_id(ib)) .ne. dom_ad(dom(ib)%knext)) then
                    call MPI_WAIT(dom(ib)%rq_p3,MPI_STATUS_IGNORE,ierr)
                 end if
                 do i=1,ni; do j=1,nj; ijk=(i-1)*nj+j
                    dom(ib)%stfcinf(6,ly+1,ijk)=dom(ib) % recvb_p3(ijk)
                 end do; end do
                 do i=ispr,iepr; do j=jspr,jepr; ijk=(i-1)*nj+j
                    fi(i,j,ke+1+ly)=dom(ib) % recvb_p3(ijk)
                 end do; end do
              end if
              end if

!--------------------------------------------------------------------------
           end do
!==========================================================================
        end do

        return
        end subroutine exchange_smlvl
!########################################################################## 
