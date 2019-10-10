!##########################################################################
        subroutine localparameters
!##########################################################################
        use vars
        use mpi
        use multidata
        implicit none
        integer :: i,j,ib,nicell,njcell,nkcell
        integer :: buffer_nmax,nmax,buffer_nemax,nemax,ni,nj,nk
        integer :: buffer_rdiv,maxrdiv,act_ngrid
        integer :: glevel,gl,mgc_i,mgc_j,mgc_k,pow2
        real    :: n1,n2,n3
        integer,allocatable,dimension(:,:) :: b_rv1,b_rv2
        nmax=0
        nemax=0

        if(LMR.eq.1) then
           pl=1+pl_ex
        else if(LMR.eq.2) then
           pl=2+pl_ex
        else
           print*,'error: wrong LMR selection!!, STOP'
           stop
        end if

        buffer_rdiv=1
        do ib=1,nbp
           buffer_rdiv = max(buffer_rdiv,rdiv(dom_id(ib)))
        end do
        call MPI_ALLREDUCE(buffer_rdiv,maxrdiv,1,MPI_INTEGER,MPI_MAX,
     &MPI_COMM_WORLD,ierr)

        act_ngrid=ngrid_input
        ngrd_gl=ngrid_input+int(log(real(maxrdiv))/log(2.0))
        allocate (rdv(0:num_domains-1,ngrd_gl))
        allocate (b_rv1(num_domains,ngrd_gl),b_rv2(num_domains,ngrd_gl))

        rdv=0

        do ib=1,nbp

        if(rdiv(dom_id(ib)).gt.1) then   
           dom(ib)%ngrid=ngrid_input+
     & int(log(real(rdiv(dom_id(ib))))/log(2.0))
        end if

        do glevel=1,ngrd_gl
           if(glevel.gt.dom(ib)%ngrid) then
              rdv(dom_id(ib),glevel)=1
           else
              if(glevel.le.act_ngrid) then
                 rdv(dom_id(ib),glevel)=rdiv(dom_id(ib))
              else
                 rdv(dom_id(ib),glevel)=int(rdv(dom_id(ib),glevel-1)/2)
              end if
           end if
        end do

        end do


        do i=0,num_domains-1
           do glevel=1,ngrd_gl
              b_rv1(i+1,glevel)= rdv(i,glevel)
           end do
        end do

        call MPI_ALLREDUCE(b_rv1,b_rv2,num_domains*ngrd_gl,
     & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)


        do i=0,num_domains-1
           do glevel=1,ngrd_gl
             rdv(i,glevel)= b_rv2(i+1,glevel)
           end do
        end do

        do ib=1,nbp

           n1=(dom(ib)%xel-dom(ib)%xsl)/g_dx
           n2=(dom(ib)%yel-dom(ib)%ysl)/g_dy
           n3=(dom(ib)%zel-dom(ib)%zsl)/g_dz

           if (abs(n1-nint(n1)).gt.1e-8) then 
              print*,'error1 in localparameters',
     & dom_id(ib),(dom(ib)%xel-dom(ib)%xsl),g_dx,n1
              stop
           end if
           if (abs(n2-nint(n2)).gt.1e-8) then 
              print*,'error2 in localparameters',
     & dom_id(ib),(dom(ib)%yel-dom(ib)%ysl),g_dy,n2
              stop
           end if
           if (abs(n3-nint(n3)).gt.1e-8) then 
              print*,'error3 in localparameters',
     & dom_id(ib),(dom(ib)%zel-dom(ib)%zsl),g_dz,n3
              stop
           end if

           nicell=rdiv(dom_id(ib))*nint(n1)
           njcell=rdiv(dom_id(ib))*nint(n2)
           nkcell=rdiv(dom_id(ib))*nint(n3)

           dom(ib)%isp=1+pl
           dom(ib)%isu=dom(ib)%isp
           dom(ib)%isv=dom(ib)%isp
           dom(ib)%isw=dom(ib)%isp

           dom(ib)%iep=nicell+dom(ib)%isp-1 
           dom(ib)%ieu=dom(ib)%iep
           dom(ib)%iev=dom(ib)%iep
           dom(ib)%iew=dom(ib)%iep

           dom(ib)%jsp=1+pl
           dom(ib)%jsu=dom(ib)%jsp
           dom(ib)%jsv=dom(ib)%jsp
           dom(ib)%jsw=dom(ib)%jsp

           dom(ib)%jep=njcell+dom(ib)%jsp-1 
           dom(ib)%jeu=dom(ib)%jep
           dom(ib)%jev=dom(ib)%jep
           dom(ib)%jew=dom(ib)%jep

           dom(ib)%ksp=1+pl
           dom(ib)%ksu=dom(ib)%ksp
           dom(ib)%ksv=dom(ib)%ksp
           dom(ib)%ksw=dom(ib)%ksp

           dom(ib)%kep=nkcell+dom(ib)%ksp-1 
           dom(ib)%keu=dom(ib)%kep
           dom(ib)%kev=dom(ib)%kep
           dom(ib)%kew=dom(ib)%kep

           ipref=dom(ib)%isp
           jpref=dom(ib)%jsp
           kpref=dom(ib)%ksp

           dom(ib)%ttc_i=nicell+2*pl
           dom(ib)%ttc_j=njcell+2*pl
           dom(ib)%ttc_k=nkcell+2*pl
           dom(ib)%ttc_ijk=dom(ib)%ttc_i*dom(ib)%ttc_j*dom(ib)%ttc_k

!=========================================================================
           if(LMR.eq.2) then

           if (dom(ib)%iprev.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%iprev)) then
                 dom(ib)%isu=dom(ib)%isp-1
              end if
           end if

           if (dom(ib)%inext.ge.0) then
              if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%inext)) then
                 dom(ib)%ieu=dom(ib)%ieu-1
              end if
           end if

           if (dom(ib)%jprev.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%jprev)) then
                 dom(ib)%jsv=dom(ib)%jsp-1
              end if
           end if

           if (dom(ib)%jnext.ge.0) then
              if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%jnext)) then
                 dom(ib)%jev=dom(ib)%jev-1
              end if
           end if

           if (dom(ib)%kprev.ge.0) then
              if(rdiv(dom_id(ib)).gt.rdiv(dom(ib)%kprev)) then
                 dom(ib)%ksw=dom(ib)%ksp-1
              end if
           end if

           if (dom(ib)%knext.ge.0) then
              if(rdiv(dom_id(ib)).lt.rdiv(dom(ib)%knext)) then
                 dom(ib)%kew=dom(ib)%kew-1
              end if
           end if
 
           end if
!=========================================================================

           if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5)
     & dom(ib)%ieu=dom(ib)%ieu-1
           if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5)
     & dom(ib)%jev=dom(ib)%jev-1 
           if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5)
     & dom(ib)%kew=dom(ib)%kew-1

           ni=dom(ib)%ttc_i
           nj=dom(ib)%ttc_j
           nk=dom(ib)%ttc_k
           nmax=max(ni*nj,ni*nk,nj*nk,nmax)
           nemax=max(ni,nj,nk,nemax)

           dom(ib)%nwork=0
           dom(ib)%nvars=0
           if(solver.eq.2) then
              pow2=2**(dom(ib)%ngrid-1)
              if(mod(nicell,pow2).ne.0 .or. mod(njcell,pow2).ne.0 
     & .or. mod(nkcell,pow2).ne.0) then
                 print*,'2:the given number not correct for mgsolver!'
                 stop
              end if

              allocate (dom(ib)%cntp(ngrd_gl))
              dom(ib)%cntp(1)=0
              do glevel=1,ngrd_gl-1
                 if(glevel.gt.dom(ib)%ngrid) then
                    gl=dom(ib)%ngrid
                 else
                    gl=glevel
                 end if
           mgc_i=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
           mgc_j=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
           mgc_k=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2
                 dom(ib)%cntp(glevel+1)=dom(ib)%cntp(glevel)+
     & mgc_i*mgc_j*mgc_k
              end do

              do glevel=1,ngrd_gl
                 if(glevel.gt.dom(ib)%ngrid) then
                    gl=dom(ib)%ngrid
                 else
                    gl=glevel
                 end if
           mgc_i=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
           mgc_j=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
           mgc_k=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2
                 dom(ib)%nwork=dom(ib)%nwork+mgc_i*mgc_j*mgc_k+
     & 8*(mgc_i-2)*(mgc_j-2)*(mgc_k-2)
                 dom(ib)%nvars=dom(ib)%nvars+mgc_i*mgc_j*mgc_k
              end do

		  allocate (dom(ib)%cof(dom(ib)%nwork)) 

           end if


           write (6,88) 'mydom#:',dom_id(ib),' i:',dom(ib)%ttc_i,
     & ' j:',dom(ib)%ttc_j,' k:',dom(ib)%ttc_k

           dom(ib)%dx=g_dx/real(rdiv(dom_id(ib)))
           dom(ib)%dy=g_dy/real(rdiv(dom_id(ib)))
           dom(ib)%dz=g_dz/real(rdiv(dom_id(ib)))
 
        end do



        buffer_nemax = nemax
        call MPI_ALLREDUCE(buffer_nemax,nemax,1,MPI_INTEGER,MPI_MAX,
     &MPI_COMM_WORLD,ierr)

        buffer_nmax = nmax
        call MPI_ALLREDUCE(buffer_nmax,nmax,1,MPI_INTEGER,MPI_MAX,
     &MPI_COMM_WORLD,ierr)

        nge=nemax*(pl+1)**2
        ngc=(pl+1)**3
        ngg=nmax

88      format (a,i5,a,i4,a,i4,a,i4)

        end subroutine localparameters
!##########################################################################
