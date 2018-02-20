!#######################################################################
        subroutine  exchangepp(g)
!#######################################################################
        use multidata
        use mpi
        use vars
        implicit none
        integer :: i,j,k,ijk,sync_dir,ib
        integer :: g,tag,ta,gl
        integer :: ni,nj,nk,ic,jc,kc,ii,jj,kk,ll
        integer :: nijk,nij,nn
        integer :: cpu_next,cpu_prev,tsend,trecv
        integer :: ns1,ns2,nsr1,nsr2,nsc1,nsc2
        integer :: nic,nif,njc,njf,nkc,nkf
        double precision    :: chc,chc1,chc2,chc3
        double precision, pointer, dimension(:) :: fi
        double precision, pointer, dimension(:) :: sbuf_m,sbuf_p
        double precision, pointer, dimension(:) :: rbuf_m,rbuf_p
        double precision, pointer, dimension(:) :: sbufc1m,sbufc1p
        double precision, pointer, dimension(:) :: rbufc1m,rbufc1p
        double precision, pointer, dimension(:) :: sbufc2m,sbufc2p
        double precision, pointer, dimension(:) :: rbufc2m,rbufc2p
        double precision, pointer, dimension(:) :: sbufc3m,sbufc3p
        double precision, pointer, dimension(:) :: rbufc3m,rbufc3p
        double precision, pointer, dimension(:) :: sbufc4m,sbufc4p
        double precision, pointer, dimension(:) :: rbufc4m,rbufc4p
        double precision, pointer, dimension(:) :: sbufe1m,sbufe1p
        double precision, pointer, dimension(:) :: rbufe1m,rbufe1p
        double precision, pointer, dimension(:) :: sbufe2m,sbufe2p
        double precision, pointer, dimension(:) :: rbufe2m,rbufe2p
        double precision, pointer, dimension(:) :: sbufe3m,sbufe3p
        double precision, pointer, dimension(:) :: rbufe3m,rbufe3p
        double precision, pointer, dimension(:) :: sbufe4m,sbufe4p
        double precision, pointer, dimension(:) :: rbufe4m,rbufe4p
        double precision, pointer, dimension(:) :: sbufe5m,sbufe5p
        double precision, pointer, dimension(:) :: rbufe5m,rbufe5p
        double precision, pointer, dimension(:) :: sbufe6m,sbufe6p
        double precision, pointer, dimension(:) :: rbufe6m,rbufe6p

        MPI_FLT   = MPI_DOUBLE_PRECISION


        do sync_dir = 1,3

!=======================================================================
!=======================================================================
           do ib=1,nbp
              if(g.gt.dom(ib)%ngrid) then
                 gl=dom(ib)%ngrid
              else
                 gl=g
              end if
        ni=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        nj=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nk=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

        fi => dom(ib)%cof
        nijk=dom(ib)%faz(g)-ni*nj*nk
        nij=ni*nj

              nif=2.0*(ni-2)+2; njf=2.0*(nj-2)+2; nkf=2.0*(nk-2)+2
              chc1=real(ni-2); chc2=real(nj-2); chc3=real(nk-2)
              if(rdiv(dom_id(ib)).eq.1) then
                 nic=ni; njc=nj; nkc=nk
              else
        if((mod(chc1,2.0).ne.0.0).and.g.lt.dom(ib)%ngrid) print*,'er-i'
        if((mod(chc2,2.0).ne.0.0).and.g.lt.dom(ib)%ngrid) print*,'er-j'
        if((mod(chc3,2.0).ne.0.0).and.g.lt.dom(ib)%ngrid) print*,'er-k'
                 nic=int((ni-2)/2)+2; njc=int((nj-2)/2)+2; 
                 nkc=int((nk-2)/2)+2
              end if


              if (sync_dir.eq.1)  then
                 cpu_prev=dom(ib)%iprev
                 cpu_next=dom(ib)%inext
                 ns1=nj  
                 ns2=nk 

                 nsr1=2.0*(ns1-2)+2
                 nsr2=2.0*(ns2-2)+2

                 chc=real(ns1-2)
                 if(rdiv(dom_id(ib)).eq.1) then
                    nsc1=ns1
                 else
                    nsc1=int((ns1-2)/2)+2
                 end if

                 chc=real(ns2-2)
                 if(rdiv(dom_id(ib)).eq.1) then
                    nsc2=ns2
                 else
                    nsc2=int((ns2-2)/2)+2
                 end if 

              else if (sync_dir.eq.2)  then
                 cpu_prev=dom(ib)%jprev
                 cpu_next=dom(ib)%jnext
                 ns1=ni 
                 ns2=nk 

                 nsr1=2.0*(ns1-2)+2
                 nsr2=2.0*(ns2-2)+2

                 chc=real(ns1-2)
                 if(rdiv(dom_id(ib)).eq.1) then
                    nsc1=ns1
                 else
                    nsc1=int((ns1-2)/2)+2
                 end if

                 chc=real(ns2-2)
                 if(rdiv(dom_id(ib)).eq.1) then
                    nsc2=ns2
                 else
                    nsc2=int((ns2-2)/2)+2
                 end if 

              else if (sync_dir.eq.3)  then
                 cpu_prev=dom(ib)%kprev
                 cpu_next=dom(ib)%knext
                 ns1=nj  
                 ns2=ni 

                 nsr1=2.0*(ns1-2)+2
                 nsr2=2.0*(ns2-2)+2

                 chc=real(ns1-2)
                 if(rdiv(dom_id(ib)).eq.1) then
                    nsc1=ns1
                 else
                    nsc1=int((ns1-2)/2)+2
                 end if

                 chc=real(ns2-2)
                 if(rdiv(dom_id(ib)).eq.1) then
                    nsc2=ns2
                 else
                    nsc2=int((ns2-2)/2)+2
                 end if 

              end if

!..........................................................................
!=== Previous Neighbor  ===> 
!..........................................................................
        if (cpu_prev.ge.0) then

           if (dom_ad(dom_id(ib)) .eq. dom_ad(cpu_prev)) then

              if (sync_dir.eq.1)  then

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                    do k=1,ns2
           do j=1,ns1
           ijk=(k-1)*ns1+j
           i=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=fi(ll)
           end do
                    end do
                 else
                    if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
           do k=2,nsc2-1
              do j=2,nsc1-1
                 ijk=(k-1)*nsc1+j
                 ii=2; jj=2*j-2; kk=2*k-2
                 ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
                 dom(dom_indid(cpu_prev))%recvb_p1(ijk)=0.125*
     & (fi(ll)+fi(ll+1)+fi(ll+nij)+fi(ll+1+nij)+
     & fi(ll+nj)+fi(ll+1+nj)+fi(ll+nj+nij)+fi(ll+1+nj+nij))
              end do
           end do
                    else
                    do kc=2,ns2-1
           do jc=2,ns1-1
           ic=2; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk

           k=2*kc-2; j=2*jc-2; ijk=(k-1)*nsr1+j
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
           k=2*kc-2; j=2*jc-1; ijk=(k-1)*nsr1+j
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
           k=2*kc-1; j=2*jc-2; ijk=(k-1)*nsr1+j
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
           k=2*kc-1; j=2*jc-1; ijk=(k-1)*nsr1+j
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
           end do
                    end do
                    end if
                 end if
              else if (sync_dir.eq.2)  then

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                    do k=1,ns2
           do i=1,ns1
           ijk=(k-1)*ns1+i
           j=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=fi(ll)
           end do
                    end do
                 else
                    if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
           do k=2,nsc2-1
              do i=2,nsc1-1
                 ijk=(k-1)*nsc1+i
                 ii=2*i-2; jj=2; kk=2*k-2
                 ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
                 dom(dom_indid(cpu_prev))%recvb_p1(ijk)=0.125*
     & (fi(ll)+fi(ll+nj)+fi(ll+nij)+fi(ll+nj+nij)+
     & fi(ll+1)+fi(ll+1+nj)+fi(ll+1+nij)+fi(ll+1+nj+nij))
              end do
           end do
                    else
                 do kc=2,ns2-1
                    do ic=2,ns1-1
           jc=2; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk

           k=2*kc-2; i=2*ic-2; ijk=(k-1)*nsr1+i
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
           k=2*kc-2; i=2*ic-1; ijk=(k-1)*nsr1+i
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
           k=2*kc-1; i=2*ic-2; ijk=(k-1)*nsr1+i
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
           k=2*kc-1; i=2*ic-1; ijk=(k-1)*nsr1+i
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
                    end do
                 end do
                    end if
                 end if
              else

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                    do i=1,ns2
           do j=1,ns1
           ijk=(i-1)*ns1+j
           k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=fi(ll)
           end do
                    end do
                 else
                    if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
           do i=2,nsc2-1
              do j=2,nsc1-1
                 ijk=(i-1)*nsc1+j
                 ii=2*i-2; jj=2*j-2; kk=2
                 ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
                 dom(dom_indid(cpu_prev))%recvb_p1(ijk)=0.125*
     & (fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))
              end do
           end do
                    else
                 do ic=2,ns2-1
                    do jc=2,ns1-1
           kc=2; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk
           i=2*ic-2; j=2*jc-2; ijk=(i-1)*nsr1+j
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
           i=2*ic-2; j=2*jc-1; ijk=(i-1)*nsr1+j
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
           i=2*ic-1; j=2*jc-2; ijk=(i-1)*nsr1+j
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
           i=2*ic-1; j=2*jc-1; ijk=(i-1)*nsr1+j
           dom(dom_indid(cpu_prev))%recvb_p1(ijk)=2.0d0*fi(ll)
                    end do
                 end do
                    end if
                 end if
!======================================================================
!======================================================================
              end if


           else !if (dom_ad(dom_id(ib)) .ne. dom_ad(cpu_prev)) then

              sbuf_m => dom(ib) % sendb_m1
              rbuf_m => dom(ib) % recvb_m1

!======================================================================
!======================================================================
              if (sync_dir.eq.1)  then

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                    tsend=ns1*ns2
                    trecv=ns1*ns2
                    do k=1,ns2
           do j=1,ns1
           ijk=(k-1)*ns1+j
           i=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
           sbuf_m(ijk)=fi(ll)
           end do
                    end do
                 else
                    if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
           tsend=nsc1*nsc2
           trecv=ns1*ns2
           do k=2,nsc2-1
              do j=2,nsc1-1
                 ijk=(k-1)*nsc1+j
                 ii=2; jj=2*j-2; kk=2*k-2
                 ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
                 sbuf_m(ijk)=0.125*
     & (fi(ll)+fi(ll+1)+fi(ll+nij)+fi(ll+1+nij)+
     & fi(ll+nj)+fi(ll+1+nj)+fi(ll+nj+nij)+fi(ll+1+nj+nij))
              end do
           end do
                    else
           tsend=nsr1*nsr2
           trecv=ns1*ns2
                    do kc=2,ns2-1
           do jc=2,ns1-1
           ic=2; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk

           k=2*kc-2; j=2*jc-2; ijk=(k-1)*nsr1+j
           sbuf_m(ijk)=2.0d0*fi(ll)
           k=2*kc-2; j=2*jc-1; ijk=(k-1)*nsr1+j
           sbuf_m(ijk)=2.0d0*fi(ll)
           k=2*kc-1; j=2*jc-2; ijk=(k-1)*nsr1+j
           sbuf_m(ijk)=2.0d0*fi(ll)
           k=2*kc-1; j=2*jc-1; ijk=(k-1)*nsr1+j
           sbuf_m(ijk)=2.0d0*fi(ll)
           end do
                    end do
                    end if
                 end if
!======================================================================
!======================================================================
              else if (sync_dir.eq.2)  then

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                    tsend=ns1*ns2
                    trecv=ns1*ns2
                    do k=1,ns2
           do i=1,ns1
           ijk=(k-1)*ns1+i
           j=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
           sbuf_m(ijk)=fi(ll)
           end do
                    end do
                 else
                    if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
           tsend=nsc1*nsc2
           trecv=ns1*ns2
           do k=2,nsc2-1
              do i=2,nsc1-1
                 ijk=(k-1)*nsc1+i
                 ii=2*i-2; jj=2; kk=2*k-2
                 ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
                 sbuf_m(ijk)=0.125*
     & (fi(ll)+fi(ll+nj)+fi(ll+nij)+fi(ll+nj+nij)+
     & fi(ll+1)+fi(ll+1+nj)+fi(ll+1+nij)+fi(ll+1+nj+nij))
              end do
           end do
                    else
           tsend=nsr1*nsr2
           trecv=ns1*ns2
                 do kc=2,ns2-1
                    do ic=2,ns1-1
           jc=2; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk

           k=2*kc-2; i=2*ic-2; ijk=(k-1)*nsr1+i
           sbuf_m(ijk)=2.0d0*fi(ll)
           k=2*kc-2; i=2*ic-1; ijk=(k-1)*nsr1+i
           sbuf_m(ijk)=2.0d0*fi(ll)
           k=2*kc-1; i=2*ic-2; ijk=(k-1)*nsr1+i
           sbuf_m(ijk)=2.0d0*fi(ll)
           k=2*kc-1; i=2*ic-1; ijk=(k-1)*nsr1+i
           sbuf_m(ijk)=2.0d0*fi(ll)
                    end do
                 end do
                    end if
                 end if
!======================================================================
!======================================================================
              else

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                    tsend=ns1*ns2
                    trecv=ns1*ns2
                    do i=1,ns2
           do j=1,ns1
           ijk=(i-1)*ns1+j
           k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
           sbuf_m(ijk)=fi(ll)
           end do
                    end do
                 else
                    if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
           tsend=nsc1*nsc2
           trecv=ns1*ns2
           do i=2,nsc2-1
              do j=2,nsc1-1
                 ijk=(i-1)*nsc1+j
                 ii=2*i-2; jj=2*j-2; kk=2
                 ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
                 sbuf_m(ijk)=0.125*
     & (fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))
              end do
           end do
                    else
           tsend=nsr1*nsr2
           trecv=ns1*ns2
                 do ic=2,ns2-1
                    do jc=2,ns1-1
           kc=2; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk
           i=2*ic-2; j=2*jc-2; ijk=(i-1)*nsr1+j
           sbuf_m(ijk)=2.0d0*fi(ll)
           i=2*ic-2; j=2*jc-1; ijk=(i-1)*nsr1+j
           sbuf_m(ijk)=2.0d0*fi(ll)
           i=2*ic-1; j=2*jc-2; ijk=(i-1)*nsr1+j
           sbuf_m(ijk)=2.0d0*fi(ll)
           i=2*ic-1; j=2*jc-1; ijk=(i-1)*nsr1+j
           sbuf_m(ijk)=2.0d0*fi(ll)
                    end do
                 end do
                    end if
                 end if
!======================================================================
!======================================================================
              end if

              tag=(2*sync_dir-1)*10**5+cpu_prev
              ta=2*sync_dir
              call MPI_IRECV  (dom(ib)%recvb_m1(1),trecv,MPI_FLT,
     &dom_ad(cpu_prev),dom(ib)%tg(ta),MPI_COMM_WORLD,dom(ib)%rq_m1,ierr)
              call MPI_SEND (dom(ib)%sendb_m1(1),tsend,MPI_FLT,
     &dom_ad(cpu_prev),tag,MPI_COMM_WORLD,ierr)


                 end if
              end if
!======================================================================
!======================================================================
        if (sync_dir.eq.3)  then

!......................................................................
!=== Previous Corner Neighbors  ===> 
!......................................................................

!=====> previous cor #1
        if (dom(ib)%corprev1.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%corprev1)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%corprev1,g)) then
              i=2; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
              dom(dom_indid(dom(ib)%corprev1))%rc1p(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%corprev1,g)) then
                 i=2; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%corprev1))%rc1p(1)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=2; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%corprev1))%rc1p(1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj-1-nij))/64.0d0
              end if
              end if
           else
              sbufc1m => dom(ib) % sc1m
              rbufc1m => dom(ib) % rc1m
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%corprev1,g)) then
              i=2; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
              sbufc1m(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%corprev1,g)) then
                 i=2; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc1m(1)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=2; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc1m(1)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj-1-nij))/64.0d0
              end if
              end if
              tag=7*10**5+dom(ib)%corprev1
              call MPI_IRECV  (dom(ib)%rc1m,1,MPI_FLT,
     &dom_ad(dom(ib)%corprev1),dom(ib)%tg(8),
     &MPI_COMM_WORLD,dom(ib)%rq_c1m,ierr)
              call MPI_SEND (dom(ib)%sc1m,1,MPI_FLT,
     &dom_ad(dom(ib)%corprev1),tag,MPI_COMM_WORLD,ierr)

           end if
        end if

!=====> previous cor #2
        if (dom(ib)%corprev2.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%corprev2)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%corprev2,g)) then
              i=2; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
              dom(dom_indid(dom(ib)%corprev2))%rc2p(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%corprev2,g)) then
                 i=2; j=nj-2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%corprev2))%rc2p(1)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=2; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%corprev2))%rc2p(1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj+1-nij))/64.0d0
              end if
              end if
           else
              sbufc2m => dom(ib) % sc2m
              rbufc2m => dom(ib) % rc2m
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%corprev2,g)) then
              i=2; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
              sbufc2m(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%corprev2,g)) then
                 i=2; j=nj-2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc2m(1)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=2; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc2m(1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj+1-nij))/64.0d0
              end if
              end if
              tag=9*10**5+dom(ib)%corprev2
              call MPI_IRECV  (dom(ib)%rc2m,1,MPI_FLT,
     &dom_ad(dom(ib)%corprev2),dom(ib)%tg(10),
     &MPI_COMM_WORLD,dom(ib)%rq_c2m,ierr)
              call MPI_SEND (dom(ib)%sc2m,1,MPI_FLT,
     &dom_ad(dom(ib)%corprev2),tag,MPI_COMM_WORLD,ierr)

           end if
        end if

!=====> previous cor #3
        if (dom(ib)%corprev3.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%corprev3)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%corprev3,g)) then
              i=ni-1; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
              dom(dom_indid(dom(ib)%corprev3))%rc3p(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%corprev3,g)) then
                 i=ni-2; j=nj-2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%corprev3))%rc3p(1)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=ni-1; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%corprev3))%rc3p(1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj+1-nij))/64.0d0
              end if
              end if
           else
              sbufc3m => dom(ib) % sc3m
              rbufc3m => dom(ib) % rc3m
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%corprev3,g)) then
              i=ni-1; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
              sbufc3m(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%corprev3,g)) then
                 i=ni-2; j=nj-2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc3m(1)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=ni-1; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc3m(1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj+1-nij))/64.0d0
              end if
              end if
              tag=11*10**5+dom(ib)%corprev3
              call MPI_IRECV  (dom(ib)%rc3m,1,MPI_FLT,
     &dom_ad(dom(ib)%corprev3),dom(ib)%tg(12),
     &MPI_COMM_WORLD,dom(ib)%rq_c3m,ierr)
              call MPI_SEND (dom(ib)%sc3m,1,MPI_FLT,
     &dom_ad(dom(ib)%corprev3),tag,MPI_COMM_WORLD,ierr)

           end if
        end if

!=====> previous cor #4
        if (dom(ib)%corprev4.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%corprev4)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%corprev4,g)) then
              i=ni-1; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
              dom(dom_indid(dom(ib)%corprev4))%rc4p(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%corprev4,g)) then
                 i=ni-2; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%corprev4))%rc4p(1)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=ni-1; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%corprev4))%rc4p(1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj-1-nij))/64.0d0
              end if
              end if
           else
              sbufc4m => dom(ib) % sc4m
              rbufc4m => dom(ib) % rc4m
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%corprev4,g)) then
              i=ni-1; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
              sbufc4m(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%corprev4,g)) then
                 i=ni-2; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc4m(1)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=ni-1; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc4m(1)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj-1-nij))/64.0d0
              end if
              end if
              tag=13*10**5+dom(ib)%corprev4
              call MPI_IRECV  (dom(ib)%rc4m,1,MPI_FLT,
     &dom_ad(dom(ib)%corprev4),dom(ib)%tg(14),
     &MPI_COMM_WORLD,dom(ib)%rq_c4m,ierr)
              call MPI_SEND (dom(ib)%sc4m,1,MPI_FLT,
     &dom_ad(dom(ib)%corprev4),tag,MPI_COMM_WORLD,ierr)

           end if
        end if


!..........................................................................
!=== Previous Edge Neighbors  ===> 
!..........................................................................

!=====> previous edge #1
        if (dom(ib)%edgprev1.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev1)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev1,g)) then
                 do nn=2,nj-1
                    i=2; j=nn; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev1))%re1p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev1,g)) then
                 do nn=2,njc-1
                    i=2; j=2*nn-2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev1))%re1p(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,nj-1
                    i=2; j=nn; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev1))%re1p(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj-1-nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgprev1))%re1p(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj+1-nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe1m => dom(ib) % se1m
              rbufe1m => dom(ib) % re1m
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev1,g)) then
                 tsend=nj
                 trecv=nj
                 do nn=2,nj-1
                    i=2; j=nn; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe1m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev1,g)) then
                 tsend=njc
                 trecv=nj
                 do nn=2,njc-1
                    i=2; j=2*nn-2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe1m(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=njf
                 trecv=nj
                 do nn=2,nj-1
                    i=2; j=nn; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe1m(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj-1-nij))/64.0d0
                    sbufe1m(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj+1-nij))/64.0d0
                 end do
              end if
              end if
              tag=15*10**5+dom(ib)%edgprev1
              call MPI_IRECV  (dom(ib)%re1m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev1),dom(ib)%tg(16),
     &MPI_COMM_WORLD,dom(ib)%rq_e1m,ierr)
              call MPI_SEND (dom(ib)%se1m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev1),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!=====> previous edge #2
        if (dom(ib)%edgprev2.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev2)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev2,g)) then
                 do nn=2,ni-1
                    i=nn; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev2))%re2p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev2,g)) then
                 do nn=2,nic-1
                    i=2*nn-2; j=nj-2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev2))%re2p(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,ni-1
                    i=nn; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev2))%re2p(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj+1-nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgprev2))%re2p(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj+1-nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe2m => dom(ib) % se2m
              rbufe2m => dom(ib) % re2m
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev2,g)) then
                 tsend=ni
                 trecv=ni
                 do nn=2,ni-1
                    i=nn; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe2m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev2,g)) then
                 tsend=nic
                 trecv=ni
                 do nn=2,nic-1
                    i=2*nn-2; j=nj-2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe2m(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=nif
                 trecv=ni
                 do nn=2,ni-1
                    i=nn; j=nj-1; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe2m(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj+1-nij))/64.0d0
                    sbufe2m(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj+1-nij))/64.0d0
                 end do
              end if
              end if
              tag=17*10**5+dom(ib)%edgprev2
              call MPI_IRECV  (dom(ib)%re2m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev2),dom(ib)%tg(18),
     &MPI_COMM_WORLD,dom(ib)%rq_e2m,ierr)
              call MPI_SEND (dom(ib)%se2m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev2),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!=====> previous edge #3
        if (dom(ib)%edgprev3.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev3)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev3,g)) then
                 do nn=2,nj-1
                    i=ni-1; j=nn; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev3))%re3p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev3,g)) then
                 do nn=2,njc-1
                    i=ni-2; j=2*nn-2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev3))%re3p(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,nj-1
                    i=ni-1; j=nn; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev3))%re3p(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj-1-nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgprev3))%re3p(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj+1-nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe3m => dom(ib) % se3m
              rbufe3m => dom(ib) % re3m
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev3,g)) then
                 tsend=nj
                 trecv=nj
                 do nn=2,nj-1
                    i=ni-1; j=nn; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe3m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev3,g)) then
                 tsend=njc
                 trecv=nj
                 do nn=2,njc-1
                    i=ni-2; j=2*nn-2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe3m(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=njf
                 trecv=nj
                 do nn=2,nj-1
                    i=ni-1; j=nn; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe3m(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj-1-nij))/64.0d0
                    sbufe3m(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj+1-nij))/64.0d0
                 end do
              end if
              end if
              tag=19*10**5+dom(ib)%edgprev3
              call MPI_IRECV  (dom(ib)%re3m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev3),dom(ib)%tg(20),
     &MPI_COMM_WORLD,dom(ib)%rq_e3m,ierr)
              call MPI_SEND (dom(ib)%se3m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev3),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!=====> previous edge #4
        if (dom(ib)%edgprev4.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev4)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev4,g)) then
                 do nn=2,ni-1
                    i=nn; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev4))%re4p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev4,g)) then
                 do nn=2,nic-1
                    i=2*nn-2; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev4))%re4p(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,ni-1
                    i=nn; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev4))%re4p(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj-1-nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgprev4))%re4p(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj-1-nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe4m => dom(ib) % se4m
              rbufe4m => dom(ib) % re4m
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev4,g)) then
                 tsend=ni
                 trecv=ni
                 do nn=2,ni-1
                    i=nn; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe4m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev4,g)) then
                 tsend=nic
                 trecv=ni
                 do nn=2,nic-1
                    i=2*nn-2; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe4m(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=nif
                 trecv=ni
                 do nn=2,ni-1
                    i=nn; j=2; k=2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe4m(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj-1-nij))/64.0d0
                    sbufe4m(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj-1-nij))/64.0d0
                 end do
              end if
              end if
              tag=21*10**5+dom(ib)%edgprev4
              call MPI_IRECV  (dom(ib)%re4m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev4),dom(ib)%tg(22),
     &MPI_COMM_WORLD,dom(ib)%rq_e4m,ierr)
              call MPI_SEND (dom(ib)%se4m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev4),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!=====> previous edge #5
        if (dom(ib)%edgprev5.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev5)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev5,g)) then
                 do nn=2,nk-1
                    i=2; j=2; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev5))%re5p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev5,g)) then
                 do nn=2,nkc-1
                    i=2; j=2; k=2*nn-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev5))%re5p(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,nk-1
                    i=2; j=2; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev5))%re5p(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj-1-nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgprev5))%re5p(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj-1+nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe5m => dom(ib) % se5m
              rbufe5m => dom(ib) % re5m
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev5,g)) then
                 tsend=nk
                 trecv=nk
                 do nn=2,nk-1
                    i=2; j=2; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe5m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev5,g)) then
                 tsend=nkc
                 trecv=nk
                 do nn=2,nkc-1
                    i=2; j=2; k=2*nn-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe5m(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=nkf
                 trecv=nk
                 do nn=2,nk-1
                    i=2; j=2; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe5m(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj-1-nij))/64.0d0
                    sbufe5m(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj-1+nij))/64.0d0
                 end do
              end if
              end if
              tag=23*10**5+dom(ib)%edgprev5
              call MPI_IRECV  (dom(ib)%re5m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev5),dom(ib)%tg(24),
     &MPI_COMM_WORLD,dom(ib)%rq_e5m,ierr)
              call MPI_SEND (dom(ib)%se5m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev5),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!=====> previous edge #6
        if (dom(ib)%edgprev6.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev6)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev6,g)) then
                 do nn=2,nk-1
                    i=2; j=nj-1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev6))%re6p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev6,g)) then
                 do nn=2,nkc-1
                    i=2; j=nj-2; k=2*nn-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev6))%re6p(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,nk-1
                    i=2; j=nj-1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgprev6))%re6p(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj+1-nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgprev6))%re6p(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj+1+nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe6m => dom(ib) % se6m
              rbufe6m => dom(ib) % re6m
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgprev6,g)) then
                 tsend=nk
                 trecv=nk
                 do nn=2,nk-1
                    i=2; j=nj-1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe6m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgprev6,g)) then
                 tsend=nkc
                 trecv=nk
                 do nn=2,nkc-1
                    i=2; j=nj-2; k=2*nn-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe6m(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=nkf
                 trecv=nk
                 do nn=2,nk-1
                    i=2; j=nj-1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe6m(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll-nj-nij)+ 1.0d0*fi(ll-nj+1-nij))/64.0d0
                    sbufe6m(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj+1+nij))/64.0d0
                 end do
              end if
              end if
              tag=25*10**5+dom(ib)%edgprev6
              call MPI_IRECV  (dom(ib)%re6m(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgprev6),dom(ib)%tg(26),
     &MPI_COMM_WORLD,dom(ib)%rq_e6m,ierr)
              call MPI_SEND (dom(ib)%se6m(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgprev6),tag,MPI_COMM_WORLD,ierr)

           end if
        end if

!======================================================================
!======================================================================
        end if
!..........................................................................
!=== Next Neighbor ===> 
!..........................................................................
        if (cpu_next.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(cpu_next)) then

              if (sync_dir.eq.1)  then

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                 do k=1,ns2
                    do j=1,ns1
           ijk=(k-1)*ns1+j
           i=ni-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=fi(ll)
                    end do
                 end do
                 else
                 if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
                    do k=2,nsc2-1
           do j=2,nsc1-1
              ijk=(k-1)*nsc1+j
              ii=ni-2; jj=2*j-2; kk=2*k-2
              ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
              dom(dom_indid(cpu_next))%recvb_m1(ijk)=0.125*
     & (fi(ll)+fi(ll+1)+fi(ll+nij)+fi(ll+1+nij)+
     & fi(ll+nj)+fi(ll+1+nj)+fi(ll+nj+nij)+fi(ll+1+nj+nij))
           end do
                    end do
                 else
                 do kc=2,ns2-1
                    do jc=2,ns1-1
           ic=ni-1; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk

           k=2*kc-2; j=2*jc-2; ijk=(k-1)*nsr1+j
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
           k=2*kc-2; j=2*jc-1; ijk=(k-1)*nsr1+j
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
           k=2*kc-1; j=2*jc-2; ijk=(k-1)*nsr1+j
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
           k=2*kc-1; j=2*jc-1; ijk=(k-1)*nsr1+j
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
                    end do
                 end do
                 end if
                 end if
              else if (sync_dir.eq.2)  then

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                    do k=1,ns2
           do i=1,ns1
           ijk=(k-1)*ns1+i
           j=nj-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=fi(ll)
           end do
                    end do
                 else
                    if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
           do k=2,nsc2-1
              do i=2,nsc1-1
                 ijk=(k-1)*nsc1+i
                 ii=2*i-2; jj=nj-2; kk=2*k-2
                 ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
                 dom(dom_indid(cpu_next))%recvb_m1(ijk)=0.125*
     & (fi(ll)+fi(ll+nj)+fi(ll+nij)+fi(ll+nj+nij)+
     & fi(ll+1)+fi(ll+1+nj)+fi(ll+1+nij)+fi(ll+1+nj+nij))
              end do
           end do
                    else
                 do kc=2,ns2-1
                    do ic=2,ns1-1
           jc=nj-1; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk

           k=2*kc-2; i=2*ic-2; ijk=(k-1)*nsr1+i
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
           k=2*kc-2; i=2*ic-1; ijk=(k-1)*nsr1+i
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
           k=2*kc-1; i=2*ic-2; ijk=(k-1)*nsr1+i
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
           k=2*kc-1; i=2*ic-1; ijk=(k-1)*nsr1+i
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
                    end do
                 end do
                    end if
                 end if
              else

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                    do i=1,ns2
           do j=1,ns1
           ijk=(i-1)*ns1+j
           k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=fi(ll)
           end do
                    end do
                 else
                    if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
           do i=2,nsc2-1
              do j=2,nsc1-1
                 ijk=(i-1)*nsc1+j
                 ii=2*i-2; jj=2*j-2; kk=nk-2
                 ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
                 dom(dom_indid(cpu_next))%recvb_m1(ijk)=0.125*
     & (fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))
              end do
           end do
                    else
                 do ic=2,ns2-1
                    do jc=2,ns1-1
           kc=nk-1; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk

           i=2*ic-2; j=2*jc-2; ijk=(i-1)*nsr1+j
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
           i=2*ic-2; j=2*jc-1; ijk=(i-1)*nsr1+j
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
           i=2*ic-1; j=2*jc-2; ijk=(i-1)*nsr1+j
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
           i=2*ic-1; j=2*jc-1; ijk=(i-1)*nsr1+j
           dom(dom_indid(cpu_next))%recvb_m1(ijk)=2.0d0*fi(ll)
                    end do
                 end do
                    end if
                 end if
!======================================================================
!======================================================================
              end if

           else !if (dom_ad(dom_id(ib)) .ne. dom_ad(cpu_next)) then

                    sbuf_p => dom(ib) % sendb_p1
                    rbuf_p => dom(ib) % recvb_p1

!======================================================================
!======================================================================
              if (sync_dir.eq.1)  then

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                 tsend=ns1*ns2
                 trecv=ns1*ns2
                 do k=1,ns2
                    do j=1,ns1
           ijk=(k-1)*ns1+j
           i=ni-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
           sbuf_p(ijk)=fi(ll)
                    end do
                 end do
                 else
                 if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
                    tsend=nsc1*nsc2
                    trecv=ns1*ns2
                    do k=2,nsc2-1
           do j=2,nsc1-1
              ijk=(k-1)*nsc1+j
              ii=ni-2; jj=2*j-2; kk=2*k-2
              ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
              sbuf_p(ijk)=0.125*
     & (fi(ll)+fi(ll+1)+fi(ll+nij)+fi(ll+1+nij)+
     & fi(ll+nj)+fi(ll+1+nj)+fi(ll+nj+nij)+fi(ll+1+nj+nij))
           end do
                    end do
                 else
                    tsend=nsr1*nsr2
                    trecv=ns1*ns2
                 do kc=2,ns2-1
                    do jc=2,ns1-1
           ic=ni-1; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk

           k=2*kc-2; j=2*jc-2; ijk=(k-1)*nsr1+j
           sbuf_p(ijk)=2.0d0*fi(ll)
           k=2*kc-2; j=2*jc-1; ijk=(k-1)*nsr1+j
           sbuf_p(ijk)=2.0d0*fi(ll)
           k=2*kc-1; j=2*jc-2; ijk=(k-1)*nsr1+j
           sbuf_p(ijk)=2.0d0*fi(ll)
           k=2*kc-1; j=2*jc-1; ijk=(k-1)*nsr1+j
           sbuf_p(ijk)=2.0d0*fi(ll)
                    end do
                 end do
                 end if
                 end if
!======================================================================
!======================================================================
              else if (sync_dir.eq.2)  then

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                    tsend=ns1*ns2
                    trecv=ns1*ns2
                    do k=1,ns2
           do i=1,ns1
           ijk=(k-1)*ns1+i
           j=nj-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
           sbuf_p(ijk)=fi(ll)
           end do
                    end do
                 else
                    if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
           tsend=nsc1*nsc2
           trecv=ns1*ns2
           do k=2,nsc2-1
              do i=2,nsc1-1
                 ijk=(k-1)*nsc1+i
                 ii=2*i-2; jj=nj-2; kk=2*k-2
                 ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
                 sbuf_p(ijk)=0.125*
     & (fi(ll)+fi(ll+nj)+fi(ll+nij)+fi(ll+nj+nij)+
     & fi(ll+1)+fi(ll+1+nj)+fi(ll+1+nij)+fi(ll+1+nj+nij))
              end do
           end do
                    else
           tsend=nsr1*nsr2
           trecv=ns1*ns2
                 do kc=2,ns2-1
                    do ic=2,ns1-1
           jc=nj-1; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk

           k=2*kc-2; i=2*ic-2; ijk=(k-1)*nsr1+i
           sbuf_p(ijk)=2.0d0*fi(ll)
           k=2*kc-2; i=2*ic-1; ijk=(k-1)*nsr1+i
           sbuf_p(ijk)=2.0d0*fi(ll)
           k=2*kc-1; i=2*ic-2; ijk=(k-1)*nsr1+i
           sbuf_p(ijk)=2.0d0*fi(ll)
           k=2*kc-1; i=2*ic-1; ijk=(k-1)*nsr1+i
           sbuf_p(ijk)=2.0d0*fi(ll)
                    end do
                 end do
                    end if
                 end if
!======================================================================
!======================================================================
              else

                 if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                    tsend=ns1*ns2
                    trecv=ns1*ns2
                    do i=1,ns2
           do j=1,ns1
           ijk=(i-1)*ns1+j
           k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
           sbuf_p(ijk)=fi(ll)
           end do
                    end do
                 else
                    if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
           tsend=nsc1*nsc2
           trecv=ns1*ns2
           do i=2,nsc2-1
              do j=2,nsc1-1
                 ijk=(i-1)*nsc1+j
                 ii=2*i-2; jj=2*j-2; kk=nk-2
                 ll=(kk-1)*nij+(ii-1)*nj+jj+nijk
                 sbuf_p(ijk)=0.125*
     & (fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))
              end do
           end do
                    else
           tsend=nsr1*nsr2
           trecv=ns1*ns2
                 do ic=2,ns2-1
                    do jc=2,ns1-1
           kc=nk-1; ll=(kc-1)*nij+(ic-1)*nj+jc+nijk

           i=2*ic-2; j=2*jc-2; ijk=(i-1)*nsr1+j
           sbuf_p(ijk)=2.0d0*fi(ll)
           i=2*ic-2; j=2*jc-1; ijk=(i-1)*nsr1+j
           sbuf_p(ijk)=2.0d0*fi(ll)
           i=2*ic-1; j=2*jc-2; ijk=(i-1)*nsr1+j
           sbuf_p(ijk)=2.0d0*fi(ll)
           i=2*ic-1; j=2*jc-1; ijk=(i-1)*nsr1+j
           sbuf_p(ijk)=2.0d0*fi(ll)
                    end do
                 end do
                    end if
                 end if
!======================================================================
!======================================================================
              end if
              tag=(2*sync_dir)*10**5+cpu_next
              ta=2*sync_dir-1
              call MPI_IRECV  (dom(ib)%recvb_p1(1),trecv,MPI_FLT,
     &dom_ad(cpu_next),dom(ib)%tg(ta),MPI_COMM_WORLD,dom(ib)%rq_p1,ierr)
              call MPI_SEND (dom(ib)%sendb_p1(1),tsend,MPI_FLT,
     &dom_ad(cpu_next),tag,MPI_COMM_WORLD,ierr)

                 end if
              end if

!======================================================================
!======================================================================
        if (sync_dir.eq.3)  then

!..........................................................................
!=== Next Corner Neighbors  ===> 
!..........................................................................

!=====> next cor #1
        if (dom(ib)%cornext1.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%cornext1)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%cornext1,g)) then
              i=ni-1; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              dom(dom_indid(dom(ib)%cornext1))%rc1m(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%cornext1,g)) then
                 i=ni-2; j=nj-2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%cornext1))%rc1m(1)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=ni-1; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%cornext1))%rc1m(1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj+1+nij))/64.0d0
              end if
              end if
           else
              sbufc1p => dom(ib) % sc1p
              rbufc1p => dom(ib) % rc1p
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%cornext1,g)) then
              i=ni-1; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              sbufc1p(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%cornext1,g)) then
                 i=ni-2; j=nj-2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc1p(1)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=ni-1; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc1p(1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj+1+nij))/64.0d0
              end if
              end if
              tag=8*10**5+dom(ib)%cornext1
              call MPI_IRECV  (dom(ib)%rc1p,1,MPI_FLT,
     &dom_ad(dom(ib)%cornext1),dom(ib)%tg(7),
     &MPI_COMM_WORLD,dom(ib)%rq_c1p,ierr)
              call MPI_SEND (dom(ib)%sc1p,1,MPI_FLT,
     &dom_ad(dom(ib)%cornext1),tag,MPI_COMM_WORLD,ierr)

           end if
        end if

!=====> next cor #2
        if (dom(ib)%cornext2.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%cornext2)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%cornext2,g)) then
              i=ni-1; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              dom(dom_indid(dom(ib)%cornext2))%rc2m(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%cornext2,g)) then
                 i=ni-2; j=2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%cornext2))%rc2m(1)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=ni-1; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%cornext2))%rc2m(1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj-1+nij))/64.0d0
              end if
              end if
           else
              sbufc2p => dom(ib) % sc2p
              rbufc2p => dom(ib) % rc2p
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%cornext2,g)) then
              i=ni-1; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              sbufc2p(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%cornext2,g)) then
                 i=ni-2; j=2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc2p(1)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=ni-1; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc2p(1)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj-1+nij))/64.0d0
              end if
              end if
              tag=10*10**5+dom(ib)%cornext2
              call MPI_IRECV  (dom(ib)%rc2p,1,MPI_FLT,
     &dom_ad(dom(ib)%cornext2),dom(ib)%tg(9),
     &MPI_COMM_WORLD,dom(ib)%rq_c2p,ierr)
              call MPI_SEND (dom(ib)%sc2p,1,MPI_FLT,
     &dom_ad(dom(ib)%cornext2),tag,MPI_COMM_WORLD,ierr)

           end if
        end if

!=====> next cor #3
        if (dom(ib)%cornext3.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%cornext3)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%cornext3,g)) then
              i=2; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              dom(dom_indid(dom(ib)%cornext3))%rc3m(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%cornext3,g)) then
                 i=2; j=2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%cornext3))%rc3m(1)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=2; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%cornext3))%rc3m(1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj-1+nij))/64.0d0
              end if
              end if
           else
              sbufc3p => dom(ib) % sc3p
              rbufc3p => dom(ib) % rc3p
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%cornext3,g)) then
              i=2; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              sbufc3p(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%cornext3,g)) then
                 i=2; j=2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc3p(1)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=2; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc3p(1)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj-1+nij))/64.0d0
              end if
              end if
              tag=12*10**5+dom(ib)%cornext3
              call MPI_IRECV  (dom(ib)%rc3p,1,MPI_FLT,
     &dom_ad(dom(ib)%cornext3),dom(ib)%tg(11),
     &MPI_COMM_WORLD,dom(ib)%rq_c3p,ierr)
              call MPI_SEND (dom(ib)%sc3p,1,MPI_FLT,
     &dom_ad(dom(ib)%cornext3),tag,MPI_COMM_WORLD,ierr)

           end if
        end if

!=====> next cor #4
        if (dom(ib)%cornext4.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%cornext4)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%cornext4,g)) then
              i=2; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              dom(dom_indid(dom(ib)%cornext4))%rc4m(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%cornext4,g)) then
                 i=2; j=nj-2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%cornext4))%rc4m(1)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=2; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 dom(dom_indid(dom(ib)%cornext4))%rc4m(1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj+1+nij))/64.0d0
              end if
              end if
           else
              sbufc4p => dom(ib) % sc4p
              rbufc4p => dom(ib) % rc4p
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%cornext4,g)) then
              i=2; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              sbufc4p(1)=fi(ll)
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%cornext4,g)) then
                 i=2; j=nj-2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc4p(1)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
              else
                 i=2; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 sbufc4p(1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj+1+nij))/64.0d0
              end if
              end if
              tag=14*10**5+dom(ib)%cornext4
              call MPI_IRECV  (dom(ib)%rc4p,1,MPI_FLT,
     &dom_ad(dom(ib)%cornext4),dom(ib)%tg(13),
     &MPI_COMM_WORLD,dom(ib)%rq_c4p,ierr)
              call MPI_SEND (dom(ib)%sc4p,1,MPI_FLT,
     &dom_ad(dom(ib)%cornext4),tag,MPI_COMM_WORLD,ierr)

           end if
        end if

!..........................................................................
!=== Next Edge Neighbors  ===> 
!..........................................................................

!=====> next edge #1
        if (dom(ib)%edgnext1.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext1)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext1,g)) then
                 do nn=2,nj-1
                    i=ni-1; j=nn; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext1))%re1m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext1,g)) then
                 do nn=2,njc-1
                 i=ni-2; j=2*nn-2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext1))%re1m(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,nj-1
                    i=ni-1; j=nn; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext1))%re1m(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj-1+nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgnext1))%re1m(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj+1+nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe1p => dom(ib) % se1p
              rbufe1p => dom(ib) % re1p
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext1,g)) then
                 tsend=nj
                 trecv=nj
                 do nn=2,nj-1
                    i=ni-1; j=nn; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe1p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext1,g)) then
                 tsend=njc
                 trecv=nj
                 do nn=2,njc-1
                 i=ni-2; j=2*nn-2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe1p(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=njf
                 trecv=nj
                 do nn=2,nj-1
                    i=ni-1; j=nn; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe1p(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj-1+nij))/64.0d0
                    sbufe1p(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj+1+nij))/64.0d0
                 end do
              end if
              end if
              tag=16*10**5+dom(ib)%edgnext1
              call MPI_IRECV  (dom(ib)%re1p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext1),dom(ib)%tg(15),
     &MPI_COMM_WORLD,dom(ib)%rq_e1p,ierr)
              call MPI_SEND (dom(ib)%se1p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext1),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!=====> next edge #2
        if (dom(ib)%edgnext2.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext2)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext2,g)) then
                 do nn=2,ni-1
                    i=nn; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext2))%re2m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext2,g)) then
                 do nn=2,nic-1
                    i=2*nn-2; j=2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext2))%re2m(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,ni-1
                    i=nn; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext2))%re2m(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj-1+nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgnext2))%re2m(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj-1+nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe2p => dom(ib) % se2p
              rbufe2p => dom(ib) % re2p
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext2,g)) then
                 tsend=ni
                 trecv=ni
                 do nn=2,ni-1
                    i=nn; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe2p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext2,g)) then
                 tsend=nic
                 trecv=ni
                 do nn=2,nic-1
                    i=2*nn-2; j=2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe2p(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=nif
                 trecv=ni
                 do nn=2,ni-1
                    i=nn; j=2; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe2p(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj-1+nij))/64.0d0
                    sbufe2p(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj-1+nij))/64.0d0
                 end do
              end if
              end if
              tag=18*10**5+dom(ib)%edgnext2
              call MPI_IRECV  (dom(ib)%re2p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext2),dom(ib)%tg(17),
     &MPI_COMM_WORLD,dom(ib)%rq_e2p,ierr)
              call MPI_SEND (dom(ib)%se2p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext2),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!=====> next edge #3
        if (dom(ib)%edgnext3.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext3)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext3,g)) then
                 do nn=2,nj-1
                    i=2; j=nn; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext3))%re3m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext3,g)) then
                 do nn=2,njc-1
                    i=2; j=2*nn-2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext3))%re3m(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,nj-1
                    i=2; j=nn; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext3))%re3m(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj-1+nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgnext3))%re3m(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj+1+nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe3p => dom(ib) % se3p
              rbufe3p => dom(ib) % re3p
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext3,g)) then
                 tsend=nj
                 trecv=nj
                 do nn=2,nj-1
                    i=2; j=nn; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe3p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext3,g)) then
                 tsend=njc
                 trecv=nj
                 do nn=2,njc-1
                    i=2; j=2*nn-2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe3p(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=njf
                 trecv=nj
                 do nn=2,nj-1
                    i=2; j=nn; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe3p(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj-1+nij))/64.0d0
                    sbufe3p(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj+1+nij))/64.0d0
                 end do
              end if
              end if
              tag=20*10**5+dom(ib)%edgnext3
              call MPI_IRECV  (dom(ib)%re3p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext3),dom(ib)%tg(19),
     &MPI_COMM_WORLD,dom(ib)%rq_e3p,ierr)
              call MPI_SEND (dom(ib)%se3p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext3),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!=====> next edge #4
        if (dom(ib)%edgnext4.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext4)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext4,g)) then
                 do nn=2,ni-1
                    i=nn; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext4))%re4m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext4,g)) then
                 do nn=2,nic-1
                 i=2*nn-2; j=nj-2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext4))%re4m(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,ni-1
                    i=nn; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext4))%re4m(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj+1+nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgnext4))%re4m(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj+1+nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe4p => dom(ib) % se4p
              rbufe4p => dom(ib) % re4p
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext4,g)) then
                 tsend=ni
                 trecv=ni
                 do nn=2,ni-1
                    i=nn; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe4p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext4,g)) then
                 tsend=nic
                 trecv=ni
                 do nn=2,nic-1
                 i=2*nn-2; j=nj-2; k=nk-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe4p(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=nif
                 trecv=ni
                 do nn=2,ni-1
                    i=nn; j=nj-1; k=nk-1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe4p(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll-nj)+ 3.0d0*fi(ll-nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll-nj+nij)+ 1.0d0*fi(ll-nj+1+nij))/64.0d0
                    sbufe4p(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj+1+nij))/64.0d0
                 end do
              end if
              end if
              tag=22*10**5+dom(ib)%edgnext4
              call MPI_IRECV  (dom(ib)%re4p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext4),dom(ib)%tg(21),
     &MPI_COMM_WORLD,dom(ib)%rq_e4p,ierr)
              call MPI_SEND (dom(ib)%se4p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext4),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!=====> next edge #5
        if (dom(ib)%edgnext5.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext5)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext5,g)) then
                 do nn=2,nk-1
                    i=ni-1; j=nj-1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext5))%re5m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext5,g)) then
                 do nn=2,nkc-1
                 i=ni-2; j=nj-2; k=2*nn-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext5))%re5m(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,nk-1
                    i=ni-1; j=nj-1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext5))%re5m(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj+1-nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgnext5))%re5m(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj+1+nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe5p => dom(ib) % se5p
              rbufe5p => dom(ib) % re5p
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext5,g)) then
                 tsend=nk
                 trecv=nk
                 do nn=2,nk-1
                    i=ni-1; j=nj-1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe5p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext5,g)) then
                 tsend=nkc
                 trecv=nk
                 do nn=2,nkc-1
                 i=ni-2; j=nj-2; k=2*nn-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe5p(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=nkf
                 trecv=nk
                 do nn=2,nk-1
                    i=ni-1; j=nj-1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe5p(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll+1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj+1-nij))/64.0d0
                    sbufe5p(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll+1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj+1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll+1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj+1+nij))/64.0d0
                 end do
              end if
              end if
              tag=24*10**5+dom(ib)%edgnext5
              call MPI_IRECV  (dom(ib)%re5p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext5),dom(ib)%tg(23),
     &MPI_COMM_WORLD,dom(ib)%rq_e5p,ierr)
              call MPI_SEND (dom(ib)%se5p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext5),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!=====> next edge #6
        if (dom(ib)%edgnext6.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext6)) then
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext6,g)) then
                 do nn=2,nk-1
                    i=ni-1; j=2; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext6))%re6m(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext6,g)) then
                 do nn=2,nkc-1
                    i=ni-2; j=2; k=2*nn-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext6))%re6m(nn)=
     &(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 do nn=2,nk-1
                    i=ni-1; j=2; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    dom(dom_indid(dom(ib)%edgnext6))%re6m(2*nn-2)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj-1-nij))/64.0d0
                    dom(dom_indid(dom(ib)%edgnext6))%re6m(2*nn-1)=
     &(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj-1+nij))/64.0d0
                 end do
              end if
              end if
           else
              sbufe6p => dom(ib) % se6p
              rbufe6p => dom(ib) % re6p
              if(rdv(dom_id(ib),g).eq.rdv(dom(ib)%edgnext6,g)) then
                 tsend=nk
                 trecv=nk
                 do nn=2,nk-1
                    i=ni-1; j=2; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe6p(nn)=fi(ll)
                 end do
              else
              if(rdv(dom_id(ib),g).gt.rdv(dom(ib)%edgnext6,g)) then
                 tsend=nkc
                 trecv=nk
                 do nn=2,nkc-1
                    i=ni-2; j=2; k=2*nn-2; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe6p(nn)=(fi(ll)+fi(ll+1)+fi(ll+nj)+fi(ll+1+nj)+
     & fi(ll+nij)+fi(ll+1+nij)+fi(ll+nj+nij)+fi(ll+1+nj+nij))/8.0
                 end do
              else
                 tsend=nkf
                 trecv=nk
                 do nn=2,nk-1
                    i=ni-1; j=2; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                    sbufe6p(2*nn-2)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll-nij)+ 3.0d0*fi(ll-1-nij)+
     & 3.0d0*fi(ll+nj-nij)+ 1.0d0*fi(ll+nj-1-nij))/64.0d0
                    sbufe6p(2*nn-1)=(27.0d0*fi(ll)+9.0d0*fi(ll-1)+
     & 9.0d0*fi(ll+nj)+ 3.0d0*fi(ll+nj-1)+
     & 9.0d0*fi(ll+nij)+ 3.0d0*fi(ll-1+nij)+
     & 3.0d0*fi(ll+nj+nij)+ 1.0d0*fi(ll+nj-1+nij))/64.0d0
                 end do
              end if
              end if
              tag=26*10**5+dom(ib)%edgnext6
              call MPI_IRECV  (dom(ib)%re6p(1),trecv,MPI_FLT,
     &dom_ad(dom(ib)%edgnext6),dom(ib)%tg(25),
     &MPI_COMM_WORLD,dom(ib)%rq_e6p,ierr)
              call MPI_SEND (dom(ib)%se6p(1),tsend,MPI_FLT,
     &dom_ad(dom(ib)%edgnext6),tag,MPI_COMM_WORLD,ierr)

           end if
        end if
!======================================================================
!======================================================================
        end if
           end do
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
           do ib=1,nbp
              if(g.gt.dom(ib)%ngrid) then
                 gl=dom(ib)%ngrid
              else
                 gl=g
              end if
        ni=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        nj=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nk=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

        fi => dom(ib)%cof
        nijk=dom(ib)%faz(g)-ni*nj*nk
        nij=ni*nj

              if (sync_dir.eq.1)  then
                 cpu_prev=dom(ib)%iprev
                 cpu_next=dom(ib)%inext
                 ns1=nj 
                 ns2=nk  
              else if (sync_dir.eq.2)  then
                 cpu_prev=dom(ib)%jprev
                 cpu_next=dom(ib)%jnext
                 ns1=ni  
                 ns2=nk  
              else if (sync_dir.eq.3)  then
                 cpu_prev=dom(ib)%kprev
                 cpu_next=dom(ib)%knext
                 ns1=nj  
                 ns2=ni 
              end if

!..............................................................................
!=== Previous Neighbor  ===> 
!..............................................................................
              if (cpu_prev.ge.0) then

                 if (dom_ad(dom_id(ib)) .eq. dom_ad(cpu_prev)) then

                    if (sync_dir.eq.1)  then

                       if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                          do k=1,ns2
                             do j=1,ns1
                                ijk=(k-1)*ns1+j
                                i=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=dom(ib)%recvb_m1(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
                             do k=1,ns2
                                do j=1,ns1
                                   ijk=(k-1)*ns1+j
                                   i=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (dom(ib)%recvb_m1(ijk)+3.0*fi(ll+nj)-fi(ll+2*nj))/4.0
                                end do
                             end do
                          else
                             do k=2,ns2-1
                                do j=2,ns1-1
                                   ijk=(k-1)*ns1+j
                                   i=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=dom(ib)%recvb_m1(ijk)
                                end do
                             end do
                          end if
                       end if

                    else if (sync_dir.eq.2)  then
                       if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                          do k=1,ns2
                             do i=1,ns1
                                ijk=(k-1)*ns1+i
                                j=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=dom(ib)%recvb_m1(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
                             do k=1,ns2
                                do i=1,ns1
                                   ijk=(k-1)*ns1+i
                                   j=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (dom(ib)%recvb_m1(ijk)+3.0*fi(ll+1)-fi(ll+2))/4.0
                                end do
                             end do
                          else
                             do k=2,ns2-1
                                do i=2,ns1-1
                                   ijk=(k-1)*ns1+i
                                   j=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=dom(ib)%recvb_m1(ijk)
                                end do
                             end do
                          end if
                       end if

                    else
                       if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                          do i=1,ns2
                             do j=1,ns1
                                ijk=(i-1)*ns1+j
                                k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=dom(ib)%recvb_m1(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
                             do i=1,ns2
                                do j=1,ns1
                                   ijk=(i-1)*ns1+j
                                   k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (dom(ib)%recvb_m1(ijk)+3.0*fi(ll+nij)-fi(ll+2*nij))/4.0
                                end do
                             end do

                          else
                             do i=2,ns2-1
                                do j=2,ns1-1
                                   ijk=(i-1)*ns1+j
                                   k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=dom(ib)%recvb_m1(ijk)
                                end do
                             end do
                          end if
                       end if

                    end if


                 else !if (dom_ad(dom_id(ib)) .ne. dom_ad(cpu_prev)) then

                 call MPI_WAIT(dom(ib)%rq_m1,MPI_STATUS_IGNORE,ierr)
                 rbuf_m => dom(ib) % recvb_m1
!======================================================================
!======================================================================
                    if (sync_dir.eq.1)  then

                       if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                          do k=1,ns2
                             do j=1,ns1
                                ijk=(k-1)*ns1+j
                                i=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=rbuf_m(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
                             do k=1,ns2
                                do j=1,ns1
                                   ijk=(k-1)*ns1+j
                                   i=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (rbuf_m(ijk)+3.0*fi(ll+nj)-fi(ll+2*nj))/4.0
                                end do
                             end do
                          else
                             do k=2,ns2-1
                                do j=2,ns1-1
                                   ijk=(k-1)*ns1+j
                                   i=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=rbuf_m(ijk)
                                end do
                             end do
                          end if
                       end if

!======================================================================
!======================================================================
                    else if (sync_dir.eq.2)  then

                       if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                          do k=1,ns2
                             do i=1,ns1
                                ijk=(k-1)*ns1+i
                                j=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=rbuf_m(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
                             do k=1,ns2
                                do i=1,ns1
                                   ijk=(k-1)*ns1+i
                                   j=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (rbuf_m(ijk)+3.0*fi(ll+1)-fi(ll+2))/4.0
                                end do
                             end do
                          else
                             do k=2,ns2-1
                                do i=2,ns1-1
                                   ijk=(k-1)*ns1+i
                                   j=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=rbuf_m(ijk)
                                end do
                             end do
                          end if
                       end if
!======================================================================
!======================================================================
                    else

                       if(rdv(dom_id(ib),g).eq.rdv(cpu_prev,g)) then
                          do i=1,ns2
                             do j=1,ns1
                                ijk=(i-1)*ns1+j
                                k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=rbuf_m(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_prev,g)) then
                             do i=1,ns2
                                do j=1,ns1
                                   ijk=(i-1)*ns1+j
                                   k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (rbuf_m(ijk)+3.0*fi(ll+nij)-fi(ll+2*nij))/4.0
                                end do
                             end do
                          else
                             do i=2,ns2-1
                                do j=2,ns1-1
                                   ijk=(i-1)*ns1+j
                                   k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=rbuf_m(ijk)
                                end do
                             end do
                          end if
                       end if

                    end if

                 end if
              end if

!..............................................................................
!=== Next Neighbor  ===> 
!..............................................................................
              if (cpu_next.ge.0) then

                 if (dom_ad(dom_id(ib)) .eq. dom_ad(cpu_next)) then

                    if (sync_dir.eq.1)  then
                       if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                          do k=1,ns2
                             do j=1,ns1
                                ijk=(k-1)*ns1+j
                                i=ni; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=dom(ib)%recvb_p1(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
                             do k=1,ns2
                                do j=1,ns1
                                   ijk=(k-1)*ns1+j
                                   i=ni; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (dom(ib)%recvb_p1(ijk)+3.0*fi(ll-nj)-fi(ll-2*nj))/4.0
                                end do
                             end do
                          else
                             do k=2,ns2-1
                                do j=2,ns1-1
                                   ijk=(k-1)*ns1+j
                                   i=ni; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=dom(ib)%recvb_p1(ijk)
                                end do
                             end do
                          end if
                       end if

                    else if (sync_dir.eq.2)  then
                       if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                          do k=1,ns2
                             do i=1,ns1
                                ijk=(k-1)*ns1+i
                                j=nj; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=dom(ib)%recvb_p1(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
                             do k=1,ns2
                                do i=1,ns1
                                   ijk=(k-1)*ns1+i
                                   j=nj; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (dom(ib)%recvb_p1(ijk)+3.0*fi(ll-1)-fi(ll-2))/4.0
                                end do
                             end do
                          else
                             do k=2,ns2-1
                                do i=2,ns1-1
                                   ijk=(k-1)*ns1+i
                                   j=nj; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=dom(ib)%recvb_p1(ijk)
                                end do
                             end do
                          end if
                       end if

                    else
                       if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                          do i=1,ns2
                             do j=1,ns1
                                ijk=(i-1)*ns1+j
                                k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=dom(ib)%recvb_p1(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
                             do i=1,ns2
                                do j=1,ns1
                                   ijk=(i-1)*ns1+j
                                   k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (dom(ib)%recvb_p1(ijk)+3.0*fi(ll-nij)-fi(ll-2*nij))/4.0
                                end do
                             end do
                          else
                             do i=2,ns2-1
                                do j=2,ns1-1
                                   ijk=(i-1)*ns1+j
                                   k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=dom(ib)%recvb_p1(ijk)
                                end do
                             end do
                          end if
                       end if
                    end if

                 else !if (dom_ad(dom_id(ib)) .ne. dom_ad(cpu_next)) then

                    call MPI_WAIT(dom(ib)%rq_p1,MPI_STATUS_IGNORE,ierr)
                    rbuf_p => dom(ib) % recvb_p1

!======================================================================
!======================================================================
                    if (sync_dir.eq.1)  then
                       if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                          do k=1,ns2
                             do j=1,ns1
                                ijk=(k-1)*ns1+j
                                i=ni; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=rbuf_p(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
                             do k=1,ns2
                                do j=1,ns1
                                   ijk=(k-1)*ns1+j
                                   i=ni; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (rbuf_p(ijk)+3.0*fi(ll-nj)-fi(ll-2*nj))/4.0
                                end do
                             end do
                          else
                             do k=2,ns2-1
                                do j=2,ns1-1
                                   ijk=(k-1)*ns1+j
                                   i=ni; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=rbuf_p(ijk)
                                end do
                             end do
                         end if
                       end if

!======================================================================
!======================================================================
                    else if (sync_dir.eq.2)  then
                       if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                          do k=1,ns2
                             do i=1,ns1
                                ijk=(k-1)*ns1+i
                                j=nj; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=rbuf_p(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
                             do k=1,ns2
                                do i=1,ns1
                                   ijk=(k-1)*ns1+i
                                   j=nj; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (rbuf_p(ijk)+3.0*fi(ll-1)-fi(ll-2))/4.0
                                end do
                             end do
                          else
                             do k=2,ns2-1
                                do i=2,ns1-1
                                   ijk=(k-1)*ns1+i
                                   j=nj; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=rbuf_p(ijk)
                                end do
                             end do
                          end if
                       end if
!======================================================================
!======================================================================
                    else
                       if(rdv(dom_id(ib),g).eq.rdv(cpu_next,g)) then
                          do i=1,ns2
                             do j=1,ns1
                                ijk=(i-1)*ns1+j
                                k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                fi(ll)=rbuf_p(ijk)
                             end do
                          end do
                       else
                          if(rdv(dom_id(ib),g).gt.rdv(cpu_next,g)) then
                             do i=1,ns2
                                do j=1,ns1
                                   ijk=(i-1)*ns1+j
                                   k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=
     & (rbuf_p(ijk)+3.0*fi(ll-nij)-fi(ll-2*nij))/4.0
                                end do
                             end do
                          else
                             do i=2,ns2-1
                                do j=2,ns1-1
                                   ijk=(i-1)*ns1+j
                                   k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                                   fi(ll)=rbuf_p(ijk)
                                end do
                             end do
                          end if
                       end if

                    end if

                 end if
              end if

!======================================================================
!======================================================================
        if (sync_dir.eq.3)  then
!..............................................................................
!=== Previous Corner Neighbors  ===> 
!..............................................................................

!=====> previous cor #1
        if (dom(ib)%corprev1.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%corprev1)) then
              i=1; j=1; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=dom(ib) % rc1m(1)
           else
              call MPI_WAIT(dom(ib)%rq_c1m,MPI_STATUS_IGNORE,ierr)
              rbufc1m => dom(ib) % rc1m
              i=1; j=1; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=rbufc1m(1)
           end if
        end if
!=====> previous cor #2
        if (dom(ib)%corprev2.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%corprev2)) then
              i=1; j=nj; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=dom(ib) % rc2m(1)
           else
              call MPI_WAIT(dom(ib)%rq_c2m,MPI_STATUS_IGNORE,ierr)
              rbufc2m => dom(ib) % rc2m
              i=1; j=nj; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=rbufc2m(1)
           end if
        end if
!=====> previous cor #3
        if (dom(ib)%corprev3.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%corprev3)) then
              i=ni; j=nj; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=dom(ib) % rc3m(1)
           else
              call MPI_WAIT(dom(ib)%rq_c3m,MPI_STATUS_IGNORE,ierr)
              rbufc3m => dom(ib) % rc3m
              i=ni; j=nj; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=rbufc3m(1)
           end if
        end if
!=====> previous cor #4
        if (dom(ib)%corprev4.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%corprev4)) then
              i=ni; j=1; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=dom(ib) % rc4m(1)
           else
              call MPI_WAIT(dom(ib)%rq_c4m,MPI_STATUS_IGNORE,ierr)
              rbufc4m => dom(ib) % rc4m
              i=ni; j=1; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=rbufc4m(1)
           end if
        end if

!..............................................................................
!=== Previous Edge Neighbors  ===> 
!..............................................................................

!=====> previous edge #1
        if (dom(ib)%edgprev1.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev1)) then
              do nn=2,nj-1
                 i=1; j=nn; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re1m(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e1m,MPI_STATUS_IGNORE,ierr)
              rbufe1m => dom(ib) % re1m
              do nn=2,nj-1
                 i=1; j=nn; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe1m(nn)
              end do
           end if
        end if
!=====> previous edge #2
        if (dom(ib)%edgprev2.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev2)) then
              do nn=2,ni-1
                 i=nn; j=nj; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re2m(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e2m,MPI_STATUS_IGNORE,ierr)
              rbufe2m => dom(ib) % re2m
              do nn=2,ni-1
                 i=nn; j=nj; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe2m(nn)
              end do
           end if
        end if
!=====> previous edge #3
        if (dom(ib)%edgprev3.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev3)) then
              do nn=2,nj-1
                 i=ni; j=nn; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re3m(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e3m,MPI_STATUS_IGNORE,ierr)
              rbufe3m => dom(ib) % re3m
              do nn=2,nj-1
                 i=ni; j=nn; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe3m(nn)
              end do
           end if
        end if
!=====> previous edge #4
        if (dom(ib)%edgprev4.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev4)) then
              do nn=2,ni-1
                 i=nn; j=1; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re4m(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e4m,MPI_STATUS_IGNORE,ierr)
              rbufe4m => dom(ib) % re4m
              do nn=2,ni-1
                 i=nn; j=1; k=1; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe4m(nn)
              end do
           end if
        end if
!=====> previous edge #5
        if (dom(ib)%edgprev5.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev5)) then
              do nn=2,nk-1
                 i=1; j=1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re5m(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e5m,MPI_STATUS_IGNORE,ierr)
              rbufe5m => dom(ib) % re5m
              do nn=2,nk-1
                 i=1; j=1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe5m(nn)
              end do
           end if
        end if
!=====> previous edge #6
        if (dom(ib)%edgprev6.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgprev6)) then
              do nn=2,nk-1
                 i=1; j=nj; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re6m(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e6m,MPI_STATUS_IGNORE,ierr)
              rbufe6m => dom(ib) % re6m
              do nn=2,nk-1
                 i=1; j=nj; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe6m(nn)
              end do
           end if
        end if
!..............................................................................
!=== Next Corner Neighbors  ===> 
!..............................................................................

!=====> next cor #1
        if (dom(ib)%cornext1.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%cornext1)) then
              i=ni; j=nj; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=dom(ib) % rc1p(1)
           else
              call MPI_WAIT(dom(ib)%rq_c1p,MPI_STATUS_IGNORE,ierr)
              rbufc1p => dom(ib) % rc1p
              i=ni; j=nj; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=rbufc1p(1)
           end if
        end if
!=====> next cor #2
        if (dom(ib)%cornext2.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%cornext2)) then
              i=ni; j=1; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=dom(ib) % rc2p(1)
           else
              call MPI_WAIT(dom(ib)%rq_c2p,MPI_STATUS_IGNORE,ierr)
              rbufc2p => dom(ib) % rc2p
              i=ni; j=1; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=rbufc2p(1)
           end if
        end if
!=====> next cor #3
        if (dom(ib)%cornext3.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%cornext3)) then
              i=1; j=1; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=dom(ib) % rc3p(1)
           else
              call MPI_WAIT(dom(ib)%rq_c3p,MPI_STATUS_IGNORE,ierr)
              rbufc3p => dom(ib) % rc3p
              i=1; j=1; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=rbufc3p(1)
           end if
        end if
!=====> next cor #4
        if (dom(ib)%cornext4.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%cornext4)) then
              i=1; j=nj; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=dom(ib) % rc4p(1)
           else
              call MPI_WAIT(dom(ib)%rq_c4p,MPI_STATUS_IGNORE,ierr)
              rbufc4p => dom(ib) % rc4p
              i=1; j=nj; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
              fi(ll)=rbufc4p(1)
           end if
        end if

!..............................................................................
!=== Next Edge Neighbors  ===> 
!..............................................................................

!=====> next edge #1
        if (dom(ib)%edgnext1.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext1)) then
              do nn=2,nj-1
                 i=ni; j=nn; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re1p(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e1p,MPI_STATUS_IGNORE,ierr)
              rbufe1p => dom(ib) % re1p
              do nn=2,nj-1
                 i=ni; j=nn; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe1p(nn)
              end do
           end if
        end if
!=====> next edge #2
        if (dom(ib)%edgnext2.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext2)) then
              do nn=2,ni-1
                 i=nn; j=1; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re2p(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e2p,MPI_STATUS_IGNORE,ierr)
              rbufe2p => dom(ib) % re2p
              do nn=2,ni-1
                 i=nn; j=1; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe2p(nn)
              end do
           end if
        end if
!=====> next edge #3
        if (dom(ib)%edgnext3.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext3)) then
              do nn=2,nj-1
                 i=1; j=nn; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re3p(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e3p,MPI_STATUS_IGNORE,ierr)
              rbufe3p => dom(ib) % re3p
              do nn=2,nj-1
                 i=1; j=nn; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe3p(nn)
              end do
           end if
        end if
!=====> next edge #4
        if (dom(ib)%edgnext4.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext4)) then
              do nn=2,ni-1
                 i=nn; j=nj; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re4p(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e4p,MPI_STATUS_IGNORE,ierr)
              rbufe4p => dom(ib) % re4p
              do nn=2,ni-1
                 i=nn; j=nj; k=nk; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe4p(nn)
              end do
           end if
        end if
!=====> next edge #5
        if (dom(ib)%edgnext5.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext5)) then
              do nn=2,nk-1
                 i=ni; j=nj; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re5p(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e5p,MPI_STATUS_IGNORE,ierr)
              rbufe5p => dom(ib) % re5p
              do nn=2,nk-1
                 i=ni; j=nj; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe5p(nn)
              end do
           end if
        end if
!=====> next edge #6
        if (dom(ib)%edgnext6.ge.0) then
           if (dom_ad(dom_id(ib)) .eq. dom_ad(dom(ib)%edgnext6)) then
              do nn=2,nk-1
                 i=ni; j=1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=dom(ib) % re6p(nn)
              end do
           else
              call MPI_WAIT(dom(ib)%rq_e6p,MPI_STATUS_IGNORE,ierr)
              rbufe6p => dom(ib) % re6p
              do nn=2,nk-1
                 i=ni; j=1; k=nn; ll=(k-1)*nij+(i-1)*nj+j+nijk
                 fi(ll)=rbufe6p(nn)
              end do
           end if
        end if
!======================================================================
!======================================================================
        end if

        end do
!--------------------------------------------------------------------------

        end do

!==========================================================================


        call mgbound(g)

        if (PERIODIC) call exchbc_mgpp(g)

        return
        end subroutine exchangepp
!##########################################################################  
        subroutine  exchbc_mgpp(g)
!##########################################################################
        use multidata
        use mpi
        use vars
        implicit none
        integer :: i,j,k,ijk,ib,totdom,ijkc,ijkn,g,gl
        integer :: ni,nj,nk,nijk
        integer :: my_cor,tsend,tag1,tag2,tag3,tag4
        double precision, pointer, dimension(:)   :: fi
        double precision, pointer, dimension(:)   :: sbuf,rbuf

        MPI_FLT   = MPI_DOUBLE_PRECISION

        totdom=num_domains

!--------------------------------------------------------------------------
        do ib=1,nbp
        if(dom(ib)%bc_west.eq.5 .or. dom(ib)%bc_east.eq.5) then
              if(g.gt.dom(ib)%ngrid) then
                 gl=dom(ib)%ngrid
              else
                 gl=g
              end if
        ni=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        nj=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nk=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

        nijk=dom(ib)%faz(g)-ni*nj*nk
        fi => dom(ib)%cof

!..........................................................................
!=== West ===> 
!..........................................................................
           if (dom(ib)%iprev.lt.0) then

              if (dom(ib)%inext.lt.0) then
                 do k=1,nk
                    do j=1,nj
                       i=1;    ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       i=ni-1; ijkn=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)= fi(ijkn)
                    end do
                 end do

              else
                 my_cor=dom(ib)%per_ip
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    do k=1,nk
                       do j=1,nj
                       ijk=(k-1)*nj+j
                       i=2; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       dom(dom_indid(my_cor))%recvb_p1(ijk)=fi(ijkc)
                       end do
                    end do

                 else
                    tsend=nj*nk

                    sbuf => dom(ib) % sendb_m1
                    rbuf => dom(ib) % recvb_m1

                    do k=1,nk
                       do j=1,nj
                       ijk=(k-1)*nj+j
                       i=2; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       sbuf(ijk)=fi(ijkc)
                       end do
                    end do

                    tag1=1*10**8+my_cor    *10**4+dom_id(ib)
                    tag2=2*10**8+dom_id(ib)*10**4+my_cor
                 call MPI_IRECV  (dom(ib)%recvb_m1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag2,MPI_COMM_WORLD,dom(ib)%rq_m1,ierr)
                 call MPI_SEND (dom(ib)%sendb_m1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag1,MPI_COMM_WORLD,ierr)

                 end if

              end if
           end if
!..........................................................................
!=== East ===> 
!..........................................................................
           if (dom(ib)%inext.lt.0) then

              if (dom(ib)%iprev.lt.0) then
                 do k=1,nk
                    do j=1,nj
                       i=ni; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       i=2;  ijkn=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)= fi(ijkn)
                    end do
                 end do

              else
                 my_cor=dom(ib)%per_in
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    do k=1,nk
                       do j=1,nj
                       ijk=(k-1)*nj+j
                       i=ni-1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       dom(dom_indid(my_cor))%recvb_m1(ijk)=fi(ijkc)
                       end do
                    end do

                 else
                    tsend=nj*nk

                    sbuf => dom(ib) % sendb_p1
                    rbuf => dom(ib) % recvb_p1

                    do k=1,nk
                       do j=1,nj
                       ijk=(k-1)*nj+j
                       i=ni-1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       sbuf(ijk)=fi(ijkc)
                       end do
                    end do

                    tag3=2*10**8+my_cor    *10**4+dom_id(ib)
                    tag4=1*10**8+dom_id(ib)*10**4+my_cor
                 call MPI_IRECV  (dom(ib)%recvb_p1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag4,MPI_COMM_WORLD,dom(ib)%rq_p1,ierr)
                 call MPI_SEND (dom(ib)%sendb_p1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag3,MPI_COMM_WORLD,ierr)

                 end if

              end if
           end if
        end if
        end do
!**************************************************************************
!**************************************************************************
        do ib=1,nbp
        if(dom(ib)%bc_west.eq.5 .or. dom(ib)%bc_east.eq.5) then
              if(g.gt.dom(ib)%ngrid) then
                 gl=dom(ib)%ngrid
              else
                 gl=g
              end if
        ni=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        nj=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nk=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

        nijk=dom(ib)%faz(g)-ni*nj*nk
        fi => dom(ib)%cof

!..........................................................................
!=== West ===> 
!..........................................................................
           if (dom(ib)%iprev.lt.0 .and. dom(ib)%inext.ge.0) then
              my_cor=dom(ib)%per_ip

              if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                 do k=1,nk
                    do j=1,nj
                       ijk=(k-1)*nj+j
                       i=1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=dom(ib) % recvb_m1(ijk)
                    end do
                 end do

              else !if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_m1,MPI_STATUS_IGNORE,ierr)
                 rbuf => dom(ib) % recvb_m1

                 do k=1,nk
                    do j=1,nj
                       ijk=(k-1)*nj+j
                       i=1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=rbuf(ijk)
                    end do
                 end do
              end if
           end if
!..........................................................................
!=== East ===> 
!..........................................................................
           if (dom(ib)%inext.lt.0 .and. dom(ib)%iprev.ge.0) then
              my_cor=dom(ib)%per_in

              if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then

                 do k=1,nk
                    do j=1,nj
                       ijk=(k-1)*nj+j
                       i=ni; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=dom(ib) % recvb_p1(ijk)
                    end do
                 end do

              else !if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_p1,MPI_STATUS_IGNORE,ierr)
                 rbuf => dom(ib) % recvb_p1

                 do k=1,nk
                    do j=1,nj
                       ijk=(k-1)*nj+j
                       i=ni; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=rbuf(ijk)
                    end do
                 end do
              end if
           end if
        end if
        end do


!==========================================================================
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!==========================================================================


!--------------------------------------------------------------------------
        do ib=1,nbp
        if(dom(ib)%bc_south.eq.5 .or. dom(ib)%bc_north.eq.5) then
              if(g.gt.dom(ib)%ngrid) then
                 gl=dom(ib)%ngrid
              else
                 gl=g
              end if
        ni=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        nj=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nk=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

        nijk=dom(ib)%faz(g)-ni*nj*nk
        fi => dom(ib)%cof

!..........................................................................
!=== South ===> 
!..........................................................................
           if (dom(ib)%jprev.lt.0) then

              if (dom(ib)%jnext.lt.0) then
                 do k=1,nk
                    do i=1,ni
                       j=1;    ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       j=nj-1; ijkn=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)= fi(ijkn)
                    end do
                 end do

              else
                 my_cor=dom(ib)%per_jp
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    do k=1,nk
                       do i=1,ni
                       ijk=(k-1)*ni+i
                       j=2; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       dom(dom_indid(my_cor))%recvb_p1(ijk)=fi(ijkc)
                       end do
                    end do

                 else
                    tsend=ni*nk

                    sbuf => dom(ib) % sendb_m1
                    rbuf => dom(ib) % recvb_m1

                    do k=1,nk
                       do i=1,ni
                       ijk=(k-1)*ni+i
                       j=2; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       sbuf(ijk)=fi(ijkc)
                       end do
                    end do

                    tag1=3*10**8+my_cor    *10**4+dom_id(ib)
                    tag2=4*10**8+dom_id(ib)*10**4+my_cor
                 call MPI_IRECV  (dom(ib)%recvb_m1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag2,MPI_COMM_WORLD,dom(ib)%rq_m1,ierr)
                 call MPI_SEND (dom(ib)%sendb_m1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag1,MPI_COMM_WORLD,ierr)

                 end if

              end if
           end if

!..........................................................................
!=== North ===> 
!..........................................................................
           if (dom(ib)%jnext.lt.0) then

              if (dom(ib)%jprev.lt.0) then
                 do k=1,nk
                    do i=1,ni
                       j=nj; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       j=2;  ijkn=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)= fi(ijkn)
                    end do
                 end do

              else
                 my_cor=dom(ib)%per_jn 
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    do k=1,nk
                       do i=1,ni
                       ijk=(k-1)*ni+i
                       j=nj-1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       dom(dom_indid(my_cor))%recvb_m1(ijk)=fi(ijkc)
                       end do
                    end do

                 else
                    tsend=ni*nk

                    sbuf => dom(ib) % sendb_p1
                    rbuf => dom(ib) % recvb_p1

                    do k=1,nk
                       do i=1,ni
                       ijk=(k-1)*ni+i
                       j=nj-1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       sbuf(ijk)=fi(ijkc)
                       end do
                    end do

                    tag3=4*10**8+my_cor    *10**4+dom_id(ib)
                    tag4=3*10**8+dom_id(ib)*10**4+my_cor
                 call MPI_IRECV  (dom(ib)%recvb_p1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag4,MPI_COMM_WORLD,dom(ib)%rq_p1,ierr)
                 call MPI_SEND (dom(ib)%sendb_p1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag3,MPI_COMM_WORLD,ierr)

                 end if

              end if
           end if
        end if
        end do
!**************************************************************************
!**************************************************************************
        do ib=1,nbp
        if(dom(ib)%bc_south.eq.5 .or. dom(ib)%bc_north.eq.5) then
              if(g.gt.dom(ib)%ngrid) then
                 gl=dom(ib)%ngrid
              else
                 gl=g
              end if
        ni=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        nj=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nk=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

        nijk=dom(ib)%faz(g)-ni*nj*nk
        fi => dom(ib)%cof

!..........................................................................
!=== South ===> 
!..........................................................................
           if (dom(ib)%jprev.lt.0 .and. dom(ib)%jnext.ge.0) then
              my_cor=dom(ib)%per_jp

              if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                 do k=1,nk
                    do i=1,ni
                       ijk=(k-1)*ni+i
                       j=1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=dom(ib) % recvb_m1(ijk)
                    end do
                 end do

              else !if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_m1,MPI_STATUS_IGNORE,ierr)
                 rbuf => dom(ib) % recvb_m1

                 do k=1,nk
                    do i=1,ni
                       ijk=(k-1)*ni+i
                       j=1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=rbuf(ijk)
                    end do
                 end do
              end if
           end if
!..........................................................................
!=== North ===> 
!..........................................................................
           if (dom(ib)%jnext.lt.0 .and. dom(ib)%jprev.ge.0) then
              my_cor=dom(ib)%per_jn 

              if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then

                 do k=1,nk
                    do i=1,ni
                       ijk=(k-1)*ni+i
                       j=nj; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=dom(ib) % recvb_p1(ijk)
                    end do
                 end do

              else !if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_p1,MPI_STATUS_IGNORE,ierr)
                 rbuf => dom(ib) % recvb_p1

                 do k=1,nk
                    do i=1,ni
                       ijk=(k-1)*ni+i
                       j=nj; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=rbuf(ijk)
                    end do
                 end do
              end if
           end if
        end if
        end do


!==========================================================================
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!==========================================================================


!--------------------------------------------------------------------------
        do ib=1,nbp
        if(dom(ib)%bc_bottom.eq.5 .or. dom(ib)%bc_top.eq.5) then
              if(g.gt.dom(ib)%ngrid) then
                 gl=dom(ib)%ngrid
              else
                 gl=g
              end if
        ni=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        nj=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nk=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

        nijk=dom(ib)%faz(g)-ni*nj*nk
        fi => dom(ib)%cof

!..........................................................................
!=== Bottom ===> 
!..........................................................................
           if (dom(ib)%kprev.lt.0) then

              if (dom(ib)%knext.lt.0) then
                 do j=1,nj
                    do i=1,ni
                       k=1;    ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       k=nk-1; ijkn=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)= fi(ijkn)
                    end do
                 end do

              else
                 my_cor=dom(ib)%per_kp
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    do j=1,nj
                       do i=1,ni
                       ijk=(j-1)*ni+i
                       k=2; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       dom(dom_indid(my_cor))%recvb_p1(ijk)=fi(ijkc)
                       end do
                    end do

                 else
                    tsend=ni*nj

                    sbuf => dom(ib) % sendb_m1
                    rbuf => dom(ib) % recvb_m1

                    do j=1,nj
                       do i=1,ni
                       ijk=(j-1)*ni+i
                       k=2; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       sbuf(ijk)=fi(ijkc)
                       end do
                    end do

                    tag1=5*10**8+my_cor    *10**4+dom_id(ib)
                    tag2=6*10**8+dom_id(ib)*10**4+my_cor
                 call MPI_IRECV  (dom(ib)%recvb_m1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag2,MPI_COMM_WORLD,dom(ib)%rq_m1,ierr)
                 call MPI_SEND (dom(ib)%sendb_m1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag1,MPI_COMM_WORLD,ierr)

                 end if

              end if
           end if


!..........................................................................
!=== Top ===> 
!..........................................................................
           if (dom(ib)%knext.lt.0) then

              if (dom(ib)%kprev.lt.0) then
                 do j=1,nj
                    do i=1,ni
                       k=nk; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       k=2;  ijkn=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)= fi(ijkn)
                    end do
                 end do

              else
                 my_cor=dom(ib)%per_kn 
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    do j=1,nj
                       do i=1,ni
                       ijk=(j-1)*ni+i
                       k=nk-1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       dom(dom_indid(my_cor))%recvb_m1(ijk)=fi(ijkc)
                       end do
                    end do

                 else
                    tsend=ni*nj

                    sbuf => dom(ib) % sendb_p1
                    rbuf => dom(ib) % recvb_p1

                    do j=1,nj
                       do i=1,ni
                       ijk=(j-1)*ni+i
                       k=nk-1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       sbuf(ijk)=fi(ijkc)
                       end do
                    end do

                    tag3=6*10**8+my_cor    *10**4+dom_id(ib)
                    tag4=5*10**8+dom_id(ib)*10**4+my_cor
                 call MPI_IRECV  (dom(ib)%recvb_p1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag4,MPI_COMM_WORLD,dom(ib)%rq_p1,ierr)
                 call MPI_SEND (dom(ib)%sendb_p1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag3,MPI_COMM_WORLD,ierr)

                 end if

              end if
           end if
        end if
        end do
!**************************************************************************
!**************************************************************************
        do ib=1,nbp
        if(dom(ib)%bc_bottom.eq.5 .or. dom(ib)%bc_top.eq.5) then
              if(g.gt.dom(ib)%ngrid) then
                 gl=dom(ib)%ngrid
              else
                 gl=g
              end if
        ni=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        nj=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nk=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

        nijk=dom(ib)%faz(g)-ni*nj*nk
        fi => dom(ib)%cof

!..........................................................................
!=== Bottom ===> 
!..........................................................................
           if (dom(ib)%kprev.lt.0 .and. dom(ib)%knext.ge.0) then
              my_cor=dom(ib)%per_kp 

              if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                 do j=1,nj
                    do i=1,ni
                       ijk=(j-1)*ni+i
                       k=1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=dom(ib) % recvb_m1(ijk)
                    end do
                 end do

              else !if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_m1,MPI_STATUS_IGNORE,ierr)
                 rbuf => dom(ib) % recvb_m1

                 do j=1,nj
                    do i=1,ni
                       ijk=(j-1)*ni+i
                       k=1; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=rbuf(ijk)
                    end do
                 end do
              end if
           end if
!..........................................................................
!=== Top ===> 
!..........................................................................
           if (dom(ib)%knext.lt.0 .and. dom(ib)%kprev.ge.0) then
              my_cor=dom(ib)%per_kn 

              if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                 do j=1,nj
                    do i=1,ni
                       ijk=(j-1)*ni+i
                       k=nk; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=dom(ib) % recvb_p1(ijk)
                    end do
                 end do

              else !if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_p1,MPI_STATUS_IGNORE,ierr)
                 rbuf => dom(ib) % recvb_p1

                 do j=1,nj
                    do i=1,ni
                       ijk=(j-1)*ni+i
                       k=nk; ijkc=(k-1)*ni*nj+(i-1)*nj+j+nijk
                       fi(ijkc)=rbuf(ijk)
                    end do
                 end do
              end if
           end if
        end if
        end do


!==========================================================================
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!==========================================================================

        return
        end subroutine exchbc_mgpp
!##########################################################################
