!##########################################################################
        subroutine  exchange_bcphi(op,nly)
!##########################################################################
        use multidata
        use mpi
        use vars
        implicit none
        integer :: i,j,k,ijk,op,ib,tsend,ly,nly
        integer :: my_cor,tag1,tag2,tag3,tag4
        integer :: is,ie,js,je,ks,ke
        double precision, pointer, dimension(:,:,:) :: fi
        double precision, pointer, dimension(:)   :: sbuf,rbuf

        MPI_FLT   = MPI_DOUBLE_PRECISION

        do ly=0,nly      
!--------------------------------------------------------------------------
        do ib=1,nbp
           select case (op)
              Case (14)  
                 fi => dom(ib)%phi
              Case (15)  
                 fi => dom(ib)%phi_reinit
              Case (16)  
                 fi => dom(ib)%phi_new
              Case (17)  
                 fi => dom(ib)%phi_init
              Case (18)  
                 fi => dom(ib)%dens
              Case (19)  
                 fi => dom(ib)%mu
           end select

                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep

!..........................................................................
!=== West ===> 
!..........................................................................
           if (dom(ib)%bc_west.eq.5 .and. dom(ib)%iprev.lt.0) then
              if (dom(ib)%inext.lt.0) then
                 do k=1,dom(ib)%ttc_k; do j=1,dom(ib)%ttc_j
                    fi(is-1-ly,j,k)= fi(ie-ly,j,k)
                 end do; end do
              else
                 my_cor=dom(ib)%per_ip
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    sbuf => dom(dom_indid(my_cor))%recvb_p1
                 else
                    sbuf => dom(ib) % sendb_m1
                 end if

                 tsend=dom(ib)%ttc_j*dom(ib)%ttc_k
                 do k=1,dom(ib)%ttc_k; do j=1,dom(ib)%ttc_j
                    ijk=(k-1)*dom(ib)%ttc_j+j
                    sbuf(ijk)=fi(is+ly,j,k)
                 end do; end do

                 if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
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
           if (dom(ib)%bc_east.eq.5 .and. dom(ib)%inext.lt.0) then
              if (dom(ib)%iprev.lt.0) then
                 do k=1,dom(ib)%ttc_k; do j=1,dom(ib)%ttc_j
                    fi(ie+1+ly,j,k)= fi(is+ly,j,k)
                 end do; end do
              else
                 my_cor=dom(ib)%per_in
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    sbuf => dom(dom_indid(my_cor))%recvb_m1
                 else
                    sbuf => dom(ib) % sendb_p1
                 end if

                 tsend=dom(ib)%ttc_j*dom(ib)%ttc_k
                 do k=1,dom(ib)%ttc_k; do j=1,dom(ib)%ttc_j
                    ijk=(k-1)*dom(ib)%ttc_j+j
                    sbuf(ijk)=fi(ie-ly,j,k)
                 end do; end do

                 if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                    tag3=2*10**8+my_cor    *10**4+dom_id(ib)
                    tag4=1*10**8+dom_id(ib)*10**4+my_cor
                    call MPI_IRECV  (dom(ib)%recvb_p1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag4,MPI_COMM_WORLD,dom(ib)%rq_p1,ierr)
                    call MPI_SEND (dom(ib)%sendb_p1(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag3,MPI_COMM_WORLD,ierr)

                 end if
              end if
           end if
!..........................................................................
!=== South ===> 
!..........................................................................
           if (dom(ib)%bc_south.eq.5 .and. dom(ib)%jprev.lt.0) then
              if (dom(ib)%jnext.lt.0) then
                 do k=1,dom(ib)%ttc_k; do i=1,dom(ib)%ttc_i
                    fi(i,js-1-ly,k)= fi(i,je-ly,k)
                 end do; end do
              else
                 my_cor=dom(ib)%per_jp
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    sbuf => dom(dom_indid(my_cor))%recvb_p2
                 else
                    sbuf => dom(ib) % sendb_m2
                 end if

                 tsend=dom(ib)%ttc_i*dom(ib)%ttc_k
                 do k=1,dom(ib)%ttc_k; do i=1,dom(ib)%ttc_i
                    ijk=(k-1)*dom(ib)%ttc_i+i
                    sbuf(ijk)=fi(i,js+ly,k)
                 end do; end do

                 if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                    tag1=3*10**8+my_cor    *10**4+dom_id(ib)
                    tag2=4*10**8+dom_id(ib)*10**4+my_cor
                    call MPI_IRECV  (dom(ib)%recvb_m2(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag2,MPI_COMM_WORLD,dom(ib)%rq_m2,ierr)
                    call MPI_SEND (dom(ib)%sendb_m2(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag1,MPI_COMM_WORLD,ierr)

                 end if
              end if
           end if
!..........................................................................
!=== North ===> 
!..........................................................................
           if (dom(ib)%bc_north.eq.5 .and. dom(ib)%jnext.lt.0) then
              if (dom(ib)%jprev.lt.0) then
                 do k=1,dom(ib)%ttc_k; do i=1,dom(ib)%ttc_i
                    fi(i,je+1+ly,k)= fi(i,js+ly,k)
                 end do; end do
              else
                 my_cor=dom(ib)%per_jn
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    sbuf => dom(dom_indid(my_cor))%recvb_m2
                 else
                    sbuf => dom(ib) % sendb_p2
                 end if

                 tsend=dom(ib)%ttc_i*dom(ib)%ttc_k
                 do k=1,dom(ib)%ttc_k; do i=1,dom(ib)%ttc_i
                    ijk=(k-1)*dom(ib)%ttc_i+i
                    sbuf(ijk)=fi(i,je-ly,k)
                 end do; end do

                 if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                    tag3=4*10**8+my_cor    *10**4+dom_id(ib)
                    tag4=3*10**8+dom_id(ib)*10**4+my_cor
                    call MPI_IRECV  (dom(ib)%recvb_p2(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag4,MPI_COMM_WORLD,dom(ib)%rq_p2,ierr)
                    call MPI_SEND (dom(ib)%sendb_p2(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag3,MPI_COMM_WORLD,ierr)

                 end if
              end if
           end if
!..........................................................................
!=== Bottom ===> 
!..........................................................................
           if (dom(ib)%bc_bottom.eq.5 .and. dom(ib)%kprev.lt.0) then
              if (dom(ib)%knext.lt.0) then
                 do j=1,dom(ib)%ttc_j; do i=1,dom(ib)%ttc_i
                    fi(i,j,ks-1-ly)= fi(i,j,ke-ly)
                 end do; end do
              else
                 my_cor=dom(ib)%per_kp
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    sbuf => dom(dom_indid(my_cor))%recvb_p3
                 else
                    sbuf => dom(ib) % sendb_m3
                 end if

                 tsend=dom(ib)%ttc_i*dom(ib)%ttc_j
                 do j=1,dom(ib)%ttc_j; do i=1,dom(ib)%ttc_i
                    ijk=(j-1)*dom(ib)%ttc_i+i
                    sbuf(ijk)=fi(i,j,ks+ly)
                 end do; end do

                 if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                    tag1=5*10**8+my_cor    *10**4+dom_id(ib)
                    tag2=6*10**8+dom_id(ib)*10**4+my_cor
                    call MPI_IRECV  (dom(ib)%recvb_m3(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag2,MPI_COMM_WORLD,dom(ib)%rq_m3,ierr)
                    call MPI_SEND (dom(ib)%sendb_m3(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag1,MPI_COMM_WORLD,ierr)

                 end if
              end if
           end if
!..........................................................................
!=== Top ===> 
!..........................................................................
           if (dom(ib)%bc_top.eq.5 .and. dom(ib)%knext.lt.0) then
              if (dom(ib)%kprev.lt.0) then
                 do j=1,dom(ib)%ttc_j; do i=1,dom(ib)%ttc_i
                    fi(i,j,ke+1+ly)= fi(i,j,ks+ly)
                 end do; end do
              else
                 my_cor=dom(ib)%per_kn
                 if (dom_ad(dom_id(ib)) .eq. dom_ad(my_cor)) then
                    sbuf => dom(dom_indid(my_cor))%recvb_m3
                 else
                    sbuf => dom(ib) % sendb_p3
                 end if

                 tsend=dom(ib)%ttc_i*dom(ib)%ttc_j
                 do j=1,dom(ib)%ttc_j; do i=1,dom(ib)%ttc_i
                    ijk=(j-1)*dom(ib)%ttc_i+i
                    sbuf(ijk)=fi(i,j,ke-ly)
                 end do; end do

                 if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                    tag3=6*10**8+my_cor    *10**4+dom_id(ib)
                    tag4=5*10**8+dom_id(ib)*10**4+my_cor
                    call MPI_IRECV  (dom(ib)%recvb_p3(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag4,MPI_COMM_WORLD,dom(ib)%rq_p3,ierr)
                    call MPI_SEND (dom(ib)%sendb_p3(1),tsend,MPI_FLT,
     &dom_ad(my_cor),tag3,MPI_COMM_WORLD,ierr)

                 end if
              end if
           end if

        end do
!**************************************************************************
!**************************************************************************
        do ib=1,nbp
           select case (op)
              Case (14)  
                 fi => dom(ib)%phi
              Case (15)  
                 fi => dom(ib)%phi_reinit
              Case (16)  
                 fi => dom(ib)%phi_new
              Case (17)  
                 fi => dom(ib)%phi_init
              Case (18)  
                 fi => dom(ib)%dens
              Case (19)  
                 fi => dom(ib)%mu
           end select
!..........................................................................
!=== West ===> 
!..........................................................................
           if (dom(ib)%bc_west.eq.5 .and.
     & dom(ib)%iprev.lt.0 .and. dom(ib)%inext.ge.0) then
              my_cor=dom(ib)%per_ip
              if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_m1,MPI_STATUS_IGNORE,ierr)
              end if
              rbuf => dom(ib) % recvb_m1
              do k=1,dom(ib)%ttc_k; do j=1,dom(ib)%ttc_j
                 ijk=(k-1)*dom(ib)%ttc_j+j
                 fi(is-1-ly,j,k)=rbuf(ijk)
              end do; end do
           end if
!..........................................................................
!=== East ===> 
!..........................................................................
           if (dom(ib)%bc_east.eq.5 .and.
     & dom(ib)%inext.lt.0 .and. dom(ib)%iprev.ge.0) then
              my_cor=dom(ib)%per_in
              if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_p1,MPI_STATUS_IGNORE,ierr)
              end if
              rbuf => dom(ib) % recvb_p1
              do k=1,dom(ib)%ttc_k; do j=1,dom(ib)%ttc_j
                 ijk=(k-1)*dom(ib)%ttc_j+j
                 fi(ie+1+ly,j,k)=rbuf(ijk)
              end do; end do
           end if
!..........................................................................
!=== South ===> 
!..........................................................................
           if (dom(ib)%bc_south.eq.5 .and.
     & dom(ib)%jprev.lt.0 .and. dom(ib)%jnext.ge.0) then
              my_cor=dom(ib)%per_jp
              if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_m2,MPI_STATUS_IGNORE,ierr)
              end if
              rbuf => dom(ib) % recvb_m2
              do k=1,dom(ib)%ttc_k; do i=1,dom(ib)%ttc_i
                 ijk=(k-1)*dom(ib)%ttc_i+i
                 fi(i,js-1-ly,k)=rbuf(ijk)
              end do; end do
           end if
!..........................................................................
!=== North ===> 
!..........................................................................
           if (dom(ib)%bc_north.eq.5 .and.
     & dom(ib)%jnext.lt.0 .and. dom(ib)%jprev.ge.0) then
              my_cor=dom(ib)%per_jn
              if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_p2,MPI_STATUS_IGNORE,ierr)
              end if
              rbuf => dom(ib) % recvb_p2
              do k=1,dom(ib)%ttc_k; do i=1,dom(ib)%ttc_i
                 ijk=(k-1)*dom(ib)%ttc_i+i
                 fi(i,je+1+ly,k)=rbuf(ijk)
              end do; end do
           end if
!..........................................................................
!=== Bottom ===> 
!..........................................................................
           if (dom(ib)%bc_bottom.eq.5 .and.
     & dom(ib)%kprev.lt.0 .and. dom(ib)%knext.ge.0) then
              my_cor=dom(ib)%per_kp
              if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_m3,MPI_STATUS_IGNORE,ierr)
              end if
              rbuf => dom(ib) % recvb_m3
              do j=1,dom(ib)%ttc_j; do i=1,dom(ib)%ttc_i
                 ijk=(j-1)*dom(ib)%ttc_i+i
                 fi(i,j,ks-1-ly)=rbuf(ijk)
              end do; end do
           end if
!..........................................................................
!=== Top ===> 
!..........................................................................
           if (dom(ib)%bc_top.eq.5 .and.
     & dom(ib)%knext.lt.0 .and. dom(ib)%kprev.ge.0) then
              my_cor=dom(ib)%per_kn 
              if (dom_ad(dom_id(ib)) .ne. dom_ad(my_cor)) then
                 call MPI_WAIT(dom(ib)%rq_p3,MPI_STATUS_IGNORE,ierr)
              end if
              rbuf => dom(ib) % recvb_p3
              do j=1,dom(ib)%ttc_j; do i=1,dom(ib)%ttc_i
                 ijk=(j-1)*dom(ib)%ttc_i+i
                 fi(i,j,ke+1+ly)=rbuf(ijk)
              end do; end do
           end if
        end do

        end do
        return
        end subroutine exchange_bcphi
!##########################################################################
