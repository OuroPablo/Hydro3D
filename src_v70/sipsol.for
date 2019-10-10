!#######################################################################
        subroutine sipsol(ifi)
!####################################################################### 
!     This routine incorporates the Stone s SIP solver, based on  
!     ILU-decomposition.
!----------------------------------------------------------------------- 
        use vars
        use multidata
        use mpi
        implicit none
        integer i,j,k,ib,ijk,l,ifi,iff,xx,maxttc_ijk
        integer is,ie,js,je,ks,ke,nj,nij,nijk
        double precision, pointer, dimension(:,:,:) :: fi
        double precision, allocatable,dimension(:,:) :: ue,un,ut,lw
        double precision, allocatable,dimension(:,:) :: ls,lb,lpr,res
        double precision :: alfa,res1,resab,small,p1,p2,p3
        double precision :: rsm,reldif,relast
        double precision :: buffer_resab 

        MPI_FLT   = MPI_DOUBLE_PRECISION

           if(ifi.eq.11) then
              iff=1
           else if(ifi.eq.22) then
              iff=2
           else if(ifi.eq.33) then
              iff=3
           else if(ifi.eq.44) then
              iff=4
           else 
              print*,'error in sipsol'
           end if

        alfa=0.92
        small=1e-20

        xx=-1; maxttc_ijk=0
        do ib=1,nbp
           if(dom(ib)%ttc_ijk.gt.maxttc_ijk) then
              maxttc_ijk=dom(ib)%ttc_ijk
              xx=ib
           end if
        end do

        allocate(ue(nbp,dom(xx)%ttc_ijk),un(nbp,dom(xx)%ttc_ijk))
        allocate(ut(nbp,dom(xx)%ttc_ijk),lw(nbp,dom(xx)%ttc_ijk))
        allocate(ls(nbp,dom(xx)%ttc_ijk),lb(nbp,dom(xx)%ttc_ijk))
        allocate(lpr(nbp,dom(xx)%ttc_ijk),res(nbp,dom(xx)%ttc_ijk))

        ue=0.0; un=0.0; ut=0.0
        lb=0.0; lw=0.0; ls=0.0; lpr=0.0
        res=0.0

        do ib=1,nbp
           if(ifi.eq.11) then
              is=dom(ib)%isu; ie=dom(ib)%ieu
              js=dom(ib)%jsu; je=dom(ib)%jeu
              ks=dom(ib)%ksu; ke=dom(ib)%keu
           else if(ifi.eq.22) then
              is=dom(ib)%isv; ie=dom(ib)%iev
              js=dom(ib)%jsv; je=dom(ib)%jev
              ks=dom(ib)%ksv; ke=dom(ib)%kev
           else if(ifi.eq.33) then
              is=dom(ib)%isw; ie=dom(ib)%iew
              js=dom(ib)%jsw; je=dom(ib)%jew
              ks=dom(ib)%ksw; ke=dom(ib)%kew
           else
              is=dom(ib)%isp; ie=dom(ib)%iep
              js=dom(ib)%jsp; je=dom(ib)%jep
              ks=dom(ib)%ksp; ke=dom(ib)%kep
           end if

           nj=dom(ib)%ttc_j
           nij=dom(ib)%ttc_i*dom(ib)%ttc_j
           nijk=dom(ib)%ttc_ijk

!.....COEFFICIENTS OF UPPER AND LOWER TRIANG. MATRICES [L] & [U] 
           do k = ks,ke
              do i = is,ie
                 do j = js,je
                    ijk=j+(i-1)*nj+(k-1)*nij
                    lb(ib,ijk)=dom(ib)%ab(i,j,k)/
     &(1.+alfa*(un(ib,ijk-nij)+ue(ib,ijk-nij)))
                    lw(ib,ijk)=dom(ib)%aw(i,j,k)/
     &(1.+alfa*(un(ib,ijk-nj)+ut(ib,ijk-nj)))
                    ls(ib,ijk)=dom(ib)%as(i,j,k)/
     &(1.+alfa*(ue(ib,ijk-1)+ut(ib,ijk-1)))
                    p1=alfa*(lb(ib,ijk)*un(ib,ijk-nij)+
     &lw(ib,ijk)*un(ib,ijk-nj))
                    p2=alfa*(lb(ib,ijk)*ue(ib,ijk-nij)+
     &ls(ib,ijk)*ue(ib,ijk-1)) 
                    p3=alfa*(lw(ib,ijk)*ut(ib,ijk-nj)+
     &ls(ib,ijk)*ut(ib,ijk-1))
                    lpr(ib,ijk)=1./(dom(ib)%ap(i,j,k)+p1+p2+p3-
     &lb(ib,ijk)*ut(ib,ijk-nij)-lw(ib,ijk)*ue(ib,ijk-nj)-
     &ls(ib,ijk)*un(ib,ijk-1)+small)
                    un(ib,ijk)=(dom(ib)%an(i,j,k)-p1)*lpr(ib,ijk)
                    ue(ib,ijk)=(dom(ib)%ae(i,j,k)-p2)*lpr(ib,ijk)
                    ut(ib,ijk)=(dom(ib)%at(i,j,k)-p3)*lpr(ib,ijk)
                 end do
              end do 
           end do 

        end do



        relast=0.0
!.....INNER ITERATIONS LOOP  
        do nsweep=1,nswp(iff)
           resab=0.d0 

           do ib=1,nbp

              if(ifi.eq.11) then
                 fi => dom(ib)%ustar
                 is=dom(ib)%isu; ie=dom(ib)%ieu
                 js=dom(ib)%jsu; je=dom(ib)%jeu
                 ks=dom(ib)%ksu; ke=dom(ib)%keu
              else if(ifi.eq.22) then
                 fi => dom(ib)%vstar
                 is=dom(ib)%isv; ie=dom(ib)%iev
                 js=dom(ib)%jsv; je=dom(ib)%jev
                 ks=dom(ib)%ksv; ke=dom(ib)%kev
              else if(ifi.eq.33) then
                 fi => dom(ib)%wstar
                 is=dom(ib)%isw; ie=dom(ib)%iew
                 js=dom(ib)%jsw; je=dom(ib)%jew
                 ks=dom(ib)%ksw; ke=dom(ib)%kew
              else
                 fi => dom(ib)%pp
                 is=dom(ib)%isp; ie=dom(ib)%iep
                 js=dom(ib)%jsp; je=dom(ib)%jep
                 ks=dom(ib)%ksp; ke=dom(ib)%kep
              end if

              nj=dom(ib)%ttc_j
              nij=dom(ib)%ttc_i*dom(ib)%ttc_j
              nijk=dom(ib)%ttc_ijk

!.....CALCULATE RESIDUAL VECTOR AND THE SUM OF ABSOLUTE VALUES 
              do k = ks,ke
                 do i = is,ie
                    do j = js,je
                       ijk=j+(i-1)*nj+(k-1)*nij
                       res(ib,ijk)=(dom(ib)%su(i,j,k)-
     &dom(ib)%ap(i,j,k)*fi(i,j,k)-dom(ib)%an(i,j,k)*fi(i,j+1,k)-
     &dom(ib)%as(i,j,k)*fi(i,j-1,k)-dom(ib)%ae(i,j,k)*fi(i+1,j,k)-
     &dom(ib)%aw(i,j,k)*fi(i-1,j,k)-dom(ib)%ab(i,j,k)*fi(i,j,k-1)-
     &dom(ib)%at(i,j,k)*fi(i,j,k+1))
                       resab=resab+abs(res(ib,ijk))
                    end do 
                 end do 
              end do 
!.....FORWARD SUBSTITUTION 
              do k = ks,ke
                 do i = is,ie
                    do j = js,je
                       ijk=j+(i-1)*nj+(k-1)*nij
                       res(ib,ijk)=(res(ib,ijk)-
     &lb(ib,ijk)*res(ib,ijk-nij)-lw(ib,ijk)*res(ib,ijk-nj)-
     &ls(ib,ijk)*res(ib,ijk-1))*lpr(ib,ijk)
                    end do
                 end do 
              end do 
!.....BACKWARD SUBSTITUTION AND CORRECTION OF VARIABLE 
              do k=ke,ks,-1
                 do i=ie,is,-1
                    do j=je,js,-1
                       ijk=j+(i-1)*nj+(k-1)*nij
                       res(ib,ijk)=res(ib,ijk)-
     &un(ib,ijk)*res(ib,ijk+1)-ue(ib,ijk)*res(ib,ijk+nj)-
     &ut(ib,ijk)*res(ib,ijk+nij)
                       fi(i,j,k)=fi(i,j,k)+res(ib,ijk)
                    end do 
                 end do 
              end do

           end do

           buffer_resab = resab
           call MPI_ALLREDUCE(buffer_resab,resab,1,MPI_FLT,MPI_SUM,
     &MPI_COMM_WORLD,ierr)

           if(nsweep.eq.1) res1=resab
           rsm=resab/(res1+small)

           call exchange(ifi) 

!.....CHECK CONVERGENCE OF INNER ITERATIONS 
! --- TERMINATION JUDGEMENT
           reldif=abs (rsm - relast)
           relast=rsm
           if (rsm.le.0.001)    goto 2100
           if (reldif.le.0.00001) goto 2100 
        end do

        nsweep = nsweep - 1
 2100   continue

        deallocate (ue,un,ut,lw,ls,lb,lpr,res)


        return 
        end 
!#######################################################################
