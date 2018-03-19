!##########################################################################
        subroutine coeff
!##########################################################################
        use vars
        use multidata
        implicit none
        integer i,j,k,ijk,ib
        integer ijkp,ijke,ijkw,ijkn,ijks,ijkb,ijkt,ijkphi,ijksu
        integer ij1,ij2,ij3,ij4,ij5,ij6,ij7,ij8,pre_nipl,pre_njpl,
     & pre_nkpl
        integer glevel,gl,mgc_i,mgc_j,mgc_k,cnt,incr,ijk_lsm
        integer prmgci,prmgcj,prmgck
        double precision :: dxx,dyy,dzz,ndx,ndy,ndz

        do ib=1,nbp

	if (L_LSM) then
         ijk_lsm=0
         dom(ib)%dens_mg=0.0
	endif
	
        dom(ib)%faz(1)=(dom(ib)%iep-dom(ib)%isp+3)*
     & (dom(ib)%jep-dom(ib)%jsp+3)*
     & (dom(ib)%kep-dom(ib)%ksp+3)
        prmgci=(dom(ib)%iep-dom(ib)%isp+3)
        prmgcj=(dom(ib)%jep-dom(ib)%jsp+3)
        prmgck=(dom(ib)%kep-dom(ib)%ksp+3)
        do glevel=2,ngrd_gl
           if(glevel.gt.dom(ib)%ngrid) then
              gl=dom(ib)%ngrid
           else
              gl=glevel
           end if
           mgc_i=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
           mgc_j=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
           mgc_k=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2
           dom(ib)%faz(glevel)=dom(ib)%faz(glevel-1)+
     & mgc_i*mgc_j*mgc_k+8*(prmgci-2)*(prmgcj-2)*(prmgck-2)
           prmgci=mgc_i; prmgcj=mgc_j; prmgck=mgc_k
        end do

        cnt=0
        do glevel=1,ngrd_gl
           if(glevel.gt.dom(ib)%ngrid) then
              gl=dom(ib)%ngrid
           else
              gl=glevel
           end if
           mgc_i=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
           mgc_j=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
           mgc_k=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

           ndx=dom(ib)%dx*2**(gl-1)
           ndy=dom(ib)%dy*2**(gl-1)
           ndz=dom(ib)%dz*2**(gl-1)

           dxx=ndx*ndx
           dyy=ndy*ndy
           dzz=ndz*ndz

           cnt=dom(ib)%faz(glevel)-mgc_i*mgc_j*mgc_k
           do  k=1,mgc_k
              do  i=1,mgc_i
                 do j=1,mgc_j
                    ijk=cnt+(k-1)*mgc_i*mgc_j+(i-1)*mgc_j+j
                    dom(ib)%cof(ijk)=0.0
                 end do
              end do
           end do

           if (L_LSM) then! .or. L_LSMbase) THEN

             do  k=2,mgc_k-1
               do  i=2,mgc_i-1
                 do j=2,mgc_j-1 
            
             if (glevel.eq.1)  then

             ijk_lsm=(k-2)*(mgc_i-2)*(mgc_j-2)+(i-2)*(mgc_j-2)+(j-1)

       	 dom(ib)%dens_mg(ijk_lsm)=dom(ib)%dens(i+pl-1,j+pl-1,k+pl-1) 
             else if (glevel.ne.1)  then
             ijk_lsm=dom(ib)%ijkp_lsm(glevel-1)+
     &(k-2)*(mgc_i-2)*(mgc_j-2)+(i-2)*(mgc_j-2)+(j-1)

             pre_nipl=(dom(ib)%ttc_i-(2*pl))/2**(glevel-2)+2 
             pre_njpl=(dom(ib)%ttc_j-(2*pl))/2**(glevel-2)+2 
             pre_nkpl=(dom(ib)%ttc_k-(2*pl))/2**(glevel-2)+2 

             ij1 = (2*k-4)*(pre_nipl-2)*(pre_njpl-2) + 
     &(2*i-4)*(pre_njpl-2) + (2*j-3) + dom(ib)%ijkp_lsm(glevel-2)
             ij2 = ij1 + 1
             ij3 = ij1 + (pre_njpl-2)
             ij4 = ij1 + (pre_njpl-2) + 1
             ij5 = ij1 + (pre_nipl-2)*(pre_njpl-2)
             ij6 = ij1 + (pre_nipl-2)*(pre_njpl-2) + 1
             ij7 = ij1 + (pre_nipl-2)*(pre_njpl-2) + (pre_njpl-2)
             ij8 = ij1 + (pre_nipl-2)*(pre_njpl-2) + (pre_njpl-2) +1

             dom(ib)%dens_mg(ijk_lsm)=0.125*(dom(ib)%dens_mg(ij1)+
     &dom(ib)%dens_mg(ij2)+dom(ib)%dens_mg(ij3)+dom(ib)%dens_mg(ij4) +
     &dom(ib)%dens_mg(ij5)+dom(ib)%dens_mg(ij6)+dom(ib)%dens_mg(ij7) +
     &dom(ib)%dens_mg(ij8))

             end if

                end do             
              end do
            end do 

           end if

           incr=(mgc_i-2)*(mgc_j-2)*(mgc_k-2)
           do  k=2,mgc_k-1
              do  i=2,mgc_i-1
                 do j=2,mgc_j-1

                    ijk=dom(ib)%faz(glevel)+
     & (k-2)*(mgc_i-2)*(mgc_j-2)+(i-2)*(mgc_j-2)+(j-1)

	     if (L_LSM) then! .or. L_LSMbase) then

             if (glevel.eq.1) ijk_lsm=(k-2)*(mgc_i-2)*(mgc_j-2)+
     & (i-2)*(mgc_j-2)+(j-1)

             if (glevel.ne.1) ijk_lsm=dom(ib)%ijkp_lsm(glevel-1)+
     & (k-2)*(mgc_i-2)*(mgc_j-2)+(i-2)*(mgc_j-2)+(j-1)
                      
             if (i.eq.2) then
               densim12=dom(ib)%dens_mg(ijk_lsm)
             else
               densim12=0.5*(dom(ib)%dens_mg(ijk_lsm)+
     & dom(ib)%dens_mg(ijk_lsm-(mgc_j-2)))
             end if
 
             if (i.eq.mgc_i-1) then
               densip12=dom(ib)%dens_mg(ijk_lsm)
             else
               densip12=0.5*(dom(ib)%dens_mg(ijk_lsm)+    
     & dom(ib)%dens_mg(ijk_lsm+(mgc_j-2)))
             end if

             if (j.eq.2) then
               densjm12=dom(ib)%dens_mg(ijk_lsm)
             else
               densjm12=0.5*(dom(ib)%dens_mg(ijk_lsm)+
     & dom(ib)%dens_mg(ijk_lsm-1))
             end if

             if (j.eq.mgc_j-1) then
               densjp12=dom(ib)%dens_mg(ijk_lsm)       
             else
               densjp12=0.5*(dom(ib)%dens_mg(ijk_lsm)+
     & dom(ib)%dens_mg(ijk_lsm+1))
             end if

             if (k.eq.2) then
               denskm12=dom(ib)%dens_mg(ijk_lsm)
             else
               denskm12=0.5*(dom(ib)%dens_mg(ijk_lsm)+
     & dom(ib)%dens_mg(ijk_lsm-(mgc_i-2)*(mgc_j-2)))
             end if

             if (k.eq.mgc_k-1) then
               denskp12=dom(ib)%dens_mg(ijk_lsm)       
             else
               denskp12=0.5*(dom(ib)%dens_mg(ijk_lsm)+
     & dom(ib)%dens_mg(ijk_lsm+(mgc_i-2)*(mgc_j-2)))
             end if

             dom(ib)%cof(ijk)=dom(ib)%dy*dom(ib)%dz/dom(ib)%dx/
     & densim12 !aw
           dom(ib)%cof(ijk+1*incr)=dom(ib)%dy*dom(ib)%dz/dom(ib)%dx/
     & densip12 !ae
           dom(ib)%cof(ijk+2*incr)=dom(ib)%dx*dom(ib)%dz/dom(ib)%dy/
     & densjm12 !as
          dom(ib)%cof(ijk+3*incr)=dom(ib)%dx*dom(ib)%dz/dom(ib)%dy/
     & densjp12 !an
          dom(ib)%cof(ijk+4*incr)=dom(ib)%dx*dom(ib)%dy/dom(ib)%dz/
     & denskm12 !ab
          dom(ib)%cof(ijk+5*incr)=dom(ib)%dx*dom(ib)%dy/dom(ib)%dz/
     & denskp12 !at

	     else
                    dom(ib)%cof(ijk)       =1.0/dxx !aw
                    dom(ib)%cof(ijk+1*incr)=1.0/dxx !ae
                    dom(ib)%cof(ijk+2*incr)=1.0/dyy !as
                    dom(ib)%cof(ijk+3*incr)=1.0/dyy !an
                    dom(ib)%cof(ijk+4*incr)=1.0/dzz !ab
                    dom(ib)%cof(ijk+5*incr)=1.0/dzz !at
	     end if

                    if (dom(ib)%iprev.lt.0 .and. dom(ib)%bc_west.ne.5
     & .and. i.eq.2)         dom(ib)%cof(ijk)       =0.0
                    if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5
     & .and. i.eq.mgc_i-1)  dom(ib)%cof(ijk+1*incr)=0.0
                    if (dom(ib)%jprev.lt.0 .and. dom(ib)%bc_south.ne.5 
     & .and. j.eq.2)         dom(ib)%cof(ijk+2*incr)=0.0
                    if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5
     & .and. j.eq.mgc_j-1)  dom(ib)%cof(ijk+3*incr)=0.0
                    if (dom(ib)%kprev.lt.0 .and. dom(ib)%bc_bottom.ne.5 
     & .and. k.eq.2)         dom(ib)%cof(ijk+4*incr)=0.0
                    if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5
     & .and. k.eq.mgc_k-1)  dom(ib)%cof(ijk+5*incr)=0.0

                    dom(ib)%cof(ijk+6*incr)=-1.0*(dom(ib)%cof(ijk)+
     & dom(ib)%cof(ijk+1*incr)+dom(ib)%cof(ijk+2*incr)+
     & dom(ib)%cof(ijk+3*incr)+dom(ib)%cof(ijk+4*incr)+
     & dom(ib)%cof(ijk+5*incr)) !ap

                    dom(ib)%cof(ijk+7*incr)=0.0 !res

                 end do
              end do
           end do

        end do

        end do


        return
        end subroutine coeff
!##########################################################################
        subroutine mgkcyc
!##########################################################################
        use vars
        use multidata
        implicit none
        integer i,j,k,ib,ijk,l
        integer glevel,mgc_i,mgc_j,mgc_k,kcycle,nrel
        integer incr,incrp
        integer, allocatable,dimension(:) :: kount

        allocate (kount(ngrd_gl))

        do ib=1,nbp
           mgc_i=(dom(ib)%iep-dom(ib)%isp+1)+2
           mgc_j=(dom(ib)%jep-dom(ib)%jsp+1)+2
           mgc_k=(dom(ib)%kep-dom(ib)%ksp+1)+2

        do  k=1,mgc_k
           do  i=1,mgc_i
              do j=1,mgc_j
                 ijk=(k-1)*mgc_i*mgc_j+(i-1)*mgc_j+j
                 dom(ib)%cof(ijk)=dom(ib)%pp(i+pl-1,j+pl-1,k+pl-1)
              end do
           end do
        end do

        incrp=7*(mgc_i-2)*(mgc_j-2)*(mgc_k-2)
        do  k=2,mgc_k-1
           do  i=2,mgc_i-1
              do j=2,mgc_j-1
                 ijk=dom(ib)%faz(1)+(k-2)*(mgc_i-2)*(mgc_j-2)+
     & (i-2)*(mgc_j-2)+(j-1)+incrp
	           if (L_LSM) then  
                   dom(ib)%cof(ijk)=dom(ib)%sup(i+pl-1,j+pl-1,k+pl-1)
		     else
                   dom(ib)%cof(ijk)=dom(ib)%su(i+pl-1,j+pl-1,k+pl-1)
		     endif
              end do
           end do
        end do

        end do



        glevel=1
        kcycle=1

c pre-relax at current finest grid level
        call mgrelax(1,irestr)

c if at coarsest level, post-relax
        if (ngrd_gl.eq.1) goto 5

c calculate residual and restrict it to ngrid-1 
c (note: the residuals are stored in the memory space used by the
c the rhs array rhsf, which is therefore erased)
        call mgrestr(glevel+1) !send coarser grid level

c set counter for grid levels to zero
        do l=1,ngrd_gl
           kount(l)=0
        end do

c set new level and continue K-cycling
        glevel=glevel+1
        nrel=irestr

c K-cycle control point
10    continue

c post-relax when kcur revisited
        if (glevel.eq.1) goto 5

c count "hit" at current level
        kount(glevel)=kount(glevel)+1

c relax at current level
        call mgrelax(glevel,nrel)


        if (kount(glevel).eq.(kcycle+1)) then

c K-cycle(iprer,ipost) complete at glevel
c inject correction to finer grid
           call mgcorr(glevel-1) !send finer grid level

c reset counter to zero at glevel
           kount(glevel)=0

c ascend to next higher level and set to post-relax there
           glevel=glevel-1
           nrel=iproln
           goto 10
        
        else

c K-cycle not complete so descend unless at coarsest
           if (glevel.lt.ngrd_gl) then
              call mgrestr(glevel+1) !send coarser grid level

c pre-relax at next coarser level
              glevel=glevel+1
              nrel=irestr
              goto 10

           else
c post-relax at coarsest level
              call mgrelax(glevel,iproln)

c inject correction to grid level 2
              call mgcorr(glevel-1) !send finer grid level

c set to post-relax at level ngrid-1
              nrel=iproln
              glevel=ngrd_gl-1
              goto 10

           end if

        end if

5     continue

c post-relax at finest level
        call mgrelax(1,iproln)

        do ib=1,nbp
           mgc_i=(dom(ib)%iep-dom(ib)%isp+1)+2
           mgc_j=(dom(ib)%jep-dom(ib)%jsp+1)+2
           mgc_k=(dom(ib)%kep-dom(ib)%ksp+1)+2

           do  k=1,mgc_k
              do  i=1,mgc_i
                 do j=1,mgc_j
                    ijk=(k-1)*mgc_i*mgc_j+(i-1)*mgc_j+j
                    dom(ib)%pp(i+pl-1,j+pl-1,k+pl-1)=dom(ib)%cof(ijk)
                 end do
              end do
           end do
        end do

        deallocate (kount)

        return
        end subroutine mgkcyc
!##########################################################################
        subroutine mgrelax(glevel,it)
!##########################################################################
        use vars
        use multidata
        implicit none
        integer i,j,k,ib
        integer ijkp,ijke,ijkw,ijkn,ijks,ijkb,ijkt,ijkphi,ijksu
        integer glevel,it,gl
        integer incr,itr,ni,nj,nk,nij,nij2,nijk
        real apr
        real, dimension(ngg) :: a,c


        do itr=1,it

        do ib=1,nbp

        if(glevel.gt.dom(ib)%ngrid) then
           gl=dom(ib)%ngrid
        else
           gl=glevel
        end if

        ni=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        nj=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nk=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2
        nij=ni*nj
        nij2=(ni-2)*(nj-2)
        nijk=ni*nj*nk

        incr=(ni-2)*(nj-2)*(nk-2)

        if(mg_itrsch.eq.1) then
           do k=2,nk-1
              do j=2,nj-1
                 do  i=2,ni-1

                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 ijkw=dom(ib)%faz(glevel)+(k-2)*nij2+(i-2)*(nj-2)+(j-1)
                    ijke=ijkw+1*incr
                    ijks=ijkw+2*incr
                    ijkn=ijkw+3*incr
                    ijkb=ijkw+4*incr
                    ijkt=ijkw+5*incr
                    ijkp=ijkw+6*incr
                    ijksu=ijkw+7*incr

                    dom(ib)%cof(ijkphi)=((dom(ib)%cof(ijksu)-
     & (dom(ib)%cof(ijkw)*dom(ib)%cof(ijkphi-nj)+
     & dom(ib)%cof(ijke)*dom(ib)%cof(ijkphi+nj)+
     & dom(ib)%cof(ijks)*dom(ib)%cof(ijkphi-1)+
     & dom(ib)%cof(ijkn)*dom(ib)%cof(ijkphi+1)+
     & dom(ib)%cof(ijkb)*dom(ib)%cof(ijkphi-nij)+
     & dom(ib)%cof(ijkt)*dom(ib)%cof(ijkphi+nij)))/dom(ib)%cof(ijkp))

                 end do
              end do
           end do

        else if(mg_itrsch.eq.2) then

           a=0.0
           c=0.0

C.....SOLVE WITH "TDMA" ALONG  J-LINES
           do k=2,nk-1
              do i=2,ni-1
                 j=1
                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 c(1)=dom(ib)%cof(ijkphi)

                 do  j=2,nj-1

                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 ijkw=dom(ib)%faz(glevel)+(k-2)*nij2+(i-2)*(nj-2)+(j-1)
                    ijke=ijkw+1*incr
                    ijks=ijkw+2*incr
                    ijkn=ijkw+3*incr
                    ijkb=ijkw+4*incr
                    ijkt=ijkw+5*incr
                    ijkp=ijkw+6*incr
                    ijksu=ijkw+7*incr


                    apr=1.0/(dom(ib)%cof(ijkp)-dom(ib)%cof(ijks)*a(j-1)) 
                    a(j)=dom(ib)%cof(ijkn)*apr
                    c(j)=(dom(ib)%cof(ijksu)-
     & dom(ib)%cof(ijke)*dom(ib)%cof(ijkphi+nj)-
     & dom(ib)%cof(ijkw)*dom(ib)%cof(ijkphi-nj)-
     & dom(ib)%cof(ijkt)*dom(ib)%cof(ijkphi+nij)- 
     & dom(ib)%cof(ijkb)*dom(ib)%cof(ijkphi-nij)-
     & dom(ib)%cof(ijks)*c(j-1))*apr
                 end do

                 do j=nj-1,2,-1
                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 dom(ib)%cof(ijkphi)=(c(j)-a(j)*dom(ib)%cof(ijkphi+1))

                 end do
              end do
           end do

C.....SOLVE WITH "TDMA" ALONG  I-LINES
           do k=2,nk-1
              do j=2,nj-1
                 i=1
                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 c(1)=dom(ib)%cof(ijkphi)

                 do  i=2,ni-1

                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 ijkw=dom(ib)%faz(glevel)+(k-2)*nij2+(i-2)*(nj-2)+(j-1)
                    ijke=ijkw+1*incr
                    ijks=ijkw+2*incr
                    ijkn=ijkw+3*incr
                    ijkb=ijkw+4*incr
                    ijkt=ijkw+5*incr
                    ijkp=ijkw+6*incr
                    ijksu=ijkw+7*incr

                    apr=1.0/(dom(ib)%cof(ijkp)-dom(ib)%cof(ijkw)*a(i-1)) 
                    a(i)=dom(ib)%cof(ijke)*apr
                    c(i)=(dom(ib)%cof(ijksu)-
     & dom(ib)%cof(ijkn)*dom(ib)%cof(ijkphi+1)-
     & dom(ib)%cof(ijks)*dom(ib)%cof(ijkphi-1)-
     & dom(ib)%cof(ijkt)*dom(ib)%cof(ijkphi+nij)-
     & dom(ib)%cof(ijkb)*dom(ib)%cof(ijkphi-nij)-
     & dom(ib)%cof(ijkw)*c(i-1))*apr
                 end do

                 do i=ni-1,2,-1
                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 dom(ib)%cof(ijkphi)=(c(i)-a(i)*dom(ib)%cof(ijkphi+nj))

                 end do
              end do
           end do

C.....SOLVE WITH "TDMA" ALONG  K-LINES
           do i=2,ni-1
              do j=2,nj-1
                 k=1
                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 c(1)=dom(ib)%cof(ijkphi)

                 do  k=2,nk-1

                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 ijkw=dom(ib)%faz(glevel)+(k-2)*nij2+(i-2)*(nj-2)+(j-1)
                    ijke=ijkw+1*incr
                    ijks=ijkw+2*incr
                    ijkn=ijkw+3*incr
                    ijkb=ijkw+4*incr
                    ijkt=ijkw+5*incr
                    ijkp=ijkw+6*incr
                    ijksu=ijkw+7*incr

                    apr=1.0/(dom(ib)%cof(ijkp)-dom(ib)%cof(ijkb)*a(k-1)) 
                    a(k)=dom(ib)%cof(ijkt)*apr
                    c(k)=(dom(ib)%cof(ijksu)-
     & dom(ib)%cof(ijkn)*dom(ib)%cof(ijkphi+1)-
     & dom(ib)%cof(ijks)*dom(ib)%cof(ijkphi-1)-
     & dom(ib)%cof(ijke)*dom(ib)%cof(ijkphi+nj)-
     & dom(ib)%cof(ijkw)*dom(ib)%cof(ijkphi-nj)-
     & dom(ib)%cof(ijkb)*c(k-1))*apr
                 end do

                 do k=nk-1,2,-1
                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 dom(ib)%cof(ijkphi)=(c(k)-a(k)*dom(ib)%cof(ijkphi+nij))

                 end do
              end do
           end do

        else if(mg_itrsch.eq.3) then

           do k=2,nk-1
              do j=2,nj-1
                 do  i=2,ni-1
                 if (MOD(i+j+k,2).eq.0) then
                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 ijkw=dom(ib)%faz(glevel)+(k-2)*nij2+(i-2)*(nj-2)+(j-1)
                    ijke=ijkw+1*incr
                    ijks=ijkw+2*incr
                    ijkn=ijkw+3*incr
                    ijkb=ijkw+4*incr
                    ijkt=ijkw+5*incr
                    ijkp=ijkw+6*incr
                    ijksu=ijkw+7*incr

                    dom(ib)%cof(ijkphi)=((dom(ib)%cof(ijksu)-
     & (dom(ib)%cof(ijkw)*dom(ib)%cof(ijkphi-nj)+
     & dom(ib)%cof(ijke)*dom(ib)%cof(ijkphi+nj)+
     & dom(ib)%cof(ijks)*dom(ib)%cof(ijkphi-1)+
     & dom(ib)%cof(ijkn)*dom(ib)%cof(ijkphi+1)+
     & dom(ib)%cof(ijkb)*dom(ib)%cof(ijkphi-nij)+
     & dom(ib)%cof(ijkt)*dom(ib)%cof(ijkphi+nij)))/dom(ib)%cof(ijkp))

                 end if
                 end do
              end do
           end do

           do k=2,nk-1
              do j=2,nj-1
                 do  i=2,ni-1
                 if (MOD(i+j+k,2).eq.1) then
                 ijkphi=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 ijkw=dom(ib)%faz(glevel)+(k-2)*nij2+(i-2)*(nj-2)+(j-1)
                    ijke=ijkw+1*incr
                    ijks=ijkw+2*incr
                    ijkn=ijkw+3*incr
                    ijkb=ijkw+4*incr
                    ijkt=ijkw+5*incr
                    ijkp=ijkw+6*incr
                    ijksu=ijkw+7*incr

                    dom(ib)%cof(ijkphi)=((dom(ib)%cof(ijksu)-
     & (dom(ib)%cof(ijkw)*dom(ib)%cof(ijkphi-nj)+
     & dom(ib)%cof(ijke)*dom(ib)%cof(ijkphi+nj)+
     & dom(ib)%cof(ijks)*dom(ib)%cof(ijkphi-1)+
     & dom(ib)%cof(ijkn)*dom(ib)%cof(ijkphi+1)+
     & dom(ib)%cof(ijkb)*dom(ib)%cof(ijkphi-nij)+
     & dom(ib)%cof(ijkt)*dom(ib)%cof(ijkphi+nij)))/dom(ib)%cof(ijkp))

                 end if
                 end do
              end do
           end do

        else 
           print*,'error about mg iteration selection'
        end if


        end do

        call mgbound(glevel)
        call exchangepp(glevel)
        end do
 
!        call mgbound(glevel)
     
        return
        end subroutine mgrelax
!##########################################################################
        subroutine mgrestr(glevel)
!##########################################################################
        use vars
        use multidata
        implicit none
        integer i,j,k,ib
        integer ijkp,ijke,ijkw,ijkn,ijks,ijkb,ijkt,ijkphi,ijksu
        integer glevel,gl
        integer nif,njf,nkf,nic,njc,nkc
        integer nijc,nijf,nijc2,nijf2,nijkc,nijkf
        integer incr_c,incr_f
        real, allocatable,dimension(:,:,:) :: resf


        do ib=1,nbp

        if(glevel.gt.dom(ib)%ngrid) then
           gl=dom(ib)%ngrid
        else
           gl=glevel
        end if

        nic=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        njc=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nkc=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

        nijc=nic*njc
        nijc2=(nic-2)*(njc-2)
        nijkc=nic*njc*nkc
        incr_c=(nic-2)*(njc-2)*(nkc-2)

        do  k=1,nkc
           do  i=1,nic
              do j=1,njc
              ijkphi=(k-1)*nijc+(i-1)*njc+j+dom(ib)%faz(glevel)-nijkc
                 dom(ib)%cof(ijkphi)=0.0
              end do
           end do
        end do

        do  k=2,nkc-1
           do  i=2,nic-1
              do j=2,njc-1
                 ijksu=dom(ib)%faz(glevel)+(k-2)*nijc2+
     & (i-2)*(njc-2)+(j-1)+7*incr_c
                 dom(ib)%cof(ijksu)=0.0
              end do
           end do
        end do

c calbculate residual for fine grid
        if(glevel.gt.dom(ib)%ngrid) then
           nif=nic
           njf=njc
           nkf=nkc
        else
           nif=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-2)+2
           njf=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-2)+2
           nkf=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-2)+2
        end if


        nijf=nif*njf
        nijf2=(nif-2)*(njf-2)
        nijkf=nif*njf*nkf
        incr_f=(nif-2)*(njf-2)*(nkf-2)

        allocate(resf(nif,njf,nkf))

        do  k=2,nkf-1
           do  i=2,nif-1
              do j=2,njf-1

           ijkphi=(k-1)*nijf+(i-1)*njf+j+dom(ib)%faz(glevel-1)-nijkf
           ijkw=dom(ib)%faz(glevel-1)+(k-2)*nijf2+(i-2)*(njf-2)+(j-1)
              ijke=ijkw+1*incr_f
              ijks=ijkw+2*incr_f
              ijkn=ijkw+3*incr_f
              ijkb=ijkw+4*incr_f
              ijkt=ijkw+5*incr_f
              ijkp=ijkw+6*incr_f
              ijksu=ijkw+7*incr_f

              resf(i,j,k)=dom(ib)%cof(ijksu)-
     & (dom(ib)%cof(ijkw)*dom(ib)%cof(ijkphi-njf)+
     & dom(ib)%cof(ijke)*dom(ib)%cof(ijkphi+njf)+
     & dom(ib)%cof(ijks)*dom(ib)%cof(ijkphi-1)+
     & dom(ib)%cof(ijkn)*dom(ib)%cof(ijkphi+1)+
     & dom(ib)%cof(ijkb)*dom(ib)%cof(ijkphi-nijf)+
     & dom(ib)%cof(ijkt)*dom(ib)%cof(ijkphi+nijf)+
     & dom(ib)%cof(ijkp)*dom(ib)%cof(ijkphi))
              end do
           end do
        end do

c calculate the right-hand side at the coarser grid
c level from the averages of the values at the 4 surrounding points;
c if there is no coarsifying in one direction, only 2 points are
c used; no exchange of boundary data is necessary

        if(glevel.gt.dom(ib)%ngrid) then

        do  k=2,nkc-1
           do  i=2,nic-1
              do j=2,njc-1

              ijksu=dom(ib)%faz(glevel)+(k-2)*nijc2+
     & (i-2)*(njc-2)+(j-1)+7*incr_c


              dom(ib)%cof(ijksu)=resf(i,j,k)        
              end do
           end do
        end do

        else

        do  k=2,nkc-1
           do  i=2,nic-1
              do j=2,njc-1

              ijksu=dom(ib)%faz(glevel)+(k-2)*nijc2+
     & (i-2)*(njc-2)+(j-1)+7*incr_c

              dom(ib)%cof(ijksu)=0.125*(resf(2*i-2,2*j-2,2*k-2)+
     & resf(2*i-1,2*j-2,2*k-2)+resf(2*i-2,2*j-1,2*k-2)+
     & resf(2*i-1,2*j-1,2*k-2)+resf(2*i-2,2*j-2,2*k-1)+
     & resf(2*i-1,2*j-2,2*k-1)+resf(2*i-2,2*j-1,2*k-1)+
     & resf(2*i-1,2*j-1,2*k-1))

              end do
           end do
        end do

        end if 

        deallocate(resf)

        end do


        return
        end subroutine mgrestr
!##########################################################################
        subroutine mgcorr(glevel)
!##########################################################################
        use vars
        use multidata
        implicit none
        integer i,j,k,ib,ijkc,ijkf
        integer glevel,gl
        integer ic,jc,kc,nic,njc,nkc,nif,njf,nkf
        integer nijf,nijkf,nijc,nijkc


        do ib=1,nbp

        if(glevel.gt.dom(ib)%ngrid) then
           gl=dom(ib)%ngrid
        else
           gl=glevel
        end if

        nif=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        njf=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nkf=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2

        nijf=nif*njf
        nijkf=nif*njf*nkf

c calculate nodes for coarser grid
        if(glevel.ge.dom(ib)%ngrid) then
           nic=nif
           njc=njf
           nkc=nkf
        else
           nic=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl)+2
           njc=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl)+2
           nkc=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl)+2
        end if

        nijc=nic*njc
        nijkc=nic*njc*nkc

c new version: the correction is the weighted average of either two 
c or four points at the coarser grid level depending on whether 
c coarsifying takes place in all directions or not

        if(glevel.ge.dom(ib)%ngrid) then
        do  kc=1,nkc
           do  ic=1,nic
              do jc=1,njc

           ijkc=(kc-1)*nijc+(ic-1)*njc+jc+dom(ib)%faz(glevel+1)-nijkc
           ijkf=(kc-1)*nijf+(ic-1)*njf+jc+dom(ib)%faz(glevel)-nijkf

              dom(ib)%cof(ijkf)=dom(ib)%cof(ijkf)+dom(ib)%cof(ijkc)

              end do
           end do
        end do

        else

        do  kc=1,nkc-1
           k=2*kc-1
           do  ic=1,nic-1
              i=2*ic-1
              do jc=1,njc-1
                 j=2*jc-1

           ijkc=(kc-1)*nijc+(ic-1)*njc+jc+dom(ib)%faz(glevel+1)-nijkc
           ijkf=(k-1)*nijf+(i-1)*njf+j+dom(ib)%faz(glevel)-nijkf

                 dom(ib)%cof(ijkf)=dom(ib)%cof(ijkf)+
     & (27.0d0*dom(ib)%cof(ijkc)+
     & 9.0d0*dom(ib)%cof(ijkc+1)+9.0d0*dom(ib)%cof(ijkc+nijc)+
     & 3.0d0*dom(ib)%cof(ijkc+nijc+1)+9.0d0*dom(ib)%cof(ijkc+njc)+
     & 3.0d0*dom(ib)%cof(ijkc+njc+1)+3.0d0*dom(ib)%cof(ijkc+njc+nijc)+
     & 1.0d0*dom(ib)%cof(ijkc+1+njc+nijc))/64.0d0

                 dom(ib)%cof(ijkf+njf)=dom(ib)%cof(ijkf+njf)+
     & (9.0d0*dom(ib)%cof(ijkc)+
     & 3.0d0*dom(ib)%cof(ijkc+1)+3.0d0*dom(ib)%cof(ijkc+nijc)+
     & 1.0d0*dom(ib)%cof(ijkc+nijc+1)+27.0d0*dom(ib)%cof(ijkc+njc)+
     & 9.0d0*dom(ib)%cof(ijkc+njc+1)+9.0d0*dom(ib)%cof(ijkc+njc+nijc)+
     & 3.0d0*dom(ib)%cof(ijkc+1+njc+nijc))/64.0d0

                 dom(ib)%cof(ijkf+1)=dom(ib)%cof(ijkf+1)+
     & (9.0d0*dom(ib)%cof(ijkc)+
     & 27.0d0*dom(ib)%cof(ijkc+1)+3.0d0*dom(ib)%cof(ijkc+nijc)+
     & 9.0d0*dom(ib)%cof(ijkc+nijc+1)+3.0d0*dom(ib)%cof(ijkc+njc)+
     & 9.0d0*dom(ib)%cof(ijkc+njc+1)+1.0d0*dom(ib)%cof(ijkc+njc+nijc)+
     & 3.0d0*dom(ib)%cof(ijkc+1+njc+nijc))/64.0d0

                 dom(ib)%cof(ijkf+nijf)=dom(ib)%cof(ijkf+nijf)+
     & (9.0d0*dom(ib)%cof(ijkc)+
     & 3.0d0*dom(ib)%cof(ijkc+1)+27.0d0*dom(ib)%cof(ijkc+nijc)+
     & 9.0d0*dom(ib)%cof(ijkc+nijc+1)+3.0d0*dom(ib)%cof(ijkc+njc)+
     & 1.0d0*dom(ib)%cof(ijkc+njc+1)+9.0d0*dom(ib)%cof(ijkc+njc+nijc)+
     & 3.0d0*dom(ib)%cof(ijkc+1+njc+nijc))/64.0d0

                 dom(ib)%cof(ijkf+njf+1)=dom(ib)%cof(ijkf+njf+1)+
     & (3.0d0*dom(ib)%cof(ijkc)+
     & 9.0d0*dom(ib)%cof(ijkc+1)+1.0d0*dom(ib)%cof(ijkc+nijc)+
     & 3.0d0*dom(ib)%cof(ijkc+nijc+1)+9.0d0*dom(ib)%cof(ijkc+njc)+
     & 27.0d0*dom(ib)%cof(ijkc+njc+1)+3.0d0*dom(ib)%cof(ijkc+njc+nijc)+
     & 9.0d0*dom(ib)%cof(ijkc+1+njc+nijc))/64.0d0

                 dom(ib)%cof(ijkf+nijf+njf)=dom(ib)%cof(ijkf+nijf+njf)+
     & (3.0d0*dom(ib)%cof(ijkc)+
     & 1.0d0*dom(ib)%cof(ijkc+1)+9.0d0*dom(ib)%cof(ijkc+nijc)+
     & 3.0d0*dom(ib)%cof(ijkc+nijc+1)+9.0d0*dom(ib)%cof(ijkc+njc)+
     & 3.0d0*dom(ib)%cof(ijkc+njc+1)+27.0d0*dom(ib)%cof(ijkc+njc+nijc)+
     & 9.0d0*dom(ib)%cof(ijkc+1+njc+nijc))/64.0d0

                 dom(ib)%cof(ijkf+nijf+1)=dom(ib)%cof(ijkf+nijf+1)+
     & (3.0d0*dom(ib)%cof(ijkc)+
     & 9.0d0*dom(ib)%cof(ijkc+1)+9.0d0*dom(ib)%cof(ijkc+nijc)+
     & 27.0d0*dom(ib)%cof(ijkc+nijc+1)+1.0d0*dom(ib)%cof(ijkc+njc)+
     & 3.0d0*dom(ib)%cof(ijkc+njc+1)+3.0d0*dom(ib)%cof(ijkc+njc+nijc)+
     & 9.0d0*dom(ib)%cof(ijkc+1+njc+nijc))/64.0d0

                 dom(ib)%cof(ijkf+nijf+njf+1)=
     & dom(ib)%cof(ijkf+nijf+njf+1)+(1.0d0*dom(ib)%cof(ijkc)+
     & 3.0d0*dom(ib)%cof(ijkc+1)+3.0d0*dom(ib)%cof(ijkc+nijc)+
     & 9.0d0*dom(ib)%cof(ijkc+nijc+1)+3.0d0*dom(ib)%cof(ijkc+njc)+
     & 9.0d0*dom(ib)%cof(ijkc+njc+1)+9.0d0*dom(ib)%cof(ijkc+njc+nijc)+
     & 27.0d0*dom(ib)%cof(ijkc+1+njc+nijc))/64.0d0

              end do
           end do
        end do

        end if

        end do

c impose Neumann and Dirichlet boundary conditions
C TEMP: periodicity is not enforced to save one call to gxch1lin;
c check whether it has an impact or not...
        

        call mgbound(glevel)

        return
        end subroutine mgcorr
!##########################################################################
        subroutine mgbound(glevel)
!##########################################################################
        use vars
        use multidata
        implicit none
        integer i,j,k,ib,ni,nj,nk,nij,nijk,ijk
        integer glevel,gl


        do ib=1,nbp

        if(glevel.gt.dom(ib)%ngrid) then
           gl=dom(ib)%ngrid
        else
           gl=glevel
        end if

        ni=(dom(ib)%iep-dom(ib)%isp+1)/2**(gl-1)+2
        nj=(dom(ib)%jep-dom(ib)%jsp+1)/2**(gl-1)+2
        nk=(dom(ib)%kep-dom(ib)%ksp+1)/2**(gl-1)+2
        nij=ni*nj
        nijk=ni*nj*nk

!...................... west and east............................
        if (dom(ib)%iprev.lt.0 .and. dom(ib)%bc_west.ne.5) then
           do k=1,nk 
              do j=1,nj
                 i=1
                 ijk=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 dom(ib)%cof(ijk)=dom(ib)%cof(ijk+nj)
              end do
           end do
        end if
        if (dom(ib)%inext.lt.0 .and. dom(ib)%bc_east.ne.5) then
           do k=1,nk
              do j=1,nj
                 i=ni
                 ijk=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 dom(ib)%cof(ijk)=dom(ib)%cof(ijk-nj)
              end do
           end do
        end if
!.....................south and north.........................
        if (dom(ib)%jprev.lt.0 .and. dom(ib)%bc_south.ne.5) then
           do k=1,nk
              do i=1,ni
                 j=1
                 ijk=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 dom(ib)%cof(ijk)=dom(ib)%cof(ijk+1)
              end do
           end do
        end if
        if (dom(ib)%jnext.lt.0 .and. dom(ib)%bc_north.ne.5) then
           do k=1,nk
              do i=1,ni
                 j=nj
                 ijk=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 dom(ib)%cof(ijk)=dom(ib)%cof(ijk-1)
              end do
           end do
        end if
!.....................bottom and top.........................
        if (dom(ib)%kprev.lt.0 .and. dom(ib)%bc_bottom.ne.5) then
           do j=1,nj
              do i=1,ni
                 k=1
                 ijk=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 dom(ib)%cof(ijk)=dom(ib)%cof(ijk+nij)
              end do
           end do
        end if
        if (dom(ib)%knext.lt.0 .and. dom(ib)%bc_top.ne.5) then
           do j=1,nj
              do i=1,ni
                 k=nk
                 ijk=(k-1)*nij+(i-1)*nj+j+dom(ib)%faz(glevel)-nijk
                 dom(ib)%cof(ijk)=dom(ib)%cof(ijk-nij)
              end do
           end do
        end if

        end do


        return
        end subroutine mgbound
!##########################################################################
