MODULE ITMM_module

  ! User derived-type modules

  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE cross_section_types
  USE geometry_types
  USE angle_types
  USE multindex_types
  USE global_variables
  USE SDD_global_variables
  USE mpi

  ! Use modules that pertain setting up problem

  USE termination_module
  USE dump_inguess_module

	IMPLICIT NONE

CONTAINS

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

	SUBROUTINE ITMM_factor

		INTEGER(kind=li) :: ind,indprime,eg,INFO

		ALLOCATE(IPVT(num_cells,egmax),IJ(num_cells,num_cells,egmax),IPVTtemp(num_cells),&
		 & Jpsitemp(4*nangle*N_side_SDbound,num_cells))
		ALLOCATE(IJtemp(num_cells,num_cells),Jphitemp(num_cells,num_cells),&
		& Kphitemp(num_cells,4*nangle*N_side_SDbound))

		!Loop over energy groups
		DO eg=1,egmax

			!Calculate I-Jphi
			DO ind=1,num_cells
				DO indprime=1,num_cells
					IJtemp(ind,indprime)=Jphi(ind,indprime,eg)
				END DO
			END DO
			DO ind=1,num_cells
				DO indprime=1,num_cells
					IF (indprime.EQ.ind) THEN
						IJtemp(ind,indprime)=1.0d0-IJtemp(ind,indprime)
					ELSE
						IJtemp(ind,indprime)=-1.0d0*IJtemp(ind,indprime)
					END IF
				END DO
			END DO
			!Factor I-Jphi
			INFO=0
			CALL dgetrf(num_cells,num_cells,IJtemp,num_cells,IPVTtemp,INFO)
   			DO ind=1,num_cells
   				IPVT(ind,eg)=IPVTtemp(ind)
       	 		DO indprime=1,num_cells

        			IJ(ind,indprime,eg)=IJtemp(ind,indprime)

        		END DO
        	END DO

        !End loop over groups
		END DO

	END SUBROUTINE ITMM_factor

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

	SUBROUTINE ITMMsolve(eg,sc_flux,q_external)

		INTEGER,INTENT(IN)::eg
		INTEGER::i,j,INFO,mpi_err
		REAL(kind=d_t),INTENT(OUT)::sc_flux(num_moments_v,namom,num_cells)
		REAL(kind=d_t),INTENT(IN)::q_external(num_moments_v,namom,num_cells)
		REAL(kind=d_t)::psiintemp(4*nangle*N_side_SDbound,1)&
		& ,sq(num_cells,1),psiouttemp(4*nangle*N_side_SDbound,1),ITMM_RHS(num_cells,1)
		

		!First, calculate cell-averaged scalar fluxes

		!Calculate RHS vector
		IF (inner.EQ.1) THEN
		    DO i=1,num_cells
			    DO j=1,num_cells
				    Jphitemp(i,j)=Jphi(i,j,eg)
			    END DO
		    END DO
	    END IF

		DO i=1,num_cells
			sq(i,1)=1.0d0/(sigma_scat(reg2mat(cells(i)%reg),1,eg,eg)%xs*&
            dens_fact(cells(i)%reg))
            sq(i,1)=sq(i,1)*q_external(1,1,i)
        END DO

        ITMM_RHS=MATMUL(Jphitemp,sq)

        DO j=1,4*nangle*N_side_SDbound
        	psiintemp(j,1)=psiin(j,eg)
        	IF (inner.EQ.1) THEN
           		DO i=1,num_cells
            		Kphitemp(i,j)=Kphi(i,j,eg)
            	END DO
        	END IF
        END DO

        ITMM_RHS=ITMM_RHS+MATMUL(Kphitemp,psiintemp)

        !Solve for scalar flux
        IF (inner.EQ.1) THEN
            DO i=1,num_cells
            	IPVTtemp(i)=IPVT(i,eg)
            	DO j=1,num_cells
            		IJtemp(i,j)=IJ(i,j,eg)
            	END DO
            END DO
        END IF
        INFO=0
        CALL dgetrs('N',num_cells,1,IJtemp,num_cells,IPVTtemp,ITMM_RHS,num_cells,INFO)
        DO i=1,num_cells
        	sc_flux(1,1,i)=ITMM_RHS(i,1)
        END DO

        !Then, calculate Subdomain edge outgoing angular fluxes

        !Calculate phi+sigma_s^-1*q
        sq(:,1)=ITMM_RHS(:,1)+sq(:,1)

        !Multiply by Jpsi
        IF (inner.EQ.1) THEN
            DO i=1,4*nangle*N_side_SDbound
            	DO j=1,num_cells
            		Jpsitemp(i,j)=Jpsi(i,j,eg)
            	END DO
            END DO
        END IF
        psiouttemp=MATMUL(Jpsitemp,sq)

        !Add in nonzero values of Kpsi*psiin
        DO i=1,nonzero
        	psiouttemp(KpsiIndexes(1,i),1)=psiouttemp(KpsiIndexes(1,i),1)&
        	+KpsiElements(i,eg)*psiintemp(KpsiIndexes(2,i),1)
        END DO

        DO i=1,4*nangle*N_side_SDbound
        	psiout(i,eg)=psiouttemp(i,1)
        END DO


	END SUBROUTINE ITMMsolve

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

    SUBROUTINE ITMM_communicate_instructions

        !This subroutine performs the initial communication between processors for PBJ
        !methods. This is a one time communication that sends the information between
        !processors, letting the proccessors known exactly what values they will be
        !receiving in what order each iteration.

        INTEGER::unique_neighs_temp(10000),i,j,temp_proc,k,new,f,octant,q,ierr=0,m
        INTEGER::mpi_status(MPI_STATUS_SIZE)=0,temp,k_dest,f_dest,refl,tpe,mate,p
        INTEGER,ALLOCATABLE :: ii(:),comm_instructions_cell_send(:,:),comm_instructions_face_send(:,:)
        INTEGER,ALLOCATABLE::comm_instructions_octant_send(:,:),comm_instructions_angle_send(:,:)

        !Determine how many neighboring processors there are
        num_neigh=1
        unique_neighs_temp=PBJrank+1
        DO i=1,num_cells

            DO j=0,3

                temp_proc=SDD_adjacency_list(SDD_cells_l2g_G(i),j)%proc
                new=1

                DO k=1,num_neigh

                    IF ((temp_proc.EQ.0).OR.(unique_neighs_temp(k).EQ.temp_proc)) new=0


                END DO

                IF (new.EQ.1) THEN

                    num_neigh=num_neigh+1
                    unique_neighs_temp(num_neigh)=temp_proc


                END IF

            END DO

        END DO

        !Store the array of neighboring processors to permanent array
        !Note: the last entry will be the total problem boundary
        ALLOCATE(unique_neighs(num_neigh))
        DO i=1,num_neigh
            unique_neighs(i)=unique_neighs_temp(i+1)
        END DO
        unique_neighs(num_neigh)=0

        !Determine how many angular fluxes will be send to each processor
        ALLOCATE(Num_ang_tosend(num_neigh))
        Num_ang_tosend=0
        DO i=1,SDD_side_cells
            k=SDDb_cells(i)%cell
            f=SDDb_cells(i)%face
            k=SDD_cells_l2g_G(k)
            temp_proc=SDD_adjacency_list(k,f)%proc
            DO j=1,num_neigh
                IF (unique_neighs(j).EQ.temp_proc) Num_ang_tosend(j)=Num_ang_tosend(j)+1
            END DO
        END DO
        Num_ang_tosend=Num_ang_tosend*4*nangle

        !Determine the maximum number of angular fluxes that will be passed
        max_ang_comm=MAXVAL(Num_ang_tosend)
        ALLOCATE(comm_instructions(4*nangle*N_side_SDbound,4,num_neigh))
        ALLOCATE(comm_instructions_cell_send(4*nangle*N_side_SDbound,num_neigh-1))
        ALLOCATE(comm_instructions_face_send(4*nangle*N_side_SDbound,num_neigh-1))
        ALLOCATE(comm_instructions_octant_send(4*nangle*N_side_SDbound,num_neigh-1))
        ALLOCATE(comm_instructions_angle_send(4*nangle*N_side_SDbound,num_neigh-1))
        ALLOCATE(pack_instructions(4*nangle*N_side_SDbound,2))
        ALLOCATE(unpack_instructions(4*nangle*N_side_SDbound,num_neigh-1))
        ALLOCATE(ii(num_neigh))
        IF (rside_cells.GE.1) THEN
            ALLOCATE(SDD_refl_BC_instructions(4*nangle*rside_cells,2))
            SDD_refl_BC_instructions=-99
        END IF
        comm_instructions=-99
        pack_instructions=-99
        unpack_instructions=-99
        ii=0

        !Loop over ITMMKindex and load comm/pack instructions
        DO i=1,4*nangle*N_side_SDbound

            !Get the psi indexes and destination proc
            k=ITMMKindex(i,1,1)
            f=ITMMKindex(i,2,1)
            octant=ITMMKindex(i,3,1)
            q=ITMMKindex(i,4,1)
            IF (k.EQ.-99) CYCLE
            k=SDD_cells_l2g_G(k)
            temp_proc=SDD_adjacency_list(k,f)%proc

            !k and f indicate the global cell/face numbers for the current proc
            !we need to change this to be the global cell for the destination proc
            k_dest=SDD_adjacency_list(k,f)%cell
            f_dest=SDD_adjacency_list(k,f)%face

            !Write the instructions
            DO j=1,num_neigh
                IF (unique_neighs(j).EQ.temp_proc) EXIT
            END DO

            ii(j)=ii(j)+1
            comm_instructions(ii(j),1,j)=k_dest
            comm_instructions(ii(j),2,j)=f_dest
            comm_instructions(ii(j),3,j)=octant
            comm_instructions(ii(j),4,j)=q
            pack_instructions(i,1)=j
            pack_instructions(i,2)=ii(j)

        END DO


        !Send/Recv communication instructions so that all procs will know what to do with the
        !angular fluxes they receive in between iterations
        ALLOCATE(SDD_send_handles_G(num_neigh-1, 4))
        ALLOCATE(SDD_recv_handles_G(num_neigh-1, 4))

        SDD_send_handles_G = 0
        SDD_recv_handles_G = 0
        comm_instructions_cell_send=0
        comm_instructions_face_send=0
        comm_instructions_octant_send=0
        comm_instructions_angle_send=0

        DO i = 1, num_neigh-1
            comm_instructions_cell_send(:,i)=comm_instructions(:,1,i)
            comm_instructions_face_send(:,i)=comm_instructions(:,2,i)
            comm_instructions_octant_send(:,i)=comm_instructions(:,3,i)
            comm_instructions_angle_send(:,i)=comm_instructions(:,4,i)
        END DO

        DO i = 1, num_neigh-1

            CALL MPI_IRECV( comm_instructions(:Num_ang_tosend(i),1,i), Num_ang_tosend(i), MPI_INTEGER, &
                            unique_neighs(i)-1, 1, MPI_COMM_WORLD, SDD_recv_handles_G(i,1), ierr)
            CALL MPI_IRECV( comm_instructions(:Num_ang_tosend(i),2,i), Num_ang_tosend(i), MPI_INTEGER, &
                            unique_neighs(i)-1, 2, MPI_COMM_WORLD, SDD_recv_handles_G(i,2), ierr)
            CALL MPI_IRECV( comm_instructions(:Num_ang_tosend(i),3,i), Num_ang_tosend(i), MPI_INTEGER, &
                            unique_neighs(i)-1, 3, MPI_COMM_WORLD, SDD_recv_handles_G(i,3), ierr)
            CALL MPI_IRECV( comm_instructions(:Num_ang_tosend(i),4,i), Num_ang_tosend(i), MPI_INTEGER, &
                            unique_neighs(i)-1, 4, MPI_COMM_WORLD, SDD_recv_handles_G(i,4), ierr)


            CALL MPI_ISEND( comm_instructions_cell_send   (:Num_ang_tosend(i),i)  , Num_ang_tosend(i), MPI_INTEGER, &
                            unique_neighs(i)-1, 1, MPI_COMM_WORLD, SDD_send_handles_G(i,1), ierr)
            CALL MPI_ISEND( comm_instructions_face_send   (:Num_ang_tosend(i),i)  , Num_ang_tosend(i), MPI_INTEGER, &
                            unique_neighs(i)-1, 2, MPI_COMM_WORLD, SDD_send_handles_G(i,2), ierr)
            CALL MPI_ISEND( comm_instructions_octant_send (:Num_ang_tosend(i),i)  , Num_ang_tosend(i), MPI_INTEGER, &
                            unique_neighs(i)-1, 3, MPI_COMM_WORLD, SDD_send_handles_G(i,3), ierr)
            CALL MPI_ISEND( comm_instructions_angle_send  (:Num_ang_tosend(i),i)  , Num_ang_tosend(i), MPI_INTEGER, &
                            unique_neighs(i)-1, 4, MPI_COMM_WORLD, SDD_send_handles_G(i,4), ierr)

        END DO
        DO i = 1, 4
            CALL MPI_WAITALL((num_neigh-1), SDD_send_handles_G(:,i), MPI_STATUSES_IGNORE, ierr)
        END DO
        DO i = 1, 4
            CALL MPI_WAITALL((num_neigh-1), SDD_recv_handles_G(:,i), MPI_STATUSES_IGNORE, ierr)
        END DO

        !Reallocate the handles
        DEALLOCATE(SDD_send_handles_G,SDD_recv_handles_G)
        ALLOCATE(SDD_send_handles_G(num_neigh-1,1),SDD_recv_handles_G(num_neigh-1,1))

        !comm_instructions now holds the global cell, face, octant, and angle associated with
        !the values that will be received from procs each iteration

        !Use the new comm_instructions obtained from adjacent processors to make unpack instructions
        DO i=1,num_neigh-1

            DO j=1,Num_ang_tosend(i)

                !From comm_instructions, what is this flux
                k=comm_instructions(j,1,i)
                f=comm_instructions(j,2,i)
                octant=comm_instructions(j,3,i)
                q=comm_instructions(j,4,i)
                IF (k.EQ.-99) CYCLE
                k=SDD_cells_g2l_G(k)

                !Search using ITMMKindex to find where this flux fits in psiin
                DO m=1,4*nangle*N_side_SDbound

                    IF ((k     .EQ.ITMMKindex(m,1,2))  .AND. &
                    &  (f      .EQ.ITMMKindex(m,2,2))  .AND. &
                    &  (octant .EQ.ITMMKindex(m,3,2))  .AND. &
                    &  (q      .EQ.ITMMKindex(m,4,2))) EXIT

                END DO

                !Set unpack instructions
                unpack_instructions(j,i)=m

            END DO

        END DO

        !pack/unpack_instructions how tell a processor how to divide its psiout vector into
        !multiple, smaller vectors that will be sent to adjacent processors, and how to reconstruct
        !their own psiin vector from those received from neighbors

        !We still need one more instruction, though, pertaining to reflective BCs
        !This requires no communication, and hence, no pack instructions
        !We will make special unpack instructions for reflective BCs that will point
        !to specific locations within the current proc's psiout vector and tell where
        !to place them within psiin
        IF (rside_cells.GE.1) THEN

            p=0

            !Loop over psiout vector
            DO i=1,4*nangle*N_side_SDbound

                !From ITMMKindex, what is this flux
                k=ITMMKindex(i,1,1)
                f=ITMMKindex(i,2,1)
                octant=ITMMKindex(i,3,1)
                q=ITMMKindex(i,4,1)
                IF (k.EQ.-99) CYCLE

                !Is this a reflective BC
                refl=0
                DO j=1,rside_cells
                    IF ((k.EQ.rb_cells(j)%cell).AND.(f.EQ.rb_cells(j)%face)) THEN
                        refl=1
                        p=p+1
                    END IF
                    IF (refl.EQ.1) EXIT
                END DO

                !If this is a reflective BC, set instructions for where it goes in psiin
                IF (refl.EQ.1) THEN

                    !Determine which axis this is reflective about to get mate octant
                    tpe=refl_face_tpe(j)
                    IF      (tpe .EQ. 1_li .OR. tpe .EQ. -1_li) THEN
                        mate = mu_mate(octant)
                    ELSE IF (tpe .EQ. 2_li .OR. tpe .EQ. -2_li) THEN
                        mate = eta_mate(octant)
                    ELSE IF (tpe .EQ. 3_li .OR. tpe .EQ. -3_li) THEN
                        mate = xi_mate(octant)
                    END IF

                    !Find the reflective mate's location in psiin and write instructions
                    refl=0
                    DO m=1,4*nangle*N_side_SDbound
                        IF ((ITMMKindex(m,1,2).EQ.k)   .AND. &
                        & (ITMMKindex(m,2,2).EQ.f)    .AND. &
                        & (ITMMKindex(m,3,2).EQ.mate) .AND. &
                        & (ITMMKindex(m,4,2).EQ.q))   THEN
                            refl=1
                        END IF
                        IF (refl.EQ.1) THEN
                            SDD_refl_BC_instructions(p,1)=i
                            SDD_refl_BC_instructions(p,2)=m
                            EXIT
                        END IF
                    END DO

                END IF

            END DO

        END IF


        !Allocate the communication psi array
        ALLOCATE(psi_comm(max_ang_comm,num_neigh-1),psi_comm_send(max_ang_comm,num_neigh-1))

    END SUBROUTINE ITMM_communicate_instructions

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

    SUBROUTINE ITMM_pack(eg)

        !This subroutine packs the psi_out vector into multiple smaller vectors
        !, one for each neighboring processor

        INTEGER::i
        INTEGER,INTENT(IN)::eg

        DO i=1,4*nangle*N_side_SDbound

            IF (pack_instructions(i,1).NE.num_neigh.AND.pack_instructions(i,2).NE.-99.AND.pack_instructions(i,1).NE.-99) &
                & psi_comm(pack_instructions(i,2),pack_instructions(i,1))=psiout(i,eg)

        END DO

    END SUBROUTINE ITMM_pack

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

    SUBROUTINE ITMM_comm(eg)

        !This subroutine sends and receives angular fluxes between processors

        INTEGER::i,ierr
        INTEGER,INTENT(IN)::eg

        psi_comm_send=psi_comm
        SDD_send_handles_G=0
        SDD_recv_handles_G=0

        DO i=1,num_neigh-1

            CALL MPI_IRECV( psi_comm(:Num_ang_tosend(i),i), Num_ang_tosend(i), MPI_DOUBLE_PRECISION, &
                            unique_neighs(i)-1, 0, MPI_COMM_WORLD, SDD_recv_handles_G(i,1), ierr)

            CALL MPI_ISEND( psi_comm_send(:Num_ang_tosend(i),i), Num_ang_tosend(i), MPI_DOUBLE_PRECISION, &
                            unique_neighs(i)-1, 0, MPI_COMM_WORLD, SDD_send_handles_G(i,1), ierr)

        END DO

        CALL MPI_WAITALL((num_neigh-1), SDD_send_handles_G(:,1), MPI_STATUSES_IGNORE, ierr)
        CALL MPI_WAITALL((num_neigh-1), SDD_recv_handles_G(:,1), MPI_STATUSES_IGNORE, ierr)

    END SUBROUTINE ITMM_comm

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

    SUBROUTINE ITMM_unpack(eg)

        !This subroutine unpacks the psiout messages received to create psiin

        INTEGER,INTENT(IN)::eg
        INTEGER::i,j,ierr

        !Zero psiin to make everything start as vacuum
        psiin=0.0d0

        !Unpack messages from neighboring procs
        DO i=1,num_neigh-1

            DO j=1,Num_ang_tosend(i)

                IF (unpack_instructions(j,i).NE.-99) psiin(unpack_instructions(j,i),eg)=psi_comm(j,i)

            END DO

        END DO

        !If there are reflective BCs, obtain fluxes from this proc's psiout
        IF (rside_cells.GE.1) THEN

            DO i=1,4*nangle*rside_cells
                IF (SDD_refl_BC_instructions(i,1).NE.-99) &
                 & psiin(SDD_refl_BC_instructions(i,2),eg)=psiout(SDD_refl_BC_instructions(i,1),eg)

            END DO

        END IF

    END SUBROUTINE ITMM_unpack

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

END MODULE ITMM_module
