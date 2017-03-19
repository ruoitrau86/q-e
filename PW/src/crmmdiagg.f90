!
! Copyright (C) 2017 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0._DP, 0._DP )
#define ONE  ( 1._DP, 0._DP )
!
!----------------------------------------------------------------------------
SUBROUTINE crmmdiagg( npwx, npw, nbnd, npol, psi, e, btype, precondition, &
                      ethr, ndiis, uspp, reorder, notconv, rmm_iter )
  !----------------------------------------------------------------------------
  !
  ! ... Iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned RMM-DIIS algorithm.
  !
  USE constants,        ONLY : pi
  USE kinds,            ONLY : DP
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npwx, npw, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
  REAL(DP),    INTENT(INOUT) :: e(nbnd)
  INTEGER,     INTENT(IN)    :: btype(nbnd)
  REAL(DP),    INTENT(IN)    :: precondition(npwx*npol)
  REAL(DP),    INTENT(IN)    :: ethr
  INTEGER,     INTENT(IN)    :: ndiis
  LOGICAL,     INTENT(IN)    :: uspp
  LOGICAL,     INTENT(IN)    :: reorder
  INTEGER,     INTENT(OUT)   :: notconv
  INTEGER,     INTENT(OUT)   :: rmm_iter
  !
  ! ... local variables
  !
  INTEGER                  :: ierr
  INTEGER                  :: idiis, ibnd
  REAL(DP)                 :: empty_ethr
  COMPLEX(DP), ALLOCATABLE :: phi(:,:,:), dphi(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dpsi(:,:), hpsi(:,:), spsi(:,:)
  !
  CALL start_clock( 'crmmdiagg' )
  !
  empty_ethr = MAX( ( ethr * 5._DP ), 1.E-5_DP )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx * npol
     kdmx = npwx * npol
     !
  END IF
  !
  ALLOCATE( phi( kdmx, nbnd, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ',' cannot allocate phi ', ABS(ierr) )
  !
  ALLOCATE( dphi( kdmx, nbnd, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ',' cannot allocate dphi ', ABS(ierr) )
  !
  ALLOCATE( dpsi( kdmx, nbnd ) )
  ALLOCATE( hpsi( kdmx, nbnd ) )
  IF ( uspp ) ALLOCATE( spsi( kdmx, nbnd ) )
  !
  rmm_iter = 0
  notconv  = 0
  !
  DO idiis = 1, ndiis
     !
     rmm_iter = rmm_iter + 1
     !
     ! ... Operate the Hamiltonian : H |psi>
     !
     CALL h_psi( npwx, npw, nbnd, psi, hpsi )
     !
     ! ... Operate the Overlap : S |psi>
     !
     IF ( uspp ) CALL s_psi( npwx, npw, nbnd, psi, spsi )
     !
     ! ... Residual vector : |R> = (H - e S) |psi>
     !
     IF ( uspp ) THEN
        !
        FORALL ( ibnd = 1:nbnd ) &
        dpsi(:,ibnd) = hpsi(:,ibnd) - e(ibnd) * spsi(:,ibnd)
        !
     ELSE
        !
        FORALL ( ibnd = 1:nbnd ) &
        dpsi(:,ibnd) = hpsi(:,ibnd) - e(ibnd) * psi(:,ibnd)
        !
     END IF
     !
     ! ... Save current wave functions and residual vectors
     !
     phi (:,:,idiis) = psi (:,:)
     dphi(:,:,idiis) = dpsi(:,:)
     !
     ! ... Perform DIIS
     !
     IF ( idiis > 1 ) THEN
        !
        ! TODO
        ! TODO phi, dphi -> psi, dpsi
        ! TODO
        !
     END IF
     !
     ! ... Line searching
     !
     DO ibnd = 1, nbnd
        !
        ! TODO
        ! TODO |psi> = |psi> + s * precondition |dpsi>
        ! TODO
        !
     END DO
     !
     ! ... Eigenvalues
     !
     ! TODO
     ! TODO
     ! TODO
     !
     ! ... Check convergence
     !
     ! TODO
     ! TODO
     ! TODO
     !
  END DO
  !
  ! ... Sort eigenvalues and eigenvectors
  !
  IF ( reorder ) THEN
     !
     ! TODO
     ! TODO
     ! TODO
     !
  END IF
  !
  DEALLOCATE( phi )
  DEALLOCATE( dphi )
  DEALLOCATE( dpsi )
  DEALLOCATE( hpsi )
  IF ( uspp ) DEALLOCATE( spsi )
  !
  CALL stop_clock( 'crmmdiagg' )
  !
  RETURN
  !
END SUBROUTINE crmmdiagg
