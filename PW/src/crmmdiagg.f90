!
! Copyright (C) 2017 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0._DP, 0._DP )
!
!----------------------------------------------------------------------------
SUBROUTINE crmmdiagg( npwx, npw, nbnd, npol, psi, e, btype, precondition, &
                      ethr, ndiis, uspp, reorder, notconv, rmm_iter )
  !----------------------------------------------------------------------------
  !
  ! ... Iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned RMM-DIIS algorithm.
  !
  USE constants, ONLY : eps14
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum, mp_bcast
  USE mp_bands,  ONLY : intra_bgrp_comm
    CALL mp_bcast( conv, root_bgrp_id, inter_bgrp_comm)
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
  REAL(DP),    ALLOCATABLE :: ew(:), hw(:), sw(:)
  LOGICAL,     ALLOCATABLE :: conv(:)
  REAL(DP),    ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
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
  ALLOCATE( ew( nbnd ) )
  ALLOCATE( hw( nbnd ) )
  ALLOCATE( sw( nbnd ) )
  ALLOCATE( conv( nbnd ) )
  !
  ALLOCATE( hc( ndiis, ndiis ) )
  ALLOCATE( sc( ndiis, ndiis ) )
  ALLOCATE( vc( ndiis, ndiis ) )
  !
  phi  = ZERO
  dphi = ZERO
  !
  dpsi = ZERO
  hpsi = ZERO
  IF ( uspp ) spsi = ZERO
  !
  ew   = 0.0_DP
  hw   = 0.0_DP
  sw   = 0.0_DP
  conv = .FALSE.
  !
  hc   = ZERO
  sc   = ZERO
  vc   = ZERO
  !
  ! ... Initial eigenvalues
  !
  CALL eigenvalues( )
  !
  ! ... RMM-DIIS's loop
  !
  rmm_iter = 0
  !
  DO idiis = 1, ndiis
     !
     rmm_iter = rmm_iter + 1
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
     ! ... Calculate eigenvalues and check convergence
     !
     CALL eigenvalues( )
     !
     IF ( notconv == 0 ) EXIT
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
  DEALLOCATE( ew )
  DEALLOCATE( hw )
  DEALLOCATE( sw )
  DEALLOCATE( conv )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  DEALLOCATE( vc )
  !
  CALL stop_clock( 'crmmdiagg' )
  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE eigenvalues( )
    !
    IMPLICIT NONE
    !
    INTEGER            :: ibnd
    REAL(DP)           :: h, s
    REAL(DP), EXTERNAL :: ZDOTC
    !
    ! ... Operate the Hamiltonian : H |psi>
    !
    CALL h_psi( npwx, npw, nbnd, psi, hpsi )
    !
    ! ... Operate the Overlap : S |psi>
    !
    IF ( uspp ) CALL s_psi( npwx, npw, nbnd, psi, spsi )
    !
    ! ... Eigenvalues
    !
    FORALL ( ibnd = 1:nbnd ) &
    hw(ibnd) = ZDOTC( kdim, psi(1,ibnd), 1, hpsi(1,ibnd), 1 )
    !
    CALL mp_sum( hw, intra_bgrp_comm )
    !
    IF ( uspp ) THEN
       !
       FORALL ( ibnd = 1:nbnd ) &
       sw(ibnd) = ZDOTC( kdim, psi(1,ibnd), 1, spsi(1,ibnd), 1 )
       !
    ELSE
       !
       FORALL ( ibnd = 1:nbnd ) &
       sw(ibnd) = ZDOTC( kdim, psi(1,ibnd), 1, psi(1,ibnd), 1 )
       !
    END IF
    !
    CALL mp_sum( sw, intra_bgrp_comm )
    !
    ew(:) = hw(:) / sw(:)
    !
    ! ... Check convergence
    !
    WHERE( btype(:) == 1 )
       !
       conv(:) = ( ( ABS( ew(:) - e(:) ) < ethr ) )
       !
    ELSEWHERE
       !
       conv(:) = ( ( ABS( ew(:) - e(:) ) < empty_ethr ) )
       !
    END WHERE
    !
    CALL mp_bcast( conv, root_bgrp_id, inter_bgrp_comm)
    !
    notconv = COUNT( .NOT. conv(:) )
    !
    e(:) = ew(:)
    !
    RETURN
    !
  END SUBROUTINE eigenvalues
  !
  !
END SUBROUTINE crmmdiagg
