!
! Copyright (C) 2017 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
!
!--------------------------------------------------------------------------
SUBROUTINE gram_schmidt_gamma( npwx, npw, nbnd, psi, uspp, gstart, nbsize )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for Gamma-only calculations.
  ! ... blocking algorithm is used.
  !
  USE constants, ONLY : eps16
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum, mp_max, mp_bcast
  USE mp_bands,  ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id, set_bgrp_indices
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx,nbnd)
  LOGICAL,     INTENT(IN)    :: uspp
  INTEGER,     INTENT(IN)    :: gstart
  INTEGER,     INTENT(IN)    :: nbsize
  !
  ! ... local variables
  !
  INTEGER                  :: iblock, nblock
  INTEGER                  :: iblock_start, iblock_end
  INTEGER                  :: jblock_start, jblock_end
  INTEGER                  :: ibnd_start, ibnd_end
  INTEGER                  :: jbnd_start, jbnd_end
  COMPLEX(DP), ALLOCATABLE :: phi(:,:), spsi(:,:), sphi(:,:)
  INTEGER,     ALLOCATABLE :: owner_bgrp_id(:)
  !
  nblock = nbnd / nbsize
  IF ( MOD( nbnd, nbsize ) /= 0 ) nblock = nblock + 1
  !
  CALL set_bgrp_indices( nblock, iblock_start, iblock_end )
  !
  IF ( my_bgrp_id >= nblock ) THEN
     !
     iblock_start = nblock + 1
     iblock_end   = nblock
     !
  END IF
  !
  ALLOCATE( phi ( npwx, nbnd ) )
  IF ( uspp ) ALLOCATE( spsi( npwx, nbnd ) )
  IF ( uspp ) ALLOCATE( sphi( npwx, nbnd ) )
  ALLOCATE( owner_bgrp_id( nblock ) )
  !
  phi = ZERO
  IF ( uspp ) spsi = ZERO
  IF ( uspp ) sphi = ZERO
  !
  ! ... Set owers of blocks
  !
  owner_bgrp_id = 0
  !
  DO iblock = 1, nblock
     !
     IF ( iblock_start <= iblock .AND. iblock <= iblock_end ) &
     owner_bgrp_id(iblock) = my_bgrp_id
     !
  END DO
  !
  CALL mp_max( owner_bgrp_id, inter_bgrp_comm )
  !
  ! ... Set Im[ psi(G=0) ] - needed for numerical stability
  !
  IF ( gstart == 2 ) psi(1,1:nbnd) = CMPLX( DBLE( psi(1,1:nbnd) ), 0._DP, kind=DP )
  !
  ! ... Operate the overlap : S |psi_j>
  !
  IF ( uspp ) CALL s_psi( npwx, npw, nbnd, psi, spsi )
  !
  ! ... Set initial : |phi_j> = |psi_j>
  !
  phi = psi
  !
  ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
  IF ( gstart == 2 ) phi(1,1:nbnd) = CMPLX( DBLE( phi(1,1:nbnd) ), 0._DP, kind=DP )
  !
  IF ( uspp ) THEN
     !
     sphi = spsi
     !
     ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
     IF ( gstart == 2 ) sphi(1,1:nbnd) = CMPLX( DBLE( sphi(1,1:nbnd) ), 0._DP, kind=DP )
     !
  END IF
  !
  ! ... Blocking loop
  !
  DO iblock = 1, nblock
     !
     ! ... Orthogonalize diagonal block by standard Gram-Schmidt
     !
     ibnd_start = ( iblock - 1 ) * nbsize + 1
     ibnd_end   = MIN( iblock * nbsize, nbnd )
     !
     IF ( owner_bgrp_id(iblock) == my_bgrp_id ) &
     CALL gram_schmidt_diag( ibnd_start, ibnd_end )
     !
     ! ... Bcast diagonal block
     !
     CALL mp_bcast( phi(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     IF ( uspp ) &
     CALL mp_bcast( sphi(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     ! ... Project off-diagonal block outside of diagonal block
     !
     jblock_start = MAX( iblock_start, iblock + 1 )
     jblock_end   = iblock_end
     !
     jbnd_start = ( jblock_start - 1 ) * nbsize + 1
     jbnd_end   = MIN( jblock_end * nbsize, nbnd )
     !
     IF ( jblock_start <= jblock_end .AND. jbnd_start <= jbnd_end ) &
     CALL project_offdiag( ibnd_start, ibnd_end, jbnd_start, jbnd_end )
     !
  END DO
  !
  ! ... Copy psi <- phi
  !
  psi = phi
  !
  DEALLOCATE( phi )
  IF ( uspp ) DEALLOCATE( spsi )
  IF ( uspp ) DEALLOCATE( sphi )
  DEALLOCATE( owner_bgrp_id )
  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE gram_schmidt_diag( ibnd_start, ibnd_end )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: ibnd_start, ibnd_end
    !
    INTEGER               :: ibnd
    REAL(DP), ALLOCATABLE :: sr(:)
    REAL(DP)              :: norm
    REAL(DP), EXTERNAL    :: DDOT
    !
    ALLOCATE( sr( ibnd_start:ibnd_end ) )
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( ibnd > ibnd_start ) THEN
          !
          ! ... <phi_j| S |psi_i>
          !
          IF ( uspp ) THEN
             !
             CALL DGEMV( 'T', 2 * npw, ibnd - ibnd_start, 2._DP, phi(1,ibnd_start), 2 * npwx, &
                         spsi(1,ibnd), 1, 0._DP, sr(ibnd_start), 1 )
             !
             IF ( gstart == 2 ) &
             CALL DAXPY( ibnd - ibnd_start, -spsi(1,ibnd), phi(1,ibnd_start), 2 * npwx, &
                         sr(ibnd_start), 1 )
             !
          ELSE
             !
             CALL DGEMV( 'T', 2 * npw, ibnd - ibnd_start, 2._DP, phi(1,ibnd_start), 2 * npwx, &
                         psi(1,ibnd), 1, 0._DP, sr(ibnd_start), 1 )
             !
             IF ( gstart == 2 ) &
             CALL DAXPY( ibnd - ibnd_start, -psi(1,ibnd), phi(1,ibnd_start), 2 * npwx, &
                         sr(ibnd_start), 1 )
             !
          END IF
          !
          CALL mp_sum( sr, intra_bgrp_comm )
          !
          ! ... phi_i = phi_i - phi_j * <phi_j| S |psi_i>
          !
          CALL DGEMV( 'N', 2 * npw, ibnd - ibnd_start, -1._DP, phi(1,ibnd_start), 2 * npwx, &
                      sr(ibnd_start), 1, 1._DP, phi(1,ibnd), 1 )
          !
          ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
          IF ( gstart == 2 ) phi(1,ibnd) = CMPLX( DBLE( phi(1,ibnd) ), 0._DP, kind=DP )
          !
          IF ( uspp ) THEN
             !
             CALL DGEMV( 'N', 2 * npw, ibnd - ibnd_start, -1._DP, sphi(1,ibnd_start), 2 * npwx, &
                         sr(ibnd_start), 1, 1._DP, sphi(1,ibnd), 1 )
             !
             ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
             IF ( gstart == 2 ) sphi(1,ibnd) = CMPLX( DBLE( sphi(1,ibnd) ), 0._DP, kind=DP )
             !
          END IF
          !
       END IF
       !
       ! ... Normalize : phi_i = phi_i / SQRT(<phi_i| S |phi_i>)
       !
       IF ( uspp ) THEN
          !
          norm = 2.0_DP * DDOT( 2 * npw, phi(1,ibnd), 1, sphi(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) norm = norm - DBLE( phi(1,ibnd) ) * DBLE ( sphi(1,ibnd) )
          !
       ELSE
          !
          norm = 2.0_DP * DDOT( 2 * npw, phi(1,ibnd), 1, phi(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) norm = norm - DBLE( phi(1,ibnd) ) * DBLE ( phi(1,ibnd) )
          !
       END IF
       !
       CALL mp_sum( norm, intra_bgrp_comm )
       !
       norm = SQRT( MAX( norm, 0.0_DP ) )
       !
       IF ( norm < eps16 ) &
       CALL errore( ' gram_schmidt_gamma ', ' vectors are linear dependent ', 1 )
       !
       phi(:,ibnd) = phi(:,ibnd) / norm
       !
       ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) phi(1,ibnd) = CMPLX( DBLE( phi(1,ibnd) ), 0._DP, kind=DP )
       !
       IF ( uspp ) THEN
          !
          sphi(:,ibnd) = sphi(:,ibnd) / norm
          !
          ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
          IF ( gstart == 2 ) sphi(1,ibnd) = CMPLX( DBLE( sphi(1,ibnd) ), 0._DP, kind=DP )
          !
       END IF
       !
    END DO
    !
    DEALLOCATE( sr )
    !
    RETURN
    !
  END SUBROUTINE gram_schmidt_diag
  !
  !
  SUBROUTINE project_offdiag( ibnd_start, ibnd_end, jbnd_start, jbnd_end )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: ibnd_start, ibnd_end
    INTEGER, INTENT(IN)  :: jbnd_start, jbnd_end
    !
    INTEGER               :: ibnd_size
    INTEGER               :: jbnd_size
    REAL(DP), ALLOCATABLE :: sr(:,:)
    !
    ibnd_size = ibnd_end - ibnd_start + 1
    jbnd_size = jbnd_end - jbnd_start + 1
    !
    ALLOCATE( sr( ibnd_start:ibnd_end, jbnd_start:jbnd_end ) )
    !
    ! ... <phi_i| S |psi_j>
    !
    IF ( uspp ) THEN
       !
       CALL DGEMM( 'T', 'N', ibnd_size, jbnd_size, 2 * npw, 2._DP, phi(1,ibnd_start), 2 * npwx, &
                   spsi(1,jbnd_start), 2 * npwx, 0._DP, sr(ibnd_start,jbnd_start), ibnd_size )
       !
       IF ( gstart == 2 ) &
       CALL DGER( ibnd_size, jbnd_size, -1._DP, psi(1,ibnd_start), 2 * npwx, &
                  spsi(1,jbnd_start), 2 * npwx, sr(ibnd_start,jbnd_start), ibnd_size )
       !
    ELSE
       !
       CALL DGEMM( 'T', 'N', ibnd_size, jbnd_size, 2 * npw, 2._DP, phi(1,ibnd_start), 2 * npwx, &
                   psi(1,jbnd_start), 2 * npwx, 0._DP, sr(ibnd_start,jbnd_start), ibnd_size )
       !
       IF ( gstart == 2 ) &
       CALL DGER( ibnd_size, jbnd_size, -1._DP, psi(1,ibnd_start), 2 * npwx, &
                  psi(1,jbnd_start), 2 * npwx, sr(ibnd_start,jbnd_start), ibnd_size )
       !
    END IF
    !
    CALL mp_sum( sr, intra_bgrp_comm )
    !
    ! ... phi_j = phi_j - phi_i * <phi_i| S |psi_j>
    !
    CALL DGEMM( 'N', 'N', 2 * npw, jbnd_size, ibnd_size, -1._DP, phi(1,ibnd_start), 2 * npwx, &
                sr(ibnd_start,jbnd_start), ibnd_size, 1._DP, phi(1,jbnd_start), 2 * npwx )
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) phi(1,jbnd_start:jbnd_end) = &
                       CMPLX( DBLE( phi(1,jbnd_start:jbnd_end) ), 0._DP, kind=DP )
    !
    IF ( uspp ) THEN
       !
       CALL DGEMM( 'N', 'N', 2 * npw, jbnd_size, ibnd_size, -1._DP, sphi(1,ibnd_start), 2 * npwx, &
                   sr(ibnd_start,jbnd_start), ibnd_size, 1._DP, sphi(1,jbnd_start), 2 * npwx )
       !
       ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) sphi(1,jbnd_start:jbnd_end) = &
                          CMPLX( DBLE( sphi(1,jbnd_start:jbnd_end) ), 0._DP, kind=DP )
       !
    END IF
    !
    DEALLOCATE( sr )
    !
    RETURN
    !
  END SUBROUTINE project_offdiag
  !
  !
END SUBROUTINE gram_schmidt_gamma
