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
#define MONE (-1._DP, 0._DP )
!
!--------------------------------------------------------------------------
SUBROUTINE gram_schmidt_k( npwx, npw, nbnd, npol, psi, uspp, nbsize )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for k-point calculations.
  ! ... blocking algorithm is used.
  !
  USE constants, ONLY : eps16
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum, mp_max
  USE mp_bands,  ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id, set_bgrp_indices
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
  LOGICAL,     INTENT(IN)    :: uspp
  INTEGER,     INTENT(IN)    :: nbsize
  !
  ! ... local variables
  !
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: iblock, nblock
  INTEGER                  :: iblock_start, iblock_end
  INTEGER                  :: jblock_start, jblock_end
  INTEGER                  :: ibnd_start, ibnd_end
  INTEGER                  :: jbnd_start, jbnd_end
  COMPLEX(DP), ALLOCATABLE :: phi(:,:), spsi(:,:), sphi(:,:)
  INTEGER,     ALLOCATABLE :: owner_bgrp_id(:)
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
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
  ALLOCATE( phi ( kdmx, nbnd ) )
  IF ( uspp ) ALLOCATE( spsi( kdmx, nbnd ) )
  IF ( uspp ) ALLOCATE( sphi( kdmx, nbnd ) )
  ALLOCATE( owner_bgrp_id( nblock ) )
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
  ! ... Operate the overlap : S |psi_j>
  !
  IF ( uspp ) CALL s_psi( npwx, npw, nbnd, psi, spsi )
  !
  ! ... Set initial : |phi_j> = |psi_j>
  !
  phi = psi
  !
  IF ( uspp ) sphi = spsi
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
  DEALLOCATE( phi  )
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
    INTEGER                  :: ibnd
    COMPLEX(DP), ALLOCATABLE :: sc(:)
    REAL(DP)                 :: norm
    REAL(DP),    EXTERNAL    :: ZDOTC
    !
    ALLOCATE( sc( ibnd_start:ibnd_end ) )
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( ibnd > ibnd_start ) THEN
          !
          ! ... <phi_j| S |psi_i>
          !
          IF ( uspp ) THEN
             !
             CALL ZGEMV( 'C', kdim, ibnd - ibnd_start, ONE, phi(1,ibnd_start), kdmx, &
                         spsi(1,ibnd), 1, ZERO, sc(ibnd_start), 1 )
             !
          ELSE
             !
             CALL ZGEMV( 'C', kdim, ibnd - ibnd_start, ONE, phi(1,ibnd_start), kdmx, &
                         psi(1,ibnd), 1, ZERO, sc(ibnd_start), 1 )
             !
          END IF
          !
          CALL mp_sum( sc, intra_bgrp_comm )
          !
          ! ... phi_i = phi_i - phi_j * <phi_j| S |psi_i>
          !
          CALL ZGEMV( 'N', kdim, ibnd - ibnd_start, MONE, phi(1,ibnd_start), kdmx, &
                      sc(ibnd_start), 1, ONE, phi(1,ibnd), 1 )
          !
          IF ( uspp ) &
          CALL ZGEMV( 'N', kdim, ibnd - ibnd_start, MONE, sphi(1,ibnd_start), kdmx, &
                      sc(ibnd_start), 1, ONE, sphi(1,ibnd), 1 )
          !
       END IF
       !
       ! ... Normalize : phi_i = phi_i / SQRT(<phi_i| S |phi_i>)
       !
       IF ( uspp ) THEN
          !
          norm = ZDOTC( kdim, phi(1,ibnd), 1, sphi(1,ibnd), 1 )
          !
       ELSE
          !
          norm = ZDOTC( kdim, phi(1,ibnd), 1, phi(1,ibnd), 1 )
          !
       END IF
       !
       CALL mp_sum( norm, intra_bgrp_comm )
       !
       norm = SQRT( MAX( norm, 0.0_DP ) )
       !
       IF ( norm < eps16 ) &
       CALL errore( ' gram_schmidt_k ', ' vectors are linear dependent ', 1 )
       !
       phi(:,ibnd) = phi(:,ibnd) / norm
       !
       IF ( uspp ) &
       sphi(:,ibnd) = sphi(:,ibnd) / norm
       !
    END DO
    !
    DEALLOCATE( sc )
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
    INTEGER,                 :: ibnd_size
    INTEGER,                 :: jbnd_size
    COMPLEX(DP), ALLOCATABLE :: sc(:,:)
    !
    ibnd_size = ibnd_end - ibnd_start + 1
    jbnd_size = jbnd_end - jbnd_start + 1
    !
    ALLOCATE( sc( ibnd_start:ibnd_end, jbnd_start:jbnd_end ) )
    !
    ! ... <phi_i| S |psi_j>
    !
    IF ( uspp ) THEN
       !
       CALL ZGEMM( 'C', 'N', ibnd_size, jbnd_size, kdim, ONE, phi(1,ibnd_start), kdmx, &
                   spsi(1,jbnd_start), kdmx, ZERO, sc(ibnd_start,jbnd_start), ibnd_size )
       !
    ELSE
       !
       CALL ZGEMM( 'C', 'N', ibnd_size, jbnd_size, kdim, ONE, phi(1,ibnd_start), kdmx, &
                   psi(1,jbnd_start), kdmx, ZERO, sc(ibnd_start,jbnd_start), ibnd_size )
       !
    END IF
    !
    CALL mp_sum( sc, intra_bgrp_comm )
    !
    ! ... phi_j = phi_j - phi_i * <phi_i| S |psi_j>
    !
    CALL ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, phi(1,ibnd_start), kdmx, &
                sc(ibnd_start,jbnd_start), ibnd_size, ONE, phi(1,jbnd_start), kdmx )
    !
    IF ( uspp ) &
    CALL ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, sphi(1,ibnd_start), kdmx, &
                sc(ibnd_start,jbnd_start), ibnd_size, ONE, sphi(1,jbnd_start), kdmx )
    !
    DEALLOCATE( sc )
    !
    RETURN
    !
  END SUBROUTINE project_offdiag
  !
  !
END SUBROUTINE gram_schmidt_k
