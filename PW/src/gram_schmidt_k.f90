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
!--------------------------------------------------------------------------
SUBROUTINE gram_schmidt_k( npwx, npw, nbnd, npol, psi, overlap, nbsize )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for k-point calculations.
  ! ... blocking algorithm is used.
  !
  USE kinds,    ONLY : DP
  USE mp,       ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm, set_bgrp_indices
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
  LOGICAL,     INTENT(IN)    :: overlap
  INTEGER,     INTENT(IN)    :: nbsize
  !
  ! ... local variables
  !
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: iblock, nblock
  INTEGER                  :: iblock_start, iblock_end
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  COMPLEX(DP), ALLOCATABLE :: sc(:,:)
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
  ALLOCATE( aux( kdmx, nbnd ) )
  ALLOCATE( sc( nbnd, nbnd ) )
  !
  ! ... Set up the Overlap matrix on the subspace :
  !
  ! ...      S_ij = <psi_i| S |psi_j>
  !
  IF ( overlap ) THEN
     !
     CALL s_psi( npwx, npw, nbnd, psi, aux )
     !
     CALL ZGEMM( 'C', 'N', nbnd, nbnd, kdim, ONE, psi, kdmx, aux, kdmx, ZERO, sc, nbnd )
     !
  ELSE
     !
     CALL ZGEMM( 'C', 'N', nbnd, nbnd, kdim, ONE, psi, kdmx, psi, kdmx, ZERO, sc, nbnd )
     !
  END IF
  !
  CALL mp_sum( sc, intra_bgrp_comm )
  !
  ! ... Blocking loop
  !
  nblock = nbnd / nbsize
  IF ( MOD( nbnd, nbsize ) /= 0 ) nblock = nblock + 1
  !
  CALL set_bgrp_indices( nblock, iblock_start, iblock_end )
  !
  DO iblock = 1, nblock
     !
     ! TODO
     ! TODO
     ! TODO
     !
  END DO
  !
  DEALLOCATE( aux )
  DEALLOCATE( sc )
  !
  RETURN
  !
END SUBROUTINE gram_schmidt_k
