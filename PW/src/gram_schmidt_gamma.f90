!
! Copyright (C) 2017 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE gram_schmidt_gamma( npwx, npw, nbnd, npol, psi, overlap, gstart, nbsize )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for Gamma-only calculations.
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
  INTEGER,     INTENT(IN)    :: gstart
  INTEGER,     INTENT(IN)    :: nbsize
  !
  ! ... local variables
  !
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: iblock, nblock
  INTEGER                  :: iblock_start, iblock_end
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  REAL(DP),    ALLOCATABLE :: sr(:,:)
  !
  ALLOCATE( aux( npwx, nbnd ) )
  ALLOCATE( sr( nbnd, nbnd ) )
  !
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  !
  IF ( gstart == 2 ) psi(1,1:nbnd) = CMPLX( DBLE( psi(1,1:nbnd) ), 0._DP, kind=DP )
  !
  ! ... Set up the Overlap matrix on the subspace :
  !
  ! ...      S_ij = <psi_i| S |psi_j>
  !
  IF ( overlap ) THEN
     !
     CALL s_psi( npwx, npw, nbnd, psi, aux )
     !
     CALL DGEMM( 'T', 'N', nbnd, nbnd, 2 * npw, 2._DP, psi, 2 * npwx, aux, 2 * npwx, 0._DP, sr, nbnd )
     !
     IF ( gstart == 2 ) CALL DGER( nbnd, nbnd, -1._DP, psi, 2 * npwx, aux, 2 * npwx, sr, nbnd )
     !
  ELSE
     !
     CALL DGEMM( 'T', 'N', nbnd, nbnd, 2 * npw, 2._DP, psi, 2 * npwx, psi, 2 * npwx, 0._DP, sr, nbnd )
     !
     IF ( gstart == 2 ) CALL DGER( nbnd, nbnd, -1._DP, psi, 2 * npwx, psi, 2 * npwx, sr, nbnd )
     !
  END IF
  !
  CALL mp_sum( sr, intra_bgrp_comm )
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
  DEALLOCATE( sr )
  !
  RETURN
  !
END SUBROUTINE gram_schmidt_gamma
