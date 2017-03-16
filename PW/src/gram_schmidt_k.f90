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
SUBROUTINE gram_schmidt_k( npwx, npw, nbnd, npol, psi, overlap, nblock )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for k-point calculations.
  ! ... blocking algorithm is used.
  !
  USE kinds,    ONLY : DP
  USE mp,       ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
  LOGICAL,     INTENT(IN)    :: overlap
  INTEGER,     INTENT(IN)    :: nblock
  !
  ! ... local variables
  !
  INTEGER                  :: kdim, kdmx
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
  CALL mp_sum( sc , intra_bgrp_comm )
  !

  !
  ! TODO
  ! TODO
  ! TODO
  !

  !
  DEALLOCATE( aux )
  DEALLOCATE( sc )
  !
  RETURN
  !
END SUBROUTINE gram_schmidt_k
