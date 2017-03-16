!
! Copyright (C) 2017 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE gram_schmidt( npwx, npw, nbnd, npol, psi, overlap, gstart, nblock )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization.
  !
  USE kinds,         ONLY : DP
  USE control_flags, ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
  LOGICAL,     INTENT(IN)    :: overlap
  INTEGER,     INTENT(IN)    :: gstart
  INTEGER,     INTENT(IN)    :: nblock
  !
  CALL start_clock( 'gsorth' )
  !
  IF ( gamma_only ) THEN
     !
     CALL gram_schmidt_gamma( npwx, npw, nbnd, npol, psi, overlap, gstart, nblock )
     !
  ELSE
     !
     CALL gram_schmidt_k( npwx, npw, nbnd, npol, psi, overlap, nblock )
     !
  END IF
  !
  CALL stop_clock( 'gsorth' )
  !
END SUBROUTINE gram_schmidt
