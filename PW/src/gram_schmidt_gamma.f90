!
! Copyright (C) 2017 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE gram_schmidt_gamma( npwx, npw, nbnd, npol, psi, overlap )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for Gamma-only calculations.
  ! ... blocking algorithm is used.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd, npol
  LOGICAL,     INTENT(IN)    :: overlap
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
  !
  ! TODO
  ! TODO
  ! TODO
  !
END SUBROUTINE gram_schmidt_gamma
