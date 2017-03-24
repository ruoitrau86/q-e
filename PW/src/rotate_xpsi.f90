!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! >>> This file is copied from rotate_wfc.f90,
! >>> and modified by Satomichi Nishihara
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_xpsi &
            ( npwx, npw, nstart, gstart, nbnd, psi, npol, overlap, evc, hevc, sevc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine (maybe it should be an interface) for
  ! ... Hamiltonian diagonalization in the subspace spanned
  ! ... by nstart states psi ( atomic or random wavefunctions ).
  ! ... Produces on output nbnd eigenvectors ( nbnd <= nstart ) in evc.
  ! ... Calls h_psi, s_psi to calculate H|psi> and S|psi>,
  ! ... which are saved in hevc and sevc.
  !
  USE kinds,         ONLY : DP
  USE control_flags, ONLY : use_para_diag, gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, gstart, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! first G with nonzero norm
    ! number of spin polarizations
  LOGICAL, INTENT(IN) :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP), INTENT(IN)  :: psi(npwx*npol,nstart)
  COMPLEX(DP), INTENT(OUT) :: evc(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT) :: hevc(npwx*npol,nbnd), sevc(npwx*npol,nbnd)
    ! H|psi> and S|psi>
  REAL(DP), INTENT(OUT) :: e(nbnd)
    ! eigenvalues
  !
  CALL start_clock( 'wfcrot' )
  !
  IF( use_para_diag ) THEN
     !
     ! use data distributed subroutine
     !
     IF ( gamma_only ) THEN
        !
        CALL protate_xpsi_gamma &
            ( npwx, npw, nstart, gstart, nbnd, psi, overlap, evc, hevc, sevc, e )
        !
     ELSE
        !
        CALL protate_xpsi_k &
            ( npwx, npw, nstart, nbnd, npol, psi, overlap, evc, hevc, sevc, e )
        !
     END IF
     !
  ELSE
     !
     ! use serial subroutines
     !
     IF ( gamma_only ) THEN
        !
        CALL rotate_xpsi_gamma &
            ( npwx, npw, nstart, gstart, nbnd, psi, overlap, evc, hevc, sevc, e )
        !
     ELSE
        !
        CALL rotate_xpsi_k &
            ( npwx, npw, nstart, nbnd, npol, psi, overlap, evc, hevc, sevc, e )
        !
     END IF
     !
  END IF
  !
  CALL stop_clock( 'wfcrot' )
  !
END SUBROUTINE rotate_xpsi
