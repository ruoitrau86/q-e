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
SUBROUTINE rrmmdiagg( npwx, npw, nbnd, psi, spsi, e, &
                      g2kin, btype, ethr, ndiis, uspp, gstart, notconv, rmm_iter )
  !----------------------------------------------------------------------------
  !
  ! ... Iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned RMM-DIIS algorithm.
  !
  USE constants, ONLY : eps14, eps16
  USE kinds,     ONLY : DP
  USE funct,     ONLY : exx_is_active
  USE mp,        ONLY : mp_sum, mp_bcast
  USE mp_bands,  ONLY : inter_bgrp_comm, intra_bgrp_comm, me_bgrp, root_bgrp, &
                        root_bgrp_id, use_bgrp_in_hpsi, set_bgrp_indices
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npwx, npw, nbnd
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx,nbnd)
  COMPLEX(DP), INTENT(OUT)   :: spsi(npwx,nbnd)
  REAL(DP),    INTENT(OUT)   :: e(nbnd)
  REAL(DP),    INTENT(IN)    :: g2kin(npwx)
  INTEGER,     INTENT(IN)    :: btype(nbnd)
  REAL(DP),    INTENT(IN)    :: ethr
  INTEGER,     INTENT(IN)    :: ndiis
  LOGICAL,     INTENT(IN)    :: uspp
  INTEGER,     INTENT(IN)    :: gstart
  INTEGER,     INTENT(OUT)   :: notconv
  INTEGER,     INTENT(OUT)   :: rmm_iter
  !
  ! ... local variables
  !
  INTEGER                  :: ierr
  INTEGER                  :: idiis
  INTEGER                  :: ibnd_start, ibnd_end, ibnd_size
  REAL(DP)                 :: empty_ethr
  COMPLEX(DP), ALLOCATABLE :: phi(:,:,:), hphi(:,:,:), sphi(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), kpsi(:,:), hkpsi(:,:), skpsi(:,:)
  REAL(DP),    ALLOCATABLE :: hr(:,:,:), sr(:,:,:)
  REAL(DP),    ALLOCATABLE :: php(:,:), psp(:,:)
  REAL(DP),    ALLOCATABLE :: ew(:), hw(:), sw(:)
  LOGICAL,     ALLOCATABLE :: conv(:)
  !
  REAL(DP),    PARAMETER   :: SREF = 0.5_DP
  REAL(DP),    PARAMETER   :: SMIN = 0.1_DP
  REAL(DP),    PARAMETER   :: SMAX = 2.0_DP
  !
  CALL start_clock( 'rrmmdiagg' )
  !
  empty_ethr = MAX( ( ethr * 5._DP ), 1.E-5_DP )
  !
  CALL set_bgrp_indices( nbnd, ibnd_start, ibnd_end )
  !
  ibnd_size = MAX( ibnd_end - ibnd_start + 1, 0 )
  !
  ALLOCATE( phi( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate phi ', ABS(ierr) )
  !
  ALLOCATE( hphi( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate hphi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( sphi( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate sphi ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( hpsi( npwx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate hpsi ', ABS(ierr) )
  !
  ALLOCATE( kpsi( npwx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate kpsi ', ABS(ierr) )
  !
  ALLOCATE( hkpsi( npwx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate hkpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( skpsi( npwx, nbnd ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate skpsi ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( hr( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate hr ', ABS(ierr) )
  !
  ALLOCATE( sr( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate sr ', ABS(ierr) )
  !
  ALLOCATE( php( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate php ', ABS(ierr) )
  !
  ALLOCATE( psp( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate psp ', ABS(ierr) )
  !
  ALLOCATE( ew( nbnd ) )
  ALLOCATE( hw( nbnd ) )
  ALLOCATE( sw( nbnd ) )
  ALLOCATE( conv( nbnd ) )
  !
  phi  = ZERO
  hphi = ZERO
  IF ( uspp ) sphi = ZERO
  !
  hpsi  = ZERO
  IF ( uspp ) spsi = ZERO
  kpsi  = ZERO
  hkpsi = ZERO
  IF ( uspp ) skpsi = ZERO
  !
  hr = 0._DP
  sr = 0._DP
  !
  php = 0._DP
  psp = 0._DP
  !
  e    = 0._DP
  ew   = 0._DP
  hw   = 0._DP
  sw   = 0._DP
  conv = .FALSE.
  !
  ! ... Set Im[ psi(G=0) ] - needed for numerical stability
  !
  IF ( gstart == 2 ) psi(1,1:nbnd) = CMPLX( DBLE( psi(1,1:nbnd) ), 0._DP, kind=DP )
  !
  ! ... Initial eigenvalues
  !
  CALL eigenvalues( .TRUE. )
  !
  ! ... RMM-DIIS's loop
  !
  rmm_iter = 0
  notconv  = nbnd
  !
  DO idiis = 1, ndiis
     !
     rmm_iter = rmm_iter + 1
     !
     ! ... Perform DIIS
     !
     CALL do_diis( idiis )
     !
     ! ... Line searching
     !
     CALL line_search( )
     !
     ! ... Calculate eigenvalues and check convergence
     !
     CALL eigenvalues( .FALSE. )
     !
     IF ( notconv == 0 ) EXIT
     !
  END DO
  !
  ! ... Merge wave functions
  !
  IF ( ibnd_start > 1 ) THEN
     !
     psi(:,1:(ibnd_start-1)) = ZERO
     IF ( uspp ) spsi(:,1:(ibnd_start-1)) = ZERO
     !
  END IF
  !
  IF ( ibnd_end < nbnd ) THEN
     !
     psi(:,(ibnd_end+1):nbnd) = ZERO
     IF ( uspp ) spsi(:,(ibnd_end+1):nbnd) = ZERO
     !
  END IF
  !
  CALL mp_sum( psi, inter_bgrp_comm )
  IF ( uspp ) CALL mp_sum( spsi, inter_bgrp_comm )
  !
  DEALLOCATE( phi )
  DEALLOCATE( hphi )
  IF ( uspp ) DEALLOCATE( sphi )
  DEALLOCATE( hpsi )
  DEALLOCATE( kpsi )
  DEALLOCATE( hkpsi )
  IF ( uspp ) DEALLOCATE( skpsi )
  DEALLOCATE( hr )
  DEALLOCATE( sr )
  DEALLOCATE( php )
  DEALLOCATE( psp )
  DEALLOCATE( ew )
  DEALLOCATE( hw )
  DEALLOCATE( sw )
  DEALLOCATE( conv )
  !
  CALL stop_clock( 'rrmmdiagg' )
  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE do_diis( idiis )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: idiis
    !
    INTEGER                  :: ibnd
    INTEGER                  :: kdiis
    REAL(DP)                 :: norm
    REAL(DP)                 :: er
    COMPLEX(DP), ALLOCATABLE :: vec1(:)
    COMPLEX(DP), ALLOCATABLE :: vec2(:,:)
    REAL(DP),    ALLOCATABLE :: vr(:)
    REAL(DP),    ALLOCATABLE :: tr(:,:)
    !
    ALLOCATE( vec1( npwx ) )
    ALLOCATE( vec2( npwx, idiis ) )
    IF ( idiis > 1 ) ALLOCATE( vr( idiis ) )
    ALLOCATE( tr( idiis, ibnd_start:ibnd_end ) )
    !
    ! ... Save current wave functions and eigenvalues
    !
    CALL DCOPY( 2 * npwx * ibnd_size, psi (1,ibnd_start), 1, phi (1,ibnd_start,idiis), 1 )
    CALL DCOPY( 2 * npwx * ibnd_size, hpsi(1,ibnd_start), 1, hphi(1,ibnd_start,idiis), 1 )
    IF ( uspp ) &
    CALL DCOPY( 2 * npwx * ibnd_size, spsi(1,ibnd_start), 1, sphi(1,ibnd_start,idiis), 1 )
    !
    CALL DCOPY( ibnd_size, hw(ibnd_start), 1, php(ibnd_start,idiis), 1 )
    CALL DCOPY( ibnd_size, sw(ibnd_start), 1, psp(ibnd_start,idiis), 1 )
    !
    ! ... <R_i|R_j>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       ! ... Residual vectors : |R> = (H - e S) |psi>
       !
       DO kdiis = 1, idiis
          !
          er = php(ibnd,kdiis)
          !
          CALL DCOPY( 2 * npw, hphi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
          !
          IF ( uspp ) THEN
             !
             CALL DAXPY( 2 * npw, -er, sphi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
             !
          ELSE
             !
             CALL DAXPY( 2 * npw, -er, phi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
             !
          END IF
          !
       END DO
       !
       er = php(ibnd,idiis)
       !
       CALL DCOPY( 2 * npw, hphi(1,ibnd,idiis), 1, vec1(1), 1 )
       !
       IF ( uspp ) THEN
          !
          CALL DAXPY( 2 * npw, -er, sphi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       ELSE
          !
          CALL DAXPY( 2 * npw, -er, phi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       END IF
       !
       CALL DGEMV( 'T', 2 * npw, idiis, 2._DP, vec2(1,1), 2 * npwx, &
                   vec1(1), 1, 0._DP, tr(1,ibnd), 1 )
       !
       IF ( gstart == 2 ) &
       tr(1:idiis,ibnd) = tr(1:idiis,ibnd) - DBLE( vec2(1,1:idiis) ) * DBLE( vec1(1) )
       !
    END DO
    !
    CALL mp_sum( tr, intra_bgrp_comm )
    hr(1:idiis,idiis,:) = tr(1:idiis,:)
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       hr(idiis,1:idiis,ibnd) = hr(1:idiis,idiis,ibnd)
       !
    END DO
    !
    ! ... <phi_i| S |phi_j>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       DO kdiis = 1, idiis
          !
          CALL DCOPY( 2 * npw, phi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
          !
       END DO
       !
       IF ( uspp ) THEN
          !
          CALL DCOPY( 2 * npw, sphi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       ELSE
          !
          CALL DCOPY( 2 * npw, phi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       END IF
       !
       CALL DGEMV( 'T', 2 * npw, idiis, 2._DP, vec2(1,1), 2 * npwx, &
                   vec1(1), 1, 0._DP, tr(1,ibnd), 1 )
       !
       IF ( gstart == 2 ) &
       tr(1:idiis,ibnd) = tr(1:idiis,ibnd) - DBLE( vec2(1,1:idiis) ) * DBLE( vec1(1) )
       !
    END DO
    !
    CALL mp_sum( tr, intra_bgrp_comm )
    sr(1:idiis,idiis,:) = tr(1:idiis,:)
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       sr(idiis,1:idiis,ibnd) = sr(1:idiis,idiis,ibnd)
       !
    END DO
    !
    ! ... Update current wave functions and residual vectors
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( idiis > 1 ) THEN
          !
          ! ... solve Rv = eSv
          !
          IF ( me_bgrp == root_bgrp ) CALL diag_diis( ibnd, idiis, vr(:) )
          CALL mp_bcast( vr, root_bgrp, intra_bgrp_comm )
          !
          psi (:,ibnd) = ZERO
          hpsi(:,ibnd) = ZERO
          IF ( uspp ) spsi(:,ibnd) = ZERO
          kpsi(:,ibnd) = ZERO
          !
          DO kdiis = 1, idiis
             !
             ! ... Wave functions
             !
             CALL DAXPY( 2 * npw, vr(kdiis), phi (1,ibnd,kdiis), 1, psi (1,ibnd), 1 )
             CALL DAXPY( 2 * npw, vr(kdiis), hphi(1,ibnd,kdiis), 1, hpsi(1,ibnd), 1 )
             IF ( uspp ) &
             CALL ZAXPY( 2 * npw, vr(kdiis), sphi(1,ibnd,kdiis), 1, spsi(1,ibnd), 1 )
             !
             ! ... Residual vectors
             !
             er = php(ibnd,kdiis)
             !
             CALL DCOPY( 2 * npw, hphi(1,ibnd,kdiis), 1, vec1(1), 1 )
             !
             IF ( uspp ) THEN
                !
                CALL DAXPY( 2 * npw, -er, sphi(1,ibnd,kdiis), 1, vec1(1), 1 )
                !
             ELSE
                !
                CALL DAXPY( 2 * npw, -er, phi(1,ibnd,kdiis), 1, vec1(1), 1 )
                !
             END IF
             !
             CALL DAXPY( 2 * npw, vr(kdiis), vec1(1), 1, kpsi(1,ibnd), 1 )
             !
          END DO
          !
       ELSE
          !
          ! ... Wave functions
          !
          norm = SQRT( sw(ibnd) )
          CALL DSCAL( 2 * npw, 1._DP / norm, psi (1,ibnd), 1 )
          CALL DSCAL( 2 * npw, 1._DP / norm, hpsi(1,ibnd), 1 )
          IF ( uspp ) &
          CALL DSCAL( 2 * npw, 1._DP / norm, spsi(1,ibnd), 1 )
          !
          ! ... Residual vectors
          !
          er = hw(ibnd)
          !
          CALL DCOPY( 2 * npw, hpsi(1,ibnd), 1, kpsi(1,ibnd), 1 )
          !
          IF ( uspp ) THEN
             !
             CALL DAXPY( 2 * npw, -er, spsi(1,ibnd), 1, kpsi(1,ibnd), 1 )
             !
          ELSE
             !
             CALL DAXPY( 2 * npw, -er, spsi(1,ibnd), 1, kpsi(1,ibnd), 1 )
             !
          END IF
          !
       END IF
       !
    END DO
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) THEN
       !
       psi (1,1:nbnd) = CMPLX( DBLE( psi (1,1:nbnd) ), 0._DP, kind=DP )
       hpsi(1,1:nbnd) = CMPLX( DBLE( hpsi(1,1:nbnd) ), 0._DP, kind=DP )
       IF ( uspp ) &
       spsi(1,1:nbnd) = CMPLX( DBLE( spsi(1,1:nbnd) ), 0._DP, kind=DP )
       kpsi(1,1:nbnd) = CMPLX( DBLE( kpsi(1,1:nbnd) ), 0._DP, kind=DP )
       !
    END IF
    !
    DEALLOCATE( vec1 )
    DEALLOCATE( vec2 )
    IF ( idiis > 1 ) DEALLOCATE( vr )
    DEALLOCATE( tr )
    !
    RETURN
    !
  END SUBROUTINE do_diis
  !
  !
  SUBROUTINE diag_diis( ibnd, idiis, vr )
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN)  :: ibnd
    INTEGER,  INTENT(IN)  :: idiis
    REAL(DP), INTENT(OUT) :: vr(idiis)
    !
    INTEGER               :: info
    INTEGER               :: ndim, kdim
    INTEGER               :: i, imin
    REAL(DP)              :: emin
    REAL(DP)              :: vnrm
    REAL(DP), ALLOCATABLE :: h1(:,:)
    REAL(DP), ALLOCATABLE :: h2(:,:)
    REAL(DP), ALLOCATABLE :: h3(:,:)
    REAL(DP), ALLOCATABLE :: s1(:,:)
    REAL(DP), ALLOCATABLE :: x1(:,:)
    REAL(DP), ALLOCATABLE :: u1(:)
    REAL(DP), ALLOCATABLE :: e1(:)
    INTEGER               :: nwork
    REAL(DP), ALLOCATABLE :: work(:)
    !
    REAL(DP), EXTERNAL    :: DDOT
    !
    ndim  = idiis
    nwork = 3 * ndim
    !
    ALLOCATE( h1( ndim, ndim ) )
    ALLOCATE( h2( ndim, ndim ) )
    ALLOCATE( h3( ndim, ndim ) )
    ALLOCATE( s1( ndim, ndim ) )
    ALLOCATE( x1( ndim, ndim ) )
    ALLOCATE( u1( ndim ) )
    ALLOCATE( e1( ndim ) )
    ALLOCATE( work( nwork ) )
    !
    h1(1:ndim,1:ndim) = hr(1:ndim,1:ndim,ibnd)
    s1(1:ndim,1:ndim) = sr(1:ndim,1:ndim,ibnd)
    !
    CALL DSYEV( 'V', 'U', ndim, s1, ndim, e1, work, nwork, info )
    !
    IF( info /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot solve diis ', ABS(info) )
    !
    kdim = 0
    !
    x1 = 0._DP
    !
    DO i = 1, ndim
       !
       IF ( e1(i) > eps14 ) THEN
          !
          kdim = kdim + 1
          !
          x1(:,kdim) = s1(:,i) / SQRT(e1(i))
          !
       END IF
       !
    END DO
    !
    IF ( kdim <= 1 ) THEN
       !
       vr        = 0._DP
       vr(idiis) = 1._DP
       !
       GOTO 10
       !
    END IF
    !
    h2 = 0._DP
    !
    CALL DGEMM( 'N', 'N', ndim, kdim, ndim, 1._DP, h1, ndim, x1, ndim, 0._DP, h2, ndim )
    !
    h3 = 0._DP
    !
    CALL DGEMM( 'T', 'N', kdim, kdim, ndim, 1._DP, x1, ndim, h2, ndim, 0._DP, h3, ndim )
    !
    e1 = 0._DP
    !
    CALL DSYEV( 'V', 'U', kdim, h3, ndim, e1, work, nwork, info )
    !
    IF( info /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot solve diis ', ABS(info) )
    !
    imin = 1
    emin = e1(1)
    !
    DO i = 2, kdim
       !
       IF ( e1(i) < emin ) imin = i
       !
    END DO
    !
    CALL DGEMV( 'N', ndim, kdim, 1._DP, x1, ndim, h3(:,imin), 1, 0._DP, vr, 1 )
    !
    s1(1:ndim,1:ndim) = sr(1:ndim,1:ndim,ibnd)
    !
    CALL DGEMV( 'N', ndim, ndim, 1._DP, s1, ndim, vr, 1, 0._DP, u1, 1 )
    !
    vnrm = SQRT( DDOT( ndim, vr, 1, u1, 1 ) )
    !
    vr = vr / vnrm
    !
10  DEALLOCATE( h1 )
    DEALLOCATE( h2 )
    DEALLOCATE( h3 )
    DEALLOCATE( s1 )
    DEALLOCATE( x1 )
    DEALLOCATE( u1 )
    DEALLOCATE( e1 )
    DEALLOCATE( work )
    !
    RETURN
    !
  END SUBROUTINE diag_diis
  !
  !
  SUBROUTINE line_search( )
    !
    IMPLICIT NONE
    !
    INTEGER               :: ibnd, ig
    REAL(DP)              :: psir, psii, psi2
    REAL(DP)              :: kdiag, k1, k2
    REAL(DP)              :: x, x2, x3, x4
    REAL(DP), ALLOCATABLE :: ekin(:)
    REAL(DP)              :: a, b
    REAL(DP)              :: ene0, ene1
    REAL(DP)              :: step, norm
    REAL(DP)              :: php, khp, khk
    REAL(DP)              :: psp, ksp, ksk
    REAL(DP), ALLOCATABLE :: hmat(:,:), smat(:,:)
    REAL(DP)              :: c1, c2
    REAL(DP), ALLOCATABLE :: coef(:,:)
    !
    REAL(DP), EXTERNAL    :: DDOT
    !
    ALLOCATE( ekin( ibnd_start:ibnd_end ) )
    ALLOCATE( hmat( 3, ibnd_start:ibnd_end ) )
    ALLOCATE( smat( 3, ibnd_start:ibnd_end ) )
    ALLOCATE( coef( 2, ibnd_start:ibnd_end ) )
    !
    ! ... Kinetic energy
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       ekin(ibnd) = 0._DP
       !
       DO ig = gstart, npw
          !
          psir = DBLE ( psi(ig,ibnd) )
          psii = AIMAG( psi(ig,ibnd) )
          psi2 = psir * psir + psii * psii
          ekin(ibnd) = ekin(ibnd) + 2._DP * g2kin(ig) * psi2
          !
       END DO
       !
       IF ( gstart == 2 ) THEN
          !
          psir = DBLE ( psi(ig,ibnd) )
          psi2 = psir * psir
          ekin(ibnd) = ekin(ibnd) + g2kin(1) * psi2
          !
       END IF
       !
    END DO
    !
    CALL mp_sum( ekin, intra_bgrp_comm )
    !
    ! ... Preconditioning vectors : K (H - e S) |psi>
    !
    ! ... G.Kresse and J.Furthmuller, PRB 54, 11169 (1996)
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       DO ig = 1, npw
          !
          x  = g2kin(ig) / ( 1.5_DP * ekin(ibnd) )
          x2 = x * x
          x3 = x * x2
          x4 = x * x3
          !
          k1 = 27._DP + 18._DP * x + 12._DP * x2 + 8._DP * x3
          k2 = k1 + 16._DP * x4
          kdiag = ( -4._DP / 3._DP / ekin(ibnd) ) * k1 / k2
          !
          kpsi(ig,ibnd) = kdiag * kpsi(ig,ibnd)
          !
       END DO
       !
    END DO
    !
    ! ... Share kpsi for all band-groups
    !
    IF ( ( .NOT. use_bgrp_in_hpsi ) .OR. exx_is_active() ) THEN
       !
       DO ibnd = 1, ( ibnd_start - 1)
          !
          kpsi(:,ibnd) = ZERO
          !
       END DO
       !
       DO ibnd = ( ibnd_end + 1 ), nbnd
          !
          kpsi(:,ibnd) = ZERO
          !
       END DO
       !
       CALL mp_sum( kpsi, inter_bgrp_comm )
       !
    END IF
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) &
    kpsi (1,1:nbnd) = CMPLX( DBLE( kpsi (1,1:nbnd) ), 0._DP, kind=DP )
    !
    ! ... Operate the Hamiltonian : H K (H - eS) |psi>
    !
    CALL h_psi( npwx, npw, nbnd, kpsi, hkpsi )
    !
    ! ... Operate the Overlap : S K (H - eS) |psi>
    !
    IF ( uspp ) CALL s_psi( npwx, npw, nbnd, kpsi, skpsi )
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) THEN
       !
       hkpsi(1,1:nbnd) = CMPLX( DBLE( hkpsi(1,1:nbnd) ), 0._DP, kind=DP )
       IF ( uspp ) &
       skpsi(1,1:nbnd) = CMPLX( DBLE( skpsi(1,1:nbnd) ), 0._DP, kind=DP )
       !
    END IF
    !
    ! ... Create 2 x 2 matrix
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       php = 2._DP * DDOT( 2 * npw, psi (1,ibnd), 1, hpsi (1,ibnd), 1 )
       khp = 2._DP * DDOT( 2 * npw, kpsi(1,ibnd), 1, hpsi (1,ibnd), 1 )
       khk = 2._DP * DDOT( 2 * npw, kpsi(1,ibnd), 1, hkpsi(1,ibnd), 1 )
       !
       IF ( gstart == 2 ) THEN
          !
          php = php - DBLE( psi (1,ibnd) ) * DBLE ( hpsi (1,ibnd) )
          khp = khp - DBLE( kpsi(1,ibnd) ) * DBLE ( hpsi (1,ibnd) )
          khk = khk - DBLE( kpsi(1,ibnd) ) * DBLE ( hkpsi(1,ibnd) )
          !
       END IF
       !
       IF ( uspp ) THEN
          !
          psp = 2._DP * DDOT( 2 * npw, psi (1,ibnd), 1, spsi (1,ibnd), 1 )
          ksp = 2._DP * DDOT( 2 * npw, kpsi(1,ibnd), 1, spsi (1,ibnd), 1 )
          ksk = 2._DP * DDOT( 2 * npw, kpsi(1,ibnd), 1, skpsi(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) THEN
             !
             psp = psp - DBLE( psi (1,ibnd) ) * DBLE ( spsi (1,ibnd) )
             ksp = ksp - DBLE( kpsi(1,ibnd) ) * DBLE ( spsi (1,ibnd) )
             ksk = ksk - DBLE( kpsi(1,ibnd) ) * DBLE ( skpsi(1,ibnd) )
             !
          END IF
          !
       ELSE
          !
          psp = 2._DP * DDOT( 2 * npw, psi (1,ibnd), 1, psi (1,ibnd), 1 )
          ksp = 2._DP * DDOT( 2 * npw, kpsi(1,ibnd), 1, psi (1,ibnd), 1 )
          ksk = 2._DP * DDOT( 2 * npw, kpsi(1,ibnd), 1, kpsi(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) THEN
             !
             psp = psp - DBLE( psi (1,ibnd) ) * DBLE ( psi (1,ibnd) )
             ksp = ksp - DBLE( kpsi(1,ibnd) ) * DBLE ( psi (1,ibnd) )
             ksk = ksk - DBLE( kpsi(1,ibnd) ) * DBLE ( kpsi(1,ibnd) )
             !
          END IF
          !
       END IF
       !
       hmat(1,ibnd) = php
       hmat(2,ibnd) = khp
       hmat(3,ibnd) = khk
       !
       smat(1,ibnd) = psp
       smat(2,ibnd) = ksp
       smat(3,ibnd) = ksk
       !
    END DO
    !
    CALL mp_sum( hmat, intra_bgrp_comm )
    CALL mp_sum( smat, intra_bgrp_comm )
    !
    ! ... Line searching for each band
    !
    IF ( me_bgrp == root_bgrp ) THEN
       !
       DO ibnd = ibnd_start, ibnd_end
          !
          php = hmat(1,ibnd)
          khp = hmat(2,ibnd)
          khk = hmat(3,ibnd)
          !
          psp = smat(1,ibnd)
          ksp = smat(2,ibnd)
          ksk = smat(3,ibnd)
          IF( psp <= eps16 ) CALL errore( ' rrmmdiagg ', ' psp <= 0 ', 1 )
          !
          norm = psp + 2._DP * ksp * SREF + ksk * SREF * SREF
          IF( norm <= eps16 ) CALL errore( ' rrmmdiagg ', ' norm <= 0 ', 1 )
          !
          ene0 = php / psp
          ene1 = ( php + 2._DP * khp * SREF + khk * SREF * SREF ) / norm
          !
          a = 2._DP * ( khp * psp - php * ksp ) / psp / psp
          b = ( ene1 - ene0 - a * SREF ) / SREF / SREF
          IF( ABS( b ) < eps16 ) CALL errore( ' rrmmdiagg ', ' b == 0 ', 1 )
          !
          step  = -0.5_DP * a / b
          step  = MAX( SMIN, step )
          step  = MIN( SMAX, step )
          norm  = psp + 2._DP * ksp * step + ksk * step * step
          IF( norm <= eps16 ) CALL errore( ' rrmmdiagg ', ' norm <= 0 ', 1 )
          norm  = SQRT( norm )
          !
          coef(1,ibnd) = 1._DP / norm
          coef(2,ibnd) = step  / norm
          !
          ! ... Update current matrix elements
          !
          c1 = coef(1,ibnd)
          c2 = coef(2,ibnd)
          !
          hw(ibnd) = php * c1 * c1 + 2._DP * khp * c1 * c2 + khk * c2 * c2
          sw(ibnd) = psp * c1 * c1 + 2._DP * ksp * c1 * c2 + ksk * c2 * c2
          !
       END DO
       !
    END IF
    !
    CALL mp_bcast( coef, root_bgrp, intra_bgrp_comm )
    CALL mp_bcast( hw(ibnd_start:ibnd_end), root_bgrp, intra_bgrp_comm )
    CALL mp_bcast( sw(ibnd_start:ibnd_end), root_bgrp, intra_bgrp_comm )
    !
    ! ... Update current wave functions
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       c1 = coef(1,ibnd)
       c2 = coef(2,ibnd)
       !
       CALL DSCAL( 2 * npw, c1, psi (1,ibnd), 1 )
       CALL DAXPY( 2 * npw, c2, kpsi(1,ibnd), 1, psi(1,ibnd), 1 )
       !
       CALL DSCAL( 2 * npw, c1, hpsi (1,ibnd), 1 )
       CALL DAXPY( 2 * npw, c2, hkpsi(1,ibnd), 1, hpsi(1,ibnd), 1 )
       !
       IF ( uspp ) THEN
          !
          CALL DSCAL( 2 * npw, c1, spsi (1,ibnd), 1 )
          CALL DAXPY( 2 * npw, c2, skpsi(1,ibnd), 1, spsi(1,ibnd), 1 )
          !
       END IF
       !
    END DO
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) THEN
       !
       psi (1,1:nbnd) = CMPLX( DBLE( psi (1,1:nbnd) ), 0._DP, kind=DP )
       hpsi(1,1:nbnd) = CMPLX( DBLE( hpsi(1,1:nbnd) ), 0._DP, kind=DP )
       IF ( uspp ) &
       spsi(1,1:nbnd) = CMPLX( DBLE( spsi(1,1:nbnd) ), 0._DP, kind=DP )
       !
    END IF
    !
    DEALLOCATE( ekin )
    DEALLOCATE( hmat )
    DEALLOCATE( smat )
    DEALLOCATE( coef )
    !
    RETURN
    !
  END SUBROUTINE line_search
  !
  !
  SUBROUTINE eigenvalues( first )
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: first
    !
    INTEGER :: ibnd
    !
    REAL(DP), EXTERNAL :: DDOT
    !
    IF ( first ) THEN
       !
       ! ... Operate the Hamiltonian : H |psi>
       !
       CALL h_psi( npwx, npw, nbnd, psi, hpsi )
       !
       ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) &
       hpsi(1,1:nbnd) = CMPLX( DBLE( hpsi(1,1:nbnd) ), 0._DP, kind=DP )
       !
       ! ... Operate the Overlap : S |psi>
       !
       IF ( uspp ) THEN
          !
          CALL s_psi( npwx, npw, nbnd, psi, spsi )
          !
          ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
          IF ( gstart == 2 ) &
          spsi(1,1:nbnd) = CMPLX( DBLE( spsi(1,1:nbnd) ), 0._DP, kind=DP )
          !
       END IF
       !
       ! ... Matrix elements
       !
       DO ibnd = ibnd_start, ibnd_end
          !
          hw(ibnd) = 2._DP * DDOT( 2 * npw, psi(1,ibnd), 1, hpsi(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) &
          hw(ibnd) = hw(ibnd) - DBLE( psi(1,ibnd) ) * DBLE ( hpsi(1,ibnd) )
          !
       END DO
       !
       CALL mp_sum( hw(ibnd_start:ibnd_end), intra_bgrp_comm )
       !
       IF ( uspp ) THEN
          !
          DO ibnd = ibnd_start, ibnd_end
             !
             sw(ibnd) = 2._DP * DDOT( 2 * npw, psi(1,ibnd), 1, spsi(1,ibnd), 1 )
             !
             IF ( gstart == 2 ) &
             sw(ibnd) = sw(ibnd) - DBLE( psi(1,ibnd) ) * DBLE ( spsi(1,ibnd) )
             !
          END DO
          !
       ELSE
          !
          DO ibnd = ibnd_start, ibnd_end
             !
             sw(ibnd) = 2._DP * DDOT( 2 * npw, psi(1,ibnd), 1, psi(1,ibnd), 1 )
             !
             IF ( gstart == 2 ) &
             sw(ibnd) = sw(ibnd) - DBLE( psi(1,ibnd) ) * DBLE ( psi(1,ibnd) )
             !
          END DO
          !
       END IF
       !
       CALL mp_sum( sw(ibnd_start:ibnd_end), intra_bgrp_comm )
       !
    END IF
    !
    ! ... Energy eigenvalues
    !
    IF( ANY( sw(ibnd_start:ibnd_end) <= eps16 ) ) &
    CALL errore( ' rrmmdiagg ', ' sw <= 0 ', 1 )
    !
    ew(1:nbnd) = 0._DP
    ew(ibnd_start:ibnd_end) = hw(ibnd_start:ibnd_end) / sw(ibnd_start:ibnd_end)
    !
    CALL mp_sum( ew, inter_bgrp_comm )
    !
    ! ... Check convergence
    !
    WHERE( btype(1:nbnd) == 1 )
       !
       conv(1:nbnd) = ( ( ABS( ew(1:nbnd) - e(1:nbnd) ) < ethr ) )
       !
    ELSEWHERE
       !
       conv(1:nbnd) = ( ( ABS( ew(1:nbnd) - e(1:nbnd) ) < empty_ethr ) )
       !
    END WHERE
    !
    CALL mp_bcast( conv, root_bgrp_id, inter_bgrp_comm )
    !
    notconv = COUNT( .NOT. conv(1:nbnd) )
    !
    e(1:nbnd) = ew(1:nbnd)
    !
    RETURN
    !
  END SUBROUTINE eigenvalues
  !
  !
END SUBROUTINE rrmmdiagg
