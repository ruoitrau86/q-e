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
  USE mp_bands,  ONLY : inter_bgrp_comm, intra_bgrp_comm, me_bgrp, &
                        root_bgrp, root_bgrp_id, set_bgrp_indices
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
  INTEGER                  :: idiis
  INTEGER                  :: ibnd_start, ibnd_end
  REAL(DP)                 :: empty_ethr
  COMPLEX(DP), ALLOCATABLE :: phi(:,:,:), hphi(:,:,:), sphi(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), spsi(:,:)
  COMPLEX(DP), ALLOCATABLE :: kpsi(:,:), hkpsi(:,:), skpsi(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc(:,:,:), sc(:,:,:)
  REAL(DP),    ALLOCATABLE :: php(:,:), psp(:,:)
  REAL(DP),    ALLOCATABLE :: ew(:), hw(:), sw(:)
  LOGICAL,     ALLOCATABLE :: conv(:)
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
  CALL set_bgrp_indices( nbnd, ibnd_start, ibnd_end )
  !
  ALLOCATE( phi( kdmx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate phi ', ABS(ierr) )
  !
  ALLOCATE( hphi( kdmx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hphi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( sphi( kdmx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate sphi ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( hpsi( kdmx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( spsi( kdmx, nbnd ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate spsi ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( kpsi( kdmx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate kpsi ', ABS(ierr) )
  !
  ALLOCATE( hkpsi( kdmx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hkpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( skpsi( kdmx, nbnd ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate skpsi ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( hc( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hc ', ABS(ierr) )
  !
  ALLOCATE( sc( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate sc ', ABS(ierr) )
  !
  ALLOCATE( php( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate php ', ABS(ierr) )
  !
  ALLOCATE( psp( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate psp ', ABS(ierr) )
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
  hpsi = ZERO
  IF ( uspp ) spsi = ZERO
  !
  kpsi  = ZERO
  hkpsi = ZERO
  IF ( uspp ) skpsi = ZERO
  !
  hc = ZERO
  sc = ZERO
  !
  php = 0.0_DP
  psp = 0.0_DP
  !
  ew   = 0.0_DP
  hw   = 0.0_DP
  sw   = 0.0_DP
  conv = .FALSE.
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
     CALL eigenvalues( )
     !
     IF ( notconv == 0 ) EXIT
     !
  END DO
  !
  ! ... Sort eigenvalues and eigenvectors
  !
  IF ( reorder ) CALL sort_vectors( )
  !
  DEALLOCATE( phi )
  DEALLOCATE( hphi )
  IF ( uspp ) DEALLOCATE( sphi )
  DEALLOCATE( hpsi )
  IF ( uspp ) DEALLOCATE( spsi )
  DEALLOCATE( kpsi )
  DEALLOCATE( hkpsi )
  IF ( uspp ) DEALLOCATE( skpsi )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  DEALLOCATE( php )
  DEALLOCATE( psp )
  DEALLOCATE( ew )
  DEALLOCATE( hw )
  DEALLOCATE( sw )
  DEALLOCATE( conv )
  !
  CALL stop_clock( 'crmmdiagg' )
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
    COMPLEX(DP), ALLOCATABLE :: res1(:,:)
    COMPLEX(DP), ALLOCATABLE :: res2(:)
    COMPLEX(DP), ALLOCATABLE :: vc(:)
    COMPLEX(DP), ALLOCATABLE :: tc(:,:)
    !
    ALLOCATE( res1( kdmx, idiis ) )
    ALLOCATE( res2( kdmx ) )
    IF ( idiis > 1 ) ALLOCATE( vc( idiis ) )
    ALLOCATE( tc( idiis, ibnd_start:ibnd_end ) )
    !
    ! ... Save current wave functions and eigenvalues
    !
    phi (:,ibnd_start:ibnd_end,idiis) = psi (:,ibnd_start:ibnd_end)
    hphi(:,ibnd_start:ibnd_end,idiis) = hpsi(:,ibnd_start:ibnd_end)
    IF ( uspp ) &
    sphi(:,ibnd_start:ibnd_end,idiis) = spsi(:,ibnd_start:ibnd_end)
    php   (ibnd_start:ibnd_end,idiis) = hw    (ibnd_start:ibnd_end)
    psp   (ibnd_start:ibnd_end,idiis) = sw    (ibnd_start:ibnd_end)
    !
    ! ... <R_i|R_j>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       ! ... Residual vectors : |R> = (H - e S) |psi>
       !
       FORALL ( kdiis = 1:idiis ) &
       res1(:,kdiis) = hphi(:,ibnd,kdiis) - php(ibnd,kdiis) * sphi(:,ibnd,kdiis)
       res2(:)       = hphi(:,ibnd,idiis) - php(ibnd,idiis) * sphi(:,ibnd,idiis)
       !
       CALL ZGEMV( 'C', kdim, idiis, ONE, res1(1,1), kdmx, &
                   res2(1), 1, ZERO, hc(1,idiis,ibnd), 1 )
       !
    END DO
    !
    tc(:,:) = hc(:,idiis,:)
    CALL mp_sum( tc, intra_bgrp_comm )
    hc(:,idiis,:) = tc(:,:)
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       hc(idiis,idiis,ibnd) = CMPLX( DBLE( hc(idiis,idiis,ibnd), 0._DP, kind=DP )
       hc(idiis,1:idiis,ibnd) = CONJG( hc(1:idiis,idiis,ibnd) )
       !
    END DO
    !
    ! ... <phi_i| S |phi_j>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( uspp ) THEN
          !
          CALL ZGEMV( 'C', kdim, idiis, ONE, phi(1,ibnd,1), kdmx, &
                      sphi(1,ibnd,idiis), 1, ZERO, sc(1,idiis,ibnd), 1 )
          !
       ELSE
          !
          CALL ZGEMV( 'C', kdim, idiis, ONE, phi(1,ibnd,1), kdmx, &
                      phi(1,ibnd,idiis), 1, ZERO, sc(1,idiis,ibnd), 1 )
          !
       END IF
       !
    END DO
    !
    tc(:,:) = sc(:,idiis,:)
    CALL mp_sum( tc, intra_bgrp_comm )
    sc(:,idiis,:) = tc(:,:)
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       sc(idiis,idiis,ibnd) = CMPLX( DBLE( sc(idiis,idiis,ibnd), 0._DP, kind=DP )
       sc(idiis,1:idiis,ibnd) = CONJG( sc(1:idiis,idiis,ibnd) )
       !
    END DO
    !
    ! ... Update current wave functions
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( idiis > 1 ) THEN
          !
          ! ... solve Rv = eSv
          !
          IF ( me_bgrp == root_bgrp ) CALL diag_diis( ibnd, idiis, vc(:) )
          CALL mp_bcast( vc, root_bgrp, intra_bgrp_comm )
          !
          psi (:,ibnd) = ZERO
          hpsi(:,ibnd) = ZERO
          IF ( uspp ) &
          spsi(:,ibnd) = ZERO
          kpsi(:,ibnd) = ZERO
          !
          DO kdiis = 1, idiis
             !
             psi (:,ibnd) = psi (:,ibnd) + vc(kdiis) * phi (:,ibnd,kdiis)
             hpsi(:,ibnd) = hpsi(:,ibnd) + vc(kdiis) * hphi(:,ibnd,kdiis)
             IF ( uspp ) &
             spsi(:,ibnd) = spsi(:,ibnd) + vc(kdiis) * sphi(:,ibnd,kdiis)
             !
             res2(:) = hphi(:,ibnd,kdiis) - php(ibnd,kdiis) * sphi(:,ibnd,kdiis)
             kpsi(:,ibnd) = kpsi(:,ibnd) + vc(kdiis) * res2(:)
             !
          END DO
          !
       ELSE
          !
          psi (:,ibnd) = phi (:,ibnd,1)
          hpsi(:,ibnd) = hphi(:,ibnd,1)
          IF ( uspp ) &
          spsi(:,ibnd) = sphi(:,ibnd,1)
          kpsi(:,ibnd) = hphi(:,ibnd,1) - php(ibnd,1) * sphi(:,ibnd,1)
          !
       END IF
       !
    END DO
    !
    DEALLOCATE( res1 )
    DEALLOCATE( res2 )
    IF ( idiis > 1 ) DEALLOCATE( vc )
    DEALLOCATE( tc )
    !
    RETURN
    !
  END SUBROUTINE do_diis
  !
  !
  SUBROUTINE diag_diis( ibnd, idiis, vc )
    !
    IMPLICIT NONE
    !
    INTEGER,     INTENT(IN)  :: ibnd
    INTEGER,     INTENT(IN)  :: idiis
    COMPLEX(DP), INTENT(OUT) :: vc(idiis)
    !
    INTEGER                  :: info
    INTEGER                  :: ndim, ndep
    INTEGER                  :: i, imin
    REAL(DP)                 :: emin
    COMPLEX(DP), ALLOCATABLE :: h1(:,:)
    COMPLEX(DP), ALLOCATABLE :: h2(:,:)
    COMPLEX(DP), ALLOCATABLE :: h3(:,:)
    COMPLEX(DP), ALLOCATABLE :: s1(:,:)
    COMPLEX(DP), ALLOCATABLE :: s2(:,:)
    COMPLEX(DP), ALLOCATABLE :: x1(:,:)
    REAL(DP),    ALLOCATABLE :: e1(:)
    INTEGER                  :: nwork
    COMPLEX(DP), ALLOCATABLE :: work(:)
    !
    ndim  = idiis
    nwork = 3 * ndim
    !
    ALLOCATE( h1( ndim, ndim ) )
    ALLOCATE( h2( ndim, ndim ) )
    ALLOCATE( h3( ndim, ndim ) )
    ALLOCATE( s1( ndim, ndim ) )
    ALLOCATE( s2( ndim, ndim ) )
    ALLOCATE( x1( ndim, ndim ) )
    ALLOCATE( e1( ndim ) )
    ALLOCATE( work( nwork ) )
    !
    h1(1:ndim,1:ndim) = hc(1:ndim,1:ndim,ibnd)
    s1(1:ndim,1:ndim) = sc(1:ndim,1:ndim,ibnd)
    !
    CALL ZHEEV( 'V', 'U', ndim, s1, ndim, e1, work, nwork, info )
    !
    IF( info /= 0 ) CALL errore( ' crmmdiagg ', ' cannot solve diis ', ABS(info) )
    !
    ndep = 0
    !
    DO i = 1, ndim
       !
       IF ( e1(i) > eps14 ) THEN
          !
          s2(:,i) = s1(:,i) / SQRT(e1(i))
          !
       ELSE
          !
          ndep = ndep + 1
          !
          s2(:,i) = ZERO
          !
       END IF
       !
    END DO
    !
    IF ( (ndim - ndep) <= 1 ) THEN
       !
       vc        = ZERO
       vc(idiis) = ONE
       !
       GOTO 10
       !
    END IF
    !
    CALL ZGEMM( 'N', 'C', ndim, ndim, ndim, ONE, s2, ndim, s1, ndim, ZERO, x1, ndim )
    !
    CALL ZGEMM( 'N', 'N', ndim, ndim, ndim, ONE, h1, ndim, x1, ndim, ZERO, h2, ndim )
    !
    CALL ZGEMM( 'N', 'N', ndim, ndim, ndim, ONE, x1, ndim, h2, ndim, ZERO, h3, ndim )
    !
    CALL ZHEEV( 'V', 'U', ndim, h3, ndim, e1, work, nwork, info )
    !
    IF( info /= 0 ) CALL errore( ' crmmdiagg ', ' cannot solve diis ', ABS(info) )
    !
    imin = 1
    emin = e1(1)
    !
    DO i = 2, ndim
       !
       IF ( e1(i) < emin ) imin = i
       !
    END DO
    !
    CALL ZGEMV( 'N', ndim, ndim, ONE, x1, ndim, h3(:,imin), 1, ZERO, vc, 1 )
    !
10  DEALLOCATE( h1 )
    DEALLOCATE( h2 )
    DEALLOCATE( h3 )
    DEALLOCATE( s1 )
    DEALLOCATE( s2 )
    DEALLOCATE( x1 )
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
    INTEGER :: ibnd
    !
    ! ... Preconditioning vectors : K (H - e S) |psi>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       ! TODO
       ! TODO
       ! TODO
       !
    END DO
    !
    ! ... Operate the Hamiltonian : H K (H - eS) |psi>
    !
    CALL h_psi( npwx, npw, nbnd, kpsi, hkpsi )
    !
    ! ... Operate the Overlap : S K (H - eS) |psi>
    !
    IF ( uspp ) CALL s_psi( npwx, npw, nbnd, kpsi, skpsi )
    !
    ! ... Line searching for each band
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       ! TODO
       ! TODO
       ! TODO
       !
    END DO
    !
    RETURN
    !
  END SUBROUTINE line_search
  !
  !
  SUBROUTINE eigenvalues( )
    !
    IMPLICIT NONE
    !
    INTEGER            :: ibnd
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
    FORALL ( ibnd = ibnd_start:ibnd_end ) &
    hw(ibnd) = ZDOTC( kdim, psi(1,ibnd), 1, hpsi(1,ibnd), 1 )
    !
    CALL mp_sum( hw(ibnd_start:ibnd_end), intra_bgrp_comm )
    !
    IF ( uspp ) THEN
       !
       FORALL ( ibnd = ibnd_start:ibnd_end ) &
       sw(ibnd) = ZDOTC( kdim, psi(1,ibnd), 1, spsi(1,ibnd), 1 )
       !
    ELSE
       !
       FORALL ( ibnd = ibnd_start:ibnd_end ) &
       sw(ibnd) = ZDOTC( kdim, psi(1,ibnd), 1, psi(1,ibnd), 1 )
       !
    END IF
    !
    CALL mp_sum( sw(ibnd_start:ibnd_end), intra_bgrp_comm )
    !
    ew(1:nbnd) = 0.0_DP
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
  SUBROUTINE sort_vectors( )
    !
    IMPLICIT NONE
    !
    INTEGER  :: ibnd
    !
    ! TODO
    ! TODO
    ! TODO
    !
    RETURN
    !
  END SUBROUTINE sort_vectors
  !
  !
END SUBROUTINE crmmdiagg
