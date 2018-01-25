  SUBROUTINE allocate_lanczos()

  USE klist,                ONLY : nks
  USE uspp,                 ONLY : okvan, nkb
  USE wvfct,                ONLY : npwx, nbnd
  USE noncollin_module,     ONLY : noncolin, npol
  USE qpoint,               ONLY : nksq
  USE control_lr,           ONLY : lgamma
  USE lrus,                 ONLY : bbg
  USE control_flags,        ONLY : gamma_only

  USE lr_global,    ONLY : evc0, sevq0, evq0, d0psi, d0psi2, size_evc1, rpert
  USE lr_lanczos,   ONLY : evc1, evc1_new, evc1_old, sevc1, beta_store, &
                           gamma_store, zeta_store, beta_store_ext,     &
                           gamma_store_ext, lanczos_steps, lanczos_steps_ext, &
                           bbk, bbnc

  IMPLICIT NONE

  INTEGER :: ncopy

  rpert=1
  IF (lgamma) rpert=3
  ncopy=2
  size_evc1= npwx*npol*nbnd*nksq*rpert

  ALLOCATE(evc0(npwx*npol,nbnd,nksq))
  IF (.NOT.lgamma) THEN
     ALLOCATE(evq0(npwx*npol,nbnd,nksq))
  ELSE
     evq0 => evc0
  ENDIF
  IF (okvan) THEN
     ALLOCATE(sevq0(npwx*npol,nbnd,nksq))
  ELSE
     sevq0 => evq0
  ENDIF

  ALLOCATE(evc1_old(npwx*npol,nbnd,nksq*rpert,ncopy))
  ALLOCATE(evc1(npwx*npol,nbnd,nksq*rpert,ncopy))
  ALLOCATE(evc1_new(npwx*npol,nbnd,nksq*rpert,ncopy))
  !
  ALLOCATE(sevc1(npwx*npol,nbnd,nksq*rpert,ncopy))
  !
  ALLOCATE(d0psi(npwx*npol,nbnd,nksq,rpert))

  evc0(:,:,:)       = (0.0d0,0.0d0)
  evq0(:,:,:)       = (0.0d0,0.0d0)
  sevq0(:,:,:)      = (0.0d0,0.0d0)
  evc1_old(:,:,:,:) = (0.0d0,0.0d0)
  evc1(:,:,:,:)     = (0.0d0,0.0d0)
  evc1_new(:,:,:,:) = (0.0d0,0.0d0)
  sevc1(:,:,:,:)    = (0.0d0,0.0d0)
  d0psi(:,:,:,:)    = (0.0d0,0.0d0)

  IF (.NOT.lgamma) THEN
     ALLOCATE(d0psi2(npwx*npol,nbnd,nksq))
     d0psi2(:,:,:) = (0.0d0,0.0d0)
  ENDIF

  ALLOCATE(beta_store(lanczos_steps))
  ALLOCATE(gamma_store(lanczos_steps))

  ALLOCATE(beta_store_ext(lanczos_steps_ext))
  ALLOCATE(gamma_store_ext(lanczos_steps_ext))
  ALLOCATE(zeta_store(rpert,rpert,lanczos_steps))

  IF (okvan) THEN
!
!  Allocate space for the coefficients of S^-1
!
     IF (gamma_only) THEN
        ALLOCATE(bbg(nkb, nkb))
     ELSE
        IF (noncolin) THEN
           ALLOCATE(bbnc(nkb*npol, nkb*npol, nksq))
        ELSE
           ALLOCATE(bbk(nkb, nkb, nksq))
        ENDIF
     ENDIF
  ENDIF

  RETURN
  END SUBROUTINE allocate_lanczos

  SUBROUTINE deallocate_lanczos()

  USE lr_global,  ONLY : evc0, sevq0, evq0, d0psi, d0psi2, size_evc1, rpert
  USE lr_lanczos, ONLY : evc1, evc1_new, evc1_old, sevc1, beta_store, &
                         gamma_store, zeta_store, beta_store_ext, &
                         gamma_store_ext, bbk, bbnc
  USE lrus,       ONLY : bbg

  IMPLICIT NONE

  IF (ASSOCIATED(sevq0,evq0)) THEN
     NULLIFY(sevq0)
  ELSEIF(ASSOCIATED(sevq0)) THEN
     DEALLOCATE(sevq0)
  ENDIF

  IF (ASSOCIATED(evq0,evc0)) THEN
     NULLIFY(evq0)
  ELSEIF(ASSOCIATED(evq0)) THEN
     DEALLOCATE(evq0)
  ENDIF

  IF (ALLOCATED(evc0)) DEALLOCATE(evc0)

  IF (ALLOCATED(evc1_old)) DEALLOCATE(evc1_old)
  IF (ALLOCATED(evc1)) DEALLOCATE(evc1)
  IF (ALLOCATED(evc1_new)) DEALLOCATE(evc1_new)
   
  IF (ALLOCATED(sevc1)) DEALLOCATE(sevc1)
   
  IF (ALLOCATED(d0psi)) DEALLOCATE(d0psi)
  IF (ALLOCATED(d0psi2)) DEALLOCATE(d0psi2)

  IF (ALLOCATED(beta_store_ext)) DEALLOCATE(beta_store_ext)
  IF (ALLOCATED(gamma_store_ext)) DEALLOCATE(gamma_store_ext)
  IF (ALLOCATED(zeta_store)) DEALLOCATE(zeta_store)

  IF (ALLOCATED(bbg)) DEALLOCATE(bbg)
  IF (ALLOCATED(bbk)) DEALLOCATE(bbk)
  IF (ALLOCATED(bbnc)) DEALLOCATE(bbnc)

  RETURN
  END SUBROUTINE deallocate_lanczos
