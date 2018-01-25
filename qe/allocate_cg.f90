!
! Copyright (C) 2018 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

  SUBROUTINE allocate_cg(npe, nc)

  USE klist,                ONLY : nks
  USE uspp,                 ONLY : okvan, nkb
  USE wvfct,                ONLY : npwx, nbnd
  USE noncollin_module,     ONLY : npol
  USE qpoint,               ONLY : nksq
  USE control_lr,           ONLY : lgamma

  USE lr_global,    ONLY : evc0, sevq0, evq0, d0psi, d0psi2, size_evc1, rpert
  USE lr_cg,        ONLY : evc1, dir, pres, res, prec_vec, dir_new
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npe, nc

  rpert=npe
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
     IF (lgamma) THEN
        sevq0 => evc0
     ELSE
        sevq0 => evq0
     ENDIF
  ENDIF

  ALLOCATE(evc1(npwx*npol,nbnd,nksq*rpert,1))
  ALLOCATE(dir(npwx*npol,nbnd,nksq*rpert,nc))
  ALLOCATE(dir_new(npwx*npol,nbnd,nksq*rpert,nc))
  ALLOCATE(res(npwx*npol,nbnd,nksq*rpert,nc))
  ALLOCATE(pres(npwx*npol,nbnd,nksq*rpert,nc))
  ALLOCATE(d0psi(npwx*npol,nbnd,nksq,rpert))
  ALLOCATE(prec_vec(npwx*npol,nbnd,nksq))

  evc0(:,:,:)       = (0.0d0,0.0d0)
  evq0(:,:,:)       = (0.0d0,0.0d0)
  sevq0(:,:,:)      = (0.0d0,0.0d0)
  evc1(:,:,:,:)     = (0.0d0,0.0d0)
  dir(:,:,:,:)      = (0.0d0,0.0d0)
  dir_new(:,:,:,:)  = (0.0d0,0.0d0)
  res(:,:,:,:)      = (0.0d0,0.0d0)
  pres(:,:,:,:)     = (0.0d0,0.0d0)
  prec_vec(:,:,:)  = 0.0d0
  d0psi(:,:,:,:)    = (0.0d0,0.0d0)

  IF (.NOT.lgamma) THEN
     ALLOCATE(d0psi2(npwx*npol,nbnd,nksq))
     d0psi2(:,:,:) = (0.0d0,0.0d0)
  ENDIF


  RETURN
  END SUBROUTINE allocate_cg

  SUBROUTINE deallocate_cg()

  USE lr_global,  ONLY : evc0, evq0, sevq0, d0psi, d0psi2
  USE lr_cg,      ONLY : evc1, dir, res, pres, prec_vec, dir_new

  IMPLICIT NONE

  IF (ASSOCIATED(sevq0, evq0)) THEN 
     NULLIFY(sevq0)
  ELSEIF (ASSOCIATED(sevq0)) THEN
     DEALLOCATE(sevq0)
  ENDIF
  IF (ASSOCIATED(evq0,evc0)) THEN
     NULLIFY(evq0)
  ELSEIF(ASSOCIATED(evq0)) THEN
     DEALLOCATE(evq0)
  ENDIF
  IF (ALLOCATED(evc0))    DEALLOCATE(evc0)

  IF (ALLOCATED(evc1))    DEALLOCATE(evc1)
  IF (ALLOCATED(res))     DEALLOCATE(res)
  IF (ALLOCATED(pres))    DEALLOCATE(pres)
  IF (ALLOCATED(dir))     DEALLOCATE(dir)
  IF (ALLOCATED(dir_new)) DEALLOCATE(dir_new)
  IF (ALLOCATED(prec_vec)) DEALLOCATE(prec_vec)
   
  IF (ALLOCATED(d0psi))   DEALLOCATE(d0psi)
  IF (ALLOCATED(d0psi2))  DEALLOCATE(d0psi2)

  RETURN
  END SUBROUTINE deallocate_cg
