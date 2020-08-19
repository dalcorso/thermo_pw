!
! Copyright (C) 2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE collect_info
!
!  This is an auxiliary module with a structure, the collect_info_type, that
!  has sufficient variables to keep in memory the information on the 
!  grid of q points and irreducible representations for a single geometry.
!  It saves also the status of the disk with the information on the
!  representations that have been already calculated.
!  Presently it has four routines:
!
!  init_collect_info copies in the structure the variables of the phonon
!  copy_collect_info copies in the phonon variables the information
!                    of the structure
!  comm_collect_info sends to all images the available information
!                    so that any image can do the task that has been
!                    calculated by the phonon for another image.
!  destroy_collect_info deallocates the memory allocated by init_collect_info. 
!
!
USE kinds, ONLY : DP

SAVE
PRIVATE

TYPE :: collect_info_type
    !
    INTEGER :: nqs                   ! number of q points
    !
    INTEGER, ALLOCATABLE :: irr_iq(:)   ! number of irrep for each q
    !
    INTEGER, ALLOCATABLE :: comp_irr_iq(:,:,:) ! irrep, iq, image
    !                                      equal one if must be computed
    INTEGER, ALLOCATABLE :: done_irr_iq(:,:,:) ! irrep, iq, image
    !                                      equal one if already computed
    INTEGER, ALLOCATABLE :: comp_iq(:,:)     ! iq, image
    !                                      equal one if must be computed 
    INTEGER, ALLOCATABLE :: done_iq(:,:)     ! iq, image
    !                                      equal one if already computed
END TYPE collect_info_type

PUBLIC collect_info_type, init_collect_info, copy_collect_info, &
       comm_collect_info, destroy_collect_info_type

CONTAINS

!----------------------------------------------------------------------------
   SUBROUTINE init_collect_info(info, nqs, nat, nima, pos, comp_irr_iq, &
                                done_irr_iq, comp_iq, done_iq, irr_iq)
!----------------------------------------------------------------------------
!
!  This routine copies the variables of the phonon into the info
!  structure allocating the necessary variables of the structure.
!
!  Every image writes in a different position (given by pos) in the
!  array, so that they can comunicate to each other who does what.
!  nima is the number of images used to set the dimension of the array.
!
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: nqs, nat, nima, pos, irr_iq(nqs)
   LOGICAL, INTENT(IN) :: comp_irr_iq(0:3*nat,nqs), done_irr_iq(0:3*nat,nqs),&
                          comp_iq(nqs), done_iq(nqs)
   TYPE(collect_info_type), INTENT(INOUT) :: info

   INTEGER :: iq, irr
!
!  these variables are the same for all images
!
   info%nqs=nqs
   ALLOCATE(info%irr_iq(nqs))
   info%irr_iq(1:nqs)= irr_iq(1:nqs)
!
!  these instead differ and are set in different positions (pos)
!  first allocate space to contain the info of each image
!
   ALLOCATE(info%comp_irr_iq(0:3*nat,nqs,nima))
   ALLOCATE(info%done_irr_iq(0:3*nat,nqs,nima))
   ALLOCATE(info%comp_iq(nqs,nima))
   ALLOCATE(info%done_iq(nqs,nima))
!
!  sets all to zero, so the variables can be broadcasted with an mp_sum
!
   info%comp_irr_iq=0
   info%done_irr_iq=0
   info%comp_iq=0
   info%done_iq=0
!
!  set to 1 the iq and irrep of this image.
!
   DO iq=1,nqs
      DO irr=0, irr_iq(iq)
         IF (comp_irr_iq(irr,iq)) &
            info%comp_irr_iq(irr,iq,pos)=1
         IF (done_irr_iq(irr,iq)) &
            info%done_irr_iq(irr,iq,pos)=1
      ENDDO
      IF (comp_iq(iq)) info%comp_iq(iq,pos)=1
      IF (done_iq(iq)) info%done_iq(iq,pos)=1
   ENDDO

   RETURN
   END SUBROUTINE init_collect_info

!----------------------------------------------------------------------------
   SUBROUTINE copy_collect_info(info, nqs, nat, nima, pos, comp_irr_iq, &
                                done_irr_iq, comp_iq, done_iq, irr_iq)
!----------------------------------------------------------------------------
!
!  This routine is the inverse of init_collect_info. It uses the information
!  saved in the info structure to set the phonon variables.
!
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: nqs, nat, nima, pos, irr_iq(nqs)
   LOGICAL, INTENT(INOUT) :: comp_irr_iq(0:3*nat,nqs), &
                    done_irr_iq(0:3*nat, nqs), comp_iq(nqs), done_iq(nqs)
   TYPE(collect_info_type), INTENT(IN) :: info

   INTEGER :: iq, irr

   comp_irr_iq=.FALSE.
   comp_iq=.FALSE.
   done_irr_iq=.FALSE.
   done_iq=.FALSE.
!
!  Note that the structure uses 0 or 1, while the phonon variables are logicals
!  This is because there is no mp_sum routines for logicals.
!
   DO iq=1,nqs
      DO irr=0, irr_iq(iq)
         IF (info%comp_irr_iq(irr,iq,pos)==1) comp_irr_iq(irr,iq)=.TRUE.
         IF (info%done_irr_iq(irr,iq,pos)==1) done_irr_iq(irr,iq)=.TRUE.
      ENDDO
      IF (info%comp_iq(iq,pos)==1) comp_iq(iq)=.TRUE.
      IF (info%done_iq(iq,pos)==1) done_iq(iq)=.TRUE.
   ENDDO
   DO iq=1,nqs
      DO irr=0, irr_iq(iq)
         IF (done_irr_iq(irr,iq)) comp_irr_iq(irr,iq)=.FALSE.
      ENDDO
      IF (done_iq(iq)) comp_iq(iq)=.FALSE.
   ENDDO

   RETURN
   END SUBROUTINE copy_collect_info

!----------------------------------------------------------------------------
   SUBROUTINE comm_collect_info(info, comm)
!----------------------------------------------------------------------------
!
!  This routine comunicates the variables of the info structure to all 
!  processors of a communicator group. 
!
   USE mp, ONLY : mp_sum

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: comm
   TYPE(collect_info_type) :: info

   CALL mp_sum(info%comp_irr_iq, comm)
   CALL mp_sum(info%done_irr_iq, comm)
   CALL mp_sum(info%comp_iq, comm)
   CALL mp_sum(info%done_iq, comm)

   RETURN
   END SUBROUTINE comm_collect_info

!----------------------------------------------------------------------------
   SUBROUTINE destroy_collect_info_type(info)
!----------------------------------------------------------------------------
!
!  This routine deallocate the space allocated by init_collect_info
!
   IMPLICIT NONE
   TYPE(collect_info_type) :: info

   IF (ALLOCATED(info%irr_iq))      DEALLOCATE(info%irr_iq)
   IF (ALLOCATED(info%comp_irr_iq)) DEALLOCATE(info%comp_irr_iq)
   IF (ALLOCATED(info%done_irr_iq)) DEALLOCATE(info%done_irr_iq)
   IF (ALLOCATED(info%comp_iq))     DEALLOCATE(info%comp_iq)
   IF (ALLOCATED(info%done_iq))     DEALLOCATE(info%done_iq)

   RETURN
   END SUBROUTINE destroy_collect_info_type

END MODULE collect_info
