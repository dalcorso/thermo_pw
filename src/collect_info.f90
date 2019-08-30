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
!  Presently it has three routines:
!  init_collect_info that copy in the structure the variables read by the
!                     phonon
!  comm_collect_info that sends to all images the available information
!  destroy_collect_info that deallocate the memory allocated by 
!                       init_collect_info. It has to be called to free the
!                       memory.
!
USE kinds, ONLY : DP

SAVE
PRIVATE

TYPE :: collect_info_type
    !
    INTEGER :: nqs
    !
    INTEGER, ALLOCATABLE :: irr_iq(:)
    !
    INTEGER, ALLOCATABLE :: comp_irr_iq(:,:,:)
    !
    INTEGER, ALLOCATABLE :: done_irr_iq(:,:,:)
    !
    INTEGER, ALLOCATABLE :: comp_iq(:,:)
    
    INTEGER, ALLOCATABLE :: done_iq(:,:)
    !
END TYPE collect_info_type

PUBLIC collect_info_type, init_collect_info, copy_collect_info, &
       comm_collect_info, destroy_collect_info_type

CONTAINS

!----------------------------------------------------------------------------
   SUBROUTINE init_collect_info(info, nqs, nat, nima, pos, comp_irr_iq, &
                                done_irr_iq, comp_iq, done_iq, irr_iq)
!----------------------------------------------------------------------------
!
!  This routine copy the variables of the phonon into the info
!  structure allocating the necessary variables of the structure.
!
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: nqs, nat, nima, pos, irr_iq(nqs)
   LOGICAL, INTENT(IN) :: comp_irr_iq(0:3*nat,nqs), done_irr_iq(0:3*nat, nqs), &
                          comp_iq(nqs), done_iq(nqs)
   TYPE(collect_info_type), INTENT(INOUT) :: info

   INTEGER :: iq, irr

   info%nqs=nqs
   ALLOCATE(info%irr_iq(nqs))
   info%irr_iq(1:nqs)= irr_iq(1:nqs)

   ALLOCATE(info%comp_irr_iq(0:3*nat,nqs,nima))
   ALLOCATE(info%done_irr_iq(0:3*nat,nqs,nima))
   ALLOCATE(info%comp_iq(nqs,nima))
   ALLOCATE(info%done_iq(nqs,nima))

   info%comp_irr_iq=0
   info%done_irr_iq=0
   info%comp_iq=0
   info%done_iq=0

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
!  This routine comunicate the variables of the info structure to all 
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

   IMPLICIT NONE
   TYPE(collect_info_type) :: info

   IF (ALLOCATED(info%irr_iq)) DEALLOCATE(info%irr_iq)
   IF (ALLOCATED(info%comp_irr_iq)) DEALLOCATE(info%comp_irr_iq)
   IF (ALLOCATED(info%done_irr_iq)) DEALLOCATE(info%done_irr_iq)
   IF (ALLOCATED(info%comp_iq)) DEALLOCATE(info%comp_iq)
   IF (ALLOCATED(info%done_iq)) DEALLOCATE(info%done_iq)

   RETURN
   END SUBROUTINE destroy_collect_info_type

END MODULE collect_info
