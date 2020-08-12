!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine init_representations_tpw()
  !-----------------------------------------------------------------------
  !
  !  This subroutine initializes the modes of all irreducible representations
  !  for all q points. It writes the files patterns.#q.xml in the outdir 
  !  directory. It is used by unrecovered  phonon runs. The small group of 
  !  q must be calculated for each q. Note that all images receives the 
  !  same modes calculated by the root processor and save them on file. 
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat
  USE cell_base,     ONLY : at, bg
  USE io_global,     ONLY : stdout
  USE symm_base,     ONLY : nsym, sr, irt, time_reversal, t_rev, s, sname
  USE control_ph,    ONLY : search_sym, current_iq, u_from_file, &
                            search_sym_save
  USE modes,         ONLY : u, npert, nirr, nmodes, name_rap_mode, &
                            num_rap_mode
  USE disp,          ONLY : x_q, nqs, lgamma_iq
  USE cryst_ph,      ONLY : magnetic_sym
  USE ph_restart,    ONLY : ph_writefile
  USE control_flags, ONLY : modenum, noinv
  USE mp,            ONLY : mp_bcast
  USE mp_world,      ONLY : root, world_comm

  USE lr_symm_base,  ONLY : gi, gimq, irotmq, minus_q, nsymq, invsymq, rtau
  USE qpoint,        ONLY : xq
  USE control_lr,    ONLY : lgamma

  implicit none

  integer ::  isym, irr, iq
  ! counters
  LOGICAL, EXTERNAL :: symmorphic_or_nzb
  integer :: ierr

  call start_clock ('init_rep')

  allocate (rtau ( 3, 48, nat))
  allocate (u ( 3 * nat, 3 * nat))
  allocate (name_rap_mode( 3 * nat))
  allocate (num_rap_mode( 3 * nat))
  allocate (npert ( 3 * nat))

  u_from_file=.FALSE.
  !
  ! allocate and calculate rtau, the rotated position of each atom
  !
  nmodes = 3 * nat
  minus_q = (modenum .eq. 0)
  IF ( .not. time_reversal ) minus_q = .false.
  ! if minus_q=.t. set_irr will search for Sq=-q+G symmetry.
  ! On output minus_q=.t. if such a symmetry has been found
  DO iq=1, nqs
     xq(1:3)  = x_q(1:3,iq)
     lgamma = lgamma_iq(iq)
!
!    search for the small group of q
!
     CALL set_small_group_of_q(nsymq,invsymq,minus_q)
!
!    calculate rtau with the new symmetry order
!
     CALL sgam_lr (at, bg, nsym, s, irt, tau, rtau, nat)
!
!    and calculate the vectors G associated to the symmetry Sq = q + G
!    if minus_q is true calculate also irotmq and the G associated to Sq=-g+G
!
     CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
!
!    Check if we can search symmetry for this q point
!
     search_sym = search_sym_save .AND. symmorphic_or_nzb()
     num_rap_mode=-1
     name_rap_mode=' '
     IF (search_sym) CALL prepare_sym_analysis_tpw(nsymq,sr,sname,&
                                                         t_rev,magnetic_sym)

     CALL find_irrep_tpw()
!
!  Only the modes calculated by node zero are sent to all images
!
     CALL mp_bcast (u, root, world_comm)
     CALL mp_bcast (nsymq, root, world_comm)
     CALL mp_bcast (npert, root, world_comm)
     CALL mp_bcast (nirr, root, world_comm)
     CALL mp_bcast (name_rap_mode, root, world_comm)
     CALL mp_bcast (num_rap_mode, root, world_comm)

     CALL ph_writefile('data_u',iq,0,ierr)
  ENDDO
  u_from_file=.TRUE.
  search_sym=search_sym_save

  DEALLOCATE (rtau)
  DEALLOCATE (u)
  DEALLOCATE (num_rap_mode)
  DEALLOCATE (name_rap_mode)
  DEALLOCATE (npert)

  CALL stop_clock ('init_rep')
  RETURN
END SUBROUTINE init_representations_tpw

