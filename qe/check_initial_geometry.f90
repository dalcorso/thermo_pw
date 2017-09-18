!
! Copyright (C) 2012-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE check_initial_geometry(auxdyn)
  !-----------------------------------------------------------------------
  !
  ! This routine checks the initial status of the phonon run and sets
  ! the variables that control the run, dealing with the image
  ! and GRID parallelization features of the phonon code.
  ! 
  ! The size of the grid is determined by the following variables:
  ! nqs : the number of q points
  ! x_q : the coordinates of the q points
  !
  ! nfs : the number of imaginary frequencies
  ! fiu : which frequencies 
  !
  ! The flags that control which tensors to calculate
  !
  ! In a recover calculation the q grid variables are already known, 
  ! read from file in phq_readin. In a calculation starting from
  ! scratch this routine sets them. The frequencies variables and the
  ! tensors flags are read from input. 
  ! The amount of work to do for each representation of each q
  ! point depends on the size of the representation and the 
  ! order of the small group of q. In a recover calculation
  ! these information are on file, when recover=.false. this
  ! routine writes the modes and their degeneration on files 
  ! and calculates the order of the small group of q. The following
  ! variables are set
  !
  ! irr_iq : for each q point how many irreducible representations
  ! npert_irr_iq : how many perturbation per representation and per q
  ! nsymq_iq : the order of the small group of q for each q
  !
  ! The following variables are set by this routine on the basis of
  ! start_irr, last_irr, start_iq, last_iq, OR of modenum, OR of ifat and 
  ! atomo:
  !
  ! comp_iq : =.TRUE. if the q point is calculated in this run
  ! comp_irr_iq : =.TRUE. if the representation is calculated in this run
  ! comp_iu : =.TRUE. if this frequency is calculated in this run
  !                   NB: start_iu, last_iu is not yet programmed
  ! 
  ! After knowing this info the routine divides the total work among
  ! the images (when nimage > 1) INDEPENDENTLY of what has been already
  ! calculated and is available on file.
  !
  ! Then, when recover=.true., the routine looks on files for pieces
  ! already calculated and sets the array
  !
  ! done_irr_iq : =.TRUE. if the representation has been already calculated
  ! done_iq : =.TRUE. if the q point has been already calculated
  ! done_iu : =.TRUE. if already calculated
  ! done_bands_iq : .TRUE. if the bands for the q point are on file.
  !
  ! If recover=.false. all these array are initialized to .false.
  !
  !
  USE io_global,       ONLY : stdout
  USE ions_base,       ONLY : nat
  USE io_files,        ONLY : tmp_dir
  USE disp,            ONLY : nqs, x_q, lgamma_iq
  USE qpoint,          ONLY : xq
  USE control_lr,      ONLY : lgamma
  USE output,          ONLY : fildyn
  USE control_ph,      ONLY : ldisp, recover, rec_code, &
                              start_q, last_q, current_iq, tmp_dir_ph, &
                              ext_recover, ext_restart, qplot, &
                              done_zeu, done_start_zstar, done_epsil, &
                              done_zue, always_run
  USE control_qe,      ONLY : use_ph_images, tcollect_all
  USE ph_restart,      ONLY : check_directory_phsave, check_available_bands,&
                              allocate_grid_variables, ph_writefile
  USE freq_ph,         ONLY : current_iu
  USE mp_images,       ONLY : nimage, intra_image_comm
  USE mp,              ONLY : mp_bcast
  USE mp_global,       ONLY : mp_global_end
  USE el_phon,         ONLY : elph_mat
  ! YAMBO >
  USE YAMBO,           ONLY : elph_yambo,dvscf_yambo
  ! YAMBO <
  !
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=256) :: auxdyn, filename
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  INTEGER :: iq, iq_start, ierr
  !
  tmp_dir=tmp_dir_ph
  !
  ! If this not a recover run, we generate the q mesh. Otherwise at this
  ! point the code has read the q mesh from the files contained in 
  ! prefix.phsave
  !
  IF (.NOT.recover.AND..NOT.ALLOCATED(x_q)) THEN
     !
     ! recover file not found or not looked for
     !
     current_iu=1
     current_iq=1
     IF (ldisp) THEN
        !
        ! ... Calculate the q-points for the dispersion
        !
        IF(elph_mat) then
           CALL q_points_wannier()
        ELSE
           IF (.NOT. qplot) CALL q_points_tpw()
        END IF
        !
        ! YAMBO >
     ELSE IF (.NOT.elph_yambo .AND. .NOT. dvscf_yambo) then
        ! YAMBO <
        !
        nqs = 1
        last_q = 1
        ALLOCATE(x_q(3,1))
        ALLOCATE(lgamma_iq(1))
        x_q(:,1)=xq(:)
        lgamma_iq(1)=lgamma
        !
     END IF
     !
     !   Save the mesh of q and the control flags on file
     !
     CALL ph_writefile('init',0,0,ierr)
     !
     !   Initialize the representations and write them on file.
     !
     CALL init_representations()
     !
  ENDIF

  IF (last_q<1.or.last_q>nqs) last_q=nqs
  IF (start_q<1.or.start_q>last_q) call errore('check_initial_status', &
     'wrong start_q',1)
!
!  now we allocate the variables needed to describe the grid
!
  CALL allocate_grid_variables()
!
!  This routine assumes that the modes are on file, either written by 
!  init_representation or written by a previous run. It takes care
!  of dealing with start_irr, last_irr flags and ifat or modenum
!  restricted  computation, moreover it sets the size of each 
!  representation and the size of the small group of q for each point.
!
  CALL initialize_grid_variables()
!
! If there are more than one image, divide the work among the images
!
  IF (nimage > 1 .AND. use_ph_images .AND. .NOT. tcollect_all ) &
                                           CALL image_q_irr_tpw()
!
  IF (recover) THEN
!
! ... Checking the status of the calculation
!
!  sets which q point and representations have been already calculated
!
     CALL check_directory_phsave()
  ELSE
     done_zeu=.FALSE.
     done_start_zstar=.FALSE.
     done_epsil=.FALSE.
     done_zue=.FALSE.
  ENDIF
  !
  auxdyn = fildyn
  RETURN
  END SUBROUTINE check_initial_geometry
