!
! Copyright (C) Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

SUBROUTINE realspace_grid_init_tpw( dfft, at, bg, gcutm, fft_fact )
    !
    ! ... Sets optimal values for dfft%nr[123] and dfft%nr[123]x
    ! ... If fft_fact is present, force nr[123] to be multiple of fft_fac([123])
    !
    USE kinds, ONLY : DP
    USE fft_support, only: good_fft_dimension, good_fft_order
    USE fft_types, ONLY : fft_type_descriptor
    USE io_global, ONLY : stdout
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: at(3,3), bg(3,3)
    REAL(DP), INTENT(IN) :: gcutm
    INTEGER, INTENT(IN)  :: fft_fact(3)
    TYPE(fft_type_descriptor), INTENT(INOUT) :: dfft
    !
    IF( dfft%nr1 == 0 .OR. dfft%nr2 == 0 .OR. dfft%nr3 == 0 ) THEN
      !
      ! ... calculate the size of the real-space dense grid for FFT
      ! ... first, an estimate of nr1,nr2,nr3, based on the max values
      ! ... of n_i indices in:   G = i*b_1 + j*b_2 + k*b_3
      ! ... We use G*a_i = n_i => n_i .le. |Gmax||a_i|
      !
      dfft%nr1 = int ( sqrt (gcutm) * &
            sqrt (at(1, 1)**2 + at(2, 1)**2 + at(3, 1)**2) ) + 1
      dfft%nr2 = int ( sqrt (gcutm) * &
            sqrt (at(1, 2)**2 + at(2, 2)**2 + at(3, 2)**2) ) + 1
      dfft%nr3 = int ( sqrt (gcutm) * &
            sqrt (at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2) ) + 1
      !
      CALL grid_set_tpw( dfft, bg, gcutm, dfft%nr1, dfft%nr2, dfft%nr3 )
      !
    ELSE
       WRITE( stdout, '( /, 3X,"Info: using nr1, nr2, nr3 values from input" )' )
    END IF

    dfft%nr1 = good_fft_order( dfft%nr1, fft_fact(1) )
    dfft%nr2 = good_fft_order( dfft%nr2, fft_fact(2) )
    dfft%nr3 = good_fft_order( dfft%nr3, fft_fact(3) )

    dfft%nr1x  = good_fft_dimension( dfft%nr1 )
    dfft%nr2x  = dfft%nr2
    dfft%nr3x  = good_fft_dimension( dfft%nr3 )

END SUBROUTINE realspace_grid_init_tpw

   SUBROUTINE grid_set_tpw( dfft, bg, gcut, nr1, nr2, nr3 )

!  this routine returns in nr1, nr2, nr3 the minimal 3D real-space FFT 
!  grid required to fit the G-vector sphere with G^2 <= gcut
!  On input, nr1,nr2,nr3 must be set to values that match or exceed
!  the largest i,j,k (Miller) indices in G(i,j,k) = i*b1 + j*b2 + k*b3
!  ----------------------------------------------

      USE kinds, ONLY : DP
      USE fft_types, ONLY : fft_type_descriptor
      IMPLICIT NONE

#if defined(__MPI)
      INCLUDE 'mpif.h'
#endif

! ... declare arguments
      TYPE(fft_type_descriptor), INTENT(IN) :: dfft
      INTEGER, INTENT(INOUT) :: nr1, nr2, nr3
      REAL(DP), INTENT(IN) :: bg(3,3), gcut

! ... declare other variables
      INTEGER :: i, j, k, nr, nb(3)
      REAL(DP) :: gsq, g(3)

!  ----------------------------------------------

      nb     = 0

! ... calculate moduli of G vectors and the range of indices where
! ... |G|^2 < gcut (in parallel whenever possible)

      DO k = -nr3, nr3
        !
        ! ... me_image = processor number, starting from 0
        !
        IF( MOD( k + nr3, dfft%nproc ) == dfft%mype ) THEN
          DO j = -nr2, nr2
            DO i = -nr1, nr1

              g( 1 ) = DBLE(i)*bg(1,1) + DBLE(j)*bg(1,2) + DBLE(k)*bg(1,3)
              g( 2 ) = DBLE(i)*bg(2,1) + DBLE(j)*bg(2,2) + DBLE(k)*bg(2,3)
              g( 3 ) = DBLE(i)*bg(3,1) + DBLE(j)*bg(3,2) + DBLE(k)*bg(3,3)

! ...         calculate modulus

              gsq =  g( 1 )**2 + g( 2 )**2 + g( 3 )**2 

              IF( gsq < gcut ) THEN

! ...           calculate maximum index
                nb(1) = MAX( nb(1), ABS( i ) )
                nb(2) = MAX( nb(2), ABS( j ) )
                nb(3) = MAX( nb(3), ABS( k ) )
              END IF

            END DO
          END DO
        END IF
      END DO

#ifdef __MPI
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, nb, 3, MPI_INTEGER, MPI_MAX, dfft%comm, i )
#endif

! ... the size of the required (3-dimensional) matrix depends on the
! ... maximum indices. Note that the following choice is slightly
! ... "small": 2*nb+2 would be needed in order to guarantee that the
! ...  sphere in G-space never overlaps its periodic image

      nr1 = 2 * nb(1) + 1
      nr2 = 2 * nb(2) + 1
      nr3 = 2 * nb(3) + 1

      RETURN
   
   END SUBROUTINE grid_set_tpw

