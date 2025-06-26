PROGRAM wurtzite
!
!  This program reads the a and c/a of a wurtzite structure, 
!  the atomic coordinates in units of a and writes on output
!  the value of u of the wurtzite structure.
!  Conversely it receives a, c/a and u and writes the cartesian
!  coordinates of the atoms.
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP)  :: a, csua, tau(3,4), u, orig
CHARACTER(LEN=6) :: atm(4)
INTEGER   :: stdin, stdout, na, nat, ichoice

stdin=5
stdout=6
WRITE(stdout,'(5x,"Choose what to do: 1 tau -> u, 2 u -> tau")') 
READ(stdin,*) ichoice
WRITE(stdout,'(5x,"insert a (in a.u.) and c/a of hexagonal lattice")')  
READ(stdin,*) a, csua  
nat=4
IF (ichoice==1) THEN
   DO na=1,nat
      WRITE(stdout,'(5x,"insert the atomic coordinates of atom ",i3)') na  
      READ(stdin,*) atm(na), tau(1,na), tau(2,na), tau(3,na)
      WRITE(stdout,'(a6,3f20.12)') atm(na), tau(1,na), tau(2,na), tau(3,na)
   ENDDO
!
!  Brings z coordinate of atom 1 in the z=0 plane
!
   orig=tau(3,1)
   DO na=1,nat
      tau(3,na)=tau(3,na)-orig
   ENDDO

   u=0.0_DP
   DO na=1,nat
      IF (u==0.0_DP.AND.ABS(tau(3,na)-csua*0.5_DP)>1.D-2.AND. &
         ABS(tau(3,na)-csua*0.5_DP)>1.D-2) THEN
         u=tau(3,na) / csua
         IF (u>0.5_DP) u=u-0.5_DP
      ENDIF
   ENDDO

   WRITE(stdout,'(5x,"celldm(1)=",f15.8)') a
   WRITE(stdout,'(5x,"celldm(3)=",f15.8)') csua
   WRITE(stdout,'(5x,"u=",f15.8)') u
ELSE
   WRITE(stdout,'(5x,"insert u")')  
   READ(stdin,*) u
   WRITE(stdout,'(f16.10)') u

   atm(1)='Zn'
   atm(2)='Zn'
   atm(3)='O '
   atm(4)='O '
   tau=0.0_DP
   tau(1,2)=0.5_DP
   tau(2,2)=-sqrt(3.0_DP) / 6.0_DP
   tau(3,2)= csua/2.0_DP
   tau(3,3)= csua * u 
   tau(1,4)=0.5_DP
   tau(2,4)=-sqrt(3.0_DP) / 6.0_DP
   tau(3,4)= csua*(1.0_DP/2.0_DP+u)

   WRITE(stdout,'(5x,"celldm(1)=",f15.8)') a
   WRITE(stdout,'(5x,"celldm(3)=",f15.8)') csua
ENDIF
!
!  Write on output the atomic coodinates
!
DO na=1,nat
   WRITE(stdout,'(a6,3f20.12)') atm(na), tau(1,na), tau(2,na), tau(3,na)
ENDDO
END PROGRAM wurtzite
