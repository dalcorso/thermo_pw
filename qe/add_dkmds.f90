! 
! Copyright (C) 2012-2017 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
subroutine add_dkmds_tpw(ik,uact,jpol,becp2,alphadk,bedp,alphapp,weight,zstar)
  !--------=========-------------------------------------------------------
  !
  ! This subroutine calculates a part of the contribution to zstar 
  ! that does not depend on the perturbed wavefunctions. The other
  ! part is calculated in compute_drhous.
  ! It assumes that the variable dpqq has been set. In the noncollinear
  ! and spin_orbit case the variable dpqq_so must be set. 
  !

  USE kinds, ONLY : DP
  USE constants, ONLY : eps12
  USE spin_orb, ONLY : lspinorb
  USE uspp, ONLY : nkb, qq_nt, qq_so
  USE cell_base, ONLY : at
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp_param, ONLY: nh
  USE phus,       ONLY: alphap
  USE lrus,   ONLY : becp1, dpqq, dpqq_so
  USE control_lr, ONLY : nbnd_occ
  USE becmod,  ONLY : bec_type

  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum


  implicit none

  integer, intent(in) :: ik, jpol
  real(DP), intent(in) :: weight
  complex(DP), intent(inout) :: zstar
  complex(DP), intent(in) :: uact (3 * nat)
  TYPE(bec_type), intent(in) :: alphadk(3), becp2, bedp, alphapp(3)

  INTEGER :: ipol, ijkb0, nt, na, ih, jh, ikb, jkb, ibnd, mu, is, &
             js, ijs, startb, lastb, nbnd_eff

  complex(DP), allocatable :: ps1(:,:,:), ps2(:,:),  ps3(:,:,:), ps4(:,:) 
  complex(DP), allocatable :: ps1_nc(:,:,:,:), ps2_nc(:,:,:),  &
                              ps3_nc(:,:,:,:), ps4_nc(:,:,:) 

  complex(DP), external :: zdotc

  call start_clock('add_dkmds')
  nbnd_eff=nbnd_occ(ik)
  if (nkb.gt.0) then 
     if (noncolin) then
        allocate (ps1_nc(nkb,npol,nbnd_eff,3))
        allocate (ps2_nc(nkb,npol,nbnd_eff))
        allocate (ps3_nc(nkb,npol,nbnd_eff,3))
        allocate (ps4_nc(nkb,npol,nbnd_eff))
     else
        allocate (ps1(nkb,nbnd_eff,3))
        allocate (ps2(nkb,nbnd_eff))
        allocate (ps3(nkb,nbnd_eff,3))
        allocate (ps4(nkb,nbnd_eff))
     end if
  end if
!
!  Parallelize on the bands
!
  call divide (intra_bgrp_comm, nbnd_eff, startb, lastb)

  if (noncolin) then
     ps1_nc = (0.d0, 0.d0)
     ps2_nc = (0.d0, 0.d0)
     ps3_nc = (0.d0, 0.d0)
     ps4_nc = (0.d0, 0.d0)
  else
     ps1 = (0.d0, 0.d0)
     ps2 = (0.d0, 0.d0)
     ps3 = (0.d0, 0.d0)
     ps4 = (0.d0, 0.d0)
  endif

  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp(na).eq.nt) then
           mu = 3 * (na - 1)
           if ( abs (uact (mu + 1) ) + &
                abs (uact (mu + 2) ) + &
                abs (uact (mu + 3) ) > eps12) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ipol = 1, 3 
                       do ibnd=startb, lastb
                          if (noncolin) then 
                             if (lspinorb) then
                                ijs=0
                                do is=1,npol
                                   do js=1,npol
                                      ijs=ijs+1
                                      ps1_nc(ikb,is,ibnd,ipol)=  &
                                                ps1_nc(ikb,is,ibnd,ipol) &
                                         + qq_so(ih,jh,ijs,nt)*   &
                                              becp1(ik)%nc(jkb,js,ibnd)*  &
                                              uact (mu + ipol)
                                      ps2_nc(ikb,is,ibnd)=  &
                                                ps2_nc(ikb,is,ibnd) &
                                         + qq_so(ih,jh,ijs,nt)* &
                                           alphap(ipol,ik)%nc(jkb,js,ibnd)*  &
                                              uact (mu + ipol)
                                      ps3_nc(ikb,is,ibnd,ipol)=  &
                                                ps3_nc(ikb,is,ibnd,ipol) &
                                         + dpqq_so(ih,jh,ijs,jpol,nt)* &
                                              becp1(ik)%nc(jkb,js,ibnd)*  &
                                              uact (mu + ipol)
                                      ps4_nc(ikb,is,ibnd)=  &
                                                ps4_nc(ikb,is,ibnd) &
                                         + dpqq_so(ih,jh,ijs,jpol,nt)* &
                                           alphap(ipol,ik)%nc(jkb,js,ibnd)*  &
                                              uact (mu + ipol)
                                   enddo
                                enddo
                             else
                                do is=1,npol
                                   ps1_nc (ikb,is,ibnd,ipol)= &
                                        ps1_nc(ikb,is,ibnd,ipol)+  &
                                      qq_nt (ih, jh, nt) *                &
                                      becp1(ik)%nc(jkb,is,ibnd)*  &
                                              uact (mu + ipol)
                                   ps2_nc(ikb,is,ibnd)=  &
                                             ps2_nc(ikb,is,ibnd) &
                                     + qq_nt(ih,jh,nt)* &
                                        alphap(ipol,ik)%nc(jkb,is,ibnd)*  &
                                          uact (mu + ipol)
                                   ps3_nc(ikb,is,ibnd,ipol)=  &
                                            ps3_nc(ikb,is,ibnd,ipol) &
                                         + dpqq(ih,jh,jpol,nt)* &
                                              becp1(ik)%nc(jkb,is,ibnd)*  &
                                              uact (mu + ipol)
                                   ps4_nc(ikb,is,ibnd)=  &
                                            ps4_nc(ikb,is,ibnd) &
                                      + dpqq(ih,jh,jpol,nt)* &
                                           alphap(ipol,ik)%nc(jkb,is,ibnd)*  &
                                           uact (mu + ipol)
                                enddo
                             endif
                          else
                             ps1 (ikb,ibnd,ipol)= &
                                ps1(ikb,ibnd,ipol)+  &
                                  qq_nt (ih, jh, nt) *         &
                                  becp1(ik)%k(jkb,ibnd) * uact (mu + ipol)
                             ps2 (ikb,ibnd)=  &
                                 ps2(ikb,ibnd) + &
                                  qq_nt(ih,jh,nt)* &
                                  alphap(ipol,ik)%k(jkb,ibnd)*uact (mu + ipol)
                             ps3(ikb,ibnd,ipol)=  &
                                 ps3(ikb,ibnd,ipol) &
                                  + dpqq(ih,jh,jpol,nt)* &
                                    becp1(ik)%k(jkb,ibnd)* uact (mu + ipol)
                             ps4(ikb,ibnd)=  &
                                    ps4(ikb,ibnd) &
                                 + dpqq(ih,jh,jpol,nt)* &
                                   alphap(ipol,ik)%k(jkb,ibnd)* uact (mu + ipol)
                          endif
                       enddo
                    enddo
                 enddo
              enddo
           endif
           ijkb0=ijkb0+nh(nt)
        endif
     enddo
  enddo
   IF (noncolin) THEN
      CALL mp_sum ( ps1_nc, intra_bgrp_comm )
      CALL mp_sum ( ps2_nc, intra_bgrp_comm )
      CALL mp_sum ( ps3_nc, intra_bgrp_comm )
      CALL mp_sum ( ps4_nc, intra_bgrp_comm )
   ELSE
      CALL mp_sum ( ps1, intra_bgrp_comm )
      CALL mp_sum ( ps2, intra_bgrp_comm )
      CALL mp_sum ( ps3, intra_bgrp_comm )
      CALL mp_sum ( ps4, intra_bgrp_comm )
   END IF 
  !
  !  Finally multiply by the becp and alphap terms
  !
  IF ( nkb > 0 ) THEN
     IF (noncolin) THEN
        zstar=zstar-weight*zdotc(nkb*nbnd_eff*npol, bedp%nc, 1, ps2_nc, 1)
        zstar=zstar+weight*zdotc(nkb*nbnd_eff*npol, becp2%nc, 1, ps2_nc, 1)
        zstar=zstar-weight*zdotc(nkb*nbnd_eff*npol, becp1(ik)%nc, 1, ps4_nc, 1)
        DO ipol=1,3
           zstar=zstar-weight*zdotc(nkb*nbnd_eff*npol, alphapp(ipol)%nc, 1, ps1_nc(1,1,1,ipol), 1)
           zstar=zstar+weight*zdotc(nkb*nbnd_eff*npol, alphadk(ipol)%nc, 1, ps1_nc(1,1,1,ipol), 1)
           zstar=zstar-weight*zdotc(nkb*nbnd_eff*npol, alphap(ipol,ik)%nc, 1, ps3_nc(1,1,1,ipol), 1)

        ENDDO
     ELSE
        zstar=zstar-weight*zdotc(nkb*nbnd_eff,bedp%k,1,ps2,1)
        zstar=zstar+weight*zdotc(nkb*nbnd_eff,becp2%k,1,ps2,1)
        zstar=zstar-weight*zdotc(nkb*nbnd_eff,becp1(ik)%k,1,ps4,1)
        DO ipol=1,3
           zstar=zstar-weight*zdotc(nkb*nbnd_eff,alphapp(ipol)%k,1,ps1(1,1,ipol),1)
           zstar=zstar+weight*zdotc(nkb*nbnd_eff,alphadk(ipol)%k,1,ps1(1,1,ipol),1)
           zstar=zstar-weight*zdotc(nkb*nbnd_eff,alphap(ipol,ik)%k,1,ps3(1,1,ipol),1)
        END DO
     END IF
  END IF

  if (noncolin) then
     if (allocated(ps1_nc)) deallocate(ps1_nc)
     if (allocated(ps2_nc)) deallocate(ps2_nc)
     if (allocated(ps3_nc)) deallocate(ps3_nc)
     if (allocated(ps4_nc)) deallocate(ps4_nc)
  else
     if (allocated(ps1)) deallocate(ps1)
     if (allocated(ps2)) deallocate(ps2)
     if (allocated(ps3)) deallocate(ps3)
     if (allocated(ps4)) deallocate(ps4)
  end if

  call stop_clock('add_dkmds')
  return

end subroutine add_dkmds_tpw
