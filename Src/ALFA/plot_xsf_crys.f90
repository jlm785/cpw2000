!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

!>  writes the crystal structure part of a xsf
!>  input file for vesta or xcrysden

subroutine plot_xsf_crys(iotape, lvesta,                                 &
   adot, ntype, natom, nameat, rat,                                      &
   mxdtyp, mxdatm)

! written July 2013. JLM
! Heavily modified to use cpw crystal description, 31 January 2021. JLM

! copyright  J.L.Martins, INESC-MN.

! version 4.99

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms

  integer, intent(in)                ::  iotape                          !<  tape number
  logical, intent(in)                ::  lvesta                          !<  indicates that the xsf file follows the vesta orientation

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in )      ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

! local:

  integer          ::  nat                                   ! total number of atoms

  real(REAL64)     ::  avec(3,3), bvec(3,3)
  real(REAL64)     ::  a1(3), a2(3), a3(3)
  real(REAL64)     ::  rcar(3)
  integer          ::  ichrg

! constants

  real(REAL64), parameter  :: BOHR = 0.5291772109_REAL64
 
! counters

  integer      :: nt, j, k


! converts stuff

  call adot_to_avec_sym(adot,avec,bvec)
  
  do k = 1,3
    a1(k) = avec(k,1) * BOHR
    a2(k) = avec(k,2) * BOHR
    a3(k) = avec(k,3) * BOHR
  enddo

! writes preamble

  write(iotape,'(a8)') ' CRYSTAL'
  write(iotape,'(a8)') ' PRIMVEC'

  if(lvesta) then

! change order for xsf so that 3 axis is along y...

    write(iotape,'(3(1x,f16.10))') a1(2),a1(3),a1(1)
    write(iotape,'(3(1x,f16.10))') a2(2),a2(3),a2(1)
    write(iotape,'(3(1x,f16.10))') a3(2),a3(3),a3(1)

  else

    write(iotape,'(3(1x,f16.10))') a1(1),a1(2),a1(3)
    write(iotape,'(3(1x,f16.10))') a2(1),a2(2),a2(3)
    write(iotape,'(3(1x,f16.10))') a3(1),a3(2),a3(3)

  endif

  nat = 0
  do nt = 1,ntype
    nat = nat + natom(nt)
  enddo

  write(iotape,'(a10)') ' PRIMCOORD'
  write(iotape,'(1x,i4,1x,i1)') nat,1

! loop over atoms

  do nt = 1,ntype
    call p_tbl_charge(nameat(nt),ichrg)
    do j = 1, natom(nt)

      do k = 1,3
        rcar(k) = avec(k,1)*rat(1,j,nt)  + avec(k,2)*rat(2,j,nt) + avec(k,3)*rat(3,j,nt)
        rcar(k) = rcar(k) * BOHR
      enddo

      if(lvesta) then
        write(iotape,'(i4,1x,3(2x,f16.10))') ichrg,rcar(2),rcar(3),rcar(1)
      else
        write(iotape,'(i4,1x,3(2x,f16.10))') ichrg, (rcar(k), k = 1,3)
      endif

    enddo
  enddo

  return
end subroutine plot_xsf_crys

