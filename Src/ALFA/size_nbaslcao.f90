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

!>  Calculates the size of the atomic orbitals basis set.
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         21 April 2014, 16 October 2026.
!>  \copyright    GNU Public License v2

subroutine size_nbaslcao(ntype, natom, n_bsets, jbset, norbat, lorb,     &
    nbaslcao,                                                            &
    mxdtyp, mxdlao, mxdset)

! written April 21, 2014. JLM
! Modified, documentation, January 2020. JLM
! Indentation. 28 September 2026. JLM+claude
! Several atomic basis sets, n_bsets, jbset, mxdset. 16 October 2026. JLM+claude


  implicit none

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type
  integer, intent(in)                ::  mxdset                          !<  array dimension for number of atomic basis sets

  integer, intent(in)                ::  n_bsets(mxdtyp)                 !<  number of basis sets for atom k
  integer, intent(in)                ::  jbset(mxdtyp)                   !<  basis set chosen for atom k (if outside 1...n_bsets the nearest is used)
  integer, intent(in)                ::  norbat(mxdset,mxdtyp)           !<  number of atomic orbitals for basis set nb and atom k
  integer, intent(in)                ::  lorb(mxdlao,mxdset,mxdtyp)      !<  angular momentum of orbital n of basis nb of atom k

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i

! output

  integer, intent(out)               ::  nbaslcao                        !<  number of vectors in atomic orbital basis

! local variables

  integer   ::  ind
  integer   ::  ib                                                       !  basis set used, consistent with atomic_orbital_c16

! counters

  integer         ::  nt, j

  ind = 0
  do nt=1,ntype
    ib = max(1,min(jbset(nt),n_bsets(nt)))
    do j=1,norbat(ib,nt)
      ind = ind + (2*lorb(j,ib,nt)+1)*natom(nt)
    enddo
  enddo
  nbaslcao = ind

  return

end subroutine size_nbaslcao
