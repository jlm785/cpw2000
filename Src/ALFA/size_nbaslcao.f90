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
!>  \date         21 April 2014.
!>  \copyright    GNU Public License v2

subroutine size_nbaslcao(ntype, natom, norbat, lorb, nbaslcao,           &
    mxdtyp, mxdlao)

! written April 21, 2014. JLM
! Modified, documentation, January 2020. JLM
! Indentation. 28 September 2026. JLM+claude


  implicit none

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i

! output

  integer, intent(out)               ::  nbaslcao                        !<  number of vectors in atomic orbital basis

! local variables

  integer   ::  ind

! counters

  integer         ::  nt, j

  ind = 0
  do nt=1,ntype
    do j=1,norbat(nt)
      ind = ind + (2*lorb(j,nt)+1)*natom(nt)
    enddo
  enddo
  nbaslcao = ind

  return

end subroutine size_nbaslcao
