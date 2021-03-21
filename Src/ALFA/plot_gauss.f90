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

!>  Represents the ionic charge density as gaussians
!>  centered on the atoms

subroutine plot_gauss(sigma, rhogau,                                     &
         adot, ntype, natom, rat, izval,                                 &
         ng, kgv,                                                        &
         mxdtyp, mxdatm, mxdgve)

! Written June 3, 2014. JLM
! Documentation, 4 February 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

  implicit none

! version 4.99

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors

  real(REAL64), intent(in)           ::  sigma                           !<  width of gaussian (in real space)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  integer, intent(in)                ::  izval(mxdtyp)                   !<  valence of atom of type i

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

! output

  complex(REAL64), intent(out)       ::  rhogau(mxdgve)                  !<  gaussian charge density for G-vector j

! local variables

  real(REAL64)     ::  vcell, bdot(3,3)
  real(REAL64)     ::  gmod, amp, fi

! constants

  real(REAL64), parameter     :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  :: C_I = cmplx(ZERO,UM,REAL64)
  real(REAL64), parameter     :: PI=3.141592653589793_REAL64

! counters

  integer      :: n, nt, j,  i,k

  call adot_to_bdot(adot,vcell,bdot)
  
  do n = 1,ng
    rhogau(n) = C_ZERO
    gmod = ZERO
    do i = 1,3
    do k = 1,3
      gmod = gmod + kgv(k,n)*(bdot(k,i)*kgv(i,n))
    enddo
    enddo
    amp = exp(-gmod*sigma*sigma/4)
    do nt = 1,ntype
      do j = 1,natom(nt)
        fi = kgv(1,n)*rat(1,j,nt) + kgv(2,n)*rat(2,j,nt) +          &
             kgv(3,n)*rat(3,j,nt)
        fi = 2*PI*fi
        rhogau(n) = rhogau(n) + izval(nt)*amp*exp(-C_I*fi)
      enddo
    enddo
  enddo
  
  return

  end subroutine plot_gauss
