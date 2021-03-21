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

!>     calculates the density of states with a gaussian method

       subroutine dos_gau(el,nrk,nband,w,ezero,lidos,np,mxdbnd)

!      written november 9 1987.jlm
!      modified November 15, 2013. jlm
!      Modified, documentation, 19 September 2020. JLM
!      copyright  J.L.Martins, INESC-MN.

!      version 4.53 of cpw

       implicit none

       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input:

       integer, intent(in)             ::  nrk                           !<  number of k-points
       integer, intent(in)             ::  mxdbnd                        !<  size of number of bands
       real(REAL64), intent(in)        ::  el(mxdbnd,nrk)                !<  eigenvalues in Hartree
       integer, intent(in)             ::  nband(nrk)                    !<  number of bands for each k-points
       real(REAL64), intent(in)        ::  w(nrk)                        !<  weight for each k-point
       integer, intent(in)             ::  np                            !<  approximate number of points for plot
       real(REAL64), intent(in)        ::  ezero                         !<  zero of energy
       logical, intent(in)             ::  lidos                         !<  true if integrated density of states is to be computed.

!      local arrays

       real(REAL64), allocatable  ::  ehist(:),dhist(:),chist(:)         !<  histogram arrays (energy, dos, idos)

!      local variables

       integer  ::  nhist                                                !  number of points in histogram
       integer  ::  nbm
       real(REAL64)  ::  emin, emax, xmax, deltae, b, xx
       real(REAL64)  ::  alfa, gam                                       !  width of gaussian

!      counters

       integer   ::  i,j,k

!      constants

       real(REAL64), parameter :: HARTREE = 27.21138386_REAL64
       real(REAL64), parameter :: SPI = 1.772453850905516_REAL64
       real(REAL64), parameter :: ZERO = 0.0_REAL64 , UM = 1.0_REAL64

!      find max of bands

       nbm = nband(1)
       do i=1,nrk
         if(nbm > nband(i)) nbm = nband(i)
       enddo

!      find emin and emax

       emin = el(1,1) - ezero
       emax = el(1,1) - ezero
       do j=1,nbm
         if(emax < el(j,1) - ezero) emax = el(j,1) - ezero
       enddo
       do i=1,nrk
         xmax = el(1,i) - ezero
         do j=1,nbm
           if(xmax < el(j,i) - ezero) xmax = el(j,i) - ezero
           if(emin > el(j,i) - ezero) emin = el(j,i) - ezero
         enddo
         if(emax > xmax) emax = xmax
       enddo

!      find scale for the plot

       deltae = HARTREE*(emax-emin) / np

       call plot_step(deltae,b)

       deltae = b/HARTREE

       emin = deltae*real(int(emin/deltae)-3)
       emax = emax + deltae
       nhist = nint((emax-emin)/deltae) + 1

       allocate(ehist(nhist))
       allocate(chist(nhist))
       allocate(dhist(nhist))

!      initialize arrays ehist, dhist and chist

       do i=1,nhist
         ehist(i) = emin + deltae*(i-1)
         dhist(i) = ZERO
         chist(i) = ZERO
       enddo

!      defines the gaussian width

       alfa = 5*(emax-emin)/real(nrk*nbm)
       alfa = max(alfa,2*deltae)
       gam = UM/(SPI*alfa)

!      loop over k points and bands

       do i=1,nrk
         do j=1,nbm
           do k=1,nhist
             xx = (ehist(k) - (el(j,i) - ezero) ) / alfa
             dhist(k) = dhist(k) + w(i)*gam*exp(-xx*xx)
             chist(k) = chist(k) + w(i)*(UM + erf(xx))/2
           enddo
         enddo
       enddo

!      printout

       write(6,'("#")')
       write(6,'("#    Gaussian width: ",f7.3)') alfa*HARTREE
       write(6,'("#")')

       call dos_print_ascii(emin,deltae,nhist,dhist,chist,lidos,1,60)

       deallocate(ehist)
       deallocate(chist)
       deallocate(dhist)

       return
       end subroutine dos_gau
