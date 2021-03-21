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

!>     This subroutine calculates the band structure along a path
!>     It uses a full diagonalization
!>     Files band.agr, band_so.agr, band.gp and band_so.gp,
!>     for later ploting with gnuplot and xmgrace are written
!>     Circuit for band structure is defined in BAND_LINES.DAT

       subroutine out_band(title, subtitle,                              &
     & emax, flgdal, flgpsd, iguess, epspsi, icmax, ztot,                &
     & adot, ntype, natom, rat,                                          &
     & ng, kgv, phase, conj,                                             &
     & ns, inds, kmax, indv, ek,                                         &
     & sfact, icmplx,                                                    &
     & veff,                                                             &
     & nqnl, delqnl, vkb, nkb,                                           &
     & latorb, norbat, nqwf, delqwf, wvfao, lorb,                        &     
     & mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)

!      version 4.42. 8 may 2004. jlm
!      modified 11 february 2008 (read from file). JLM
!      modified September 5, 2012 (f90, circuit). jlm
!      modified October 18, 2013 (band index permutations). jlm
!      modified January 8, 2013 (new interface). jlm
!      modified, dimensions vkb, March 31, 2014. jlm
!      modified (title,subtitle) xmgrace plot. 4 August 2012. JLM
!      modified 4.7X November 2015. JLM
!      modified 7 November 2018. Bug, use ljump in match state. JLM
!      modified 1 August 2019. Reduce memory use, documentation. JLM
!      modified February 2020, out_band_match, lkpg. JLM
!      modified, ifail, 18 February 2020. JLM
!      Modified h_kb_dia_all, 7 June 2020. JLM
!      Modified, vmax, vmin, 27 November 2020. JLM

!      copyright  Jose Luis Martins/INESC-MN

!      version 4.99

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type
       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars
       integer, intent(in)                ::  mxdlqp                     !<  array dimension for local potential
       integer, intent(in)                ::  mxdcub                     !<  array dimension for 3-index g-space
       integer, intent(in)                ::  mxdlao                     !<  array dimension of orbital per atom type

       character(len=50), intent(in)      ::  title                      !<  title for plots
       character(len=140), intent(in)     ::  subtitle                   !<  subtitle for plots

       real(REAL64), intent(in)           ::  emax                       !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
       character(len=4), intent(in)       ::  flgdal                     !<  dual approximation if equal to 'DUAL'
       character(len=6), intent(in)       ::  flgpsd                     !<  type of pseudopotential
       real(REAL64), intent(in)           ::  epspsi                     !<  requested precision of the eigenvectors
       integer, intent(in)                ::  icmax                      !<  maximum number of iterations for diagonalization
       real(REAL64), intent(in)           ::  ztot                       !<  total charge density (electrons/cell)

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase

       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
       integer, intent(in)                ::  inds(mxdgve)               !<  star to which g-vector n belongs
       integer, intent(in)                ::  kmax(3)                    !<  max value of |kgv(i,n)|
       integer, intent(in)                ::  indv(mxdcub)               !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax
       real(REAL64), intent(in)           ::  ek(mxdnst)                 !<  kinetic energy (hartree) of g-vectors in star j

       complex(REAL64), intent(in)        ::  sfact(mxdtyp,mxdnst)       !<  structure factor
       integer, intent(in)                ::  icmplx                     !<  indicates if the structure factor is complex

       complex(REAL64), intent(in)        ::  veff(mxdnst)               !<  ionic potential (local+Hartree+XC) for the prototype g-vector in star j

       integer, intent(in)                ::  nqnl(mxdtyp)               !<  number of points for the non-local pseudopotential interpolation
       real(REAL64), intent(in)           ::  delqnl(mxdtyp)             !<  step used in the interpolation
       real(REAL64), intent(in)  ::   vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)     !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. (not normalized to vcell, hartree)
       integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)       !<   KB pseudo.  normalization for atom k, ang. mom. l

       logical, intent(in)                ::  latorb                     !<  indicates if all atoms have information about atomic orbitals
       integer, intent(in)                ::  norbat(mxdtyp)             !<  number of atomic orbitals for atom k
       integer, intent(in)                ::  nqwf(mxdtyp)               !<  number of points for wavefunction interpolation for atom k
       real(REAL64), intent(in)           ::  delqwf(mxdtyp)             !<  step used in the wavefunction interpolation for atom k
       real(REAL64), intent(in)   ::   wvfao(-2:mxdlqp,mxdlao,mxdtyp)    !<  wavefunction for atom k, ang. mom. l 
       integer, intent(in)                ::  lorb(mxdlao,mxdtyp)        !<  angular momentum of orbital n of atom k

!      input and output

       integer, intent(inout)             ::  iguess                     !<  LOCAL USE, keep for compatibility  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0
       
       
!      allocatable arrays for Brillouin zone path

       integer                            ::  nlines                     !  number of lines in reciprocal space
       integer, allocatable               ::  nkstep(:)                  !  number of steps in line
       logical, allocatable               ::  ljump(:)                   !  indicates if the new line contains a jump from the preceeding
       integer                            ::  nvert                      !  number of vertical lines in plot
       real(REAL64), allocatable          ::  xcvert(:)                  !  x coordinate of vertical line
       real(REAL64), allocatable          ::  xk(:)                      !  x coordinate of k-point in plot
       real(REAL64), allocatable          ::  rk(:,:)                    !  x coordinate of k-point in plot
       real(REAL64), allocatable          ::  e_of_k(:,:)                !  band energies of k-point in plot
       real(REAL64), allocatable          ::  e_of_k_so(:,:)             !  spin-orbit band energies of k-point in plot
       character(len=6), allocatable      ::  label(:)                   !  label of symmetry k-points
       real(REAL64), allocatable          ::  xklab(:)                   !  x coordinate of label

!      allocatable arrays with larger scope

       real(REAL64), allocatable          ::  ei(:)                      !  eigenvalue no. i. (hartree)
       real(REAL64), allocatable          ::  ei_so(:)                   !  spin-orbit eigenvalue (hartree)
       real(REAL64), allocatable          ::  hdiag(:)                   !  hamiltonian diagonal
       integer, allocatable               ::  isort(:)                   !  g-vector associated with row/column i of hamiltonian
       integer, allocatable               ::  isort_so(:)                !  g-vector associated with row/column i of spin-orbit hamiltonian
       real(REAL64), allocatable          ::  qmod(:)                    !  length of k+g-vector of row/column i
       real(REAL64), allocatable          ::  ekpg(:)                    !  kinetic energy (hartree) of k+g-vector of row/column i
       complex(REAL64), allocatable       ::  psi(:,:)                   !  |psi> component j of eigenvector i (guess on input)
       complex(REAL64), allocatable       ::  hpsi(:,:)                  !  H | psi> 
       real(REAL64), allocatable          ::  ekpsi(:)                   !  kinetic energy of eigenvector i. (hartree)
       real(REAL64), allocatable          ::  vscr(:)                    !  screened potential in the fft real space mesh
       complex(REAL64), allocatable       ::  psi_so(:,:)                !  component j of eigenvector i (guess on input)
       real(REAL64), allocatable          ::  ekpsi_so(:)                !  kinetic energy of eigenvector i. (hartree)

!      local variables

       integer           ::  mxdscr         !  array dimension for screening potential

       integer           ::  ifail          !  if ifail=0 the ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

       integer           ::  mxddim         !  array dimension for the hamiltonian
       integer           ::  mxdbnd         !  array dimension for the number of bands
       integer           ::  mxdwrk         !  array dimension for fft transform workspace

       integer           ::  mtxd           !  dimension of the hamiltonian
       integer           ::  neig           !  number of eigenvectors required (maybe modified on output)
       real(REAL64)      ::  rkpt(3)        !  j-th component in lattice coordinates of the k-point
       integer           ::  kmscr(7)       !  max value of kgv(i,n) used for the potential fft mesh
       integer           ::  idshift        !  shift of the fft mesh, used /= 0 only in highly banked memory.

       real(REAL64)      ::  vmax, vmin     !  maximum and minimum values of vscr

       real(REAL64)      ::  eref           !  reference energy for plot
       integer           ::  nocc           !  number of occupied states (different color) or recycled
       integer           ::  nstyle         !  choice of plot style

       integer           ::  irk,nrka
       integer           ::  iotape
       character(len=5)  ::  labelk
       integer           ::  ipr,nrk2,nd
       integer           ::  nsfft(3)

       integer           ::  nextline       !  indicates next line on band circuit
       integer           ::  nl             !  keeps track of next lines

!      constants

       real(REAL64), parameter  :: XSC = 1.000_REAL64                     !  criteria for reduced vector size 
!       real(REAL64), parameter  :: XSC = 0.99_REAL64                     !  Use this value or even smaller if you have problems with memory 

!      counters

       integer    ::  j


!      calculates local potential in fft mesh


       if(flgdal == 'DUAL') then
         kmscr(1) = kmax(1)/2 + 2
         kmscr(2) = kmax(2)/2 + 2
         kmscr(3) = kmax(3)/2 + 2
       else
         kmscr(1) = kmax(1)
         kmscr(2) = kmax(2)
         kmscr(3) = kmax(3)
       endif

       call size_fft(kmscr,nsfft,mxdscr,mxdwrk)

       allocate(vscr(mxdscr))

       ipr = 1

       idshift = 0
       call pot_local(ipr, vscr, vmax, vmin, veff, kmscr, idshift,       &
     & ng, kgv, phase, conj, ns, inds,                                   &
     & mxdscr, mxdgve, mxdnst)


       iotape = 13
       call out_band_circuit_size('BAND_LINES.DAT',iotape,1,adot,        &
     &                  neig,nrk2,nlines,nvert)

       allocate(xk(nrk2))
       allocate(rk(3,nrk2))
       allocate(xcvert(nvert))
       allocate(ljump(nlines))
       allocate(nkstep(nlines))
       allocate(label(nvert+nlines))
       allocate(xklab(nvert+nlines))

       call out_band_get_circuit('BAND_LINES.DAT',iotape,1,adot,         &
     &                  xk,rk,xcvert,ljump,nkstep,label,xklab,           &
     &                  neig,nrk2,nlines,nvert)


       allocate(e_of_k(neig,nrk2))
       allocate(e_of_k_so(2*neig,nrk2))

!      finds mxddim, mxdbnd

       mxdbnd = neig

       mxddim = 1
       do irk=1,nrk2

!        loop over k-points

         do j=1,3
           rkpt(j) = rk(j,irk)
         enddo         

         call size_mtxd(emax,rkpt,adot,ng,kgv,nd)

         if(nd > mxddim) mxddim = nd

       enddo

!      allocates arrays

       allocate(ei(mxdbnd))
       allocate(ei_so(2*mxdbnd))
       allocate(hdiag(mxddim))
       allocate(isort(mxddim))
       allocate(qmod(mxddim))
       allocate(ekpg(mxddim))
       allocate(psi(mxddim,mxdbnd))
       allocate(hpsi(mxddim,mxdbnd))
       allocate(ekpsi(mxdbnd))

       allocate(isort_so(2*mxddim))
       allocate(psi_so(2*mxddim,2*mxdbnd))
       allocate(ekpsi_so(2*mxdbnd))



       iguess = 0
       nextline = 1
       nl = 1

       do irk=1,nrk2

!        loop over k-points

         do j=1,3
           rkpt(j) = rk(j,irk)
         enddo
 
         write(6,*)
         write(6,'("  Working on k-point ",i7,"  out of ", i7)')         &
     &               irk,nrk2
         write(6,*)

         nocc = neig

         call h_kb_dia_all('pw  ', emax, rkpt, neig, nocc,               &
     &   flgpsd, ipr, ifail, icmax, iguess, epspsi,                      &
     &   ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                 &
     &   sfact, veff, icmplx,                                            &
     &   nqnl, delqnl, vkb, nkb,                                         &
     &   ntype, natom, rat, adot,                                        &
     &   mtxd, hdiag, isort, qmod, ekpg, .FALSE.,                        &
     &   psi, hpsi, ei,                                                  &
     &   vscr, kmscr,                                                    &
     &   latorb, norbat, nqwf, delqwf, wvfao, lorb,                      &     
     &   mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,         &
     &   mxdbnd, mxdscr, mxdlao)


         call kinetic_energy(neig,mtxd,ekpg,psi,ekpsi,                   &
     &   mxddim,mxdbnd)

         ipr = 1
         nrka = -1
         call print_eig(ipr,irk,labelk,nrka,rkpt,                        &
     &   mtxd,icmplx,neig,psi,                                           &
     &   adot,ei,ekpsi,isort,kgv,                                        &
     &   mxddim,mxdbnd,mxdgve)


         call spin_orbit_perturb(rkpt,mtxd,isort,                        &
     &   neig,psi,ei,ei_so,psi_so,.TRUE.,                                &
     &   ng,kgv,                                                         &
     &   nqnl,delqnl,vkb,nkb,                                            &
     &   ntype,natom,rat,adot,                                           &
     &   mxdtyp,mxdatm,mxdlqp,mxddim,mxdbnd,mxdgve)

         ipr = 1
         nrka = -1

         if(ipr == 2) then
           call kinetic_energy_so(neig,mtxd,ekpg,psi_so,ekpsi_so,        &
     &     mxddim,mxdbnd)
         endif

         call print_eig_so(ipr,irk,labelk,nrka,rkpt,                     &
     &   mtxd,neig,psi_so,                                               &
     &   adot,ei_so,ekpsi_so,isort,kgv,                                  &
     &   mxddim,mxdbnd,mxdgve)



         call out_band_match(irk,nlines,ljump,nkstep,                    &
     &        neig,icmplx,XSC,                                           &
     &        mtxd,psi,isort,psi_so,ei,ei_so,                            &
     &        e_of_k,e_of_k_so,nrk2,                                     &
     &        mxddim,mxdbnd)

       enddo

!      writes the output files for xmgrace and gnuplot

       iotape = 15
       nstyle = 2

       call out_band_eref(neig,nrk2,ztot,2,1,e_of_k,eref,nocc)

       call out_band_gnuplot('band.gp',iotape,                           &
     &        neig,nrk2,xk,e_of_k,eref,                                  &
     &        nvert,xcvert,nlines,ljump,nkstep,label,xklab)

       call out_band_xmgrace('band.agr',iotape,                          &
     &        title,subtitle,nstyle,                                     &
     &        neig,nrk2,xk,e_of_k,eref,nocc,                             &
     &        nvert,xcvert,nlines,ljump,nkstep,label,xklab)

       call out_band_eref(neig,nrk2,ztot,1,1,e_of_k_so,eref,nocc)

       call out_band_gnuplot('band_so.gp',iotape,                        &
     &    2*neig,nrk2,xk,e_of_k_so,eref,                                 &
     &    nvert,xcvert,nlines,ljump,nkstep,label,xklab)

       call out_band_xmgrace('band_so.agr',iotape,                       &
     &    title,subtitle,nstyle,                                         &
     &    2*neig,nrk2,xk,e_of_k_so,eref,nocc,                            &
     &    nvert,xcvert,nlines,ljump,nkstep,label,xklab)


       deallocate(nkstep)
       deallocate(ljump)

       deallocate(xcvert)
       deallocate(xk)
       deallocate(rk)
       deallocate(e_of_k)
       deallocate(label)
       deallocate(xklab)

       deallocate(vscr)

       deallocate(ei)
       deallocate(hdiag)
       deallocate(isort)
       deallocate(qmod)
       deallocate(ekpg)
       deallocate(psi)
       deallocate(ekpsi)

       deallocate(psi_so)
       deallocate(ekpsi_so)

       return
       end subroutine out_band
