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
!>     using a SVD based Generalized Luttinger-Kohn interpolation
!>     Files band.agr, band_so.agr, band.gp and band_so.gp,
!>     for later ploting with gnuplot and xmgrace are written
!>     Circuit for band structure is defined in BAND_LINES.DAT
!>
!>  \author       Jose Luis Martins
!>  \version      5.03
!>  \date         8 may 2004, 29 November 2021.
!>  \copyright    GNU Public License v2

       subroutine out_band_glk(title, subtitle,                         &
     & ninterp, nignore, xsvd, csvd,                                    &
     & emax, flgdal, flgpsd, iguess, epspsi, icmax, ztot, efermi,       &
     & adot, ntype, natom, rat,                                         &
     & ng, kgv, phase, conj,                                            &
     & ns, inds, kmax, indv, ek,                                        &
     & sfact, icmplx,                                                   &
     & veff,                                                            &
     & nqnl, delqnl, vkb, nkb,                                          &
     & latorb, norbat, nqwf, delqwf, wvfao, lorb,                       &
     & mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)


!      version 4.42. 8 may 2004. jlm
!      modified 11 february 2008 (read from file). JLM
!      modified September 5, 2012 (f90, circuit). jlm
!      modified October 18, 2013 (band index permutations). jlm
!      modified January 8, 2013 (new interface). jlm
!      modified, vkb dimensions, March 31, 2014. jlm
!      modified, July 9, 2014.  JLM
!      modified (title,subtitle) xmgrace plot. 4 August 2012. JLM
!      Modified 8 November 2015. Compatibility new libpw. JLM
!      Modified, 9 June 2020, from 4.93CJA h_kb_dia_all,
!      spin-orbit, match-state. JLM
!      Modified, 14 June 2020, icmax. JLM
!      Modified, csvd, xsvd, interpolation_glk, 19-24 August 2020. JLM
!      Modified, vmax, vmin, 27 November 2020. JLM
!      Modified, efermi, 29 November 2021. JLM

!      copyright  Jose Luis Martins/INESC-MN

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type
       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars
       integer, intent(in)                ::  mxdlqp                     !<  array dimension for local potential
       integer, intent(in)                ::  mxdcub                     !<  array dimension for 3-index g-space
       integer, intent(in)                ::  mxdlao                     !<  array dimension for 3-index g-space

       character(len=50), intent(in)      ::  title                      !<  title for plots
       character(len=140), intent(in)     ::  subtitle                   !<  subtitle for plots

       integer, intent(in)                ::  ninterp                    !<  number of interpolation points
       integer, intent(in)                ::  nignore                    !<  ignore the first nignore bands
       real(REAL64), intent(in)           ::  xsvd                       !<  Ignore states with SVD singular values smaller than xsvd.  0 < xsvd < 1.
       real(REAL64), intent(in)           ::  csvd                       !<  Use csvd*neig in SVD procedure. 1 < csvd < 2.

       real(REAL64), intent(in)           ::  emax                       !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
       character(len=6), intent(in)       ::  flgpsd                     !<  type of pseudopotential
       character(len=4), intent(in)       ::  flgdal                     !<  dual approximation if equal to 'DUAL'
       integer, intent(in)                ::  iguess                     !<  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0
       real(REAL64), intent(in)           ::  epspsi                     !<  requested precision of the eigenvectors
       integer, intent(in)                ::  icmax                      !<  maximum number of iterations for diagonalization
       real(REAL64), intent(in)           ::  ztot                       !<  total charge density (electrons/cell)
       real(REAL64), intent(in)           ::  efermi                     !  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is 1 if one must take the complex conjugate of x*phase

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
       real(REAL64), intent(in)  ::   vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)     !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
       integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)       !<  KB pseudo.  normalization for atom k, ang. mom. l

       logical, intent(in)                ::  latorb                     !<  indicates if all atoms have information about atomic orbitals
       integer, intent(in)                ::  norbat(mxdtyp)             !<  number of atomic orbitals for atom k
       integer, intent(in)                ::  nqwf(mxdtyp)               !<  number of points for wavefunction interpolation for atom k
       real(REAL64), intent(in)           ::  delqwf(mxdtyp)             !<  step used in the wavefunction interpolation for atom k
       real(REAL64), intent(in)   ::   wvfao(-2:mxdlqp,mxdlao,mxdtyp)    !<  wavefunction for atom k, ang. mom. l
       integer, intent(in)                ::  lorb(mxdlao,mxdtyp)        !<  angular momentum of orbital n of atom k

!      allocatable arrays for Brillouin zone path

       integer                            ::  nlines                     !  number of lines in reciprocal space
       integer, allocatable               ::  nkstep(:)                  !  number of steps in line
       logical, allocatable               ::  ljump(:)                   !  indicates if the new line contains a jump from the preceeding
       integer                            ::  nvert                      !  number of vertical lines in plot
       real(REAL64), allocatable          ::  xcvert(:)                  !  x coordinate of vertical line
       real(REAL64), allocatable          ::  xk(:)                      !  x coordinate of k-point in plot
       real(REAL64), allocatable          ::  rk(:,:)                    !  k-point in plot
       character(len=6), allocatable      ::  label(:)                   !  label of symmetry k-points
       real(REAL64), allocatable          ::  xklab(:)                   !  x coordinate of label

       real(REAL64), allocatable          ::  e_of_k(:,:)                !  band energies of k-point in plot
       real(REAL64), allocatable          ::  e_of_k_so(:,:)             !  spin-orbit band energies of k-point in plot

       real(REAL64), allocatable          ::  e_of_k_int(:,:)            !  interpolated band energies of k-point in plot
       real(REAL64), allocatable          ::  e_of_k_so_int(:,:)         !  interpolated band energies of k-point in plot with spin-orbit
       real(REAL64), allocatable          ::  xk_int(:)                  !  x coordinate of interpolated k-point in plot
       integer, allocatable               ::  nkstep_int(:)              !  number of steps in interpolated line

!      variables for match_state

       complex(REAL64),  allocatable      ::  psiold(:,:)                !  old eigenvectors
       integer, allocatable               ::  imatch(:)                  !  new matching
       integer, allocatable               ::  isold(:)                   !  old isort
       integer, allocatable               ::  iperm(:),ipermold(:)       !  permutation of the eigenvalues

       complex(REAL64),  allocatable      ::  psiold_so(:,:)             !  old eigenvectors
       integer, allocatable               ::  imatch_so(:)               !  new matching
       integer, allocatable               ::  iperm_so(:),ipermold_so(:) !  permutation of the eigenvalues

!      allocatable arrays with larger scope

       real(REAL64), allocatable          ::  ei(:)                      !  eigenvalue no. i. (hartree)
       real(REAL64), allocatable          ::  ei_so(:)                   !  spin-orbit eigenvalue (hartree)
       real(REAL64), allocatable          ::  hdiag(:)                   !  hamiltonian diagonal
       integer, allocatable               ::  isort(:)                   !  g-vector associated with row/column i of hamiltonian
       integer, allocatable               ::  isort_so(:)                !  g-vector associated with row/column i of spin-orbit hamiltonian
       real(REAL64), allocatable          ::  qmod(:)                    !  length of k+g-vector of row/column i
       real(REAL64), allocatable          ::  ekpg(:)                    !  kinetic energy (hartree) of k+g-vector of row/column i
       complex(REAL64), allocatable       ::  psi(:,:)                   !  component j of eigenvector i
       real(REAL64), allocatable          ::  ekpsi(:)                   !  kinetic energy of eigenvector i. (hartree)

       complex(REAL64), allocatable       ::  psi_all(:,:,:)             !  wavevector of points used in the interpolation
       integer,  allocatable              ::  isort_all(:,:)             !  g-vector associated with row/column i of hamiltonian

       real(REAL64), allocatable          ::  ei0(:)                     !  eigenvalue no. i. (hartree)
       real(REAL64), allocatable          ::  hdiag0(:)                  !  hamiltonian diagonal
       integer, allocatable               ::  isort0(:)                  !  g-vector associated with row/column i of hamiltonian
       real(REAL64), allocatable          ::  qmod0(:)                   !  length of k+g-vector of row/column i
       real(REAL64), allocatable          ::  ekpg0(:)                   !  kinetic energy (hartree) of k+g-vector of row/column i
       complex(REAL64), allocatable       ::  psi0(:,:)                  !  component j of eigenvector i
       complex(REAL64), allocatable       ::  hpsi0(:,:)                 !  H | psi > for component j of eigenvector i
       real(REAL64), allocatable          ::  ekpsi0(:)                  !  kinetic energy of eigenvector i. (hartree)

       real(REAL64), allocatable          ::  ei1(:)                     !  eigenvalue no. i. (hartree)
       real(REAL64), allocatable          ::  hdiag1(:)                  !  hamiltonian diagonal
       integer, allocatable               ::  isort1(:)                  !  g-vector associated with row/column i of hamiltonian
       real(REAL64), allocatable          ::  qmod1(:)                   !  length of k+g-vector of row/column i
       real(REAL64), allocatable          ::  ekpg1(:)                   !  kinetic energy (hartree) of k+g-vector of row/column i
       complex(REAL64), allocatable       ::  psi1(:,:)                  !  component j of eigenvector i
       real(REAL64), allocatable          ::  ekpsi1(:)                  !  kinetic energy of eigenvector i. (hartree)

       real(REAL64), allocatable          ::  vscr(:)                    !  screened potential in the fft real space mesh
       complex(REAL64), allocatable       ::  psi_so(:,:)                !  component j of eigenvector i (guess on input)

!      local variables

       integer           ::  mxdscr         !  array dimension for screening potential

       real(REAL64)      ::  xw_all(2)      !  relative coordinates of the k-point sum xk_all = 1
       integer           ::  mtxd_all(2)    !  wavefunction dimension (basis size) for k-point
       integer           ::  neig_all(2)  !  number of wavefunctions for k-point
       real(REAL64)      ::  rkpt_all(3,2)

       integer           ::  mxddim         !  array dimension for the hamiltonian
       integer           ::  mxdbnd         !  array dimension for the number of bands
       integer           ::  mxdwrk         !  array dimension for fft transform workspace

       integer           ::  mtxd           !  dimension of the hamiltonian
       integer           ::  mtxd0          !  dimension of the hamiltonian
       integer           ::  mtxd1          !  dimension of the hamiltonian
       integer           ::  neig           !  number of eigenvectors required (maybe modified on output)
       real(REAL64)      ::  rkpt(3)        !  j-th component in lattice coordinates of the k-point
       integer           ::  kmscr(7)       !  max value of kgv(i,n) used for the potential fft mesh
       integer           ::  idshift        !  shift of the fft mesh, used /= 0 only in highly banked memory.

       real(REAL64)      ::  vmax, vmin     !  maximum and minimum values of vscr

       real(REAL64)      ::  eref           !  reference energy for plot
       integer           ::  nocc           !  number of occupied states (different color) or recycled
       integer           ::  nstyle         !  choice of plot style

       integer           ::  irk,nrka
       integer           ::  irkint
       integer           ::  iotape
       character(len=5)  ::  labelk
       integer           ::  nrk2, nd
       integer           ::  nrk3           !  number of poins in interpolated grid
       integer           ::  nsfft(3)
       real(REAL64)      ::  rk0(3),rk1(3)

       real(REAL64)      ::  xw

       integer           ::  ipr, ifail

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter  :: XSC = 1.000_REAL64                     !  criteria for reduced vector size
!       real(REAL64), parameter  :: XSC = 0.99_REAL64                     !  Use this value or even smaller if you have problems with memory

!      counters

       integer    ::  i, j, n, m, nk, ni


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
       call out_band_circuit_size('BAND_LINES.DAT',iotape,ninterp,adot,  &
     &                  neig,nrk2,nlines,nvert)

       allocate(xk(nrk2))
       allocate(rk(3,nrk2))
       allocate(xcvert(nvert))
       allocate(ljump(nlines))
       allocate(nkstep(nlines))
       allocate(label(nvert+nlines))
       allocate(xklab(nvert+nlines))

       call out_band_get_circuit('BAND_LINES.DAT',iotape,ninterp,adot,   &
     &                  xk,rk,xcvert,ljump,nkstep,label,xklab,           &
     &                  neig,nrk2,nlines,nvert)


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

       irk = 0
       do n = 1,nlines
         if(ljump(n)) then

           irk = irk + 1
           do j=1,3
             rk0(j) = rk(j,irk)
           enddo

         endif

         do nk = 1,nkstep(n)

           irk = irk + 1
           do j=1,3
             rk1(j) = rk(j,irk)
           enddo

           do ni = 1,ninterp-1

             xw = (UM*ni) / (UM*ninterp)

             do j=1,3
               rkpt(j) = (UM-xw)*rk0(j) + xw*rk1(j)
             enddo

             call size_mtxd(emax,rkpt,adot,ng,kgv,nd)

             if(nd > mxddim) mxddim = nd

           enddo

           do j = 1,3
             rk0(j) = rk1(j)
           enddo

         enddo

       enddo

!      finds the total number of interpolated points so match_state deallocates correctly

       irkint = 0
       do n = 1,nlines
         if(ljump(n)) then
           irkint = irkint + 1
         endif
         do nk = 1,nkstep(n)
           do ni = 1,ninterp-1
             irkint = irkint + 1
           enddo
           irkint = irkint + 1
         enddo
       enddo

       nrk3 = irkint

        if(nignore >= neig) then
          write(6,*)
          write(6,*) '  out_band_lk:  asked to ignore ',nignore,'  bands'
          write(6,*) '  number of bands is: ',neig
          write(6,*)

          stop

        endif

!      allocates arrays


       allocate(ei(mxdbnd))
       allocate(ei_so(2*mxdbnd))
       allocate(hdiag(mxddim))
       allocate(isort(mxddim))
       allocate(qmod(mxddim))
       allocate(ekpg(mxddim))
       allocate(psi(mxddim,2*mxdbnd))
       allocate(ekpsi(mxdbnd))

       allocate(isort_so(2*mxddim))
       allocate(psi_so(2*mxddim,2*mxdbnd))

       allocate(psiold(mxddim,mxdbnd))
       allocate(imatch(mxdbnd))
       allocate(isold(mxddim))
       allocate(iperm(mxdbnd),ipermold(mxdbnd))

       allocate(psiold_so(2*mxddim,2*mxdbnd))
       allocate(imatch_so(2*mxdbnd))
       allocate(iperm_so(2*mxdbnd),ipermold_so(2*mxdbnd))


       allocate(ei0(mxddim))
       allocate(hdiag0(mxddim))
       allocate(isort0(mxddim))
       allocate(qmod0(mxddim))
       allocate(ekpg0(mxddim))
       allocate(ekpsi0(mxdbnd))

       allocate(psi0(mxddim,mxdbnd))
       allocate(hpsi0(mxddim,2*mxdbnd))

       allocate(ei1(mxddim))
       allocate(hdiag1(mxddim))
       allocate(isort1(mxddim))
       allocate(qmod1(mxddim))
       allocate(ekpg1(mxddim))
       allocate(ekpsi1(mxdbnd))

       allocate(psi1(mxddim,mxdbnd))

       allocate(e_of_k(neig,nrk2))
       allocate(e_of_k_so(2*neig,nrk2))

       allocate(e_of_k_int(neig-nignore,nrk3))
       allocate(e_of_k_so_int(2*(neig-nignore),nrk3))
       allocate(xk_int(nrk3))
       allocate(nkstep_int(nlines))

       do n = 1,nlines
         nkstep_int(n) = nkstep(n)*ninterp
       enddo


!       allocate(psi_ref(mxddim,2*mxdbnd))
       allocate(psi_all(mxddim,mxdbnd,2))
       allocate(isort_all(mxddim,2))

!       allocate(singval(2*mxdbnd))

       irk = 0
       irkint = 0

       do n = 1,nlines

!        loop over lines

         if(ljump(n)) then

           irk = irk + 1
           irkint = irkint + 1
           do j=1,3
             rk0(j) = rk(j,irk)
           enddo

           ipr = 1

           nocc = neig

           call h_kb_dia_all('pw  ', emax, rk0, neig, nocc,              &
     &     flgpsd, ipr, ifail, icmax, iguess, epspsi,                    &
     &     ng, kgv, phase, conj, ns, inds, kmax, indv, ek,               &
     &     sfact, veff, icmplx,                                          &
     &     nqnl, delqnl, vkb, nkb,                                       &
     &     ntype, natom, rat, adot,                                      &
     &     mtxd0, hdiag0, isort0, qmod0, ekpg0, .FALSE.,                 &
     &     psi0, hpsi0, ei0,                                             &
     &     vscr, kmscr,                                                  &
     &     latorb, norbat, nqwf, delqwf, wvfao, lorb,                    &
     &     mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,       &
     &     mxdbnd, mxdscr, mxdlao)

           do j = 1,neig
             e_of_k(j,irk) = ei0(j)
           enddo
           do j = 1,neig-nignore
             e_of_k_int(j,irkint) = ei0(j+nignore)
           enddo
           xk_int(irkint) = xk(irk)

           call kinetic_energy(neig,mtxd0,ekpg0,psi0,ekpsi0,             &
     &     mxddim,mxdbnd)

           ipr = 1
           nrka = -1
           call print_eig(ipr,irk,labelk,nrka,rk0,                       &
     &     mtxd0,icmplx,neig,psi0,                                       &
     &     adot,ei0,ekpsi0,isort0,kgv,                                   &
     &     mxddim,mxdbnd,mxdgve)

           call spin_orbit_perturb(rk0,mtxd0,isort0,                     &
     &     neig,psi0,ei0,ei_so,psi_so,.TRUE.,                            &
     &     ng,kgv,                                                       &
     &     nqnl,delqnl,vkb,nkb,                                          &
     &     ntype,natom,rat,adot,                                         &
     &     mxdtyp,mxdatm,mxdlqp,mxddim,mxdbnd,mxdgve)

           do j = 1,2*neig
             e_of_k_so(j,irk) = ei_so(j)
           enddo
           do j = 1,2*neig-2*nignore
             e_of_k_so_int(j,irkint) = ei_so(j+2*nignore)
           enddo

           call out_band_match(irkint, nlines, ljump, nkstep_int,        &
     &        neig-nignore, icmplx,XSC,                                  &
     &        mtxd0, psi0(:,nignore+1), isort0, psi_so(:,2*nignore+1),   &
     &        ei0(nignore+1), ei_so(2*nignore+1),                        &
     &        e_of_k_int, e_of_k_so_int, nrk3,                           &
     &        mxddim, mxdbnd-nignore)

         endif

!        end case of jump in band lines

         do nk = 1,nkstep(n)

           irk = irk + 1
           do j=1,3
             rk1(j) = rk(j,irk)
           enddo

!          reuses hpsi0
           nocc = neig

           call h_kb_dia_all('pw  ', emax, rk1, neig, nocc,              &
     &     flgpsd, ipr, ifail, icmax, iguess, epspsi,                    &
     &     ng, kgv, phase, conj, ns, inds, kmax, indv, ek,               &
     &     sfact, veff, icmplx,                                          &
     &     nqnl, delqnl, vkb, nkb,                                       &
     &     ntype, natom, rat, adot,                                      &
     &     mtxd1, hdiag1, isort1, qmod1, ekpg1, .FALSE.,                 &
     &     psi1, hpsi0, ei1,                                             &
     &     vscr, kmscr,                                                  &
     &     latorb, norbat, nqwf, delqwf, wvfao, lorb,                    &
     &     mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,       &
     &     mxdbnd, mxdscr, mxdlao)


           do j=1,neig
             e_of_k(j,irk) = ei1(j)
           enddo

           call kinetic_energy(neig,mtxd1,ekpg1,psi1,ekpsi1,             &
     &     mxddim,mxdbnd)

           ipr = 1
           nrka = -1
           call print_eig(ipr,irk,labelk,nrka,rk1,                       &
     &     mtxd1,icmplx,neig,psi1,                                       &
     &     adot,ei1,ekpsi1,isort1,kgv,                                   &
     &     mxddim,mxdbnd,mxdgve)

           mtxd_all(1) = mtxd0
           mtxd_all(2) = mtxd1

           neig_all(1) = neig-nignore
           neig_all(2) = neig-nignore

           do i=1,mtxd0
             isort_all(i,1) =  isort0(i)
           enddo
           do i=1,mtxd1
             isort_all(i,2) =  isort1(i)
           enddo

           do m=1,neig-nignore
           do i=1,mtxd0
             psi_all(i,m,1) =  psi0(i,m+nignore)
           enddo
           enddo
           do m=1,neig-nignore
           do i=1,mtxd1
             psi_all(i,m,2) = psi1(i,m+nignore)
           enddo
           enddo


           do ni = 1,ninterp-1
             irkint = irkint + 1
             xw = (UM*ni) / (UM*ninterp)
             xk_int(irkint) = (UM-xw)*xk(irk-1) + xw*xk(irk)


             xw_all(1) = UM-xw
             xw_all(2) = xw

             do j=1,3
               rkpt(j) = (UM-xw)*rk0(j) + xw*rk1(j)
               rkpt_all(j,1) = rk0(j)
               rkpt_all(j,2) = rk1(j)
             enddo


             call interpolation_glk(2, emax, neig-nignore, xsvd, csvd,   &
     &         xw_all, rkpt_all, mtxd_all, neig_all,                     &
     &         isort_all, psi_all,                                       &
     &         ei, psi, hpsi0, mtxd, isort, qmod, ekpg,                  &
     &         ng, kgv,                                                  &
     &         ntype, natom, rat, adot,                                  &
     &         nqnl, delqnl, vkb, nkb,                                   &
     &         vscr, kmscr,                                              &
     &         mxdtyp, mxdatm, mxddim, mxdlqp, mxdbnd, mxdgve, mxdscr)


             do j=1,neig-nignore
               e_of_k_int(j,irkint) = ei(j)
             enddo

             call spin_orbit_perturb(rkpt, mtxd, isort, neig-nignore,    &
     &       psi, ei, ei_so, psi_so, .TRUE.,                             &
     &       ng, kgv,                                                    &
     &       nqnl, delqnl, vkb, nkb,                                     &
     &       ntype, natom, rat, adot,                                    &
     &       mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

             do j=1,2*neig-2*nignore
               e_of_k_so_int(j,irkint) = ei_so(j)
             enddo

             call out_band_match(irkint, nlines, ljump, nkstep_int,      &
     &          neig-nignore, icmplx, XSC,                               &
     &          mtxd, psi, isort, psi_so, ei, ei_so,                     &
     &          e_of_k_int, e_of_k_so_int, nrk3,                         &
     &          mxddim, mxdbnd)

           enddo

!          stores the end point

           irkint = irkint + 1
           do j = 1,neig-nignore
             e_of_k_int(j,irkint) = ei1(j+nignore)
           enddo
           xk_int(irkint) = xk(irk)

           call spin_orbit_perturb(rk1,mtxd1,isort1,                     &
     &     neig,psi1,ei1,ei_so,psi_so,.TRUE.,                            &
     &     ng,kgv,                                                       &
     &     nqnl,delqnl,vkb,nkb,                                          &
     &     ntype,natom,rat,adot,                                         &
     &     mxdtyp,mxdatm,mxdlqp,mxddim,mxdbnd,mxdgve)

           do j = 1,2*neig
             e_of_k_so(j,irk) = ei_so(j)
           enddo
           do j=1,2*neig-2*nignore
             e_of_k_so_int(j,irkint) = ei_so(j+2*nignore)
           enddo

           call out_band_match(irkint, nlines, ljump, nkstep_int,        &
     &        neig-nignore, icmplx, XSC,                                 &
     &        mtxd1, psi1(:,nignore+1), isort1, psi_so(:,2*nignore+1),   &
     &        ei1(nignore+1), ei_so(2*nignore+1),                        &
     &        e_of_k_int, e_of_k_so_int, nrk3,                           &
     &        mxddim, mxdbnd-nignore)

!          keeps the result

           mtxd0 = mtxd1
           do j = 1,3
             rk0(j) = rk1(j)
           enddo
           do j = 1,mtxd1
             hdiag0(j) = hdiag1(j)
             isort0(j) = isort1(j)
             qmod0(j) = qmod1(j)
             ekpg0(j) = ekpg1(j)
           enddo
           do i = 1,neig
             ekpsi0(i) = ekpsi1(i)
             ei0(j) = ei1(j)
             do j=1,mtxd1
               psi0(j,i) = psi1(j,i)
             enddo
           enddo

         enddo

       enddo

!      writes the output files for xmgrace

       iotape = 15
       nstyle = 2

       call out_band_eref(neig,nrk2,ztot,efermi,2,1,e_of_k,eref,nocc)

       if(nignore /= 0) eref = ZERO


!        call out_band_gnuplot('band_lk.gp',iotape,                        &
!      &        neig,nrk2,xk,e_of_k,eref,                                  &
!      &        nvert,xcvert,nlines,ljump,nkstep,label,xklab)

       nstyle = 1

       call out_band_xmgrace('band_lk_ref.agr',iotape,                   &
     &        title,subtitle,nstyle,                                     &
     &        neig,nrk2,xk,e_of_k,eref,nocc,                             &
     &        nvert,xcvert,nlines,ljump,nkstep,label,xklab)

       nstyle = 2

!        use the same reference as e_of_k
!
!        call out_band_eref(neig-nignore,irkint,ztot,2,1,e_of_k_int,eref,nocc)
!
!        if(nignore /= 0) eref = ZERO

       call out_band_xmgrace('band_lk_int.agr', iotape,                  &
     &        title, subtitle, nstyle,                                   &
     &        neig-nignore, irkint, xk_int, e_of_k_int, eref, nocc,      &
     &        nvert, xcvert, nlines, ljump, nkstep_int, label, xklab)

!      uses the reference of e_of_k_so

       call out_band_eref(neig,nrk2,ztot,efermi,1,1,e_of_k_so,eref,nocc)

!        call out_band_eref(neig-nignore,irkint,ztot,1,1,e_of_k_so_int,eref,nocc)

       if(nignore /= 0) eref = ZERO

       call out_band_xmgrace('band_lk_int_so.agr', iotape,               &
     &    title, subtitle, nstyle,                                       &
     &    2*neig-2*nignore, irkint, xk_int, e_of_k_so_int, eref, nocc,   &
     &    nvert, xcvert, nlines, ljump, nkstep_int, label, xklab)



       deallocate(psiold)
       deallocate(imatch)
       deallocate(isold)
       deallocate(iperm,ipermold)

       deallocate(psiold_so)
       deallocate(imatch_so)
       deallocate(iperm_so,ipermold_so)

       deallocate(nkstep)
       deallocate(ljump)

       deallocate(xcvert)
       deallocate(xk)
       deallocate(rk)
       deallocate(label)
       deallocate(xklab)

       deallocate(e_of_k)
       deallocate(e_of_k_so)
       deallocate(e_of_k_int)
       deallocate(e_of_k_so_int)

       deallocate(vscr)

       deallocate(ei)
       deallocate(hdiag)
       deallocate(isort)
       deallocate(qmod)
       deallocate(ekpg)
       deallocate(psi)
       deallocate(ekpsi)

       deallocate(ei0)
       deallocate(hdiag0)
       deallocate(isort0)
       deallocate(qmod0)
       deallocate(ekpg0)
       deallocate(psi0)
       deallocate(ekpsi0)

       deallocate(ei1)
       deallocate(hdiag1)
       deallocate(isort1)
       deallocate(qmod1)
       deallocate(ekpg1)
       deallocate(psi1)
       deallocate(ekpsi1)

       deallocate(psi_so)

       return
       end subroutine out_band_glk
