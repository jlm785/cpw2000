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

!>     Calculates the self consistent potential and charge
!>     for kb pseudopotential with either an atomic orbital 
!>     or plane wave basis set.

       subroutine cpw_scf(flgaopw, iprglob, icmax, iguess, kmscr,        &
     &     efermi, elects, exc, strxc, ealpha, lkpg, lsafescf,           &
     &     dims_, crys_, flags_, pwexp_, recip_, acc_, xc_, strfac_,     &
     &     vcomp_, pseudo_, atorb_, kpoint_, hamallk_, psiallk_,         &
     &     total_, ewald_, chdens_) 

!        subroutine scf_kb_c16(flgaopw,pwexp_,acc_,flags_,                 &
!      & xc_,iprglob,iguess,                                               &
!      & kmscr,                                                            &
!      & recip_,                                                           &
!      & strfac_, vcomp_,                                                  &
!      & pseudo_,                                                          &
!      & atorb_,                                                           &
!      & crys_,                                                            &
!      & kpoint_,                                                          &
!      & hamallk_,psiallk_,                                                &
!      & efermi,elects,exc,strxc,total_,ealpha,ewald_,                     &
!      & chdens_,                                                          &
!      & dims_)

       
!      version 4.0. 19 october 1993. jlm
!      version 4.2  25 february 1999. jlm
!      modified 17 april 1999
!      modified (4f) 21 february 2000. jlm
!      modified (chd) 24 july 2002. jlm
!      modified (mixer_bfgs) 17 april 2006. jlm
!      modified (average potential) 28 september 2013. jlm
!      modified 19 november 2013. veffr(1). jlm
!      Modified c19, f90, etc. October 2015. JLM
!      Modified kmscr, 28 October 2015. JLM
!      Modified lkplusg, October 2018. JLM
!      Modified, documentations, 5 August 2019. JLM
!      Modified, types, December 2019.  JLM
!      Modified, dimension of ei, 18 February 2020. JLM
!      Modified h_kb_diag_all, 7 June 2020. JLM
!      Modified, new test of mixer failure. 28 November 2020. JLM

!      copyright inesc-mn/Jose Luis Martins

!      version 4.99
     
       use cpw_variables

       implicit none

!       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       type(dims_t)                       ::  dims_                      !<  array dimensions

!       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
!       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type
!       integer, intent(in)                ::  mxdlqp                     !<  array dimension for local potential
!       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
!       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
!       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars
!       integer, intent(in)                ::  mxdcub                     !<  array dimension for 3-index g-space
!       integer, intent(in)                ::  mxdnrk                     !<  size of k-points
!       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands
!       integer, intent(in)                ::  mxdlao                     !<  array dimension of orbital per atom type

       type(crys_t)                       ::  crys_                      !<  crystal structure

!       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in real space
!       integer, intent(in)                ::  ntype                      !<  number of types of atoms
!       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
!       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  lattice coordinates of atom j of type i

       type(flags_t)                      ::  flags_                     !<  computational flags

!       character(len=6), intent(in)       ::  flgpsd                     !<  type of pseudopotential
!       character(len=6), intent(in)       ::  flgscf                     !<  type of self consistent calculation '    PW','AO    ','AOJC  ','AOJCPW'
!NOT USED       character(len=6), intent(in)       ::  flgmix                     !<  type of scf mixing.  Unused at the moment

       type(pwexp_t)                      ::  pwexp_                     !<  plane-wave expansion choices

!       real(REAL64), intent(in)           ::  emax                       !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
!       real(REAL64), intent(in)           ::  teleck                     !<  electronic temperature (in Kelvin)
!       logical, intent(in)                ::  lkplusg                    !<  If true use the previous G-vectors (same mtxd and isort)

       type(recip_t)                      ::  recip_                     !<  reciprocal space information

!       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
!       integer, intent(in)                ::  kgv(3,dims_%mxdgve)              !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
!       complex(REAL64), intent(in)        ::  phase(dims_%mxdgve)              !<  real part of the phase factor of G-vector n
!       real(REAL64), intent(in)           ::  conj(dims_%mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase

!       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
!       integer, intent(in)                ::  inds(dims_%mxdgve)               !<  star to which g-vector n belongs
!       integer, intent(in)                ::  kmax(3)                    !<  max value of kgv(i,n)
!       integer, intent(in)                ::  indv(dims_%mxdcub)               !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax
!       integer, intent(in)                ::  mstar(dims_%mxdnst)              !<  number of G-vectors in the j-th star
!       real(REAL64), intent(in)           ::  ek(dims_%mxdnst)                 !<  kinetic energy (hartree) of g-vectors in star j

       type(xc_t)                         ::  xc_                        !<  exchange and correlation choice

!       character(len=*), intent(in)       ::  author                     !<  type of xc wanted (ca=pz , pbe, tbl)
!       real(REAL64), intent(in)           ::  tblaha                     !<  Tran-Blaha constant

       type(acc_t)                        ::  acc_                       !<  accuracy parameters

!       integer, intent(in)                ::  itmax                      !<  maximum number of iterations
!       real(REAL64), intent(in)           ::  epscv                      !<  convergence criteria for PW calculations
!       real(REAL64), intent(in)           ::  epscvao                    !<  convergence criteria for AO calculations
!       real(REAL64), intent(in)           ::  epspsi                     !<  requested precision of the eigenvectors

       type(kpoint_t)                     ::  kpoint_                    !<  k-point data

!       integer, intent(in)                ::  nrk                        !<  number of k-points for integration in the irreducible wedge of the brillouin zone
!       real(REAL64), intent(in)           ::  rk(3,dims_%mxdnrk)               !<  component in lattice coordinates of the k-point in the mesh
!       real(REAL64), intent(in)           ::  wgk(dims_%mxdnrk)                !<  weight in the integration of k-point

!       integer, intent(inout)             ::  nband(dims_%mxdnrk)              !<  number of bands for each k-points

       type(strfac_t)                     ::  strfac_                    !<  structure factors

!       complex(REAL64), intent(in)        ::  sfact(dims_%mxdtyp,dims_%mxdnst)       !<  real part of the structure factor
!       integer, intent(in)                ::  icmplx                     !<  indicates if the structure factor is complex

       type(pseudo_t)                     ::  pseudo_                    !<  pseudo-potential (Kleinman-Bylander)

!       CHECK THE CONSISTENCY   NQ <----> NQNL,  DELQ <----> DELQNL

!       integer, intent(in)                ::  nqnl(dims_%mxdtyp)               !<  number of points for the non-local pseudopotential interpolation
!       real(REAL64), intent(in)           ::  delqnl(dims_%mxdtyp)             !<  step used in the interpolation
!       real(REAL64), intent(in)  ::   vkb(-2:dims_%mxdlqp,0:3,-1:1,dims_%mxdtyp)     !<  (1/q**l) * kb nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
!       integer, intent(in)                ::  nkb(0:3,-1:1,dims_%mxdtyp)       !<   kb pseudo.  normalization for atom k, ang. mom. l

!       real(REAL64), intent(in)           ::  ztot                       !<  total charge density (electrons/cell)

       type(atorb_t)                      ::  atorb_                     !<  atomic orbitals in G-space

!       integer, intent(in)                ::  norbat(dims_%mxdtyp)             !<  number of atomic orbitals for atom k
!       integer, intent(in)                ::  nqwf(dims_%mxdtyp)               !<  number of points for wavefunction interpolation for atom k
!       real(REAL64), intent(in)           ::  delqwf(dims_%mxdtyp)             !<  step used in the wavefunction interpolation for atom k
!       real(REAL64), intent(in)      ::  wvfao(-2:dims_%mxdlqp,dims_%mxdlao,dims_%mxdtyp)  !<  wavefunction for atom k, ang. mom. l (normalized to vcell)
!       integer, intent(in)                ::  lorb(dims_%mxdlao,dims_%mxdtyp)        !<  angular momentum of orbital n of atom k

       type(enfrst_t)                     ::  total_                     !<  Total energy force stress

!       real(REAL64), intent(out)          ::  energy                     !<  total electronic(+ewald) energy.

       type(enfrst_t)                     ::  ewald_                     !<  Ewald energy force stress

!       real(REAL64), intent(in)           ::  enerew                     !<  Ewald energy
       
       type(vcomp_t)                      ::  vcomp_                     !<  Componemts of local potential

!       complex(REAL64), intent(in)        :: vion(dims_%mxdnst)                !<  ionic potential for the prototype G-vector in star j

!       complex(REAL64), intent(inout)     ::  vhar(dims_%mxdnst)               !<  Hartree potential for the prototype G-vector
!       complex(REAL64), intent(inout)     ::  vxc(dims_%mxdnst)                !<  exchange+correlation potential for the prototype G-vector
       
!       complex(REAL64), intent(out)       ::  veff(dims_%mxdnst)               !<  effective potential (local+Hartree+XC) for the prototype g-vector in star j

       type(chdens_t)                     ::  chdens_                    !<  charge densities    

!       complex(REAL64), intent(inout)     ::  den(dims_%mxdnst)                !<  total charge density for the prototype G-vector
!       complex(REAL64), intent(in)        ::  denc(dims_%mxdnst)               !<  core charge density for the prototype G-vector
       
       type(hamallk_t)                    ::  hamallk_                   !<  hamiltonian size and indexation for all k-points

!       integer, intent(inout)             ::  mtxd_allk(dims_%mxdnrk)          !<  dimension of the hamiltonian for k-point n
!       integer,  intent(inout)            ::  isort_allk(dims_%mxddim,dims_%mxdnrk)  !<  G-vector associated with k+G vector i of hamiltonian for k-point n

       type(psiallk_t)                    ::  psiallk_                   !<  psi for all k-points

!       complex(REAL64), intent(inout)     :: psi_allk(dims_%mxddim,dims_%mxdbnd,dims_%mxdnrk)  !<  eigenvectors for all k-points

!       real(REAL64), intent(out)          ::  eig_allk(dims_%mxdnrk*dims_%mxdbnd)    !<  eigenvalue j, for all the k-points
!       real(REAL64), intent(out)          ::  occ_allk(dims_%mxdnrk*dims_%mxdbnd)    !<  fractional ocupation of level j, for all the k-points
      
       
       real(REAL64), intent(in)           ::  ealpha                     !<  alpha term. (G=0)


       integer, intent(in)                ::  icmax                      !<  maximum value of outer iteration
       integer, intent(in)                ::  iprglob                    !<  Global printing level

       character(len=2), intent(in)       ::  flgaopw                    !<  type of calculation 'AO' or 'PW

!      input and output

       integer, intent(inout)             ::  kmscr(7)                   !<  max value of kgv(i,n) used for the potential FFT mesh and fft mesh size

       integer, intent(inout)             ::  iguess                     !<  tells if guess eigenvectors are available

       logical, intent(in)                ::  lkpg                       !<  If true use the previous G-vectors (same mtxd and isort)



!      output

       real(REAL64), intent(out)          ::  efermi                     !<  Fermi energy (or highest occupied state)
       real(REAL64), intent(out)          ::  elects                     !<  electronic temperature*entropy (hartree)
       real(REAL64), intent(out)          ::  exc                        !<  total XC energy (int rho * e_XC in Hartree).
       real(REAL64), intent(out)          ::  strxc(3,3)                 !<  contribution of xc to the stress tensor (contravariant,Hartree)

       logical,  intent(out)              ::  lsafescf                   !<  is true if self-consistency was reached.

!      local allocatable arrays

!       complex(REAL64), allocatable       ::  psi(:,:)
       complex(REAL64), allocatable       ::  hpsi(:,:)
       complex(REAL64), allocatable       ::  denk(:)

       complex(REAL64), allocatable       ::  vhxc(:)                    !  Hartre+xc potential that enters the diagonalization or after mixing
       complex(REAL64), allocatable       ::  vhxcout(:)                 !  Hartre+xc potential that comes out of the diagonalization
       complex(REAL64), allocatable       ::  delvhxc(:)                 !  After mixing vhxc = vhxc + delvhxc for next iteration
 
       complex(REAL64), allocatable       ::  rholap(:)
       complex(REAL64), allocatable       ::  tauk(:)
       complex(REAL64), allocatable       ::  tau(:)

       real(REAL64), allocatable          ::  hdiag(:)                   !  Hamiltonian diagonal for k+G-vector i
       real(REAL64), allocatable          ::  qmod(:)                    !  length of k+G-vector i
       real(REAL64), allocatable          ::  ekpg(:)                    !  kinetic energy (Hartree) of k+G-vector  i

       real(REAL64), allocatable          ::  ekl(:)                     !  kinetic energy of wave-function j, for all the k-points
 
       real(REAL64), allocatable          ::  ei(:)                      !  eigenvalue j
       real(REAL64), allocatable          ::  ekn(:)                     !  kinetic energy of wave-function j
       real(REAL64), allocatable          ::  occp(:)                    !  ocupation*weight*spin deg. of eigenvector j

       real(REAL64), allocatable          ::  vscr(:)                    !  screened potential in the FFT real space mesh

       complex(REAL64), allocatable       ::  vhxclow(:)                 !  Hartre+xc potential that entered the diagonalization that yelded the lowest energy
       complex(REAL64), allocatable       ::  vhxcoutlow(:)              !  Hartre+xc potential that came out of the diagonalization that yelded the lowest energy


!      other variables

       integer                ::  mxdscr                                 !  array dimension for vscr

       integer                ::  mxdupd                                 !  array dimension for Jacobian update of scf potential
       integer                ::  mxdscf                                 !  array dimension of old scf-iteration information

       real(REAL64)           ::  epsconv                                !  convergence criteria PW/AO

       integer                ::  ifail                                  !  if ifail=0 the subroutine was successfull. Otherwise ifail indicates the number of correct digits.
       integer                ::  minifail                               !  minimum value of non-zero ifail, or 100 if all ifail is 0
      
       real(REAL64)           ::  vmax, vmin                             !  maximum and minimum values of vscr
       real(REAL64)           ::  vmaxold, vminold                       !  old values of vmax,vmin
       integer                ::  nfailmix                               !  number of mixer failures detected
       integer                ::  maxnfailmix                            !  maximum number of corrections to mixer failures 

       real(REAL64)           ::  eband                                  !  band energy (Hartree)
       real(REAL64)           ::  ektot                                  !  electron kinetic energy (Hartree)
       real(REAL64)           ::  bandwid                                !  occupied band width estimate (Hartree)
       real(REAL64)           ::  penngap                                !  Penn gap estimate (Hartree)
       real(REAL64)           ::  rhovxc                                 !  correction for XC energy (int rho * v_XC in Hartree).
       real(REAL64)           ::  errvhxc                                !  maximum value of abs(vhxcout(i) - vhxc(i))

       integer                ::  neig, mtxd
       real(REAL64)           ::  rkpt(3)

       integer                ::  ipr, iconv, nrka
       integer                ::  ickin, idshift
       integer                ::  iel
       character(len=5)       ::  labelk

       real(REAL64)           ::  veffr1
       
       integer                ::  nsfft(3)
       integer                ::  mxdwrk

       real(REAL64)           ::  oldenergy                              !  energy of the previous iteration
       real(REAL64)           ::  enerlow                                !  lowest energy in previous scf iterations
       real(REAL64)           ::  eharrfou                               !  total energy of the Harris-Weinert-Foulkes functional
       
       integer                ::  itmix                                  !  iteration number for mixer, itmix = 1 resets.

       real(REAL64)           ::  tin, tout

       character(len=4)       ::  diag_type                              !  selects diagonalization, 'pw  ','ao  ','aojc'
       integer                ::  nocc

!      counters

       integer       ::  i, j
       integer       ::  irk, iter

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       real(REAL64), parameter :: DELTA = 0.2_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

!      You can reduce mxdupd,mxdscf in the unlikey case this takes too much
!      memory (linear in number of atoms, only relevant for very old machines)
!      mxdupd = max(dims_%mxdnst/10,min(dims_%mxdnst,100))

       mxdupd = dims_%mxdnst
       mxdscf = acc_%itmax

!      mixing maybe restarted if problems are detected, but not too many times

       nfailmix = 0
       maxnfailmix = max(1,mxdscf/10)

       allocate(ekl(dims_%mxdnrk*dims_%mxdbnd))

       allocate(hpsi(dims_%mxddim,dims_%mxdbnd))

       allocate(denk(dims_%mxdnst))

       allocate(vhxc(dims_%mxdnst))
       allocate(vhxcout(dims_%mxdnst))
       allocate(delvhxc(dims_%mxdnst))

       allocate(vhxclow(dims_%mxdnst))
       allocate(vhxcoutlow(dims_%mxdnst))

       allocate(ei(dims_%mxdbnd))

       allocate(hdiag(dims_%mxddim))
       allocate(qmod(dims_%mxddim))
       allocate(ekpg(dims_%mxddim))

       allocate(ekn(dims_%mxdbnd))
       allocate(occp(dims_%mxdbnd))

       call size_fft(kmscr,nsfft,mxdscr,mxdwrk)

       allocate(vscr(mxdscr))

       do i=1,recip_%ns
         vhxc(i) = vcomp_%vxc(i) + vcomp_%vhar(i)
       enddo

       if(flags_%flgscf /= '    PW' .and. flags_%flgscf /= 'AO    ' .and.              &
     &    flags_%flgscf /= 'AOJC  ' .and. flags_%flgscf /= 'AOJCPW') then

         write(6,*)
         write(6,'("    STOPPED in scf_kb_c16:   unknown type ",         &
     &       " of scf flag  ",a6)') flags_%flgscf

         stop

       endif


       if(flgaopw == 'PW') then
         epsconv = acc_%epscv
       elseif(flgaopw == 'AO') then
         epsconv = acc_%epscvao
       endif


       call zesec(tin)

!      start of self-consistent cycles
!      *******************************

       lsafescf = .FALSE.
       itmix = 1

       do iter = 1,acc_%itmax

!        do until convergence

         iel = 0

!        calculates local potential in fft mesh


         do i=1,recip_%ns
           vcomp_%veff(i) = vcomp_%vion(i) + vhxc(i)
         enddo

         ipr = 0
         if(iprglob > 1) ipr = 1
         if(iprglob == 4) ipr = 2

         idshift = 0

         call pot_local(ipr, vscr, vmax, vmin, vcomp_%veff, kmscr, idshift,        &
     &   recip_%ng, recip_%kgv, recip_%phase, recip_%conj,                         &
     &   recip_%ns, recip_%inds,                                                   &
     &   mxdscr, dims_%mxdgve, dims_%mxdnst)

!        checks if there is a problem, restart the convergence procedure

         if(itmix > 2) then
           if( (abs(vmin-vminold) > DELTA .or. abs(vmax-vmaxold) > DELTA)   &
     &            .and. nfailmix < maxnfailmix) then

             write(6,*)
             write(6,'("   Restarted in cpw_scf:  cycle is diverging",    &
     &            " local potential jumped by",g14.6)')                   &
     &            max(abs(vmin-vminold),abs(vmax-vmaxold))
             itmix = 1
             nfailmix = nfailmix + 1

             do i=1,recip_%ns
               vhxc(i) = vhxclow(i)
             enddo

             do i=1,recip_%ns
               vcomp_%veff(i) = vcomp_%vion(i) + vhxc(i)
             enddo

             call pot_local(ipr, vscr, vmax, vmin, vcomp_%veff, kmscr, idshift,      &
     &           recip_%ng, recip_%kgv, recip_%phase, recip_%conj,                   &
     &           recip_%ns, recip_%inds,                                             &
     &           mxdscr, dims_%mxdgve, dims_%mxdnst)

           endif
         endif

         vmaxold = vmax
         vminold = vmin

         minifail = 100

         do irk = 1,kpoint_%nrk

!          loop over k-points

           rkpt(1) = kpoint_%rk(1,irk)
           rkpt(2) = kpoint_%rk(2,irk)
           rkpt(3) = kpoint_%rk(3,irk)
!
           neig = kpoint_%nband(irk)

!          calculates hamiltonian and diagonalizes
!          ***************************************

!          past first iteration one has a guess of the eigenvalues

           if(iter > 1) iguess = 1

           mtxd = hamallk_%mtxd_allk(irk)

           if(flgaopw == 'PW') then

             ipr = 0
             if(iprglob > 2) ipr = 1

             diag_type = 'pw  '
             nocc = neig

             call h_kb_dia_all(diag_type, pwexp_%emax, rkpt, neig, nocc, &
     &         flags_%flgpsd, ipr, ifail, icmax, iguess, acc_%epspsi,    &
     &         recip_%ng, recip_%kgv, recip_%phase, recip_%conj,         &
     &         recip_%ns, recip_%inds, recip_%kmax,                      &
     &         recip_%indv, recip_%ek,                                   &
     &         strfac_%sfact, vcomp_%veff, strfac_%icmplx,               &
     &         pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,       &
     &         crys_%ntype,crys_%natom,crys_%rat,crys_%adot,             &
     &         mtxd, hdiag, hamallk_%isort_allk(:,irk),                  &
     &         qmod, ekpg, lkpg,                                         &
     &         psiallk_%psi_allk(:,:,irk), hpsi, ei,                     &
     &         vscr, kmscr,                                              &
     &         atorb_%latorb, atorb_%norbat, atorb_%nqwf,                &
     &         atorb_%delqwf, atorb_%wvfao, atorb_%lorb,                 &
     &         dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,   &
     &         dims_%mxdcub, dims_%mxdlqp, dims_%mxddim, dims_%mxdbnd,   &
     &         mxdscr, dims_%mxdlao)

             if(ifail /= 0) minifail = min(ifail,minifail)
      
             if(ifail < -3) then
               write(6,*)
               write(6,'("   Stopped in cpw_scf:  cycle is diverging",    &
     &            " negative number of accuracy digits",i5)') ifail

               stop

             endif

           elseif(flgaopw == 'AO') then
 
             veffr1 = real(vcomp_%veff(1),REAL64)
             nocc = neig

             if(flags_%flgscf == 'AOJCPW') diag_type = 'aojc'
             if(flags_%flgscf == 'AOJC  ') diag_type = 'aojc'
             if(flags_%flgscf == 'AO    ') diag_type = 'ao  '

             call h_kb_dia_all(diag_type, pwexp_%emax, rkpt, neig, nocc, &
     &         flags_%flgpsd, ipr, ifail, icmax, iguess, acc_%epspsi,    &
     &         recip_%ng, recip_%kgv, recip_%phase, recip_%conj,         &
     &         recip_%ns, recip_%inds, recip_%kmax,                      &
     &         recip_%indv, recip_%ek,                                   &
     &         strfac_%sfact, vcomp_%veff, strfac_%icmplx,               &
     &         pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,       &
     &         crys_%ntype,crys_%natom,crys_%rat,crys_%adot,             &
     &         mtxd, hdiag, hamallk_%isort_allk(:,irk),                  &
     &         qmod, ekpg, lkpg,                                         &
     &         psiallk_%psi_allk(:,:,irk), hpsi, ei,                     &
     &         vscr, kmscr,                                              &
     &         atorb_%latorb, atorb_%norbat, atorb_%nqwf,                &
     &         atorb_%delqwf, atorb_%wvfao, atorb_%lorb,                 &
     &         dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,   &
     &         dims_%mxdcub, dims_%mxdlqp, dims_%mxddim, dims_%mxdbnd,   &
     &         mxdscr, dims_%mxdlao)

           else

             write(6,*)
             write(6,'("     STOPPED in scf_kb_c16:   unknown type",     &
     &                  " of basis   ",a2)') flgaopw

             stop

           endif


           hamallk_%mtxd_allk(irk) = mtxd
           kpoint_%nband(irk) = neig


!          end of diagonalization


!          calculates the kinetic energy

           call kinetic_energy(neig,mtxd,ekpg,psiallk_%psi_allk(:,:,irk),EKN,     &
     &     dims_%mxddim,dims_%mxdbnd)

!          prints the eigensolutions

           ipr = 0
           if(iprglob == 3) ipr = 1
           if(iprglob == 4) ipr = 2
           
           nrka = -1

           ickin = 1

           call print_eig(ipr,irk,labelk,nrka,rkpt,                      &
     &     mtxd,ickin,neig,psiallk_%psi_allk(:,:,irk),                   &
     &     crys_%adot,ei,EKN,hamallk_%isort_allk(:,irk),recip_%kgv,      &
     &     dims_%mxddim,dims_%mxdbnd,dims_%mxdgve)

!          stores eigenvalues and kinetic energies

           do j=1,neig
             iel = iel + 1
             psiallk_%eig_allk(iel) = ei(j)
             ekl(iel) = ekn(j)
           enddo

!          end of loop over k-points

         enddo

!        finds the fermi level

         ipr = iprglob
         
         call fermi_level(psiallk_%eig_allk,pseudo_%ztot,pwexp_%teleck,        &
     &   kpoint_%nrk,kpoint_%wgk,kpoint_%nband,                                &
     &   psiallk_%occ_allk,efermi,eband,elects,bandwid,penngap,                &
     &   dims_%mxdnrk,dims_%mxdbnd) 

         call print_fermi_occup(ipr,psiallk_%eig_allk,pwexp_%teleck,           &
     &   kpoint_%nrk,kpoint_%wgk,kpoint_%nband,                                &
     &   psiallk_%occ_allk,efermi,eband,elects,bandwid,penngap,                &
     &   dims_%mxdnrk,dims_%mxdbnd) 

!       calculates the Harris-Foulkes functional energy

        if(itmix /= 1) then
          call harris_weinert_foulkes(eharrfou, eband, exc, ewald_%energy,    &
              recip_%ns, recip_%mstar, recip_%ek,                             &
              vhxc, chdens_%den,                                              &
              crys_%adot,                                                     &
              dims_%mxdnst)
        endif

!        initializes charge density

         do i=1,recip_%ns
           chdens_%den(i) = C_ZERO
         enddo

         allocate(tau(dims_%mxdnst))

        if(xc_%author == "TBL") then
          do i = 1,recip_%ns
            tau(i) = C_ZERO
          enddo
        endif         

!       second loop over k-points


         ektot = zero
         iel = 0
         do irk = 1,kpoint_%nrk

           rkpt(1) = kpoint_%rk(1,irk)
           rkpt(2) = kpoint_%rk(2,irk)
           rkpt(3) = kpoint_%rk(3,irk)
         
           neig = kpoint_%nband(irk)
           mtxd = hamallk_%mtxd_allk(irk)

!          adds to sum of occupied eigenvalues and kinetic energy

           do j = 1,neig
             iel = iel + 1
             occp(j) = 2*kpoint_%wgk(irk)*psiallk_%occ_allk(iel)
             ektot = ektot + occp(j)*ekl(iel)
           enddo

!          adds to total charge density


           call charge_by_fft(mtxd, neig, occp,                            &
     &     hamallk_%isort_allk(:,irk), psiallk_%psi_allk(:,:,irk), denk,   &
     &     recip_%ng, recip_%kgv, recip_%phase , recip_%conj, recip_%ns,   &
     &     recip_%inds, recip_%kmax, recip_%mstar,                         &
     &     dims_%mxddim, dims_%mxdbnd, dims_%mxdgve, dims_%mxdnst)

           if(xc_%author == "TBL") then

             allocate(tauk(dims_%mxdnst))

             call tau_by_fft(tauk, mtxd, neig, occp,                     &
     &       hamallk_%isort_allk(:,irk),psiallk_%psi_allk(:,:,irk),      &
     &       rkpt,crys_%adot,                                            &
     &       recip_%ng,recip_%kgv,recip_%phase,recip_%conj,recip_%ns,    &
     &       recip_%inds,recip_%kmax,recip_%mstar,                       &
     &       dims_%mxddim,dims_%mxdbnd,dims_%mxdgve,dims_%mxdnst)

             do i = 1,recip_%ns
               tau(i) = tau(i) + tauk(i)
             enddo
     
             deallocate(tauk)
     
           endif

           do i=1,recip_%ns
             chdens_%den(i) = chdens_%den(i) + denk(i)
           enddo

         enddo

!        end of loop over k-points

         allocate(rholap(dims_%mxdnst))

         if(xc_%author == "TBL") then
         
           call lap_rho(chdens_%den, rholap, crys_%adot,                 &
     &     recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns,  &
     &     recip_%inds, recip_%mstar,                                    &
     &     dims_%mxdgve, dims_%mxdnst)
     
          
         endif

!        calculates the screening potential

         ipr = 0
         if(iprglob > 1) ipr = 1
         if(iprglob == 4) ipr = 2
         
         call v_hartree_xc(ipr, xc_%author, xc_%tblaha,                  &
     &   crys_%adot, exc, strxc, rhovxc,                                 &
     &   vcomp_%vhar, vcomp_%vxc, chdens_%den, chdens_%denc, rholap, tau,   &
     &   recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns,    &
     &   recip_%inds, recip_%kmax, recip_%mstar, recip_%ek,              &
     &   dims_%mxdgve, dims_%mxdnst)
     
         deallocate(tau)
         deallocate(rholap)
        
         do i=1,recip_%ns
           vhxcout(i) = vcomp_%vxc(i) + vcomp_%vhar(i)
         enddo

!        calculates the total energy and checks for convergence

         ipr = 0
         if(iprglob == 1) ipr = -5
         if(iprglob == 2) ipr = -20
         if(iprglob > 2) ipr = 20

         if(iter /= 1) oldenergy = total_%energy

         call total_ks_energy(ipr, strfac_%icmplx, iter, itmix, eharrfou,   &
     &   iconv, errvhxc, epsconv,                                           &
     &   total_%energy, eband, ektot, exc, ealpha, ewald_%energy,           &
     &   recip_%ng, recip_%kgv, recip_%ns, recip_%mstar, recip_%ek,         &
     &   vcomp_%vion, vhxc, vhxcout, chdens_%den,                           &
     &   pseudo_%ztot, crys_%adot,                                          &
     &   dims_%mxdgve, dims_%mxdnst)

         if(itmix == 1) then

           enerlow = total_%energy
           do i = 1,recip_%ns
             vhxclow(i) = vhxc(i)
           enddo
           do i = 1,recip_%ns
             vhxcoutlow(i) = vhxcout(i)
           enddo

         else

           if(total_%energy < enerlow) then
             enerlow = total_%energy
             do i = 1,recip_%ns
               vhxclow(i) = vhxc(i)
             enddo
             do i = 1,recip_%ns
               vhxcoutlow(i) = vhxcout(i)
             enddo
           else
             if(oldenergy < total_%energy - 1.0) then
               write(6,*)
               write(6,'("   Restarted in cpw_scf:  cycle is diverging",   &
     &         " energy increased by",g14.6)') total_%energy - oldenergy

               itmix = 1
               do i = 1,recip_%ns
                 vhxc(i) = vhxclow(i)
               enddo
               do i = 1,recip_%ns
                 vhxcout(i) = vhxcoutlow(i)
               enddo

             endif
           endif

         endif

!        accelarates self-consistency

         ipr = 0
         if(iprglob == 4) ipr = 1

!         init = 0
!         if(iter == 1) init = 2

!        call mixer_broyden1(iter,adot,ztot,

         call mixer_bfgs_c16(itmix, crys_%adot, pseudo_%ztot,            &
     &   bandwid, penngap, total_%energy,                                &
     &   recip_%ng, recip_%phase, recip_%conj, recip_%ns,                &
     &   recip_%mstar, recip_%ek,                                        &
     &   vhxc, vhxcout, delvhxc,                                         &
     &   dims_%mxdgve, dims_%mxdnst, mxdupd, mxdscf)
     
         do i = 2,recip_%ns
           vhxc(i) = vhxc(i) + delvhxc(i)
         enddo


!        average potential

         vhxc(1) = vhxcout(1)


         call zesec(tout)
         if(iprglob > 1) then
           write(6,'("  computing time for iteration ",i5,3x,f12.3)')    &
     &          iter,tout-tin
         endif
         tin = tout

         if(iconv == 1 .and. iter > 1) exit

!        end of self consistent loop

         itmix = itmix + 1

       enddo
       
       deallocate(ekl)

       deallocate(hpsi)
       
       deallocate(denk)
       
       deallocate(vhxc)
       deallocate(vhxcout)
       deallocate(delvhxc)
       
       deallocate(vhxclow)
       deallocate(vhxcoutlow)

       deallocate(ei)

       deallocate(hdiag)
       deallocate(qmod)
       deallocate(ekpg)

       deallocate(ekn)
       deallocate(occp)

       deallocate(vscr)

!      if self-consistency not reached

       if(iconv == 0) then
         write(6,*)
         write(6,'("     scf_kb_c16   maximum number of iterations ",    &
     &        "exceeded: ",i6)') acc_%itmax
         write(6,*)
       else
         if(minifail == 100) then
           lsafescf = .TRUE.
         elseif(minifail > 3) then
           lsafescf = .TRUE.
           write(6,*)
           write(6,'("   WARNING:  diagonalization had only ",i3,        &
     &         " accurate digits")') minifail 
           write(6,*)
         endif
       endif
       
       return
       end subroutine cpw_scf

