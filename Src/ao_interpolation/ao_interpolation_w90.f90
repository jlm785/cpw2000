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

!>  Writes the files for wannier90 (wan.win, wan.amn, wan.mmn, wan.eig)
!>
!>  \author       Carlos Loia Reis
!>  \version      5.11
!>  \date         May 2020, 6 October 2024.
!>  \copyright    GNU Public License v2

subroutine ao_interpolation_w90(mtb,                                     &
      emax, ztot,                                                        &
      adot, ntype, natom, nameat, rat,                                   &
      ng, kgv,                                                           &
      icmplx,                                                            &
      norbat, nqwf, delqwf, wvfao, lorb,                                 &
      mxdtyp, mxdatm, mxdgve, mxdlqp, mxdlao)

! Written by Carlos Loia Reis, May 2020, based on previous code,
! ao_interpolation_out_band_fold.
! Modified, documentation, May 2020. JLM


  use NonOrthoInterp


  implicit none

  integer, parameter                 :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  type(noiData_t)                    ::  mtb                             !<  see NonOrthoInterp

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2)                   ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  icmplx                          !<  indicates if the structure factor is complex

  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  (1/q**l) * wavefunction for atom k, ang. mom. l (unnormalized to vcell)

! allocatable arrays

  complex(REAL64), allocatable       ::  psi_ao(:,:)
  complex(REAL64), allocatable       ::  ao_basis(:,:)
  complex(REAL64), allocatable       ::  ao_basis_so(:,:)

  complex(REAL64), allocatable       ::  ao_basis_ortho(:,:)


  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i

  integer, allocatable               ::  infolcao(:,:)                   !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)

  complex(REAL64), allocatable       ::  S(:,:)
  complex(REAL64), allocatable       ::  S12(:,:)
  complex(REAL64), allocatable       ::  S12_inv(:,:)
  complex(REAL64), allocatable       ::  Swrk(:,:)
  real(REAL64),    allocatable       ::  ev_wrk(:)

  complex(REAL64), allocatable       ::  prod(:,:)                       !  <bas|psi>

  integer, allocatable               :: nnlist(:, :)                          !
  integer, allocatable               :: nncell(:, :, :)                       !*

  real(REAL64),allocatable           :: kpt_lattice(:,:)

  integer, allocatable               :: isort2(:)

  complex(REAL64), allocatable       :: mmnk(:,:)

  complex(REAL64), allocatable       :: psi_ao2(:,:)

! allocatable arrays for Brillouin zone path

  integer                            ::  nlines                          !  number of lines in reciprocal space
  real(REAL64), allocatable          ::  e_of_k(:,:)                     !  band energies of k-point in plot
  real(REAL64), allocatable          ::  e_of_k_so(:,:)                  !  spin-orbit band energies of k-point in plot

! local variables

  integer                            ::  mxddim                          !  array dimension for the hamiltonian
  real(REAL64)                       ::  bdot(3,3),vcell
  integer                            ::  mtxd                            !  dimension of the hamiltonian
  integer                            ::  nd, mxdorb, lmax, k
  logical                            ::  lnewanl

  logical                            ::  true_pkn

  integer           ::  neig                                             !  number of eigenvectors required (maybe modified on output)
  real(REAL64)      ::  rkpt(3)                                          !  j-th component in lattice coordinates of the k-point

  integer           ::  nstyle                                           !  choice of plot style

  integer           ::  irk
  integer           ::  iotape
  integer           ::  nrk2

  real(REAL64)      ::  avec(3,3),bvec(3,3)

  real(REAL64), allocatable :: ev_interp(:)

  logical       ::  lkplusg                                         !  If true use the previous G-vectors (same mtxd and isort)

  integer          ::  iorb, jorb, iband, jband

  integer          ::  ir_size_wf


  integer irk2
  integer mtxd2

  real(REAL64)   :: rkpt2(3)


  integer           :: soc
  real(REAL64)      :: e_fermi

  integer :: numoccupied, num_atoms, num_kpts, norb
  integer :: ikpt

  integer ndum
  real(REAL64) k_start(3), k_end(3)
  character(len=6) label_start, label_mid, label_end

  character(len=80) sys_str
  integer, parameter            :: num_nnmax = 12

  integer                       :: nntot                                 !

! counters

  integer    ::   i, j

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
  real(REAL64), parameter     ::  EV = 27.21138505_REAL64


! Start A) First things first. Produce a Wannier90 input file called wan.win

  open (unit=22, file="wan.win", form="formatted")


  write(22,*) "#----------------------------------"
  write(22,*) "#restart=plot"
  write(22,*) "num_iter = 500"
  write(22,*) "bands_plot = true"
  write(22,*) "bands_plot_format = xmgrace"
  write(22,*) "dos=true"
!  write(22,*) "dos_plot_format = xmgrace"
  write(22,*) "kmesh = 25"
  write(22,*) "adpt_smr_max=0.1 eV"
  write(22,*) "#----------------------------------"


! input needed here

  if(mtb%lso== 1) then
    soc = 1
    numoccupied = nint(ztot)
  else
    soc = 0
    numoccupied = nint(ztot)/2
  endif

  e_fermi = mtb %eref * EV

  norb = mtb%nband


  write(22,'("num_wann", 4x, i5)') mtb%nband
  write(22,'("num_bands", 4x, i5)') mtb%nband


  write(22,'("mp_grid", 4x, 3i5)') mtb%fiData%mp_grid(:)


  call adot_to_avec_sym(adot,avec,bvec)

  write(22,'("begin unit_cell_cart")')
  write(22,'("bohr")')
  write(22,'(3f14.8)') (avec(k,1),k=1,3)
  write(22,'(3f14.8)') (avec(k,2),k=1,3)
  write(22,'(3f14.8)') (avec(k,3),k=1,3)
  write(22,'("end unit_cell_cart")')
  write(22,*)

  num_atoms = 0
  do i=1, ntype
  do j=1, natom(i)
    num_atoms = num_atoms + 1
  enddo
  enddo

  write(22,'("begin atoms_frac")')
  do i=1, ntype
  do j=1, natom(i)
      write(22,'(a2, 3f14.8)') nameat(i), rat(:,j,i)
  enddo
  enddo
  write(22,'("end atoms_frac")')

  write(22,*)


  label_start='      '
  label_mid  ='      '
  label_end  ='      '

! band lines
  open(unit=44, file="BAND_LINES.DAT", form="formatted")
  read(44,*) nlines, ndum

  write(22,*)
  write(22,'("begin kpoint_path")')
  do i=1, nlines
    read(44,*) k_start(1),k_start(2),k_start(3), k_end(1),k_end(2),k_end(3), ndum, label_start, label_mid, label_end
    write(22,'(a6,2x,3f14.8,2x,a6,2x, 3f14.8)') label_start, k_start, label_end, k_end
  enddo
  write(22,'("end kpoint_path")')

  close(unit=44)


  num_kpts = mtb%fiData%mp_grid(1)*mtb%fiData%mp_grid(2)*mtb%fiData%mp_grid(3)

  mtb%fiData%num_kpts = num_kpts

  allocate(kpt_lattice(3,num_kpts))


  write(22,*)
  write(22,'("begin kpoints")')

  ikpt=1
  do i=1,mtb%fiData%mp_grid(1)
  do j=1,mtb%fiData%mp_grid(2)
  do k=1,mtb%fiData%mp_grid(3)

    rkpt(1) = (i-1)*1.0D0/mtb%fiData%mp_grid(1)
    rkpt(2) = (j-1)*1.0D0/mtb%fiData%mp_grid(2)
    rkpt(3) = (k-1)*1.0D0/mtb%fiData%mp_grid(3)

    kpt_lattice(:,ikpt) = rkpt(:)
    mtb%fiData%kpt_latt = kpt_lattice

    write(22,'(3f14.8)')  rkpt(:)
    ikpt=ikpt+1

  enddo
  enddo
  enddo

  write(22,'("end kpoints")')
  write(22,*)

  close(unit=22)

! wan.win was written, now call wannier90.x to produce the wan.nnkp file and read the relevant information

! RELEVANT CODE SHOULD BE EXTRACTED FROM WANNIER90 TO SIMPLIFY THIS...

  call execute_command_line("rm ./wan.nnkp  ; wannier90.x -pp wan")

  call execute_command_line('grep -A1 "begin nnkpts" wan.nnkp | tail -n 1 > inf.txt')

  open(unit=33, file="inf.txt", form="formatted")

  read(33, *) nntot

  close(33)

  write(sys_str,'("grep -A", i5,  " ""begin nnkpts"" wan.nnkp > inf.txt")') (nntot*num_kpts +2)

  write(6,*) "sys_str = " , sys_str
  write(6,*)

  call execute_command_line(sys_str)

  open(unit=33, file="inf.txt", form="formatted")

  read(33, *) sys_str
  read(33, *) nntot

  write(*,*) "nntot is" , nntot
  write(6,*)

  allocate(nnlist(num_kpts,num_nnmax))
  allocate(nncell(3,num_kpts,num_nnmax))

  do ikpt=1, num_kpts
    do i=1, nntot
      read(33,*) irk, nnlist(ikpt,i), nncell(1,ikpt,i), nncell(2,ikpt,i), nncell(3,ikpt,i)
    enddo
  enddo

  close(33)

! Now we have the k-point info desired by wannier90


  neig = norb

  true_pkn = .true.

  nrk2=num_kpts

  open(unit = 100, file = "wan.mmn", form="formatted")
  open(unit = 101, file = "wan.amn", form="formatted")

  write(100,*) "autogenerated by mtb"
  write(101,*) "autogenerated by mtb"

  write(100,'(3i5)') neig, num_kpts, nntot
  write(101,'(3i5)') neig, num_kpts, norb

  open(unit=102,file="wan.eig", form="formatted")

  allocate(e_of_k(neig,nrk2))
  allocate(e_of_k_so(2*neig,nrk2))

! finds mxddim, mxdbnd

!  mxdbnd = neig
!------------------------------------------------------------------
  allocate(ev_interp(mtb%nband))
!------------------------------------------------------------------

!
!---------- atomic orbital preparation section---------------------

  call adot_to_bdot(adot,vcell,bdot)

! finds mxdddim

  mxddim = 1
  do irk=1,nrk2
    do j=1,3
      rkpt(j) = kpt_lattice(j,irk)                                       ! computed in folded k-point
    enddo
    call size_mtxd(emax,rkpt,adot,ng,kgv,nd)
    if(nd > mxddim) mxddim = nd
  enddo

! WHY IS THIS HERE?  THESE ARE VERY BIG ARRAYS!

  mxddim = 4*mxddim

! finds mxdorb

  mxdorb = 0
  lmax = 0
  do k=1,ntype
    do j=1,norbat(k)
      mxdorb = mxdorb + (2*lorb(j,k)+1)*natom(k)
      lmax = max(lorb(j,k),lmax)
    enddo
  enddo


  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(isort2(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))


! allocate wavefunctions needded for pkn calculation

  allocate(ao_basis(mxddim,mxdorb))

  if(mtb%lso==1) then

    allocate(ao_basis_so(2*mxddim,2*mxdorb))
    allocate(ao_basis_ortho(2*mxddim,2*mxdorb))
    allocate(psi_ao(2*mxddim,2*mxdorb))
    allocate(psi_ao2(2*mxddim,2*mxdorb))

    allocate(prod(2*mxdorb,2*mxdorb))

    allocate(S(2*mxdorb,2*mxdorb))
    allocate(S12(2*mxdorb,2*mxdorb))
    allocate(S12_inv(2*mxdorb,2*mxdorb))
    allocate(Swrk(2*mxdorb,2*mxdorb))
    allocate(ev_wrk(2*mxdorb))
    allocate(mmnk(neig,neig))

    allocate(infolcao(5,mxdorb))

  else

    allocate(ao_basis_ortho(mxddim,mxdorb))
    allocate(psi_ao(mxddim,mxdorb))
    allocate(psi_ao2(mxddim,mxdorb))

    allocate(prod(mxdorb,mxdorb))

    allocate(S(mxdorb,mxdorb))
    allocate(S12(mxdorb,mxdorb))
    allocate(S12_inv(mxdorb,mxdorb))
    allocate(Swrk(mxdorb,mxdorb))
    allocate(ev_wrk(mxdorb))

    allocate(mmnk(neig,neig))

    allocate(infolcao(5,mxdorb))

  endif


  inquire(iolength = ir_size_wf)  irk, rkpt,  psi_ao(:,:) , mtxd
  open(unit = 88, file ="wf_bnd_rec.dat", access="direct", recl=ir_size_wf)

  write(6,*)
  write(6,*) '   loop over k-points'
  write(6,*)

  do irk=1,nrk2

!   loop over k-points

    do j=1,3
      rkpt(j) = kpt_lattice(j,irk)
    enddo

    write(6,'(i5," of ", i5, 3f8.5)') irk,nrk2, rkpt(1),rkpt(2),rkpt(3)


    call NonOrthoInterpRun(mtb,rkpt,ev_interp)

    lkplusg = .FALSE.

    call hamilt_struct(emax, rkpt, mtxd, isort, qmod, ekpg,  lkplusg,    &
        ng, kgv, adot,                                                   &
        mxdgve, mxddim)

    lnewanl = .TRUE.

    call atomic_orbital_c16(rkpt, mtxd, isort, icmplx,                   &
        mxdorb, ao_basis, infolcao,                                      &
        ng, kgv,                                                         &
        norbat, nqwf, delqwf, wvfao, lorb,                               &
        ntype, natom, rat, adot,                                         &
        mxdtyp, mxdatm, mxdlqp, mxddim, mxdorb, mxdgve, mxdlao)


    if (irk==1) then
       do iorb = 1, mxdorb
         write(*,'("infolcao",6i6)') iorb, (infolcao(j,iorb),j=1,5)
       enddo
    endif

    if(mtb%lso==1) then

      do i=1,mxdorb
      do j=1,mtxd
        ao_basis_so(2*j-1,2*i  ) = ao_basis(j,i)
        ao_basis_so(2*j  ,2*i  ) = C_ZERO
        ao_basis_so(2*j-1,2*i-1) = C_ZERO
        ao_basis_so(2*j  ,2*i-1) = ao_basis(j,i)
      enddo
      enddo

! mtb%wrk has the eigenvectors corresponding to ev_interp

      call zgemm('N', 'N', 2*mtxd, norb, norb, C_UM, ao_basis_so,        &
           2*mxddim, mtb%wrk, 2*mxdorb, C_ZERO, psi_ao, 2*mxddim)

      call ao_SortToGammaSo(psi_ao, isort, mtxd, neig, psi_ao2, 2*mxddim)

      write(88,rec=irk)  irk, rkpt, psi_ao2(:,:), 2*mtxd

      call zgemm('C', 'N', norb, norb, 2*mtxd, C_UM, ao_basis_so,        &
          2*mxddim, ao_basis_so, 2*mxddim, C_ZERO, S, 2*mxdorb)

      call  GetS12(S, S12, S12_inv, Swrk, ev_wrk, norb)

      call zgemm('n', 'n', 2*mtxd, norb, norb, C_UM, ao_basis_so,        &
          2*mxddim,  S12, 2*mxdorb, C_ZERO, ao_basis_ortho, 2*mxddim)

      call zgemm('C', 'N', norb, neig, 2*mtxd, C_UM, ao_basis_ortho,     &
          2*mxddim, psi_ao, 2*mxddim, C_ZERO, prod, 2*mxdorb)

      psi_ao(:,:) = psi_ao2(:,:)

      do iband=1, neig
      do iorb =1, norb
        write(101,'(3i5,2f14.8)') iband, iorb, irk, prod(iorb,iband)
      enddo
      enddo

      do iband=1, neig
        write(102,'(2i5,f14.8)') iband, irk, ev_interp(iband)*27.21138505
      enddo

    else  ! NON SPIN ORBIT CASE

      call zgemm('N', 'N', mtxd, norb, norb, C_UM, ao_basis, mxddim,     &
           mtb%wrk, mxdorb, C_ZERO, psi_ao, mxddim)


      call ao_SortToGamma(psi_ao, isort, mtxd, neig, psi_ao2, mxddim)


      write(88,rec=irk)  irk, rkpt, psi_ao2(:,:), mtxd


      call zgemm('C', 'N', norb, norb, mtxd, C_UM, ao_basis,             &
          mxddim, ao_basis, mxddim, C_ZERO, S, mxdorb)

      call GetS12(S, S12, S12_inv, Swrk, ev_wrk, norb)


      call zgemm('n', 'n', mtxd, norb, norb, C_UM, ao_basis,             &
          mxddim,  S12, mxdorb, C_ZERO, ao_basis_ortho, mxddim)


      call zgemm('C', 'N', norb, neig, mtxd, C_UM, ao_basis_ortho,       &
          mxddim, psi_ao, mxddim, C_ZERO, prod, mxdorb)

      psi_ao(:,:) = psi_ao2(:,:)


      do iband=1, neig
      jorb=1
      do iorb =1, norb
!        if (infolcao(4,iorb) < 2) then
          write(101,'(3i5,2f14.8)') iband, jorb, irk, prod(iorb,iband)
!          write(101,'(3i5,2f14.8)') iband, jorb, irk, U(iorb,iband)
          jorb=jorb+1
!        endif
      enddo
      enddo

      do iband=1, neig
        write(102,'(2i5,f14.8)') iband, irk, ev_interp(iband)*27.21138505
      enddo

    endif  ! SpinOrbit selection


!   start loop over neighbours

    do i = 1,nntot

      irk2 = nnlist(irk, i)

      rkpt2(1)  = kpt_lattice(1,irk2)
      rkpt2(2)  = kpt_lattice(2,irk2)
      rkpt2(3)  = kpt_lattice(3,irk2)

      call NonOrthoInterpRun(mtb,rkpt2,ev_interp)

      rkpt2(1)  = kpt_lattice(1,irk2)  + nncell(1,irk,i)
      rkpt2(2)  = kpt_lattice(2,irk2)  + nncell(2,irk,i)
      rkpt2(3)  = kpt_lattice(3,irk2)  + nncell(3,irk,i)

      lkplusg = .FALSE.

      call hamilt_struct(emax, rkpt2,                                    &
          mtxd2, isort2, qmod, ekpg,  lkplusg,                           &
          ng, kgv, adot,                                                 &
          mxdgve, mxddim)

      lnewanl = .TRUE.

      call atomic_orbital_c16(rkpt2, mtxd2, isort2, icmplx,              &
          mxdorb, ao_basis, infolcao,                                    &
          ng, kgv,                                                       &
          norbat, nqwf, delqwf, wvfao, lorb,                             &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdorb, mxdgve, mxdlao)

      if(mtb%lso==1) then

        do iorb=1,mxdorb
        do j=1,mtxd
          ao_basis_so(2*j-1,2*iorb  ) = ao_basis(j,iorb)
          ao_basis_so(2*j  ,2*iorb  ) = C_ZERO
          ao_basis_so(2*j-1,2*iorb-1) = C_ZERO
          ao_basis_so(2*j  ,2*iorb-1) = ao_basis(j,iorb)
        enddo
        enddo

        call zgemm('N', 'N', 2*mtxd2, norb, norb, C_UM, ao_basis_so,     &
            2*mxddim, mtb%wrk, 2*mxdorb, C_ZERO, ao_basis_ortho, 2*mxddim)

        call ao_SortToGammaSO(ao_basis_ortho, isort2,  mtxd2,  neig,  psi_ao2,  2*mxddim)

        call zgemm('C', 'N', norb, norb, 2*mtxd2, C_UM, psi_ao,          &
            2*mxddim, psi_ao2, 2*mxddim, C_ZERO, mmnk, 2*mxdorb)

      else

        call zgemm('N', 'N', mtxd2, norb, norb, C_UM, ao_basis, mxddim,  &
             mtb%wrk, mxdorb, C_ZERO, ao_basis_ortho, mxddim)


        call ao_SortToGamma(ao_basis_ortho, isort2, mtxd2, neig, psi_ao2, mxddim)

        call zgemm('C', 'N', norb, norb, mtxd2, C_UM, psi_ao,            &
            mxddim, psi_ao2, mxddim, C_ZERO, mmnk, mxdorb)

      endif


!dbg  write(*,'("neighbour", 6f8.5, 2i5)') rkpt, rkpt2, mtxd, norb

    write(100,'(5i5)') irk, irk2, nncell(1,irk,i), nncell(2,irk,i) ,nncell(3,irk,i)

     do jband=1, neig
     do iband=1, neig
       write(100,'(2f22.12)') mmnk(iband,jband)
     enddo
     enddo

!dbg write(99) mmnk


    enddo ! end loop over neighbours

  enddo  ! end loop over kpoints

! close files

  close(unit=101)
  close(unit=102)


!dbg call ao_interpolation_write_chk(mtb%fiData, neig)

  iotape = 15
  nstyle = 2


  deallocate(e_of_k)
  deallocate(e_of_k_so)
  deallocate(ev_interp)

  deallocate(ao_basis_ortho)
  deallocate(psi_ao2)

  deallocate(prod)

  deallocate(S)
  deallocate(S12)
  deallocate(S12_inv)
  deallocate(Swrk)
  deallocate(ev_wrk)

  deallocate(mmnk)

  deallocate(infolcao)



  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)

  deallocate(ao_basis)

  if(mtb%lso==1) then
    deallocate(psi_ao)
    deallocate(ao_basis_so)
  else
    deallocate(psi_ao)
  endif

  write(6,*)
  write(6,*) 'W90 interface calculations Done.'
  write(6,*)

  return

end subroutine ao_interpolation_w90



! THESE TWO SUBROUTINES SHOULD BE INLINED


subroutine ao_SortToGamma(psi_in,isort_in, mtxd_in, neig, psi_out, mxddim)
  implicit none

  integer, parameter                 ::  REAL64 = selected_real_kind(12)

  integer      :: mtxd_in
  integer      :: neig
  integer      :: mxddim
  integer      :: isort_in(mxddim)
  complex(REAL64) :: psi_in(mxddim,neig)
  complex(REAL64) :: psi_out(mxddim,neig)


  integer ig

  psi_out(:,:) = cmplx(0.0_REAL64,0.0_REAL64,REAL64)

  do ig=1, mtxd_in
    psi_out(isort_in(ig),:) = psi_in(ig,:)
  enddo

end subroutine ao_SortToGamma


subroutine ao_SortToGammaSo(psi_in,isort_in, mtxd_in, neig, psi_out, mxddim)
  implicit none

  integer, parameter                 ::  REAL64 = selected_real_kind(12)

  integer      :: mtxd_in
  integer      :: neig
  integer      :: mxddim
  integer      :: isort_in(mxddim)
  complex(REAL64) :: psi_in(mxddim,neig)
  complex(REAL64) :: psi_out(mxddim,neig)


  integer ig

  psi_out(:,:) = cmplx(0.0_REAL64,0.0_REAL64,REAL64)

  do ig=1, mtxd_in
    psi_out(2*(isort_in(ig))-1,:) = psi_in(2*ig-1,:)
    psi_out(2*(isort_in(ig)),:) = psi_in(2*ig,:)
  enddo

end subroutine ao_SortToGammaSo



