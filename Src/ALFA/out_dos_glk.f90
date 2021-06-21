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

!>  Calculates the bands on an uniform grid
!>  for later processing using the GLK interpolation

  subroutine out_dos_glk(diag_type, lworkers, xsvd, csvd,                &
  title, subtitle,                                                       &
  emax, flgdal, flgpsd,                                                  &
  epspsi, icmax, ztot,                                                   &
  adot, ntype, natom, rat, ntrans, mtrx,                                 &
  ng, kgv, phase, conj,                                                  &
  ns, inds, kmax, indv, ek,                                              &
  sfact, icmplx,                                                         &
  veff,                                                                  &
  nqnl, delqnl, vkb, nkb,                                                &
  latorb, norbat, nqwf, delqwf, wvfao, lorb,                             &
  mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)
     

! version 4.42. 8 may 2004. jlm
! modified 11 february 2008 (read from file). JLM
! modified December 18, 2013 (f90, new interface). JLM
! modified, dimensions vkb, May 23, 2014. JLM
! workers and restart capabilities for version 4.93. CLR
! adpated for version 4.96 March 2020. CLR
! modified, latorb, 7 June 2020. JLM
! modified, diag_type, icmax, 8-14 June 2020. JLM
! modified, nc(3), cube2tetra, x/csvd, 15 September 2020. JLM
! Modified, vmax, vmin, 27 November 2020. JLM

! copyright  Jose Luis Martins/Carlos Loia Reis/INESC-MN


! The reference states are calculated in an inefficient way....

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  character(len=4), intent(in)       ::  diag_type                       !<  selects diagonalization, 'pw  ','ao  ','aojc'
  logical, intent(in)                ::  lworkers                        !<  use lworkers
  real(REAL64), intent(in)           ::  xsvd                            !<  Ignore states with SVD singular values smaller than xsvd.  0 < xsvd < 1.
  real(REAL64), intent(in)           ::  csvd                            !<  Use csvd*neig in SVD procedure. 1 < csvd < 2.

  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4)                   ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  integer, intent(in)                ::  icmax                           !<  maximum number of iterations for diagonalization
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase

  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(in)                ::  indv(mxdcub)                    !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j

  complex(REAL64), intent(in)        ::  sfact(mxdtyp,mxdnst)            !<  structure factor
  integer, intent(in)                ::  icmplx                          !<  indicates if the structure factor is complex

  complex(REAL64), intent(in)        ::  veff(mxdnst)                    !<  ionic potential (local+Hartree+XC) for the prototype g-vector in star j

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<   KB pseudo.  normalization for atom k, ang. mom. l

  logical, intent(in)                ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  wavefunction for atom k, ang. mom. l (normalized to vcell)
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k

! local allocatable arrays for integration points

  integer, allocatable               ::  kmap(:,:,:,:)                   !  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
  integer, allocatable               ::  nband(:)                        !  number of bands for each k-points
  integer, allocatable               ::  indk(:,:)                       !  index of the six k-points neighbouring k-point i
  real(REAL64),allocatable           ::  rk(:,:)                         !  component in lattice coordinates of the k-point in the mesh
  real(REAL64),allocatable           ::  wgk(:)                          !  weight in the integration of k-point

! allocatable arrays for eigensolutions

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !  component j of eigenvector i (guess on input)
  complex(REAL64), allocatable       ::  hpsi(:,:)                       !  H | psi> 
  real(REAL64), allocatable          ::  ekpsi(:)                        !  kinetic energy of eigenvector i. (hartree)
  real(REAL64), allocatable          ::  ekpsi_so(:)                     !  kinetic energy of eigenvector i. (hartree)
  real(REAL64), allocatable          ::  vscr(:)                         !  screened potential in the fft real space mesh

  real(REAL64), allocatable          ::  ei_so(:)                        !  spin-orbit eigenvalue (hartree)
  complex(REAL64), allocatable       ::  psi_so(:,:)                     !  component j of eigenvector i (guess on input)

  real(REAL64), allocatable          ::  e_of_k(:,:)                     !  band energies of k-point in plot
  real(REAL64), allocatable          ::  e_of_k_so(:,:)                  !  spin-orbit band energies of k-point in plot

  complex(REAL64), allocatable       ::  psi_all(:,:,:)                  !  psi for the 4 reference points
  integer,allocatable                ::  isort_all(:,:)                  !  isort for the 4 reference points

  real(REAL64), allocatable          ::  rk_ref(:,:)

! local variables

  integer                            ::  mxdscr                          !  array dimension for screening potential
  
  logical                            ::  lpsiso                          !  If true calculates the spin-orbit perturbed wave-functions 
  
  integer                            ::  mxddim                          !  array dimension for the hamiltonian
  integer                            ::  mxdbnd                          !  array dimension for the number of bands
  integer                            ::  mxdwrk                          !  array dimension for fft transform workspace
  integer                            ::  mxdnrk                          !  size of k-points

  integer                            ::  mtxd                            !  dimension of the hamiltonian
  integer                            ::  neig                            !  number of eigenvectors required (maybe modified on output)
  real(REAL64)                       ::  rkpt(3)                         !  j-th component in lattice coordinates of the k-point
  integer                            ::  kmscr(7)                        !  max value of kgv(i,n) used for the potential fft mesh
  integer                            ::  idshift                         !  shift of the fft mesh, used /= 0 only in highly banked memory.

  real(REAL64)                       ::  vmax, vmin                      !  maximum and minimum values of vscr

  integer                            ::  irk,nrka
  character(len=5)                   ::  labelk
  integer                            ::  ipr,nrk,nd
  logical                            ::  lfile
  integer                            ::  nbandi, nx,ny,nz
  real(REAL64)                       ::  sx,sy,sz
  integer                            ::  nsfft(3)

  integer                            ::  io, io11, io66, io67            !  tape numbers
  integer                            ::  ioerr

  integer                            ::  identif                         !  identifier

  character(len=12)                  ::  filemesh
  character(len=14)                  ::  file_wf
  character(len=15)                  ::  file_eig
  character(len=12)                  ::  filedos

! workers and restart  
  
  integer                            ::  nc(3)                           !  size of reference mesh

  integer                            ::  irec_err
  integer                            ::  irk_rd
  integer                            ::  iworker, nworker, cworker
  integer                            ::  ir_size, pp_flag
  logical                            ::  lmyjob

! reference states

  integer                            ::  ic1,ic2,ic3
  integer, allocatable               ::  irk_from_ic(:,:,:)
  integer, allocatable               ::  ic_from_irk(:,:)
  integer                            ::  m1, m2, m3

  integer                            ::  nrk3

  real(REAL64)                       ::  zk(4)
  real(REAL64)                       ::  y(3)
  integer                            ::  iq(4,3)
  
  integer                            ::  ix1, iy1, iz1 
  integer                            ::  ix2, iy2, iz2
  integer                            ::  ix3, iy3, iz3
  integer                            ::  ix4, iy4, iz4
  
  real(REAL64)                       ::  rkpt_all(3,4)
  real(REAL64)                       ::  rki

  integer                            ::  mtxd_all(4)                     !  mtxd for the 4 reference points
  integer                            ::  neig_all(4)                     !  neig for the 4 reference points      
  
  logical                            ::  lex

  real(REAL64)                       ::  t0, t1 , t2                     !  timings

! parameters

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer                            ::  i


  call zeelap(t0)

  io11 = 11
  io66 = 66
  io67 = 67
  io = 12

  filemesh = 'DOS_MESH.DAT'
  file_wf = 'wf_dos_rec.dat'
  file_eig = 'eig_dos_rec.dat'
  filedos = 'dos_file.dat'
  
! calculates local potential in fft mesh


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
  
  call pot_local(ipr, vscr, vmax, vmin, veff, kmscr, idshift,            &
  ng, kgv, phase, conj, ns, inds,                                        &
  mxdscr, mxdgve, mxdnst)
  
  ipr = 2

  lfile = .false.
  open(unit=io11, file=filemesh, status='old', IOSTAT=ioerr, form = 'formatted')
  if(ioerr == 0) lfile = .true.

  if(lfile) then
    read(io11,*) nbandi, nx,ny,nz, sx,sy,sz
    read(io11,*,IOSTAT=ioerr) nc(1),nc(2),nc(3)
    if(ioerr /= 0) then
      write(6,*)
      write(6,*) '    Missing second line in DOS_MESH.DAT'
      write(6,*) '    Using default nc = 2'
      write(6,*)
      nc(1) = 2
      nc(2) = 2
      nc(3) = 2
    endif
             
    close(unit=io11)       
    mxdbnd = nbandi
  else       
    mxdbnd = nint(ztot) + 4
    nbandi = mxdbnd
    nx = 8
    ny = 8
    nz = 8
    sx = ZERO
    sy = ZERO
    sz = ZERO
    nc(1) = 2
    nc(2) = 2
    nc(3) = 2
  endif
  
  write(6,*)
  write(6,*) "  nc subdivision is", nc(1),nc(2),nc(3)
  write(6,*)
  
  call size_mxdnrk(nx,ny,nz, sx,sy,sz, adot, ntrans, mtrx,               &
  mxdnrk)

  allocate(kmap(3,nx,ny,nz))
  allocate(nband(mxdnrk))
  allocate(indk(6,mxdnrk))
  allocate(rk(3,mxdnrk))
  allocate(wgk(mxdnrk))
  
  call int_pnt(nbandi, nx,ny,nz, sx,sy,sz, ipr,                          &
  adot,                                                                  &
  ntrans, mtrx,                                                          &
  nrk, rk, wgk, nband, indk, kmap,                                       &
  mxdnrk, mxdbnd)

  neig = nbandi
  
! finds mxddim

  if(lworkers) then
    write(6,*) "Using Workers and Restart Capablities"       

    write(6,*) 'how many workers ?'
    read(5,*) nworker
  
    write(6,*) 'current worker ?'
    read(5,*) iworker
  else
    nworker = 1
    iworker = 1
  endif

  mxddim = 1
  do irk = 1,nrk
    cworker = mod(irk-1,nworker) + 1
    if(cworker== iworker) then
      rkpt(1) = rk(1,irk)
      rkpt(2) = rk(2,irk)
      rkpt(3) = rk(3,irk)
      call size_mtxd(emax,rkpt,adot,ng,kgv,nd)         
      if(nd > mxddim) mxddim = nd
    endif
  enddo

  do ic1 = 0,nc(1)
  do ic2 = 0,nc(2)
  do ic3 = 0,nc(3)
    rkpt(1) = ic1*UM/(nc(1)*UM)
    rkpt(2) = ic2*UM/(nc(2)*UM)
    rkpt(3) = ic3*UM/(nc(3)*UM)
    call size_mtxd(emax,rkpt,adot,ng,kgv,nd)         
    if(nd > mxddim) mxddim = nd
  enddo
  enddo
  enddo

! mxddim = int(1.01*mxddim)

! allocates arrays

  allocate(ei(mxdbnd))
  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdbnd))
  allocate(hpsi(mxddim,mxdbnd))
  allocate(ekpsi(mxdbnd))
  allocate(ekpsi_so(2*mxdbnd))
  
  lpsiso = .FALSE.

  allocate(psi_so(2*mxddim,2*mxdbnd))
  allocate(ei_so(2*mxdbnd))

  nrk3 = (nc(1)+1)*(nc(2)+1)*(nc(3)+1)
  
  allocate(irk_from_ic(0:nc(1),0:nc(2),0:nc(3)))
  allocate(ic_from_irk(3,nrk3))
  
  allocate(isort_all(mxddim,4))      
  allocate(psi_all(mxddim,mxdbnd,4))

  
! open files

  irk=1
  inquire(iolength = ir_size) irk_rd, ei(:), ei_so(:)

! if run by a human do not restart

  if(.not. lworkers) then
    inquire(unit = io67, exist = lex)
    if(lex) close(unit = io67)
    call execute_command_line("rm " // file_eig // " 2> /dev/null ")
  endif


  open(unit = io67, file = file_eig, access="direct", recl=ir_size)

  inquire(iolength = ir_size) irk, mtxd, rkpt, psi(:,:), isort(:)

! if run by a human do not restart

  if(.not. lworkers) then
    inquire(unit = io66, exist = lex)
    if(lex) close(unit = io66)
    call execute_command_line("rm " // file_wf // " 2> /dev/null ")
  endif

  open(unit = io66, file = file_wf, access="direct", recl=ir_size)

  pp_flag = 0


! Preparatory calculation.

  allocate(rk_ref(3,nrk3))

  irk = 1
  do ic1 = 0,nc(1)
  do ic2 = 0,nc(2)
  do ic3 = 0,nc(3)

    rk_ref(1,irk) = ic1*UM/(nc(1)*UM)
    rk_ref(2,irk) = ic2*UM/(nc(2)*UM)
    rk_ref(3,irk) = ic3*UM/(nc(3)*UM)

    irk_from_ic(ic1,ic2,ic3) = irk
    ic_from_irk(1,irk) = ic1
    ic_from_irk(2,irk) = ic2
    ic_from_irk(3,irk) = ic3
    irk=irk+1  

  enddo
  enddo
  enddo



  call interp_glk_prep(diag_type, io66,                                  &
    nrk3, rk_ref,                                                        &
    emax, neig, flgpsd,                                                  &
    epspsi, icmax,                                                       &
    adot, ntype, natom, rat,                                             &
    ng, kgv, phase, conj,                                                &
    ns, inds, kmax, indv, ek,                                            &
    sfact, icmplx,                                                       &
    veff,                                                                &
    nqnl, delqnl, vkb, nkb,                                              &
    vscr, kmscr,                                                         &
    latorb, norbat, nqwf, delqwf, wvfao, lorb,                           &
    mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxddim,              &
    mxdbnd, mxdscr, mxdlao)



  call zeelap(t1)

! loop over k-points


  neig_all(1) = neig
  neig_all(2) = neig
  neig_all(3) = neig
  neig_all(4) = neig

  do irk=1,nrk

    cworker = mod(irk-1,nworker) + 1
    lmyjob = .TRUE.
    if(cworker == iworker) then
      irk_rd = -10        
      read(io67,rec=irk, iostat=irec_err) irk_rd
      if(irk_rd == irk) then
        write(6,'("not computing k-point #",i5)') irk
        lmyjob = .FALSE.                                
      endif            
    else
      write(6,'("not computing k-point #",i5)') irk
      lmyjob = .FALSE.                   
    endif 

    if(lmyjob) then

      rkpt(1) = rk(1,irk)
      rkpt(2) = rk(2,irk)
      rkpt(3) = rk(3,irk)

      m1 = mod(floor(nc(1)*rkpt(1)), nc(1))
      m2 = mod(floor(nc(2)*rkpt(2)), nc(2))
      m3 = mod(floor(nc(3)*rkpt(3)), nc(3))

      y(1) = nc(1)*rkpt(1) - m1*UM
      y(2) = nc(2)*rkpt(2) - m2*UM
      y(3) = nc(3)*rkpt(3) - m3*UM

      call cube2tetra(y,zk,iq)

      ix1 = m1 + iq(1,1)
      iy1 = m2 + iq(1,2)
      iz1 = m3 + iq(1,3)

      ix2 = m1 + iq(2,1)
      iy2 = m2 + iq(2,2)
      iz2 = m3 + iq(2,3)

      ix3 = m1 + iq(3,1)  
      iy3 = m2 + iq(3,2)
      iz3 = m3 + iq(3,3)

      ix4 = m1 + iq(4,1) 
      iy4 = m2 + iq(4,2)
      iz4 = m3 + iq(4,3)

      read(io66,rec=irk_from_ic(ix1,iy1,iz1))  irk_rd, mtxd_all(1),      &
      rkpt_all(:,1), psi_all(:,:,1), isort_all(:,1)          
       
      read(io66,rec=irk_from_ic(ix2,iy2,iz2))  irk_rd, mtxd_all(2),      &
      rkpt_all(:,2), psi_all(:,:,2), isort_all(:,2)          
       
      read(io66,rec=irk_from_ic(ix3,iy3,iz3))  irk_rd, mtxd_all(3),      &
      rkpt_all(:,3), psi_all(:,:,3), isort_all(:,3)                
       
      read(io66,rec=irk_from_ic(ix4,iy4,iz4))  irk_rd,  mtxd_all(4),     &
      rkpt_all(:,4), psi_all(:,:,4), isort_all(:,4)          

!     paranoid check

      do i = 1,3
        rki = zk(1)*rkpt_all(i,1) + zk(2)*rkpt_all(i,2)                  &
            + zk(3)*rkpt_all(i,3) + zk(4)*rkpt_all(i,4)
        if(abs(rkpt(i) - rki) > 0.000001) then
          write(6,*)
          write(6,'("    stopped in out_dos_lk:  rk = ",2 f12.4)')       &
               rkpt(i), rki

          stop

        endif
      enddo


      call interpolation_glk(4, emax, neig, xsvd, csvd,                  &
        ZK, rkpt_all, mtxd_all, neig_all,                                &
        isort_all, psi_all,                                              &
        ei, psi, hpsi, mtxd, isort, qmod, ekpg,                          &
        ng, kgv,                                                         &
        ntype, natom, rat, adot,                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        vscr, kmscr,                                                     &
        mxdtyp, mxdatm, mxddim, mxdlqp, mxdbnd, mxdgve, mxdscr)


      ipr = 1
      nrka = -1
      call print_eig(ipr,irk,labelk,nrka,rkpt,                           &
      mtxd,icmplx,neig,psi,                                              &
      adot,ei,ekpsi,isort,kgv,                                           &
      mxddim,mxdbnd,mxdgve)


      call spin_orbit_perturb(rkpt,mtxd,isort,                           &
      neig,psi,ei,ei_so,psi_so,lpsiso,                                   &
      ng,kgv,                                                            &
      nqnl,delqnl,vkb,nkb,                                               &
      ntype,natom,rat,adot,                                              &
      mxdtyp,mxdatm,mxdlqp,mxddim,mxdbnd,mxdgve)

      if(ipr == 2) then
        call kinetic_energy_so(neig,mtxd,ekpg,psi_so,ekpsi_so,           &
        mxddim,mxdbnd)
      endif

      call print_eig_so(ipr,irk,labelk,nrka,rkpt,                        &
      mtxd,neig,psi_so,                                                  &
      adot,ei_so,ekpsi_so,isort,kgv,                                     &
      mxddim,mxdbnd,mxdgve)


      write(6,'( "iworker #",i5, "   writing in irk # "                  &
                & ,i5, "   of ", i5)') iworker, irk,nrk
      write(6,*)

      write(io67, rec = irk) irk, ei(:), ei_so(:)                        

    endif

  enddo

! end loop over k-points

  do irk = 1, nrk
    read(io67, rec = irk, iostat = irec_err) irk_rd
    if(irk_rd /= irk) then 
       pp_flag = 0

       exit

    endif            
    pp_flag = 1            
  enddo

  call zeelap(t2)

  write(6,*)
  write(6,'("  Elapsed time for reference states: ",f12.3)') t1-t0
  write(6,'("  Elapsed time for DOS calculation:  ",f12.3)') t2-t1
  write(6,*)

  deallocate(vscr)

  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)
  deallocate(psi)
  deallocate(hpsi)
  deallocate(ekpsi)
  deallocate(ekpsi_so)

  deallocate(psi_so)
  
  deallocate(irk_from_ic)
  deallocate(ic_from_irk)
  deallocate(isort_all)      
  deallocate(psi_all)

  allocate(e_of_k(neig,nrk))
  allocate(e_of_k_so(2*neig,nrk))

  if (pp_flag ==1) then
    do irk = 1,nrk

      read(io67,rec=irk, iostat=irec_err) irk_rd, ei(:), ei_so(:)

      e_of_k(1:neig,irk) = ei(1:neig)
      e_of_k_so(1:2*neig,irk) = ei_so(1:2*neig)

    enddo
  endif

  deallocate(ei)
  deallocate(ei_so)

  if(pp_flag == 1) then        

    write(6,*)
    write(6,*) '  Generating DOS files!'
    write(6,*)

    identif = 0

    call out_dos_write(filedos, io, title, subtitle,                     &
      .TRUE. , .TRUE. ,identif,                                          &
      nrk, nx, ny, nz, ztot, adot, ntrans, mtrx,                         &
      nband, rk, wgk, indk, kmap, e_of_k, e_of_k_so,                     &
      mxdnrk, mxdbnd)

  else

    close(io66)
    close(io67)
    write(6,*) "this worker is done. run again when all workers are done to see results."          

  endif
 
! if run by a human clean up temporary files

  if(.not. lworkers) then

    call execute_command_line("rm " // file_wf // " 2> /dev/null ")
    call execute_command_line("rm " // file_eig // " 2> /dev/null ")

  endif

  deallocate(kmap)
  deallocate(nband)
  deallocate(indk)
  deallocate(rk)
  deallocate(wgk)

  deallocate(e_of_k)
  deallocate(e_of_k_so)

  return
  end subroutine out_dos_glk
