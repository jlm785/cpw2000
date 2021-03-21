!>     Define main variables and types for main cpw program

       module cpw_variables

!      copyright Jose Luis Martins/inesc-mn
!      May be used for research purposes only


!      version 4.95    November   2019

       implicit none
       
       integer, parameter          :: REAL64 = selected_real_kind(12)


!      dimensions

       type  ::  dims_t

         integer                            ::  mxdtyp                   !<  array dimension of types of atoms
         integer                            ::  mxdatm                   !<  array dimension of number of atoms of a given type

         integer                            ::  mxdgve                   !<  array dimension for g-space vectors
         integer                            ::  mxdnst                   !<  array dimension for g-space stars
         integer                            ::  mxdcub                   !<  array dimension for 3-index g-space

         integer                            ::  mxdlqp                   !<  array dimension for local potential
         integer                            ::  mxdlao                   !<  array dimension of orbital per atom type

         integer                            ::  mxdnrk                   !<  array dimension for number of k-points    
        
         integer                            ::  mxdbnd                   !<  array dimension for the number of bands
         integer                            ::  mxddim                   !<  array dimension for the hamiltonian
         integer                            ::  mxdsml                   !<  array dimension for the small hamiltonian

       end type dims_t


!      computational flags

       type  ::  flags_t

         character(len=4)                   ::  flgdal                   !<  controls use of dual approximation
         character(len=6)                   ::  flgscf                   !<  type of self consistent field and diagonalization
         character(len=6)                   ::  flgpsd                   !<  type of pseudopotential
         character(len=6)                   ::  flgcal                   !<  type of calculation
         character(len=6)                   ::  flgkeat                  !<  adds keating force field correction to energy
         character(len=6)                   ::  flgmix                   !<  choice of potential mixing

       end type flags_t


!      crystal structure

       type  ::  crys_t

         real(REAL64)                       ::  alatt                    !<  lattice constant

         real(REAL64)                       ::  adot(3,3)                !<  metric in direct space in atomic units (Bohr radius)

         integer                            ::  ntype                    !<  number of types of atoms
         integer, allocatable               ::  natom(:)                 !<  number of atoms of type i
         real(REAL64), allocatable          ::  rat(:,:,:)               !<  lattice coordinates of atom j of type i
         real(REAL64), allocatable          ::  atmass(:)                !<  atomic mass of atoms of type i
         character(len=2), allocatable      ::  nameat(:)                !<  chemical symbol for the type i

       end type crys_t


!      exchange and correlation choice

       type  ::  xc_t

         character(len=4)                   ::  author                   !<  type of xc wanted (ca=pz , pw92 , pbe,...)
         real(REAL64)                       ::  tblaha                   !<  Tran-Blaha constant

       end type xc_t 

!      plane-wave expansion

       type  ::  pwexp_t

         real(REAL64)                       ::  emax                     !<  kinetic energy cutoff of plane wave expansion (Hartree).
         real(REAL64)                       ::  teleck                   !<  Temperature of electrons (kelvin)

         integer                            ::  nbandin                  !<  target for number of bands      
       
         logical                            ::  lkplusg                  !<  If true use the previous G-vectors (same mtxd and isort) when epskplusg > error
         real(REAL64)                       ::  epskplusg                !<  criteria for switching to fixed k+G

       end type pwexp_t


!      accuracy parameters

       type  ::  acc_t

         integer                            ::  itmax                    !<  maximum number of self-consistent iterations
         real(REAL64)                       ::  epscv                    !<  convergence criterium for self-consistency
         real(REAL64)                       ::  epscvao                  !<  convergence criterium for self-consistency in atomic orbitals
         real(REAL64)                       ::  epspsi                   !<  convergence criterium for diagonalization

       end type acc_t 


!      space group information

       type  ::  spaceg_t  

         integer                            ::  ntrans                   !<  number of symmetry operations in the factor group
         integer                            ::  mtrx(3,3,48)             !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation
         real(REAL64)                       ::  tnp(3,48)                !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector of the k-th symmetry operation

       end type spaceg_t


!      reciprocal space

       type  ::  recip_t  

!         real(REAL64)                       ::  emax                     !<  kinetic energy cutoff of plane wave expansion (Hartree).
         integer                            ::  ng                       !<  total number of G-vectors with length less than gmax
         integer, allocatable               ::  kgv(:,:)                 !<  i-th component (reciprocal lattice coordinates) of the n-th G-vector ordered by stars of increasing length
         complex(REAL64), allocatable       ::  phase(:)                 !<  phase factor of G-vector n
         real(REAL64), allocatable          ::  conj(:)                  !<  is -1 if one must take the complex conjugate of x*phase
         integer                            ::  ns                       !<  number os stars with length less than gmax
         integer, allocatable               ::  inds(:)                  !<  star to which G-vector n belongs
         integer                            ::  kmax(3)                  !<  max value of |kgv(i,n)|
         integer, allocatable               ::  indv(:)                  !<  kgv(i,indv(jadd)) is the G-vector associated with jadd. jadd is defined by the g-vector components and kmax
         real(REAL64), allocatable          ::  ek(:)                    !<  kinetic energy (hartree) of g-vectors in star j
         integer, allocatable               ::  mstar(:)                 !<  number of G-vectors in the j-th star
         integer, allocatable               ::  izstar(:)                !<  is 0 if the phase=0    NOT USED

       end type recip_t


!      pseudo-potential

       type  ::  pseudo_t  

!        pseudo-potential (Kleinman-Bylander)

         real(REAL64)                       ::  ealraw                   !<  G=0 contrib. to the total energy. (non norm. to vcell, Hartree)
         real(REAL64), allocatable          ::  vloc(:,:)                !<  local pseudopotential for atom k (Hartree)
         real(REAL64), allocatable          ::  dcor(:,:)                !<  core charge density for atom k
         real(REAL64), allocatable          ::  dval(:,:)                !<  valence charge density for atom k

         real(REAL64), allocatable          ::  zv(:)                    !<  Valence of atom with type i

         real(REAL64)                       ::  ztot                     !<  total charge density (electrons/cell)
         integer, allocatable               ::  nq(:)                    !<  number of points for pseudo interpolation for atom k
         real(REAL64), allocatable          ::  delq(:)                  !<  step used in the interpolation for atom k
         integer, allocatable               ::  nkb(:,:,:)               !<  KB pseudo.  normalization for atom k, ang. mom. l-1

!        local pseudopotential in G-space

         real(REAL64), allocatable          ::  vql(:,:)                 !<  local pseudopotential for atom type i and prototype g-vector in star j
         real(REAL64), allocatable          ::  dnc(:,:)                 !<  core charge for atom type i and prototype g-vector in star j
         complex(REAL64), allocatable       ::  dvql(:)                  !<  derivative of the local pseudopotential for the prototype g-vector in star j
         complex(REAL64), allocatable       ::  ddc(:)                   !<  derivative of the core charge for the prototype g-vector in star j

!        KB non-local pseudo-potential in g-space

         real(REAL64), allocatable          ::  vkb(:,:,:,:)             !<  KB nonlocal pseudo. for atom k, ang. mom. l-1 (not normalized to vcell, Hartree)

       end type pseudo_t


!      structure factors

       type  ::  strfac_t

         integer                            ::  icmplx                   !<  indicates if the structure factor is complex
         complex(REAL64), allocatable       ::  sfact(:,:)               !<  structure factor for atom k and star i

       end type strfac_t
     

!      charge densities

       type  ::  chdens_t

         complex(REAL64), allocatable       ::  den(:)                   !<  total charge density for the prototype G-vector in star j
         complex(REAL64), allocatable       ::  denc(:)                  !<  core charge density for the prototype G-vector in star j
         complex(REAL64), allocatable       ::  dens(:)                  !<  spherical atomic valence charge density for the prototype G-vector in star j
         complex(REAL64), allocatable       ::  dend(:)                  !<  bonding charge density from previous md step
         complex(REAL64), allocatable       ::  dend1(:)                 !<  bonding charge density from second previous md step

       end type chdens_t


!      local potential contributions

       type  ::  vcomp_t

         complex(REAL64), allocatable       ::  vion(:)                  !<  ionic potential for the prototype G-vector in star j
         complex(REAL64), allocatable       ::  vhar(:)                  !<  Hartree potential for the prototype G-vector in star j
         complex(REAL64), allocatable       ::  vxc(:)                   !<  Hartree+exchange+correlation potential for the prototype G-vector in star j
         complex(REAL64), allocatable       ::  veff(:)                  !<  Effective potential for the prototype G-vector in star j

       end type vcomp_t


!      atomic orbitals in G-space

       type  ::  atorb_t

         logical                            ::  latorb                   !<  indicates if all atoms have information about atomic orbitals

         integer, allocatable               ::  norbat(:)                !<  number of atomic orbitals for atom k
         integer, allocatable               ::  lorb(:,:)                !<  angular momentum of orbital n of atom k
         real(REAL64), allocatable          ::  wvfao(:,:,:)             !<  wavefunction for atom k, ang. mom. l 

         integer, allocatable               ::  nqwf(:)                  !<  number of points for wavefunction interpolation for atom k
         real(REAL64), allocatable          ::  delqwf(:)                !<  step used in the wavefunction interpolation for atom k

       end type atorb_t


!      k-point data

       type  ::  kpoint_t

         integer                            ::  nx, ny, nz               !<  size of the integration mesh in k-space (nx*ny*nz)
         real(REAL64)                       ::  sx, sy, sz               !<  offset of the integration mesh (usually 0.5)

         integer                            ::  nrk                      !<  number of k-points for integration in the irreducible wedge of the brillouin zone
         real(REAL64), allocatable          ::  rk(:,:)                  !<  j-th component (reciprocal lattice coordinates) of the i-th k-point in the integration mesh
         real(REAL64), allocatable          ::  wgk(:)                   !<  weight in the integration of k-point i (sum_i w(i) = 1)
         integer, allocatable               ::  nband(:)                 !<  number of bands considered for k-point i
         integer, allocatable               ::  indk(:,:)                !<  index of the six k-points neighbouring k-point i
         integer, allocatable               ::  kmap(:,:,:,:)            !<  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation

       end type kpoint_t 


!      hamiltonian size and indexation for all k-points

       type  ::  hamallk_t

         integer, allocatable               ::  mtxd_allk(:)             !<  dimension of the hamiltonian for k-point n
         integer, allocatable               ::  isort_allk(:,:)          !<  G-vector associated with k+G vector i of hamiltonian for k-point n

       end type hamallk_t


!      psi for all k-points

       type  ::  psiallk_t

         real(REAL64), allocatable          ::  eig_allk(:)              !<  eigenvalue j, for all the k-points
         real(REAL64), allocatable          ::  occ_allk(:)              !<  fractional ocupation of level j, for all the k-points
         complex(REAL64), allocatable       ::  psi_allk(:,:,:)          !<  eigenvectors for all k-points

       end type psiallk_t


!      molecular dynamics variables

       type  ::  moldyn_t

         integer                            ::  nstep                    !<  number of molecular dynamics or maximum minimization steps
         integer                            ::  iseed                    !<  seed for random number generator
         real(REAL64)                       ::  tstep                    !<  timestep in atomic units (2.4x10^17 s)
         real(REAL64)                       ::  tempinik                 !<  initial temperature in k (Maxwell distribution)
         real(REAL64)                       ::  ekin                     !<  kinetic energy of the nuclei (atomic units)
         real(REAL64), allocatable          ::  vat(:,:,:)               !<  d rat / d t  velocity in lattice coordinates of atom j of type i
         real(REAL64), allocatable          ::  rat1(:,:,:)              !<  previous value of lattice coordinates of atom j of type i
         real(REAL64), allocatable          ::  frc1(:,:,:)              !<  previous value of force on atom j of type i
         real(REAL64)                       ::  beta                     !<  friction coefficient/mass (in a.u.)
         real(REAL64)                       ::  tempk                    !<  ionic temperature (in Kelvin)

         real(REAL64)                       ::  pgtol                    !<  criteria for force and stress convergence
         real(REAL64)                       ::  dxmax                    !<  maximum step size in forces and stresses

       end type moldyn_t 


!      variational cell shape variables

       type  ::  vcsdyn_t

         real(REAL64)                       ::  celmas                   !<  fictitious cell mass (atomic units)
         real(REAL64)                       ::  ekcell                   !<  fictitious cell kinetic energy
         real(REAL64)                       ::  vadot(3,3)               !<  d adot / d t  rate of change of metric
         real(REAL64)                       ::  adot1(3,3)               !<  previous metric
         real(REAL64)                       ::  frcel1(3,3)              !<  previous fictitious cell force

         real(REAL64)                       ::  press                    !<  external pressure (Hartree/Bohr**3)
         real(REAL64)                       ::  strext(3,3)              !<  external stress tensor in covariant lattice coordinates (Hartree/Bohr^3)

       end type vcsdyn_t 


!      energy force stress

       type  ::  enfrst_t

         real(REAL64)                       ::  energy                   !<  Eenergy (Hartree)
         real(REAL64), allocatable          ::  force(:,:,:)             !<  d energy / d rat,  force on the n-th atom of type i (contravarian components, hartree/bohr)
         real(REAL64)                       ::  stress(3,3)              !<  d energy / d adot,  stress tensor (contravariant)

       end type enfrst_t


       end module cpw_variables
