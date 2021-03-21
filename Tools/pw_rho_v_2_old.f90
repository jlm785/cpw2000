       program pw_rho_v_2_old

!      rewrites the PW_RHO_V file for older versions (old versions do not read new files)


!      Written May 25, 2020 from previous code. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.98

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer                            ::  mxdtyp                     !<  array dimension of types of atoms
       integer                            ::  mxdatm                     !<  array dimension of types of atoms
       integer                            ::  mxdgve                     !<  array dimension for g-space vectors
       integer                            ::  mxdnst                     !<  array dimension for g-space stars
       integer                            ::  mxdlqp                     !<  array dimension for local potential
       integer                            ::  mxdlao                     !<  array dimension of orbital per atom type

       integer                            ::  io                         !<  number of tape to which the pseudo is added.

       integer                            ::  ipr                        !<  printing level

!      output

       character(len=60)                  ::  pwline                     !<  identification of the calculation
       character(len=50)                  ::  title                      !<  title for plots
       character(len=140)                 ::  subtitle                   !<  subtitle for plots
       character(len=250)                 ::  meta_cpw2000               !<  metadata from cpw2000

       character(len=3)                   ::  author                     !<  type of xc wanted (CA=PZ , PW92 , PBE)
       character(len=6)                   ::  flgscf                     !<  type of self consistent field and diagonalization
       character(len=4)                   ::  flgdal                     !<  whether the dual approximation is used
       real(REAL64)                       ::  emax                       !<  kinetic energy cutoff of plane wave expansion (Hartree).
       real(REAL64)                       ::  teleck                     !<  electronic temperature (in Kelvin)

       integer                            ::  nband                      !<  target for number of bands
       integer                            ::  nx,ny,nz                   !<  divisions of Brillouin zone for integration (Monkhorst-Pack)
       real(REAL64)                       ::  sx,sy,sz                   !<  shift of points in division of Brillouin zone for integration (Monkhorst-Pack)
       real(REAL64)                       ::  alatt                      !<  lattice constant

       real(REAL64)                       ::  adot(3,3)                  !<  metric in direct space
       integer                            ::  ntype                      !<  number of types of atoms
       integer,allocatable                ::  natom(:)                   !  number of atoms of type i
       character(len=2),allocatable       ::  nameat(:)                  !  chemical symbol for the type i
       real(REAL64),allocatable           ::  rat(:,:,:)                 !  k-th component (in lattice coordinates) of the position of the n-th atom of type i

       integer                            ::  ng                         !<  total number of g-vectors with length less than gmax
       integer                            ::  kmax(3)                    !<  max value of kgv(i,n)
       integer, allocatable               ::  kgv(:,:)                   !  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
       complex(REAL64), allocatable       ::  phase(:)                   !  phase factor of G-vector n
       real(REAL64), allocatable          ::  conj(:)                    !  is -1 if one must take the complex conjugate of x*phase
       integer, allocatable               ::  inds(:)                    !  star to which g-vector n belongs
       integer, allocatable               ::  indv(:)                    !  kgv(i,indv(jadd)) is the g-vector associated with jadd jadd is defined by the g-vector components and kmax
       integer                            ::  ns                         !<  number os stars with length less than gmax
       integer, allocatable               ::  mstar(:)                   !  number of g-vectors in the j-th star
       real(REAL64), allocatable          ::  ek(:)                      !  kinetic energy (Hartree) of g-vectors in star j

       complex(REAL64), allocatable       ::  veff(:)                    !<  effective potential (local+Hartree+Xc) for the prototype g-vector in star j
       complex(REAL64), allocatable       ::  den(:)                     !<  valence charge density for the prototype g-vector in star j
       complex(REAL64), allocatable       ::  dend(:)                    !<  bonding charge density for the prototype g-vector in star j

       character(len=3), allocatable      ::  irel(:)                    !<  type of calculation relativistic/spin
       character(len=4), allocatable      ::  icore(:)                   !<  type of partial core correction
       character(len=2), allocatable      ::  icorr(:)                   !<  type of correlation
       character(len=60), allocatable     ::  iray(:)                    !<  information about pseudopotential
       character(len=70), allocatable     ::  ititle(:)                  !<  further information about pseudopotential

       real(REAL64)                       ::  ealraw                     !<  G=0 contrib. to the total energy. (non norm. to vcell,hartree)
       integer, allocatable               ::  nq(:)                      !  number of points for pseudo interpolation for atom k
       real(REAL64), allocatable          ::  delq(:)                    !  step used in the pseudo interpolation for atom k
       real(REAL64), allocatable          ::  vkb(:,:,:,:)               !  kb nonlocal pseudo. for atom k, ang. mom. l. (not normalized to vcell, hartree)
       integer, allocatable               ::  nkb(:,:,:)                 !   kb pseudo.  normalization for atom k, ang. mom. l
       real(REAL64), allocatable          ::  vloc  (:,:)                !  local pseudopotential for atom k (hartree)
       real(REAL64), allocatable          ::  dcor(:,:)                  !  core charge density for atom k
       real(REAL64), allocatable          ::  dval (:,:)                 !  valence charge density for atom k
       integer, allocatable               ::  norbat(:)                  !  number of atomic orbitals for atom k
       integer, allocatable               ::  nqwf(:)                    !  number of points for wavefunction interpolation for atom k
       real(REAL64), allocatable          ::  delqwf(:)                  !  step used in the wavefunction interpolation for atom k
       integer, allocatable               ::  lorb(:,:)                  !  angular momentum of orbital n of atom k
       real(REAL64), allocatable          ::  wvfao(:,:,:)               !  wavefunction for atom k, ang. mom. l

       logical                            ::  latorb                     !<  indicates if all atoms have information about atomic orbitals
       real(REAL64), allocatable          ::  zv(:)                      !  valence of atom with type i
       real(REAL64)                       ::  ztot                       !<  total charge density (electrons/cell)

       character(len=1)                   ::  yesno
       character(len=140)                 ::  line140
       character(len=9)                  ::  bdate
       character(len=8)                   ::  btime

       integer                            ::  mxdl

       integer                            ::  iflag

!      counters

       integer    ::  i, j, k

       
       write(6,*)
       write(6,'("  program reads the existing PW_RHO_V.DAT file ")')
       write(6,'("  and writes an old style PW_RHO_V_OLD.DAT file ")')
       write(6,*)

       write(6,*)
       write(6,'("  The default old style is the one between versions 4.54 and 4.75 ")')
       write(6,*)
       write(6,'("  Do you want an even older style? (y/n) ")')
       write(6,*)
       
       read(5,*) yesno
       
       if(yesno == 'y' .or. yesno == 'Y') then
         write(6,*)
         write(6,'("  Using style before versions 4.54 ")')
         write(6,*)
         iflag = 0
       else
         write(6,*)
         write(6,'("  Using style between versions 4.54 and 4.75 ")')
         write(6,*)
         iflag = 1
       endif

       io = 52

       call pw_rho_v_in_size('PW_RHO_V.DAT', io,                               &
     &    mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdlao)


       allocate(natom(mxdtyp))
       allocate(nameat(mxdtyp))
       allocate(rat(3,mxdatm,mxdtyp))

       allocate(kgv(3,mxdgve))
       allocate(phase(mxdgve))
       allocate(conj(mxdgve))
       allocate(mstar(mxdnst))
       allocate(den(mxdnst))
       allocate(dend(mxdnst))
       allocate(veff(mxdnst))


       allocate(nq(mxdtyp))
       allocate(delq(mxdtyp))
       allocate(vkb(-2:mxdlqp,0:3,-1:1,mxdtyp))
       allocate(nkb(0:3,-1:1,mxdtyp))
       allocate(vloc(-1:mxdlqp,mxdtyp))
       allocate(dcor(-1:mxdlqp,mxdtyp))
       allocate(dval(-1:mxdlqp,mxdtyp))
       allocate(zv(mxdtyp))
       allocate(norbat(mxdtyp))
       allocate(lorb(mxdlao,mxdtyp))
       allocate(wvfao(-2:mxdlqp,mxdlao,mxdtyp))
       allocate(nqwf(mxdtyp))
       allocate(delqwf(mxdtyp))

       allocate(irel(mxdtyp))
       allocate(icore(mxdtyp))
       allocate(icorr(mxdtyp))
       allocate(iray(mxdtyp))
       allocate(ititle(mxdtyp))


       ipr = 1

       call pw_rho_v_in('PW_RHO_V.DAT', io, ipr,                         &
     &        pwline, title, subtitle, meta_cpw2000,                     &
     &        author, flgscf, flgdal, emax, teleck,                      &
     &        nx, ny, nz, sx, sy, sz, nband,                             &
     &        alatt, adot, ntype, natom, nameat, rat,                    &
     &        ng, kmax, kgv, phase, conj, ns, mstar,                     &
     &        veff, den, dend,                                           &
     &        irel, icore, icorr, iray, ititle,                          &
     &        ealraw, zv, ztot,                                          &
     &        nq, delq, vkb, nkb, vloc, dcor, dval,                      &
     &        norbat, nqwf, delqwf, wvfao, lorb, latorb,                 &
     &        mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdlao)


       call size_mxdlqp_lao(ntype,nameat,                                &
     &        mxdtyp,mxdl,mxdlao)

       io = 52
       
       open(unit=io,file='PW_RHO_V_OLD.DAT',status='UNKNOWN',form='UNFORMATTED')

!      quantities relevant for array size

       write(io) ntype,ng,ns,mxdl
       write(io) (natom(i),i=1,ntype)

!      date and calculation identification and parameters

       call zedate(bdate)
       call zetime(btime)
       write(io) bdate,btime
       
       write(io) author,flgscf,flgdal
       write(io) emax,teleck,nx,ny,nz,sx,sy,sz,nband,alatt

       if(iflag == 0) then
         write(io) pwline
       else
         line140(1:60) = pwline
         line140(61:110) = title
         line140(111:140) = subtitle(1:40)
         write(io) line140
       endif

!      crystal structure (other than quantities relevant for array size)

       write(io) ((adot(i,j),j=1,3),i=1,3)
       write(io) (nameat(i),i=1,ntype)

       do i=1,ntype
         write(io) ((rat(j,k,i),j=1,3),k=1,natom(i))
       enddo

!      G- space structure

       write(io) (kmax(j),j=1,3)
       write(io) ((kgv(j,k),j=1,3),k=1,ng)
       write(io) (phase(i),conj(i),i=1,ng)
       write(io) (mstar(i),i=1,ns)

!      charge density, bonding charge density, effective (Hartree+XC) potential

       write(io) (den(i),i=1,ns)
       write(io) (den(i)-dend(i),i=1,ns)      
       write(io) (veff(i),i=1,ns)

       call read_write_pseudo(io,ntype,nameat,                           &
     & mxdtyp)
       
       close(unit = io)


       end  program pw_rho_v_2_old
