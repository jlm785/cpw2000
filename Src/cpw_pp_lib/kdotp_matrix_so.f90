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

!>     Calculates the matrices for the k.p method
!>     OLD SUBROUTINE keep for testing purposes

       subroutine kdotp_matrix_so(mtxd,neig,psi,ei,rkpt,isort,           &
     & hso0,dhso0drk,d2hso0drk2,                                         &
     & ng,kgv,                                                           &
     & ntype,natom,rat,adot,                                             &
     & nqnl,delqnl,vkb,nkb,                                              &
     & mxdtyp,mxdatm,mxdlqp,mxddim,mxdbnd,mxdgve)

!      Written April 14, 2014, from previous code. JLM
!      Modified 4 March 2020. Documentation, nder. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.95


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type
       integer, intent(in)                ::  mxdlqp                     !<  array dimension for local potential
       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands
       integer, intent(in)                ::  mxdgve                     !<  array dimension of G-space vectors

       integer, intent(in)                ::  mtxd                       !<  wavefunction dimension
       integer, intent(in)                ::  neig                       !<  number of wavefunctions
       real(REAL64), intent(in)           ::  rkpt(3)                    !<  k-point reciprocal lattice coordinates 
       integer, intent(in)                ::  isort(mxddim)              !<  G-vector corresponding to coefficient i of wavefunction 
       complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)         !<  wavevector
       real(REAL64), intent(in)           ::  ei(mxdbnd)                 !<  eigenvalues (Hartree) for rk0

       integer, intent(in)                ::  ng                         !<  size of g-space
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  G-vectors in reciprocal lattice coordinates 

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  lattice coordinates of atom j of type i
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in real space

       integer, intent(in)                ::  nqnl(mxdtyp)               !<  number of points for pseudo interpolation for atom k
       real(REAL64), intent(in)           ::  delqnl(mxdtyp)             !<  step used in the pseudo interpolation for atom k
       real(REAL64), intent(in)    ::    vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
       integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)       !<  KB pseudo.  normalization for atom k, ang. mom. l

!      output

       complex(REAL64), intent(out)       ::  hso0(2*mxdbnd,2*mxdbnd)          !<  <Psi|H|Psi> without spin-orbit
       complex(REAL64), intent(out)       :: dhso0drk(2*mxdbnd,2*mxdbnd,3)     !<  d <Psi|H|Psi> d k
       complex(REAL64), intent(out)       :: d2hso0drk2(2*mxdbnd,2*mxdbnd,3,3) !<  d^2 <Psi|H|Psi> d k^2

!      local allocatable arrays

       complex(REAL64), allocatable       :: pmat(:,:,:)                 !  <Psi|p_j|Psi>

       complex(REAL64), allocatable       ::  anlso(:,:)                 !  KB projectors with spin-orbit
       real(REAL64), allocatable          ::  xnlkbso(:)                 !  KB normalization without spin-orbit
       complex(REAL64), allocatable       ::  anlspin(:,:)               !  KB projectors without spin-orbit but spin representation
       real(REAL64), allocatable          ::  xnlkbspin(:)               !  KB normalization without spin-orbit but spin representation

       complex(REAL64), allocatable       :: danlspindrk(:,:,:)          !  d anlspin / d rkpt
       complex(REAL64), allocatable       :: d2anlspindrk2(:,:,:,:)      !  d^2 anlspin / d rkpt^2
       complex(REAL64), allocatable       :: danlsodrk(:,:,:)            !  d anlso / d rkpt
       complex(REAL64), allocatable       :: d2anlsodrk2(:,:,:,:)        !  d^2 anlso / d rkpt^2

       complex(REAL64), allocatable       ::  vnlso0(:,:)                !  <Psi|V_NL|Psi>
       complex(REAL64), allocatable       ::  dvnlso0drk(:,:,:)          !  d <Psi|V_NL|Psi> d k
       complex(REAL64), allocatable       ::  d2vnlso0drk2(:,:,:,:)      !  d^2 <Psi|V_NL|Psi> d k^2

       complex(REAL64), allocatable       ::  psi_so0(:,:)               !  component j of eigenvector i (guess on input)

!      local variables

       integer           ::  mxdanl         !  array dimension of number of projectors
       integer           ::  mxdaso         !  array dimension of number of projectors with spin-orbit
       integer           ::  nanl, nanlso   !  number of KB projectors
       real(REAL64)      ::  vcell, bdot(3,3)

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64
       complex(REAL64), parameter  :: C_ZERO = cmplx(ZERO,ZERO,REAL64)

!      counters

       integer    ::  i, j, m, n


       call adot_to_bdot(adot,vcell,bdot)

       allocate(pmat(3,mxdbnd,mxdbnd))

       call psi_p_psi(mtxd,neig,psi,pmat,rkpt,isort,ng,kgv,.FALSE.,      &
     & mxddim,mxdbnd,mxdgve)

       call size_proj_nl_kb(ntype,natom,nkb,nanl,nanlso,mxdtyp)
         
       mxdanl = nanl
       mxdaso = nanlso

       allocate(anlspin(2*mxddim,2*mxdanl))
       allocate(xnlkbspin(2*mxdanl))
       allocate(anlso(2*mxddim,mxdaso))
       allocate(xnlkbso(mxdaso))
         
       allocate(danlspindrk(2*mxddim,2*mxdanl,3))
       allocate(d2anlspindrk2(2*mxddim,2*mxdanl,3,3))
       allocate(danlsodrk(2*mxddim,mxdaso,3))
       allocate(d2anlsodrk2(2*mxddim,mxdaso,3,3))


       call proj_nl_kb_so_der(rkpt,mtxd,isort,nanl,nanlso,               &
     & ng,kgv,                                                           &
     & nqnl,delqnl,vkb,nkb,                                              &
     & ntype,natom,rat,adot,                                             &
     & anlspin,anlso,xnlkbspin,xnlkbso,                                  &
     & danlspindrk,danlsodrk,d2anlspindrk2,d2anlsodrk2,                  &
     & mxdtyp,mxdatm,mxdlqp,mxddim,mxdanl,mxdaso,mxdgve)

       allocate(vnlso0(2*mxdbnd,2*mxdbnd))
       allocate(dvnlso0drk(2*mxdbnd,2*mxdbnd,3))
       allocate(d2vnlso0drk2(2*mxdbnd,2*mxdbnd,3,3))

       allocate(psi_so0(2*mxddim,2*mxdbnd))
         
       do n = 1,neig
       do m = 1,mxddim
         psi_so0(2*m-1,2*n-1) = psi(m,n)
         psi_so0(2*m  ,2*n-1) = C_ZERO
         psi_so0(2*m-1,2*n  ) = C_ZERO
         psi_so0(2*m  ,2*n  ) = psi(m,n)
       enddo
       enddo

!      spin-orbit 

       call psi_vnl_psi_der(2*mtxd,2*neig,nanlso,psi_so0,2,              &
     & vnlso0,dvnlso0drk,d2vnlso0drk2,                                   &
     & anlso,xnlkbso,danlsodrk,d2anlsodrk2,                              &
     & 2*mxddim,2*mxdbnd,mxdaso)

       do i=1,2*neig
       do j=1,2*neig
         hso0(j,i) = vnlso0(j,i)
       enddo
       enddo

       call psi_vnl_psi(2*mtxd,2*neig,psi_so0,vnlso0,                    &
     & anlspin,xnlkbspin,2*nanl,                                         &
     & 2*mxddim,2*mxdbnd,2*mxdanl)

       do i=1,2*neig
       do j=1,2*neig
         hso0(j,i) = hso0(j,i) - vnlso0(j,i)
       enddo
       enddo

!      eigenvalues

       do j=1,neig
         hso0(2*j-1,2*j-1) = hso0(2*j-1,2*j-1) + cmplx(ei(j),ZERO,REAL64)
         hso0(2*j  ,2*j  ) = hso0(2*j  ,2*j  ) + cmplx(ei(j),ZERO,REAL64)
       enddo
       
!      kinetic and non-local energy

       do i=1,2*neig
       do j=1,2*neig
         do m=1,3
           dhso0drk(j,i,m) = dvnlso0drk(j,i,m)
         enddo
       enddo
       enddo

       do i=1,neig
       do j=1,neig
         do m=1,3
         do n=1,3
           dhso0drk(2*j-1,2*i-1,m)  = dhso0drk(2*j-1,2*i-1,m)            &
     &                                    + bdot(m,n)*pmat(n,j,i)
           dhso0drk(2*j  ,2*i  ,m)  = dhso0drk(2*j  ,2*i  ,m)            &
     &                                    + bdot(m,n)*pmat(n,j,i)
         enddo
         enddo
       enddo
       enddo


       do i=1,2*neig
       do j=1,2*neig
         do n=1,3
         do m=1,3
           d2hso0drk2(j,i,m,n) = d2vnlso0drk2(j,i,m,n)
         enddo
         enddo
       enddo
       enddo

       do j=1,2*neig
         do n=1,3
         do m=1,3
           d2hso0drk2(j,j,m,n) = d2hso0drk2(j,j,m,n) + bdot(m,n)
         enddo
         enddo
       enddo

       deallocate(pmat)

       deallocate(anlspin)
       deallocate(xnlkbspin)
       deallocate(anlso)
       deallocate(xnlkbso)

       deallocate(danlspindrk)
       deallocate(d2anlspindrk2)
       deallocate(danlsodrk)
       deallocate(d2anlsodrk2)

       deallocate(vnlso0)
       deallocate(dvnlso0drk)
       deallocate(d2vnlso0drk2)

       deallocate(psi_so0)

       return

       end subroutine kdotp_matrix_so

