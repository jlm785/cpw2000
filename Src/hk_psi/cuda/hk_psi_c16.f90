       subroutine hk_psi_c16(mtxd,neig,psi,hpsi,lnewanl,                 &
     & ng,kgv,                                                           &
     & qmod,isort,vscr,kmscr,                                            &
     & anlga,xnlkb,nanl,                                                 &
     & mxddim,mxdbnd,mxdanl,mxdgve,mxdscr)

!      calculates the product of the hamiltonian times
!      neig wavevectors. the non-local pseudopotential
!      is separable. the local potential is dealt with
!      fast fourier transforms.

!      INTERFACE TO CUDA LIBRARY

!      copyright INESC-MN/Jose Luis Martins


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxddim                     !  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !  array dimension for number of bands
       integer, intent(in)                ::  mxdanl                     !  array dimension of number of projectors
       integer, intent(in)                ::  mxdgve                     !  array dimension for g-space vectors
       integer, intent(in)                ::  mxdscr                     !  array dimension of vscr

       integer, intent(in)                ::  mtxd                       !  wavefunction dimension (basis size)
       integer, intent(in)                ::  neig                       !  number of wavefunctions
       real(REAL64), intent(in)           ::  qmod(mxddim)               !  length of k+g-vector of row/column i
       integer, intent(in)                ::  isort(mxddim)              !  g-vector associated with row/column i of hamiltonian

       real(REAL64), intent(in)           ::  vscr(mxdscr)               !  screened potential in the fft real space mesh
       integer, intent(in)                ::  kmscr(3)                   !  max value of kgv(i,n) used for the potential fft mesh

       integer, intent(in)                ::  ng                         !  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,mxdgve)              !  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

       integer, intent(in)                ::  nanl                       !  number of projectors
       complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)       !  Kleinman-Bylander projectors
       real(REAL64), intent(in)           ::  xnlkb(mxdanl)              !  Kleinman-Bylander normalization
       complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)         !  wavevector

!      input and output

       logical, intent(inout)             ::  lnewanl                    !  indicates that anlga has been recalculated (not used in default implementation)

!      output

       complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)        !  |hpsi> =  V_NL |psi>


       call compact_expand(psi,mtxd,neig,mxddim,'C')
       call compact_expand(anlga,mtxd,nanl,mxddim,'C')


        call hk_psi_c16_cuda(mtxd,neig,psi,hpsi, lnewanl,                &
     & ng,kgv,                                                           &
     & qmod,isort,vscr,kmscr,                                            &
     & anlga,xnlkb,nanl,                                                 &
     & mtxd,mxdbnd,mxdanl,mxdgve,mxdscr)

       call compact_expand(psi,mtxd,neig,mxddim,'E')
       call compact_expand(hpsi,mtxd,neig,mxddim,'E')
       call compact_expand(anlga,mtxd,nanl,mxddim,'E')

      

       return

       end subroutine hk_psi_c16
