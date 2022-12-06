subroutine pw_rho_v_in_crystal_sym(filename, io,                         &
         adot,ntype,natom,nameat,rat,                                    &
         ntrans, mtrx, tnp,                                              &
         mxdtyp, mxdatm)

! This subroutines reads the first part of file io,
! and returns the symmetry operations

! Written November 11 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  character(len=*), intent(in)       ::  filename                        !<  name of file
  integer, intent(in)                ::  io                              !<  number of tape to which the pseudo is added.

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms

! output

  real(REAL64), intent(out)          ::  adot(3,3)                       !<  metric in direct space
  integer, intent(out)               ::  ntype                           !<  number of types of atoms
  integer, intent(out)               ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(out)      ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(out)               ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(out)               ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(out)          ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group


! other local variables

  integer             ::  ioerr

  integer               ::  ng                              !  total number of g-vectors with length less than gmax
  integer               ::  ns                              !  number os stars with length less than gmax

  integer               ::  mxdl
! counters

  integer    ::  i, j, k


! reads the first record from the file

  open(unit=io,file=trim(filename),status='old',form='UNFORMATTED')

  read(io,iostat = ioerr) ntype, ng, ns, mxdl, ntrans

  if(ioerr == 0) then
    if(ntrans > 48 .or. ntrans < 1) then
      write(6,*)
      write(6,'("  STOPPED in pw_rho_v_in_crystal_sym")')
      write(6,'("  ntrans = ",i8)') ntrans

      stop

    endif
  else
    write(6,*)
    write(6,'("  WARNING in pw_rho_v_in_crystal_sym")')
    write(6,'("  old file style does not have the symmetry information ")')

    ntrans = 0

  endif


  read(io)    (natom(i),i=1,ntype)

  read(io)    ! bdate,btime

  read(io)    ! author,flgscf,flgdal

! backcompatibilty blues

  read(io)    ! emax, teleck, nx,ny,nz, sx,sy,sz, nband, alatt, efermi

  read(io)    ! line250, meta_cpw2000

! old files do not have spacegroup information

  if(ntrans == 0) then
    read(io) ( (adot(i,j),i=1,3),j=1,3 ),                               &
             ( ((mtrx(i,j,k),j=1,3),i=1,3), k=1,ntrans ),               &
             ( (tnp(i,k),i=1,3), k=1,ntrans )
  else
    read(io) ( (adot(i,j),i=1,3),j=1,3 )
  endif

  read(io)  (nameat(i),i=1,ntype)

  do i=1,ntype
    read(io) ((rat(j,k,i),j=1,3),k=1,natom(i))
  enddo

  close(unit=io)

  return

end subroutine pw_rho_v_in_crystal_sym

