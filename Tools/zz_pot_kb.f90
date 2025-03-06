program zz_potkb

! writes a Coulomb potential or a null potential or any other
! bespoke atomic potential in the KB format for the cpw2000 code.

! Change izv for Coulomb and null.  Or write whatever you want following
! the commented examples.

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! parameters

  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  integer, parameter          ::  mxdnqp = 101, mxdlqp = 4005
  integer, parameter          ::  mxdnga = 30

! strings

  character(len=2)            ::  namel,icorr,icorrt
  character(len=3)            ::  irel
  character(len=4)            ::  icore
  character(len=10)           ::  iray(6)
  character(len=10)           ::  psdtitle(20)

  character(len=30)           ::  filename

! pseudo arrays

  real(REAL64)                ::  vqnl(mxdnqp,mxdnqp)
  integer                     ::  ngaunl(3),xxnl(mxdnga),wgnl(mxdnga)
  real(REAL64)                ::  vda(mxdlqp)
  integer                     ::  nkb(3)

! other variables

  integer                     ::  it                                     !  tape number
  integer                     ::  izv                                    !  nuclear charge

  integer                     ::  nql                                    !  number of mesh points
  real(REAL64)                ::  delql, vql0                            !  mesh step, limit r->0

  integer                     ::  norb                                   !  number of KB projectors
  real(REAL64)                ::  eig                                    !  eigenvalue of associated orbital

  real(REAL64)                ::  q                                      !  length reciprocal vector

  integer                     ::  nqwf                                   !  number of mesh points
  real(REAL64)                ::  deltaq                                 !  mesh step, limit r->0
  integer                     ::  norbat                                 !  number of atomic orbitals

  real(REAL64)                ::  q0sq                                   !  square of the extension of wave-function (reciprocal space)


! counters

  integer      ::  i, j, lp1

! tape number

  it = 10

! "nuclear charge"

! Use izv=0 for potential without electrons, izv=1 for hydrogen, izv=2 for helium.

  izv = 0

  filename = 'ZZ_POTKB_F.DAT'
  if(izv == 1) filename = 'H_POTKB_F.DAT'
  if(izv == 2) filename = 'He_POTKB_F.DAT'

  open (unit=it, file=trim(filename), status='new', form='formatted')

!   title

  namel = 'ZZ'
  if(izv == 1) namel = 'H '
  if(izv == 2) namel = 'He'

  icorrt = 'CA'
  irel = '   '
  icore = '    '

  do j=1,6
    iray(j) = '          '
  enddo
  if(izv == 0) then
    iray(1) = 'Bespoke Po'
    iray(2) = 'tential   '
  else
    iray(1) = 'True Poten'
    iray(2) = 'tial      '
  endif

  do j=1,20
    psdtitle(j) = '          '
  enddo

  write(it,'(1x,a2,1x,a2,1x,a3,1x,a4,1x,6a10,1x,7a10)')                  &
       namel, icorrt, irel, icore, (iray(i),i=1,6), (psdtitle(i),i=1,20)

! valence and mesh sizes

  nql = 4000
  delql = 0.01_REAL64
  vql0 = 0.0
  write(it,*) izv, nql, delql, vql0

  norb = 1
  nkb(1) = 0
  eig = (izv*izv*UM)/2
  write(it,*) norb
  write(it,*) 0
  write(it,*) (nkb(i),i=1,norb)
  write(it,*) eig

! local potential

!  Example of a repulsive gaussian potential
!
!  do j=1,nql
!    q = real(j)*delql
!    vda(j) = 1.185*exp(-0.1*q*q)
!  enddo


  do j = 1,nql
    q = j*delql
    vda(j) = -izv*8*PI/(q*q)
  enddo
  do j = 1,nql
    write(it,*) vda(j)
  enddo

! non-local potential

  do lp1 = 1,norb
    do j=1,nql+1
      vda(j) = ZERO
    enddo
    do j = 1,nql+1
      write(it,*) vda(j)
    enddo
  enddo

!   core charge

  do j = 1,nql
    vda(j) = ZERO
  enddo
  do j = 1,nql
    write(it,*) vda(j)
  enddo

! valence charge

  do j = 1,nql
    q = j*delql
    vda(j) = (16*izv*UM)/((4*UM+q*q)*(4*UM+q*q))
  enddo
  do j = 1,nql
    write(it,*) vda(j)
  enddo

  nqwf = nql
  deltaq = delql
  norbat = 1
  if(izv == 0) norbat = 0
  write(it,*) nqwf,deltaq,norbat

  if(norbat > 0) then

    q0sq = UM
    if(izv == 2) q0sq = 2.73
    do j=1,nql
      q = j*delql
      vda(j) = (16*UM)/((UM+q*q/q0sq)*(UM+q*q/q0sq))
    enddo

    write(it,*) 0, eig
    do j=1,nqwf
      write(it,*) vda(j)
    enddo

  endif

  close(unit=it)


  stop

end program zz_potkb
