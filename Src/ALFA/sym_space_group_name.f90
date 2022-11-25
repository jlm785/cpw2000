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

!>  Prints the possible space groups
!>
!>  \author       Jose Luis Martins
!>  \version      5.01
!>  \date         23 April 2021
!>  \copyright    GNU Public License v2

subroutine sym_space_group_name(ibravais, code_group, ntrans, mtrx, tnp)

! information taken from http://img.chem.ucl.ac.uk/sgp/large/sgp.htm

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  ibravais                         !<  bravais lattice (1=sc,2=fcc,3=bcc,4=st,5=ct,6=so,7=bco,8=fco,9=bco,10=hex,11=romb,12=sm,13=cm,14=tri)
  integer, intent(in)               ::  code_group                       !<  point group
  !! 1  & C_1 & 12 & C_2v & 23 & D_6h \\
  !! \hline
  !! 2  & C_i & 13 & C_3v & 24 & D_2d \\
  !! \hline
  !! 3  & C_s & 14 & C_4v & 25 & D_3d \\
  !! \hline
  !! 4  & C_2 & 15 & C_6v & 26 & S_4  \\
  !! \hline
  !! 5  & C_3 & 16 & C_2h & 27 & S_6  \\
  !! \hline
  !! 6  & C_4 & 17 & C_3h & 28 & T    \\
  !! \hline
  !! 7  & C_6 & 18 & C_4h & 29 & T_h  \\
  !! \hline
  !! 8  & D_2 & 19 & C_6h & 30 & T_d  \\
  !! \hline
  !! 9  & D_3 & 20 & D_2h & 31 & O    \\
  !! \hline
  !! 10 & D_4 & 21 & D_3h & 32 & O_h  \\
  !! \hline
  !! 11 & D_6 & 22 & D_4h &    &      \\

  integer, intent(in)               ::  ntrans                           !<  number of symmetry operations in the factor group
  integer, intent(in)               ::  mtrx(3,3,48)                     !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)          ::  tnp(3,48)                        !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! local variables

  logical        ::  linversion
  logical        ::  lfrac

  integer        ::  inv
  integer        ::  mtmp(3,3)

  integer        ::  io

! parameters

!   real(REAL64), parameter    ::  UM = 1.0_REAL64
!   real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-6_REAL64

! counters

  integer        ::  i, j

! checks if inversion is present

  io = 6

  linversion = .FALSE.
  if(ntrans > 1) then
    linversion = .TRUE.
    inv = ntrans/2+1
    do i = 1,3
      do j = 1,3
        mtmp(i,j) = mtrx(i,j,inv)
      enddo
      mtmp(i,i) = mtmp(i,i) + 1
    enddo
    do i = 1,3
    do j = 1,3
      if(mtmp(i,j) /= 0) linversion = .FALSE.
    enddo
    enddo
  endif

! checks if fractional translations are present ( beware: origin may be badly picked)

  lfrac = .FALSE.

  do j = 1,ntrans
    do i = 1,3
      if(abs(tnp(i,j)) > EPS) lfrac = .TRUE.
    enddo
  enddo

  write(io,*)
  if(ibravais == 1) then

    write(io,*) '  The Bravais lattice is  simple cubic'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  Pn-3 (201) Th2, Pa-3 (205) Th6'
          write(io,*) '  Pm-3 (200) Th1, if origin is badly chosen'
        else
          write(io,*) '  Pm-3 (200) Th1'
        endif
        if(code_group /= 29) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 48) then
        if(lfrac) then
          write(io,*) '  Pn-3n (222) Oh2, Pm-3n (223) Oh3, Pn-3m (224) Oh4'
          write(io,*) '  Pm-3m (221) Oh1, if origin is badly chosen'
        else
          write(io,*) '  Pm-3m (221) Oh1'
        endif
        if(code_group /= 32) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 12) then
        if(lfrac) then
          write(io,*) '  P2_13 (198) T4'
          write(io,*) '  P23 (195) T1, if origin is badly chosen'
        else
          write(io,*) '  P23 (195) T1'
        endif
        if(code_group /= 28) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 24) then
        if(code_group == 30) then
          if(lfrac) then
            write(io,*) '  P-43n (218) Td4'
            write(io,*) '  P-43m (215) Td1, if origin is badly chosen'
          else
            write(io,*) '  P-43m (215) Td1'
          endif
        elseif(code_group == 31) then
          if(lfrac) then
            write(io,*) '  P4_232 (208) O2, P4_332 (212) O6'
            write(io,*) '  P4_132 (213) O7'
            write(io,*) '  P432 (207) O1,  if origin is badly chosen'
          else
            write(io,*) '  P432 (207) O1'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 207, 208, 212, 213, 215, 218'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 2) then

    write(io,*) '  The Bravais lattice is fcc  face centered cubic'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  Fd-3 (203) Th4'
          write(io,*) '  Fm-3 (202) Th3, if origin is badly chosen'
        else
         write(io,*) '  Fm-3 (202) Th3'
        endif
        if(code_group /= 29) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 48) then
        if(lfrac) then
          write(io,*) '  Fm-3c (226) Oh6, Fd-3m (227) Oh7, Fd-3c (228) Oh8'
          write(io,*) '  Fm-3m (225) Oh5, if origin is badly chosen'
        else
          write(io,*) '  Fm-3m (225) Oh5'
        endif
        if(code_group /= 32) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 12) then
        if(lfrac) then
          write(io,*) '  F23 (196) T2, if origin is badly chosen'
        else
          write(io,*) '  F23 (196) T2'
        endif
        if(code_group /= 28) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 24) then
        if(code_group == 30) then
          if(lfrac) then
            write(io,*) '  F-43c (219) Td5'
            write(io,*) '  F-432 (216) Td2, if origin is badly chosen'
          else
            write(io,*) '  F-432 (216) Td2'
          endif
        elseif(code_group == 31) then
          if(lfrac) then
            write(io,*) '  F4_132 (210) O4'
            write(io,*) '  F-432 (216) O3, if origin is badly chosen'
          else
            write(io,*) '  F-432 (216) O3'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 209, 210, 216, 219'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 3) then

    write(io,*) '  The Bravais lattice is bcc  body centered cubic'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  Ia-3 (206) Th7'
          write(io,*) '  Im-3 (204) Th5, if origin is badly chosen'
        else
          write(io,*) '  The space group could be:'
          write(io,*) '  Im-3 (204) Th5'
        endif
        if(code_group /= 29) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 48) then
        if(lfrac) then
          write(io,*) '  Ia-3d (230) Oh10'
          write(io,*) '  Im-3m (229) Oh9, if origin is badly chosen'
        else
          write(io,*) '  Im-3m (229) Oh9'
        endif
        if(code_group /= 32) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 12) then
        if(lfrac) then
          write(io,*) '  I2_13 (199) T5'
          write(io,*) '  I23 (197) T3, if origin is badly chosen'
        else
          write(io,*) '  I23 (197) T3'
        endif
        if(code_group /= 28) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 24) then
        if(code_group == 30) then
          if(lfrac) then
            write(io,*) '  I-43d (220) Td6'
            write(io,*) '  I-43m (217) Td3, if origin is badly chosen'
          else
            write(io,*) '  I-43m (217) Td3'
          endif
        elseif(code_group == 31) then
          if(lfrac) then
            write(io,*) '  I4_132 (214) O8'
            write(io,*) '  I432 (211) O5, if origin is badly chosen'
          else
            write(io,*) '  I432 (211) O5'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 211, 214, 217, 220'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 4) then

    write(io,*) '  The Bravais lattice is simple tetragonal'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  P4_2/m (84) C4h2, P4/n (85) C4h3, P4_2/n (86) C4h5'
          write(io,*) '  P4/m (83) C4h1, if origin is badly chosen'
        else
          write(io,*) '  P4/m (83) C4h1'
        endif
        if(code_group /= 18) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 16) then
        if(lfrac) then
          write(io,*) '  P4/mcc (124) D4h2, P4/nbm (125) D4h3, P4/nnc (126) D4h4'
          write(io,*) '  P4/mbm (127) D4h5, P4/mnc (128) D4h6, P4/nmm (129) D4h7'
          write(io,*) '  P4/ncc (130) D4h8, P4_2/mmc (131) D4h9, P4_2/mcm (132) D4h10'
          write(io,*) '  P4_2/nbc (133) D4h11, P4_2/nnm (134) D4h12, P4_2/mbc (135) D4h13'
          write(io,*) '  P4_2/mnm (136) D4h14, P4_2/nmc (137) D4h15, P4_2/ncm (138) D4h16'
          write(io,*) '  P4/mmm (123) D4h1, if origin is badly chosen'
        else
          write(io,*) '  P4/mmm (123) D4h1'
        endif
        if(code_group /= 22) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 4) then
        if(code_group == 6) then
          if(lfrac) then
            write(io,*) '  P4_1 (76) C42, P4_2 (77) C43, P4_3 (78) C44'
            write(io,*) '  P4 (75) C41, if origin is badly chosen'
          else
            write(io,*) '  P4 (75) C41'
          endif
        elseif(code_group == 26) then
          write(io,*) '  P-4 (81) S41'
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 75, 76, 77, 78, 81'
          write(6,*)
        endif
      elseif(ntrans == 8) then
        if(code_group == 10) then
          if(lfrac) then
            write(io,*) '  P42_12 (90) D42, P4_122 (91) D43, P4_2_12 (92) D44'
            write(io,*) '  P4_222 (93) D45, P42_22_12 (94) D46'
            write(io,*) '  P4_322 (95) D47, P4_32_12 (96) D48'
            write(io,*) '  P422 (89) D41,  if origin is badly chosen'
          else
            write(io,*) '  P422 (89) D41'
          endif
        elseif(code_group == 14) then
          if(lfrac) then
            write(io,*) '  P4bm (100) C4v2, P4-2cm (101) C4v3'
            write(io,*) '  P4_2nm (102) C4v4, P4cc (103) C4v5, P4nc (104) C4v6'
            write(io,*) '  P4_2mc(105) C4v7, P4_2bc(106) C4v8'
            write(io,*) '  P4mm (99) C4v1, if origin is badly chosen'
          else
            write(io,*) '  P4mm (99) C4v1'
          endif
        elseif(code_group == 24) then
          if(lfrac) then
            write(io,*) '  P-42c (112) D2d2, P-42_1m (113) D2d3, P-42_1c (114) D2d4'
            write(io,*) '  P-4c2 (116) D2d6, P-4b2 (117) D2d7, P-4n2 (118) D2d8'
            write(io,*) '  P-42m (111) D2d1, P-4m2 (115) D2d5, if origin is badly chosen'
          else
            write(io,*) '  P-42m (111) D2d1, P-4m2 (115) D2d5'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 89->96, 99->106, 111->118'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 5) then

    write(io,*) '  The Bravais lattice is centered tetragonal'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  I4_1/a (88) C4h6'
          write(io,*) '  I4/m (87) C4h5, if origin is badly chosen'
        else
          write(io,*) '  I4/m (87) C4h5'
        endif
        if(code_group /= 18) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 16) then
        if(lfrac) then
          write(io,*) '  I4/mcm (140) D4h18, I4_1/amd (141) D4h19, I4_1/acd (142) D4h20'
          write(io,*) '  I4/mmm (139) D4h17, if origin is badly chosen'
        else
          write(io,*) '  I4/mmm (139) D4h17'
        endif
        if(code_group /= 22) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 4) then
        if(code_group == 6) then
          if(lfrac) then
            write(io,*) '  I4_1 (80) C46'
            write(io,*) '  I4 (79) C45, if origin is badly chosen'
          else
            write(io,*) '  I4 (79) C45'
          endif
        elseif(code_group == 26) then
          write(io,*) '  I-4 (82) S42'
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 79, 80, 82'
          write(6,*)
        endif
      elseif(ntrans == 8) then
        if(code_group == 10) then
          if(lfrac) then
            write(io,*) '  I4_122 (98) D410'
            write(io,*) '  I422 (97) D49, if origin is badly chosen'
          else
            write(io,*) '  I422 (97) D49'
          endif
        elseif(code_group == 14) then
          if(lfrac) then
            write(io,*) '  I4cm (108) C4v10, I4_1md(109) C4v11, I4_1cd (110) C4v12'
            write(io,*) '  I4mm (107) C4v9, if origin is badly chosen'
          else
            write(io,*) '  I4mm (107) C4v9'
          endif
        elseif(code_group == 24) then
          if(lfrac) then
            write(io,*) '  I-4c2 (120) D2d10, I-42m (121) D2d11, I-42d (122) D2d12'
            write(io,*) '  I-4m2 (119) D2d9, if origin is badly chosen'
          else
            write(io,*) '  I-4m2 (119) D2d9'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 97, 98, 107->110, 119->122'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 6) then

    write(io,*) '  The Bravais lattice is simple orthorhombic'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  Pnnn (48) D2h2, Pccm (49) D2h3, Pban (50) D2h4'
          write(io,*) '  Pmma (51) D2h5, Pnna (52) D2h6, Pmna (53) D2h7'
          write(io,*) '  Pcca (54) D2h8, Pbam (55) D2h9, Pccn (56) D2h10'
          write(io,*) '  Pbcm (57) D2h11, Pnnm (58) D2h12, Pmmn (59) D2h13'
          write(io,*) '  Pbcn (60) D2h14, Pbca (61) D2h15, Pnma (62) D2h16'
          write(io,*) '  Pmmm (47) D2h1, if origin is badly chosen'
        else
          write(io,*) '  Pmmm (47) D2h1'
        endif
        if(code_group /= 20) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 4) then
        if(code_group == 8) then
          if(lfrac) then
            write(io,*) '  P222_1 (17) D22, P2_12_12 (18) D23, P2_12_12_1 (19) D24'
            write(io,*) '  P222 (16) D21, if origin is badly chosen'
          else
            write(io,*) '  P222 (16) D21'
          endif
        elseif(code_group == 12) then
          if(lfrac) then
            write(io,*) '  Pmc2_1 (26) C2v2, Pcc2 (27) C2v3, Pma2 (28) C2v4'
            write(io,*) '  Pca2_1 (29) C2v5, Pnc2 (30) C2v6, Pmn2_1 (31) C2v7'
            write(io,*) '  Pba2 (32) C2v8, Pna2_1 (33) C2v9, Pnn2 (34) C2v10'
            write(io,*) '  Pmm2 (25) C2v1, if origin is badly chosen'
          else
            write(io,*) '  Pmm2 (25) C2v1'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 16->19, 25->34'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 7) then

    write(io,*) '  The Bravais lattice is base centered orthorhombic'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  Cmcm (63) D2h17, Cmca (64) D2h18, Cccm (66) D2h20'
          write(io,*) '  Cmma (67) D2h21, Ccca (68) D2h22'
          write(io,*) '  Cmmm (65) D2h19, if origin is badly chosen'
        else
          write(io,*) '  Cmmm (65) D2h19'
        endif
        if(code_group /= 20) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 4) then
        if(code_group == 8) then
          if(lfrac) then
            write(io,*) '  C222_1 (20) D25'
            write(io,*) '  C222 (21) D26, if origin is badly chosen'
          else
            write(io,*) '  C222 (21) D26'
          endif
        elseif(code_group == 12) then
          if(lfrac) then
            write(io,*) '  Cmc2_1 (36) C2v12, Ccc2 (37) C2v13'
            write(io,*) '  Abm2 (39) C2v15, Ama2 (40) C2v16, Aba2 (41) C2v17'
            write(io,*) '  Cmm2 (35) C2v11, Amm2 (38) C2v14, if origin is badly chosen'
          else
            write(io,*) '  Cmm2 (35) C2v11, Amm2 (38) C2v14'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 20, 21, 35->41'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 8) then

    write(io,*) '  The Bravais lattice is face centered orthorhombic'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  Fddd (70) D2h24'
          write(io,*) '  Fmmm (69) D2h23, if origin is badly chosen'
        else
          write(io,*) '  Fmmm (69) D2h23'
        endif
        if(code_group /= 20) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 4) then
        if(code_group == 8) then
          write(io,*) '  F222 (22) D27'
        elseif(code_group == 12) then
          if(lfrac) then
            write(io,*) '  Fdd2 (43) C2v19'
            write(io,*) '  Fmm2 (42) C2v18, if origin is badly chosen'
          else
            write(io,*) '  Fmm2 (42) C2v18'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 22, 42, 43'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 9) then

    write(io,*) '  The Bravais lattice is body centered orthorhombic'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 8) then
        if(lfrac) then
          write(io,*) '  Ibam (72) D2h26, Ibca (73) D2h27, Imma (74) D2h28'
          write(io,*) '  Immm (71) D2h25, if origin is badly chosen'
        else
          write(io,*) '  Immm (71) D2h25'
        endif
        if(code_group /= 20) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 4) then
        if(code_group == 8) then
          if(lfrac) then
            write(io,*) '  I2_12_12_1 (24) D29'
            write(io,*) '  I222 (23) D28, if origin is badly chosen'
          else
            write(io,*) '  I222 (23) D28'
          endif
        elseif(code_group == 12) then
          if(lfrac) then
            write(io,*) '  Iba2 (45) C2v21, Ima2 (46) C2v22'
            write(io,*) '  Imm2 (44) C2v20, if origin is badly chosen'
          else
            write(io,*) '  Imm2 (44) C2v20'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 23, 24, 44->46'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 10) then

    write(io,*) '  The Bravais lattice is hexagonal'
    write(io,*) '  including non-rhombohedral trigonal'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 6) then
        write(io,*) '  P-3 (147) if origin is badly chosen'
        if(code_group /= 27) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 12) then
        if(code_group == 25) then
          if(lfrac) then
            write(io,*) '  P-31c (163) D3d2, P-3c1 (165) D3d4'
            write(io,*) '  P-31m (162) D3d1, P-3m1 (164) D3d3,  if origin is badly chosen'
          else
            write(io,*) '  P-31m (162) D3d1, P-3m1 (164) D3d3'
          endif
        elseif(code_group == 19) then
          if(lfrac) then
            write(io,*) '  P6_3/m (176) C6h2'
            write(io,*) '  P6/m (175) C6h1, if origin is badly chosen'
          else
            write(io,*) '  P6/m (175) C6h1'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 162->165, 175, 176'
          write(6,*)
        endif
      elseif(ntrans == 24) then
        if(lfrac) then
          write(io,*) '  P6/mmc (192) D6h2, P6_3/mcm (193) D6h3, P6_3/mmc (194) D6h4'
          write(io,*) '  P6/mmm (191) D6h1, if origin is badly chosen'
        else
          write(io,*) '  P6/mmm (191) D6h1'
        endif
        if(code_group /= 23) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 3) then
        if(lfrac) then
          write(io,*) '  P3_1 (144) C32, P3_2 (145) C33'
          write(io,*) '  P3 (143) C31, if origin is badly chosen'
        else
          write(io,*) '  P3 (143) C31'
        endif
        if(code_group /= 5) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 6) then
        if(code_group == 9) then
          if(lfrac) then
            write(io,*) '  P3_112 (151) D33, P3_121 (152) D34'
            write(io,*) '  P3_212 (153) D35, P3_221 (154) D36'
            write(io,*) '  P312 (149) D31, P321 (150) D32, if origin is badly chosen'
          else
            write(io,*) '  P312 (149) D31, P321 (150) D32'
          endif
        elseif(code_group == 13) then
          if(lfrac) then
            write(io,*) '  P3c1 (158) C3v3, P31c (159) C3v4'
            write(io,*) '  P3m1 (156) C3v1, P31m (157) C3v2, if origin is badly chosen'
          else
            write(io,*) '  P3m1 (156) C3v1, P31m (157) C3v2'
          endif
        elseif(code_group == 7) then
          if(lfrac) then
            write(io,*) '  P6_1 (169) C62, P6_5 (170) C63, P6_2 (171) C64'
            write(io,*) '  P6_4 (172) C65, P6_3 (173) C66'
            write(io,*) '  P6 (168) C61, if origin is badly chosen'
          else
            write(io,*) '  P6 (168) C61'
          endif
        elseif(code_group == 17) then
            write(io,*) '  P-6 (174) C3h1'
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 149-154, 156->159, 168->174'
          write(6,*)
        endif
      elseif(ntrans == 12) then
        if(code_group == 11) then
          if(lfrac) then
            write(io,*) '  P6_122 (178) D62, P6_522 (179) D63, P6_222 (180) D64'
            write(io,*) '  P6_422 (181) D65, P6_322 (182) D66'
            write(io,*) '  P622 (177) D61, if origin is badly chosen'
          else
            write(io,*) '  P622 (177) D61'
          endif
        elseif(code_group == 15) then
          if(lfrac) then
            write(io,*) '  P6cc (184) C6v2, P6_3cm (185) C6v3, P6_3mc (186) C6v4'
            write(io,*) '  P6mm (183) C6v1, if origin is badly chosen'
          else
            write(io,*) '  P6mm (183) C6v1'
          endif
        elseif(code_group == 21) then
          if(lfrac) then
            write(io,*) '  P-6c2 (188) D3h2, P-62c (190) D3h4'
            write(io,*) '  P-6m2 (187) D3h1, P-62m (189) D3h3, if origin is badly chosen'
          else
            write(io,*) '  P-6m2 (187) D3h1, P-62m (189) D3h3'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 177->190'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 11) then

    write(io,*) '  The Bravais lattice is rhombohedral'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 6) then
        write(io,*) '  R-3 (148) S62'
        if(code_group /= 27) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 12) then
        if(lfrac) then
          write(io,*) '  R-3c (167) D3d6'
          write(io,*) '  R-3m (166) D3d5, if origin is badly chosen'
        else
          write(io,*) '  R-3m (166) D3d5'
        endif
        if(code_group /= 25) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 3) then
        write(io,*) '  R3 (146) C34'
        if(code_group /= 5) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      elseif(ntrans == 6) then
        if(code_group == 9) then
          write(io,*) '  R32 (155) D37'
        elseif(code_group == 13) then
          if(lfrac) then
            write(io,*) '  R3c (161) C3v6'
            write(io,*) '  R3m (160) C3v5, if origin is badly chosen'
          else
            write(io,*) '  R3m (160) C3v5'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 155, 160, 161'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 12) then

    write(io,*) '  The Bravais lattice is simple monoclinic'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 4) then
        if(lfrac) then
          write(io,*) '  P2_1/m (11) C2h2, P2/c (13) C2h4, P2_1/c (14) C2h5'
          write(io,*) '  P2/m (10) C2h1, if origin is badly chosen'
        else
          write(io,*) '  P2/m (10) C2h1'
        endif
        if(code_group /= 16) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 2) then
        if(code_group == 4) then
          if(lfrac) then
            write(io,*) '  P2_1 (4) C22'
            write(io,*) '  P2 (3) C21, if origin is badly chosen'
          else
            write(io,*) '  P2 (3) C21'
          endif
        elseif(code_group == 3) then
          if(lfrac) then
            write(io,*) '  Pc (7) Cs2'
            write(io,*) '  Pm (6) Cs1, if origin is badly chosen'
          else
            write(io,*) '  Pm (6) Cs1'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 3, 4, 6, 7'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 13) then

    write(io,*) '  The Bravais lattice is centered monoclinic'
    write(io,*)
    write(io,*) '  The space group could be:'
    write(io,*)

    if(linversion) then

      if(ntrans == 4) then
        if(lfrac) then
          write(io,*) '  C2/c (15) C2h6'
          write(io,*) '  C2/m (12) C2h3, if origin is badly chosen'
        else
          write(io,*) '  C2/m (12) C2h3'
        endif
        if(code_group /= 16) then
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*)
        endif
      endif

    else

      if(ntrans == 2) then
        if(code_group == 4) then
            write(io,*) '  C2 (5) C23'
        elseif(code_group == 3) then
          if(lfrac) then
            write(io,*) '  Cc (9) Cs4'
            write(io,*) '  Cm (8) Cs3, if origin is badly chosen'
          else
            write(io,*) '  Cm (8) Cs3'
          endif
        else
          write(6,*)
          write(6,*) '  Inconsistent information. Inform author'
          write(6,*) '  Check 5, 8, 9'
          write(6,*)
        endif
      endif

    endif

  elseif(ibravais == 14) then

    write(io,*) '  The Bravais lattice is triclinic'
    write(io,*)
    write(io,*) '  The space group is:'
    write(io,*)

    if(linversion) then

      write(io,*) '  P-1 (2) Ci'
      if(code_group /= 2) then
        write(6,*)
        write(6,*) '  Inconsistent information. Inform author'
        write(6,*)
      endif

    else

      write(io,*) '  P1 (1) C1'
      if(code_group /= 1) then
        write(6,*)
        write(6,*) '  Inconsistent information. Inform author'
        write(6,*)
      endif

    endif

  endif

  write(io,*)
  write(io,*) '  As not all 230 possible space groups were tested,'
  write(io,*) '  check this space group information carefuly.'
  write(io,*)

  return
end subroutine sym_space_group_name
