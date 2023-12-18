!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_enefunc_amber_mod
!> @brief   define potential energy functions
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_enefunc_amber_mod

  use at_enefunc_restraints_mod
  use at_enefunc_table_mod
  use at_enefunc_pme_mod
  use at_energy_mod
  use at_boundary_str_mod
  use at_restraints_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use math_libs_mod
  use molecules_str_mod
  use fileio_prmtop_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! constants
  integer, private, parameter :: MAXWRK   = 12
  integer, private, parameter :: MAX_EXCL = 16 ! = 4(1-2) + 4x3 (1-3)
  integer, private, parameter :: MAX_NB14 = 36 ! = 3x4x3 (1-4) 

  ! subroutines
  public  :: define_enefunc_amber
  private :: setup_enefunc_bond
  private :: setup_enefunc_angl
  private :: setup_enefunc_dihe
  private :: setup_enefunc_impr
  private :: setup_enefunc_nonb
  private :: count_nonb_excl

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    define_enefunc_amber
  !> @brief        a driver subroutine for defining potential energy
  !! @authors      NT
  !! @param[in]    ene_info   : ENERGY section control parameters information
  !! @param[in]    boundary   : boundary condition
  !! @param[in]    prmtop     : AMBER parameter topology information
  !! @param[in]    molecule   : molecule information
  !! @param[in]    restraints : restraints information
  !! @param[inout] enefunc    : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine define_enefunc_amber(ene_info, boundary, prmtop, &
                                  molecule, restraints, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info 
    type(s_boundary),        intent(in)    :: boundary
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: alloc_stat
    integer                  :: i, found, found2
    integer                  :: nwork


    ! bond
    !
    call setup_enefunc_bond(prmtop, enefunc)

    ! angle
    !
    call setup_enefunc_angl(prmtop, enefunc)

    ! dihedral
    !
    call setup_enefunc_dihe(prmtop, enefunc)

    ! improper
    !
    call setup_enefunc_impr(prmtop, enefunc)

    ! nonbonded
    !
    call setup_enefunc_nonb(ene_info, prmtop, molecule, enefunc)

    ! lookup table
    !
    call setup_enefunc_table(ene_info, molecule, enefunc)

    ! PME
    !
    call define_enefunc_pme(ene_info, boundary, molecule, enefunc)

    ! restraints
    !
    call setup_enefunc_restraints(molecule, restraints, enefunc)

    ! allocate working array
    !
    nwork = max(molecule%num_atoms, &
                enefunc%num_bonds,  &
                enefunc%num_angles, &
                enefunc%num_dihedrals, &
                enefunc%num_impropers, &
                enefunc%num_restraintfuncs, &
                enefunc%num_restraintgroups)

    allocate(enefunc%work(MAXWRK,nwork), stat=alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! write summary of energy function
    !
    if (main_rank) then

      if (.not. ene_info%table) then

        found  = 0
        found2 = 0

        do i = 1, molecule%num_atoms
          found  = found  + enefunc%num_nonb_excl(i)
          found2 = found2 + enefunc%num_nb14_calc(i)
        end do

      end if

      write(MsgOut,'(A)') &
           'Define_Enefunc_AMBER> Number of Interactions in Each Term'
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  bond_ene        = ', enefunc%num_bonds,           &
           '  angle_ene       = ', enefunc%num_angles
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  torsion_ene     = ', enefunc%num_dihedrals,       &
           '  improper_ene    = ', enefunc%num_impropers
      if (.not. ene_info%table) then
        write(MsgOut,'(A20,I10,A20,I10)')                       &
           '  nb_exclusions   = ', found,                       &
           '  nb14_calc       = ', found2
      end if
      write(MsgOut,'(A20,I10,A20,I10)')                         &
           '  restraint_groups= ', enefunc%num_restraintgroups, &
           '  restraint_funcs = ', enefunc%num_restraintfuncs
      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine define_enefunc_amber

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_bond
  !> @brief        define BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_bond(prmtop, enefunc)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: nbond, i, icnt
    integer                  :: istart, iend


    nbond = prmtop%num_bondh + prmtop%num_mbonda

    call alloc_enefunc(enefunc, EneFuncBond, nbond)


    icnt = 0
    do i = 1, prmtop%num_bondh

      icnt = icnt + 1

      enefunc%bond_list(1,icnt) = prmtop%bond_inc_hy(1,i) / 3 + 1
      enefunc%bond_list(2,icnt) = prmtop%bond_inc_hy(2,i) / 3 + 1

      enefunc%bond_force_const(icnt) = &
             prmtop%bond_fcons_uniq(prmtop%bond_inc_hy(3,i))
      enefunc%bond_dist_min(icnt)    = &
             prmtop%bond_equil_uniq(prmtop%bond_inc_hy(3,i))
    end do

    do i = 1, prmtop%num_mbonda

      icnt = icnt + 1

      enefunc%bond_list(1,icnt) = prmtop%bond_wo_hy(1,i) / 3 + 1
      enefunc%bond_list(2,icnt) = prmtop%bond_wo_hy(2,i) / 3 + 1

      enefunc%bond_force_const(icnt) = &
             prmtop%bond_fcons_uniq(prmtop%bond_wo_hy(3,i))
      enefunc%bond_dist_min(icnt)    = &
             prmtop%bond_equil_uniq(prmtop%bond_wo_hy(3,i))

    end do

    enefunc%num_bonds = nbond


    call get_loop_index(enefunc%num_bonds, istart, iend)

    enefunc%istart_bond = istart
    enefunc%iend_bond   = iend

    return

  end subroutine setup_enefunc_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_angl
  !> @brief        define ANGLE term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_angl(prmtop, enefunc)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: nangl, i, icnt
    integer                  :: istart, iend


    nangl = prmtop%num_anglh + prmtop%num_mangla

    call alloc_enefunc(enefunc, EneFuncAngl, nangl)


    icnt = 0
    do i = 1, prmtop%num_anglh

      icnt = icnt + 1

      enefunc%angl_list(1,icnt) = prmtop%angl_inc_hy(1,i) / 3 + 1
      enefunc%angl_list(2,icnt) = prmtop%angl_inc_hy(2,i) / 3 + 1
      enefunc%angl_list(3,icnt) = prmtop%angl_inc_hy(3,i) / 3 + 1

      enefunc%angl_force_const(icnt) = &
             prmtop%angl_fcons_uniq(prmtop%angl_inc_hy(4,i))
      enefunc%angl_theta_min(icnt)   = &
             prmtop%angl_equil_uniq(prmtop%angl_inc_hy(4,i))

    end do

    do i = 1, prmtop%num_mangla

      icnt = icnt + 1

      enefunc%angl_list(1,icnt) = prmtop%angl_wo_hy(1,i) / 3 + 1
      enefunc%angl_list(2,icnt) = prmtop%angl_wo_hy(2,i) / 3 + 1
      enefunc%angl_list(3,icnt) = prmtop%angl_wo_hy(3,i) / 3 + 1

      enefunc%angl_force_const(icnt) = &
             prmtop%angl_fcons_uniq(prmtop%angl_wo_hy(4,i))
      enefunc%angl_theta_min(icnt)   = &
             prmtop%angl_equil_uniq(prmtop%angl_wo_hy(4,i))

    end do

    enefunc%num_angles = nangl


    call get_loop_index(enefunc%num_angles, istart, iend)

    enefunc%istart_angle = istart
    enefunc%iend_angle   = iend

    return

  end subroutine setup_enefunc_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_dihe
  !> @brief        define DIHEDRAL term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_dihe(prmtop, enefunc)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, idih
    integer                  :: istart, iend


    enefunc%num_dihedrals = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) >= 0) then
        enefunc%num_dihedrals = enefunc%num_dihedrals + 1
      end if
    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) >= 0) then
        enefunc%num_dihedrals = enefunc%num_dihedrals + 1
      end if
    end do

    call alloc_enefunc(enefunc, EneFuncDihe, enefunc%num_dihedrals)

    idih = 0

    do i = 1, prmtop%num_diheh
      
      if (prmtop%dihe_inc_hy(4,i) >= 0) then

        idih = idih + 1
        enefunc%dihe_list(1,idih) =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
        enefunc%dihe_list(2,idih) =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
        enefunc%dihe_list(3,idih) = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
        enefunc%dihe_list(4,idih) =      prmtop%dihe_inc_hy(4,i)  / 3 + 1

        enefunc%dihe_force_const(idih) = &
             prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
        enefunc%dihe_periodicity(idih) = &
             prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i))
        enefunc%dihe_phase(idih)       = &
             prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))

        if (prmtop%lscee_scale_factor) then
          enefunc%dihe_scee(idih) = &
             1.0_wp / prmtop%scee_scale_fact(prmtop%dihe_inc_hy(5,i))
        else
          enefunc%dihe_scee(idih) = 1.0_wp / 1.2_wp
        end if

        if (prmtop%lscnb_scale_factor) then
          enefunc%dihe_scnb(idih) = &
             1.0_wp / prmtop%scnb_scale_fact(prmtop%dihe_inc_hy(5,i))
        else
          enefunc%dihe_scnb(idih) = 1.0_wp / 2.0_wp
        end if

      end if

    end do

    ! check dihe_wo_hy array
    !
    do i = 1, prmtop%num_mdihea

      if (prmtop%dihe_wo_hy(4,i) >= 0) then

        idih = idih + 1
        enefunc%dihe_list(1,idih) =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
        enefunc%dihe_list(2,idih) =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
        enefunc%dihe_list(3,idih) = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
        enefunc%dihe_list(4,idih) =      prmtop%dihe_wo_hy(4,i)  / 3 + 1

        enefunc%dihe_force_const(idih) = &
             prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
        enefunc%dihe_periodicity(idih) = &
             prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i))
        enefunc%dihe_phase(idih)       = &
             prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))

        if (prmtop%lscee_scale_factor) then
          enefunc%dihe_scee(idih) = &
             1.0_wp / prmtop%scee_scale_fact(prmtop%dihe_wo_hy(5,i))
        else
          enefunc%dihe_scee(idih) = 1.0_wp / 1.2_wp
        end if

        if (prmtop%lscnb_scale_factor) then
          enefunc%dihe_scnb(idih) = &
             1.0_wp / prmtop%scnb_scale_fact(prmtop%dihe_wo_hy(5,i))
        else
          enefunc%dihe_scnb(idih) = 1.0_wp / 2.0_wp
        end if

      end if

    end do


    call get_loop_index(enefunc%num_dihedrals, istart, iend)

    enefunc%istart_dihedral = istart
    enefunc%iend_dihedral   = iend

    return

  end subroutine setup_enefunc_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_impr
  !> @brief        define IMPROPER dihedral term in potential energy function
  !! @authors      NT
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_impr(prmtop, enefunc)

    ! formal arguments
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, iimp
    integer                  :: istart, iend


    enefunc%num_impropers = 0

    do i = 1, prmtop%num_diheh
      if (prmtop%dihe_inc_hy(4,i) < 0) then
        enefunc%num_impropers = enefunc%num_impropers + 1
      end if
    end do

    do i = 1, prmtop%num_mdihea
      if (prmtop%dihe_wo_hy(4,i) < 0) then
        enefunc%num_impropers = enefunc%num_impropers + 1
      end if
    end do

    call alloc_enefunc(enefunc, EneFuncImpr, enefunc%num_impropers)

    iimp = 0

    do i = 1, prmtop%num_diheh
      
      if (prmtop%dihe_inc_hy(4,i) < 0) then

        iimp = iimp + 1
        enefunc%impr_list(1,iimp) =      prmtop%dihe_inc_hy(1,i)  / 3 + 1
        enefunc%impr_list(2,iimp) =      prmtop%dihe_inc_hy(2,i)  / 3 + 1
        enefunc%impr_list(3,iimp) = iabs(prmtop%dihe_inc_hy(3,i)) / 3 + 1
        enefunc%impr_list(4,iimp) = iabs(prmtop%dihe_inc_hy(4,i)) / 3 + 1

        enefunc%impr_force_const(iimp) = &
             prmtop%dihe_fcons_uniq(prmtop%dihe_inc_hy(5,i))
        enefunc%impr_periodicity(iimp) = &
             prmtop%dihe_perio_uniq(prmtop%dihe_inc_hy(5,i))
        enefunc%impr_phase(iimp)       = &
             prmtop%dihe_phase_uniq(prmtop%dihe_inc_hy(5,i))

      end if

    end do

    ! check dihe_wo_hy array
    !
    do i = 1, prmtop%num_mdihea

      if (prmtop%dihe_wo_hy(4,i) < 0) then

        iimp = iimp + 1
        enefunc%impr_list(1,iimp) =      prmtop%dihe_wo_hy(1,i)  / 3 + 1
        enefunc%impr_list(2,iimp) =      prmtop%dihe_wo_hy(2,i)  / 3 + 1
        enefunc%impr_list(3,iimp) = iabs(prmtop%dihe_wo_hy(3,i)) / 3 + 1
        enefunc%impr_list(4,iimp) = iabs(prmtop%dihe_wo_hy(4,i)) / 3 + 1

        enefunc%impr_force_const(iimp) = &
             prmtop%dihe_fcons_uniq(prmtop%dihe_wo_hy(5,i))
        enefunc%impr_periodicity(iimp) = &
             prmtop%dihe_perio_uniq(prmtop%dihe_wo_hy(5,i))
        enefunc%impr_phase(iimp)       = &
             prmtop%dihe_phase_uniq(prmtop%dihe_wo_hy(5,i))

      end if

    end do


    call get_loop_index(enefunc%num_impropers, istart, iend)

    enefunc%istart_improper = istart
    enefunc%iend_improper   = iend

    return

  end subroutine setup_enefunc_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_enefunc_nonb
  !> @brief        define NON-BOND term in potential energy function
  !! @authors      NT
  !! @param[in]    ene_info : ENERGY section control parameters information
  !! @param[in]    prmtop   : AMBER parameter topology information
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_nonb(ene_info, prmtop, molecule, enefunc)

    ! formal arguments
    type(s_ene_info),        intent(in)    :: ene_info
    type(s_prmtop),          intent(in)    :: prmtop
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, j, k, nnonb, numi, numl, nb_par_idx


    nnonb = prmtop%num_types

    ELECOEF          = ELECOEF_AMBER

    call alloc_enefunc(enefunc, EneFuncNbon, nnonb)

    do i = 1, nnonb
      do j = 1, nnonb

        nb_par_idx = prmtop%nb_par_idx(nnonb*(i-1)+j)

        if (nb_par_idx < 0) then
          enefunc%nb14_lj12(i,j) = 0.0_wp
          enefunc%nb14_lj6 (i,j) = 0.0_wp
          enefunc%nonb_lj12(i,j) = 0.0_wp
          enefunc%nonb_lj6 (i,j) = 0.0_wp
        else
          enefunc%nb14_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          enefunc%nb14_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
          enefunc%nonb_lj12(i,j) = prmtop%lennarda(nb_par_idx)
          enefunc%nonb_lj6 (i,j) = prmtop%lennardb(nb_par_idx)
        end if

      end do
    end do

    enefunc%num_atom_cls = nnonb


    ! treatment for 1-2, 1-3, 1-4 interactions
    !
    if (.not.ene_info%table) &
      call count_nonb_excl(molecule, enefunc)

    return

  end subroutine setup_enefunc_nonb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    count_nonb_excl
  !> @brief        exclude 1-2, 1-3 interactions
  !! @authors      NT
  !! @param[in]    molecule : molecule information
  !! @param[inout] enefunc  : potential energy functions information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine count_nonb_excl(molecule, enefunc)

    ! formal arguments
    type(s_molecule),        intent(in)    :: molecule
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer                  :: i, k, k1, k2, k3, k4
    integer                  :: num_excl, num_excl_total
    integer                  :: num_nb14, num_nb14_total
    integer                  :: max_nb14_num
    integer                  :: alloc_stat, dealloc_stat
    logical                  :: duplicate

    integer,     allocatable :: nonb_excl_list(:,:)
    integer,     allocatable :: nb14_calc_list(:,:)
    real(wp),    allocatable :: nb14_qq_scale (:,:)
    real(wp),    allocatable :: nb14_lj_scale (:,:)


    ! allocate nonbonded exclusion list and 1-4 interaction list
    !
    allocate(enefunc%num_nonb_excl(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(nonb_excl_list(MAX_EXCL,molecule%num_atoms),stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(enefunc%num_nb14_calc(molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    allocate(nb14_calc_list(MAX_NB14,molecule%num_atoms), &
             nb14_qq_scale (MAX_NB14,molecule%num_atoms), &
             nb14_lj_scale (MAX_NB14,molecule%num_atoms), stat=alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    enefunc%num_nonb_excl(1:molecule%num_atoms) = 0
    enefunc%num_nb14_calc(1:molecule%num_atoms) = 0

    ! exclude 1-2 interaction
    !
    num_excl_total = 0
    do k = 1, molecule%num_bonds

      k1 = molecule%bond_list(1,k)
      k2 = molecule%bond_list(2,k)

      if (k1 < k2) then
        num_excl = enefunc%num_nonb_excl(k1) + 1
        enefunc%num_nonb_excl(k1) = num_excl
        nonb_excl_list(num_excl,k1) = k2
        num_excl_total = num_excl_total + 1
      else
        num_excl = enefunc%num_nonb_excl(k2) + 1
        enefunc%num_nonb_excl(k2) = num_excl
        nonb_excl_list(num_excl,k2) = k1
        num_excl_total = num_excl_total + 1
      end if

    end do

    ! exclude 1-3 interaction
    !
    do k = 1, molecule%num_angles

      k1 = molecule%angl_list(1,k)
      k3 = molecule%angl_list(3,k)

      if (k1 < k3) then

        num_excl = enefunc%num_nonb_excl(k1)
        duplicate = .false.

        do i = 1, enefunc%num_nonb_excl(k1)
          if (k3 == nonb_excl_list(i,k1)) duplicate = .true. 
        end do 

        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%num_nonb_excl(k1) = num_excl
          nonb_excl_list(num_excl,k1) = k3
          num_excl_total = num_excl_total + 1
        end if

      else

        num_excl = enefunc%num_nonb_excl(k3)
        duplicate = .false.

        do i = 1, enefunc%num_nonb_excl(k3)
          if (k1 == nonb_excl_list(i,k3)) duplicate = .true. 
        end do 

        if (.not. duplicate) then
          num_excl = num_excl + 1
          enefunc%num_nonb_excl(k3) = num_excl
          nonb_excl_list(num_excl,k3) = k1
          num_excl_total = num_excl_total + 1
        end if

      end if
    end do

    ! count 1-4 interaction
    !
    num_nb14_total = 0
    do k = 1, molecule%num_dihedrals

      k1 = molecule%dihe_list(1,k)
      k4 = molecule%dihe_list(4,k)

      if (k1 < k4) then

        num_nb14 = enefunc%num_nb14_calc(k1) 
        duplicate = .false. 

        do i = 1, enefunc%num_nonb_excl(k1)
          if (k4 == nonb_excl_list(i,k1)) duplicate = .true. 
        end do 

        do i = 1, num_nb14
          if (k4 == nb14_calc_list(i,k1)) duplicate = .true. 
        end do

        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k1) = num_nb14
          nb14_calc_list(num_nb14,k1) = k4
          nb14_qq_scale (num_nb14,k1) = enefunc%dihe_scee(k)
          nb14_lj_scale (num_nb14,k1) = enefunc%dihe_scnb(k)
          num_nb14_total = num_nb14_total + 1
        end if

      else

        num_nb14 = enefunc%num_nb14_calc(k4) 
        duplicate = .false. 

        do i = 1, enefunc%num_nonb_excl(k4)
          if (k1 == nonb_excl_list(i,k4)) duplicate = .true. 
        end do 

        do i = 1, num_nb14
          if (k1 == nb14_calc_list(i,k4)) duplicate = .true. 
        end do

        if (.not. duplicate) then
          num_nb14 = num_nb14 + 1
          enefunc%num_nb14_calc(k4) = num_nb14
          nb14_calc_list(num_nb14,k4) = k1
          nb14_qq_scale (num_nb14,k4) = enefunc%dihe_scee(k)
          nb14_lj_scale (num_nb14,k4) = enefunc%dihe_scnb(k)
          num_nb14_total = num_nb14_total + 1
        end if

      end if
    end do
    
    ! pack 2D-array into 1D-array
    !
    allocate(enefunc%nonb_excl_list(num_excl_total), stat = alloc_stat) 
    if (alloc_stat /=0) &
      call error_msg_alloc

    call pack_array_i(molecule%num_atoms, enefunc%num_nonb_excl,  &
                      nonb_excl_list, enefunc%nonb_excl_list)

    deallocate(nonb_excl_list, stat = dealloc_stat)


    max_nb14_num = max(1,maxval(enefunc%num_nb14_calc(1:molecule%num_atoms)))
    allocate(enefunc%nb14_calc_list(max_nb14_num,molecule%num_atoms), &
             enefunc%nb14_qq_scale (max_nb14_num,molecule%num_atoms), &
             enefunc%nb14_lj_scale (max_nb14_num,molecule%num_atoms), &
             stat = alloc_stat)
    if (alloc_stat /=0) &
      call error_msg_alloc

    enefunc%nb14_calc_list(1:max_nb14_num,1:molecule%num_atoms) = 0
    enefunc%nb14_qq_scale (1:max_nb14_num,1:molecule%num_atoms) = 0.0_wp
    enefunc%nb14_lj_scale (1:max_nb14_num,1:molecule%num_atoms) = 0.0_wp

    do i = 1, molecule%num_atoms
      enefunc%nb14_calc_list(1:enefunc%num_nb14_calc(i),i) = &
                       nb14_calc_list(1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_lj_scale (1:enefunc%num_nb14_calc(i),i) = &
                       nb14_qq_scale (1:enefunc%num_nb14_calc(i),i)
      enefunc%nb14_lj_scale (1:enefunc%num_nb14_calc(i),i) = &
                       nb14_lj_scale (1:enefunc%num_nb14_calc(i),i)
    end do


    deallocate(nb14_calc_list, &
               nb14_qq_scale,  &
               nb14_lj_scale, stat = dealloc_stat)

    return

  end subroutine count_nonb_excl

end module at_enefunc_amber_mod
