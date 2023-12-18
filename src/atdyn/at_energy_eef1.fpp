!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_energy_eef1_mod
!> @brief   calculate nonbonded energy with EEF1/IMM1
!! @authors Takaharu Mori (TM)
! 
!  (c) Copyright 2016 RIKEN. All rights reserved.
! 
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_energy_eef1_mod

  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use molecules_str_mod
  use timers_mod
  use mpi_parallel_mod
  use messages_mod
  use constants_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  integer,               save :: istart, iend
  real(wp), allocatable, save :: FF(:)
  real(wp), allocatable, save :: dmfac(:,:)
  real(wp), allocatable, save :: coord_min(:,:), coord_tmp(:,:)

  real(wp), parameter         :: EEF1_CUTOFF_DIST = 9.0_wp

  ! subroutines
  public   :: compute_energy_nonbond_eef1
  private  :: compute_energy_reference
  private  :: compute_energy_nonbond14_eef1
  private  :: compute_energy_nonbond15_eef1
  private  :: compute_energy_nonbond14_imm1
  private  :: compute_energy_nonbond15_imm1
  private  :: compute_ellipsoid_depth
  private  :: g

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_eef1
  !> @brief        Calculate nonbond energy in no boundary condition
  !! @authors      TM
  !! @note         Shift function is not used in the ELECT term here.
  !!               Therefore, the ELECT energy differs from
  !!               that in the case of implicit_solvent = NONE.
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    nonb_ene : flag for calculate nonbonded energy
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @param[inout] esolv    : solvation free energy of target systems
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond_eef1(enefunc, molecule, pairlist, nonb_ene,&
                                         coord, force, virial, eelec, evdw,    &
                                         esolv)

    ! formal arguments
    type(s_enefunc),         intent(in)    :: enefunc
    type(s_molecule),        intent(in)    :: molecule
    type(s_pairlist),        intent(in)    :: pairlist
    logical,                 intent(in)    :: nonb_ene
    real(wp),                intent(in)    :: coord(:,:)
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eelec
    real(wp),                intent(inout) :: evdw
    real(wp),                intent(inout) :: esolv

    return

  end subroutine compute_energy_nonbond_eef1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_reference
  !> @brief        calculate self energy term
  !  @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         T.Lazaridis & M.Karplus, Proteins, 35, 133−152 (1999)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_reference(enefunc, molecule, coord, force, esolv)

    ! formal argumens
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: esolv

    ! local variables
    integer                   :: i, natom, n
    integer                   :: istart, iend, isd
    real(wp)                  :: ihalf_thick, zt, fz, hr
    real(wp)                  :: force_i, ddg, tmp, x2y2, pore_rad
    real(wp)                  :: rt, fr, x2y2z2, x1, y1, z1, x, y, z
    real(wp)                  :: a, b, c, m, s, m1, m2, fact, depth
    real(wp)                  :: ia, ib, ic, x2, y2, z2
    real(wp)                  :: nx, ny, nz, nn, d, x0, y0, z0
    real(wp)                  :: gg, grad_theta, grad_phi, diff, prev_d
    real(wp)                  :: r, theta, phi, theta1, theta2, f1, f2, phi1, phi2
    real(wp)                  :: s1, s2, t, s1_1, s2_1, t_1
    real(wp)                  :: fact_x, fact_y, fact_z, fact_xy
    logical                   :: make_pore
    integer,          pointer :: atmcls(:)
    real(wp),         pointer :: gref(:,:)

    return

  end subroutine compute_energy_reference

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_eef1
  !> @brief        calculate nonbonded14 energy
  !  @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         T.Lazaridis & M.Karplus, Proteins, 35, 133−152 (1999)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_eef1(enefunc, molecule, coord, force, &
                                           virial, eelec, evdw, esolv)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw
    real(wp),                 intent(inout) :: esolv

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, inv_rij2
    real(wp)                  :: inv_rij, inv_r3, inv_r6, inv_r10, inv_r12
    real(wp)                  :: lj6, lj12, E14FAC
    real(wp)                  :: term_elec, shift, dshift, switch, dswitch
    real(wp)                  :: cutoff2, switchdist2
    real(wp)                  :: c1_switch, c2_switch, c4_switch, coef
    real(wp)                  :: c12, c6, term_lj12, term_lj6
    real(wp)                  :: factor_i, factor_j, alpha_4pi_i
    real(wp)                  :: inv_lambda_i, inv_lambda_j, x_i, x_j
    real(wp)                  :: work(1:3), force_local(1:3)
    real(wp)                  :: elec_i, esolv_tmp, rtmp(1:3), rvdw_i, vol_i
    integer                   :: i, j, k, l, natom, id, my_id
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: atmcls_i, atmcls_j

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: dielec_const
    real(wp),         pointer :: cutoff, switchdist
    real(wp),         pointer :: alpha_4pi(:,:)
    real(wp),         pointer :: rvdw(:), vol(:,:), inv_lambda(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb14_calc(:), nb14_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    return

  end subroutine compute_energy_nonbond14_eef1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond15_eef1
  !> @brief        calculate nonbonded energy with pairlist (NOBC) with EEF1
  !! @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         T.Lazaridis & M.Karplus, Proteins, 35, 133−152 (1999)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond15_eef1(enefunc, molecule, pairlist, coord, &
                                           force, virial, eelec, evdw, esolv)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw
    real(wp),                 intent(inout) :: esolv

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, inv_rij2
    real(wp)                  :: inv_rij, inv_r6, inv_r12, inv_r10, dr126
    real(wp)                  :: lj6, lj12
    real(wp)                  :: term_elec, shift, dshift, switch, dswitch
    real(wp)                  :: cutoff2, switchdist2
    real(wp)                  :: c1_switch, c2_switch, c4_switch, coef
    real(wp)                  :: inv_lambda_i, inv_lambda_j, x_i, x_j
    real(wp)                  :: work(1:3), force_local(1:3), esolv_tmp
    real(wp)                  :: alpha_4pi_i, factor_i, factor_j
    real(wp)                  :: rtmp(1:3), elec_i, rvdw_i, vol_i
    integer                   :: i, j, k, l, natom, id, my_id
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: atmcls_i, atmcls_j

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: dielec_const
    real(wp),         pointer :: cutoff, switchdist
    real(wp),         pointer :: alpha_4pi(:,:)
    real(wp),         pointer :: rvdw(:), vol(:,:), inv_lambda(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb15_calc(:,:), nb15_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    return

  end subroutine compute_energy_nonbond15_eef1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond14_imm1
  !> @brief        calculate nonbonded14 energy with IMM1
  !  @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         T.Lazaridis & M.Karplus, Proteins, 35, 133−152 (1999)
  !!               T.Lazaridis, Proteins, 52, 176-192 (2003)
  !!               A. Rahaman & T. Lazaridis, BBA, 1838, 1440-1447 (2014)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond14_imm1(enefunc, molecule, coord, force, &
                                           virial, eelec, evdw, esolv)

    ! formal arguments
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_molecule), target, intent(in)    :: molecule
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw
    real(wp),                 intent(inout) :: esolv

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2, inv_rij2
    real(wp)                  :: inv_rij, inv_r3, inv_r6, inv_r10, inv_r12
    real(wp)                  :: lj6, lj12, elec_i, E14FAC
    real(wp)                  :: term_elec, shift, dshift, switch, dswitch
    real(wp)                  :: cutoff2, switchdist2
    real(wp)                  :: c1_switch, c2_switch, c4_switch, coef
    real(wp)                  :: c12, c6, term_lj12, term_lj6
    real(wp)                  :: factor_i, factor_j, alpha_4pi_i
    real(wp)                  :: inv_lambda_i, inv_lambda_j, x_i, x_j
    real(wp)                  :: work(1:3), force_local(1:3), rtmp(1:3)
    real(wp)                  :: esolv_tmp, esolv_i, esolv_j
    real(wp)                  :: fij, z, ihalf_thick, a
    real(wp)                  :: zt_i, zt_j, FF_i, FF_j, ddg_i, ddg_j
    real(wp)                  :: gfree_i, gfree_j, rvdw_i, vol_i
    real(wp)                  :: dqfac, sqfifj
    integer                   :: i, j, k, l, n, natom
    integer                   :: my_id, id
    integer                   :: num_nb14, ini_nb14, fin_nb14
    integer                   :: atmcls_i, atmcls_j

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: dielec_const
    real(wp),         pointer :: cutoff, switchdist
    real(wp),         pointer :: alpha_4pi(:,:), gfree_t(:,:)
    real(wp),         pointer :: rvdw(:), vol(:,:), inv_lambda(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb14_calc(:), nb14_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    return

  end subroutine compute_energy_nonbond14_imm1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_nonbond15_imm1
  !> @brief        calculate nonbonded energy with pairlist (NOBC) with IMM1
  !! @authors      TM
  !! @param[in]    enefunc  : potential energy functions information
  !! @param[in]    molecule : molecule information
  !! @param[in]    pairlist : pairlist information
  !! @param[in]    coord    : coordinates of target systems
  !! @param[inout] force    : forces of target systems
  !! @param[inout] virial   : virial of target systems
  !! @param[inout] eelec    : electrostatic energy of target systems
  !! @param[inout] evdw     : lennard-jones energy of target systems
  !! @note         T.Lazaridis & M.Karplus, Proteins, 35, 133−152 (1999)
  !!               T.Lazaridis, Proteins, 52, 176-192 (2003)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_nonbond15_imm1(enefunc, molecule, pairlist, coord, &
                                           force, virial, eelec, evdw, esolv)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(in)    :: enefunc
    type(s_pairlist), target, intent(in)    :: pairlist
    real(wp),                 intent(in)    :: coord(:,:)
    real(wp),                 intent(inout) :: force(:,:)
    real(wp),                 intent(inout) :: virial(3,3)
    real(wp),                 intent(inout) :: eelec
    real(wp),                 intent(inout) :: evdw
    real(wp),                 intent(inout) :: esolv

    ! local variables
    real(wp)                  :: dij(1:3), rij, rij2
    real(wp)                  :: inv_rij, inv_r6, inv_r12, inv_r10, dr126
    real(wp)                  :: lj6, lj12, elec_i
    real(wp)                  :: term_elec, shift, dshift, switch, dswitch
    real(wp)                  :: cutoff2, switchdist2
    real(wp)                  :: c1_switch, c2_switch, c4_switch, coef
    real(wp)                  :: inv_lambda_i, inv_lambda_j, x_i, x_j
    real(wp)                  :: rtmp(1:3), force_local(1:3)
    real(wp)                  :: work(1:3), esolv_tmp, dqfac, sqfifj
    real(wp)                  :: factor_i, factor_j, alpha_4pi_i
    real(wp)                  :: fij, z, ihalf_thick, a
    real(wp)                  :: rvdw_i, vol_i, inv_rij2
    real(wp)                  :: gfree_i, gfree_j, zt_i, zt_j, ci, cj
    real(wp)                  :: esolv_i, esolv_j, FF_i, FF_j, ddg_i, ddg_j
    integer                   :: i, j, k, l, n, natom
    integer                   :: num_nb15, ini_nb15, fin_nb15
    integer                   :: id, my_id
    integer                   :: atmcls_i, atmcls_j

    real(wp),         pointer :: charge(:)
    real(wp),         pointer :: dielec_const
    real(wp),         pointer :: cutoff, switchdist
    real(wp),         pointer :: alpha_4pi(:,:), gfree_t(:,:)
    real(wp),         pointer :: rvdw(:), vol(:,:), inv_lambda(:,:)
    integer,          pointer :: atmcls(:)
    integer,          pointer :: num_nb15_calc(:,:), nb15_calc_list(:,:)

#ifdef OMP
    integer                   :: omp_get_thread_num
#endif

    return

  end subroutine compute_energy_nonbond15_imm1

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_distance
  !> @brief        calculate minimum distance between a point and ellipsoid
  !! @authors      TM
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_ellipsoid_depth(a, b, c, m1, m2, natom, coord)

    ! formal arguments
    real(wp), intent(in)  :: a, b, c, m1, m2
    integer,  intent(in)  :: natom
    real(wp), intent(in)  :: coord(:,:)

    ! local variables
    integer               :: isd, ih, ig, il, is, i
    integer               :: ncycle, icycle, nlen, ixx
    real(wp)              :: x0, y0, z0, x, y, z, r
    real(wp)              :: x1, y1, z1, x2, y2, z2, d1, d2, ini_f
    real(wp)              :: theta, phi, dtheta, dphi
    real(wp)              :: tmp_theta(0:2), tmp_phi(0:2), f(0:2)
    real(wp)              :: com_theta, com_phi, hlamd1, hlamd2
    real(wp)              :: exp_theta, exp_phi, exp_f
    real(wp)              :: red_theta, red_phi, red_f
    real(wp)              :: ref_theta, ref_phi, ref_f
    real(wp)              :: s1, s2, t
    logical               :: middle(0:2)
    integer,  parameter   :: nsteps    = 10000
    real(wp), parameter   :: tolerance = 0.00001_wp
    real(wp), parameter   :: lamd1     = 0.33_wp
    real(wp), parameter   :: lamd2     = 1.70_wp

    return

  end subroutine compute_ellipsoid_depth

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      g
  !> @brief        calculate distance between a point and ellipsoid
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function g(theta, phi, a, b, c, m1, m2, x, y, z, x0, y0, z0)

    ! formal arguments
    real(wp), intent(in)  :: theta, phi, a, b, c, m1, m2, x, y, z
    real(wp), intent(out) :: x0, y0, z0

    ! local variables
    real(wp)              :: g, f, r, sc, ss
    real(wp)              :: sin_theta, cos_theta, sin_phi, cos_phi

    return

  end function g

end module at_energy_eef1_mod
