!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_experiments_mod
!> @brief   experimental data
!! @authors Osamu Miyashita (OM), Takaharu Mori (TM)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
  
module at_experiments_mod

  use molecules_str_mod
  use at_restraints_str_mod
  use at_experiments_str_mod
  use at_enefunc_str_mod
  use fileio_sit_mod
  use fileio_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use string_mod
  use timers_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_exp_info
    logical                         :: emfit           = .false.
    integer                         :: emfit_type      = ExperimentsEmfit
    character(MaxFilename)          :: emfit_target    = ''
    real(wp)                        :: emfit_sigma     = 2.5_wp
    real(wp)                        :: emfit_tolerance = 0.001_wp
    integer                         :: emfit_period    = 1

    real(wp)                        :: emfit_pixel_size   = 1.0
    real(wp)                        :: emfit_roll_angle   = 0.0
    real(wp)                        :: emfit_tilt_angle   = 0.0
    real(wp)                        :: emfit_yaw_angle    = 0.0
    real(wp)                        :: emfit_shift_x      = 0.0
    real(wp)                        :: emfit_shift_y      = 0.0
    real(wp)                        :: emfit_shift_z      = 0.0
    real(wp)                        :: emfit_tilt_series_angle = 60.0
    real(wp)                        :: emfit_tilt_series_step = 1.0

  end type s_exp_info

  type(s_experiments), target, save :: experiments
  integer,                     save :: natom_local, emfit_icycle
  real(wp),                    save :: corrcoeff_save
  integer,        allocatable, save :: list_local(:)
  integer,        allocatable, save :: icount(:), domain_index(:)
  integer,        allocatable, save :: ig_min_local(:,:), ig_max_local(:,:)
  integer,        allocatable, save :: num_atoms_local(:)
  integer,        allocatable, save :: ig_lower(:,:), ig_upper(:,:)
  real(wp),       allocatable, save :: dot_exp_drhodx(:), dot_exp_drhody(:), dot_exp_drhodz(:)
  real(wp),       allocatable, save :: erfa(:,:), expa(:,:)
  real(wp),       allocatable, save :: derfa_x(:,:), derfa_y(:,:), derfa_z(:,:)
  real(wp),       allocatable, save :: dexpa_x(:,:), dexpa_y(:,:), dexpa_z(:,:)
  real(wp),       allocatable, save :: dot_sim_drhodx(:), dot_sim_drhody(:), dot_sim_drhodz(:)
  real(wp),       allocatable, save :: emfit_force(:,:)
  character(MaxFilename)          :: emfit_target_test    = ''

  ! subroutines
  public  :: show_ctrl_experiments
  public  :: read_ctrl_experiments
  public  :: setup_experiments
  private :: setup_experiments_emfit
  public  :: compute_energy_experimental_restraint
  private :: compute_energy_experimental_restraint_emfit
  private :: setup_experiments_emfit_img
  private :: compute_energy_experimental_restraint_emfit_img
  private :: Euler_angle2Matrix
  private :: read_header_spi
  private :: read_data_spi

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_experiments
  !> @brief        show EXPERIMENTS section usage
  !! @authors      TM
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "remd", "min", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_experiments(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'remd')

        write(MsgOut,'(A)') '[EXPERIMENTS]'
        write(MsgOut,'(A)') 'emfit             = YES / NO '
        write(MsgOut,'(A)') 'emfit_type        = VOLUME / IMAGE '
        write(MsgOut,'(A)') ' '
        write(MsgOut,'(A)') '# emfit_target    = sample.sit  # EM data file'
        write(MsgOut,'(A)') '# emfit_sigma     = 2.5         # resolution parameter of the simulated map'
        write(MsgOut,'(A)') '# emfit_tolerance = 0.001       # Tolerance for error'
        write(MsgOut,'(A)') '# emfit_period    = 1           # emfit force update period'
        write(MsgOut,'(A)') ' '
        write(MsgOut,'(A)') '# emfit_pixel_size    = 0.0  # Pixel Size in Angstrom'
        write(MsgOut,'(A)') '# emfit_roll_angle    = 0.0  # Euleur roll angle in degrees'
        write(MsgOut,'(A)') '# emfit_tilt_angle    = 0.0  # Euleur tilt angle in degrees'
        write(MsgOut,'(A)') '# emfit_yaw_angle     = 0.0  # Euleur yaw angle in degrees'
        write(MsgOut,'(A)') '# emfit_shift_x       = 0.0  # Shift in x direction in pixel'
        write(MsgOut,'(A)') '# emfit_shift_y       = 0.0  # Shift in y direction in pixel'
        write(MsgOut,'(A)') '# emfit_shift_z       = 0.0  # Shift in z direction in pixel'
        write(MsgOut,'(A)') '# emfit_tilt_series_angle       = 60.0  # Maximum tilt series angle in degrees'
        write(MsgOut,'(A)') '# emfit_tilt_series_step        = 1.0  # Step angle of the tilt series in degrees'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_experiments
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_ctrl_experiments
  !> @brief        read EXPERIMENTS section in the control file
  !! @authors      TM
  !! @param[in]    handle   : unit number of control files
  !! @param[out]   exp_info : EXPERIMENTS section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_experiments(handle, exp_info) 

    ! parameters
    character(*),            parameter     :: Section = 'Experiments'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_exp_info),        intent(inout) :: exp_info 

    ! local variables


    ! read parameters from control file
    ! 
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'emfit',           &
                                exp_info%emfit)
    call read_ctrlfile_type(handle, Section, 'emfit_type',         &
                               exp_info%emfit_type, ExperimentsTypes)

    call read_ctrlfile_string (handle, Section, 'emfit_target',    &
                              exp_info%emfit_target)

    call read_ctrlfile_real   (handle, Section, 'emfit_sigma',     &
                              exp_info%emfit_sigma)

    call read_ctrlfile_real   (handle, Section, 'emfit_tolerance', &
                              exp_info%emfit_tolerance)

    call read_ctrlfile_integer(handle, Section, 'emfit_period',    &
                              exp_info%emfit_period)

    call read_ctrlfile_real   (handle, Section, 'emfit_roll_angle', &
                            exp_info%emfit_roll_angle)

    call read_ctrlfile_real   (handle, Section, 'emfit_tilt_angle', &
                            exp_info%emfit_tilt_angle)

    call read_ctrlfile_real   (handle, Section, 'emfit_yaw_angle',	&
                            exp_info%emfit_yaw_angle)

    call read_ctrlfile_real   (handle, Section, 'emfit_shift_x',	&
                            exp_info%emfit_shift_x)

    call read_ctrlfile_real   (handle, Section, 'emfit_shift_y',	&
                            exp_info%emfit_shift_y)

    call read_ctrlfile_real   (handle, Section, 'emfit_shift_z',	&
                            exp_info%emfit_shift_z)

    call read_ctrlfile_real(handle, Section, 'emfit_pixel_size',		&
                            exp_info%emfit_pixel_size)

    call read_ctrlfile_real(handle, Section, 'emfit_tilt_series_angle',		&
                            exp_info%emfit_tilt_series_angle)

    call read_ctrlfile_real(handle, Section, 'emfit_tilt_series_step',		&
                            exp_info%emfit_tilt_series_step)


    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then

      if (exp_info%emfit) then
        if (exp_info%emfit_type == ExperimentsEmfit) then

          write(MsgOut,'(a)') 'Read_Ctrl_Experiments > Parameters for experimental data fitting'
          write(MsgOut,'(a20)') '  emfit           = VOLUME'
          write(MsgOut,'(a20,a)') '  emfit_target    = ', trim(exp_info%emfit_target)
          write(MsgOut,'(a20,F10.4,a20,F10.4)')                      &
                '  emfit_sigma     = ', exp_info%emfit_sigma,        &
                '  emfit_tolerance = ', exp_info%emfit_tolerance
          write(MsgOut,'(a20,I10)')                                  &
                '  emfit_period    = ', exp_info%emfit_period
          write(MsgOut,'(a)') ' '

        else if (exp_info%emfit_type == ExperimentsEmfitImg) then

          write(MsgOut,'(a)') 'Read_Ctrl_Experiments > Parameters for experimental data fitting of images'
          write(MsgOut,'(a20)') '  emfit           = IMAGE'
          write(MsgOut,'(a20,a)') '  emfit_target    = ', trim(exp_info%emfit_target)
          write(MsgOut,'(a20,F10.4,a20,F10.4)')                      &
                '  emfit_sigma     = ', exp_info%emfit_sigma,        &
                '  emfit_tolerance = ', exp_info%emfit_tolerance
          write(MsgOut,'(a20,I10, a20,F10.4)')                       &
                '  emfit_period    = ', exp_info%emfit_period, &
                '  emfit_yaw_angle = ', exp_info%emfit_yaw_angle

          write(MsgOut,'(a20,F10.4,a20,F10.4)')&
                '  emfit_roll_angle = ', exp_info%emfit_roll_angle, &
                '  emfit_tilt_angle = ', exp_info%emfit_tilt_angle
          write(MsgOut,'(a20,F10.4,a20,F10.4)')&      
                '  emfit_shift_x    = ', exp_info%emfit_shift_x,  &
                '  emfit_shift_y    = ', exp_info%emfit_shift_y
          write(MsgOut,'(a)') ' '
        end if
      end if

    end if

    return

  end subroutine read_ctrl_experiments

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_experiments
  !> @brief        setup experiments information
  !! @authors      TM
  !! @param[in]    exp_info    : EXPERIMENTS section control parameters
  !! @param[inout] experiments : experiments information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_experiments(exp_info, molecule, restraints, enefunc)

    ! formal arguments
    type(s_exp_info),        intent(in)    :: exp_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer    :: i
    logical    :: do_emfit_res

    experiments%do_emfit = .false.
    experiments%emfit_type = exp_info%emfit_type

    if (exp_info%emfit) then

      do_emfit_res = .false.
      do i = 1, restraints%nfunctions
        if (enefunc%restraint_kind(i) == RestraintsFuncEM) do_emfit_res = .true.
      end do

      if (.not. do_emfit_res) then
        call error_msg('Setup_Experiments> EM is not defined in [RESTRAINTS]')
      end if

      experiments%do_emfit = .true.

      if (exp_info%emfit_type == ExperimentsEmfit) then
        call setup_experiments_emfit(exp_info, molecule)
      else if (exp_info%emfit_type == ExperimentsEmfitImg) then
        call setup_experiments_emfit_img(exp_info, molecule)
      endif

    end if

    return

  end subroutine setup_experiments

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_experiments_emfit
  !> @brief        setup experiments information
  !! @authors      OM, TM
  !! @param[in]    exp_info    : EXPERIMENTS section control parameter
  !! @param[inout] experiments : experiments information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_experiments_emfit(exp_info, molecule)

    ! formal arguments
    type(s_exp_info),        intent(in)    :: exp_info
    type(s_molecule),        intent(in)    :: molecule

    ! local variables
    type(s_sit)             :: sit
    integer                 :: i, j, k, n(3), max_nxyz
    real(wp)                :: pi, r, y


    ! read sitfile
    !
    call input_sit(exp_info%emfit_target, sit)

    ! allocate array
    !
    n(1) = sit%nx
    n(2) = sit%ny
    n(3) = sit%nz
    max_nxyz = maxval(n)
    call alloc_experiments(experiments, ExperimentsEmfit, &
                           sit%nx, sit%ny, sit%nz,        &
                           molecule%num_atoms, max_nxyz)

    experiments%emfit%sigma        = exp_info%emfit_sigma
    experiments%emfit%tolerance    = exp_info%emfit_tolerance
    experiments%emfit%emfit_period = exp_info%emfit_period

    experiments%emfit%target_map(:,:,:) = sit%map_data(:,:,:)
    experiments%emfit%nx = sit%nx
    experiments%emfit%ny = sit%ny
    experiments%emfit%nz = sit%nz
    experiments%emfit%x0 = sit%x0
    experiments%emfit%y0 = sit%y0
    experiments%emfit%z0 = sit%z0
    experiments%emfit%dx = sit%dx
    experiments%emfit%dy = sit%dy
    experiments%emfit%dz = sit%dz

    ! need to calculate norm for force calculation
    !
    experiments%emfit%norm_exp = 0.0_wp
    do i = 0, sit%nx - 1
      do j = 0, sit%ny - 1
        do k = 0, sit%nz - 1
          experiments%emfit%norm_exp = experiments%emfit%norm_exp &
                                     + experiments%emfit%target_map(i,j,k)**2
        end do
      end do
    end do
    experiments%emfit%norm_exp = sqrt(experiments%emfit%norm_exp)

    ! boundary of voxel i is [bound_x(i) bound_x(i+1)]
    ! voxel from (0,0,0) to (nx-1, ny-1, nz-1)
    !
    do i = 0, sit%nx
      experiments%emfit%bound_x(i) = sit%x0 + sit%dx * (dble(i) - 0.5_wp)
    end do
    do i = 0, sit%ny
      experiments%emfit%bound_y(i) = sit%y0 + sit%dy * (dble(i) - 0.5_wp)
    end do
    do i = 0, sit%nz
      experiments%emfit%bound_z(i) = sit%z0 + sit%dz * (dble(i) - 0.5_wp)
    end do

    ! determine_cutoff
    !
    pi = acos(-1.0_wp)
    r  = 0.0_wp
    y  = 0.0_wp

    do while (1.0_wp - y > exp_info%emfit_tolerance)
      y = erf(sqrt(3.0_wp/2.0_wp)*r) - sqrt(6.0_wp/pi)*r*exp(-3.0_wp*r*r/2.0_wp)
      r = r + 0.01_wp
    end do

    experiments%emfit%n_grid_cut_x = ceiling(r*exp_info%emfit_sigma/sit%dx)
    experiments%emfit%n_grid_cut_y = ceiling(r*exp_info%emfit_sigma/sit%dy)
    experiments%emfit%n_grid_cut_z = ceiling(r*exp_info%emfit_sigma/sit%dz)

    call dealloc_sit(sit)

    ! for MPI parallelization
    !
    allocate(icount(nproc_city))
    allocate(ig_min_local(nproc_city,3))
    allocate(ig_max_local(nproc_city,3))
    allocate(domain_index(molecule%num_atoms))
    allocate(num_atoms_local(nproc_city))

    ! force update
    !
    allocate(list_local (  molecule%num_atoms))
    allocate(emfit_force(3,molecule%num_atoms))
    emfit_icycle = -1

    ! write summary
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Experiments_Emfit> Setup variables for EMFIT'
      write(MsgOut,'(A20,F10.3)') '  radius/sigma    = ', r
      write(MsgOut,'(A20,F10.3)') '  radius          = ', r*experiments%emfit%sigma
      write(MsgOut,'(A20,F10.3)') '  dx              = ', experiments%emfit%dx
      write(MsgOut,'(A20,F10.3)') '  dy              = ', experiments%emfit%dy
      write(MsgOut,'(A20,F10.3)') '  dz              = ', experiments%emfit%dz
      write(MsgOut,'(A,I10)') '  adjacent grids to calculate density along x = ', &
                              experiments%emfit%n_grid_cut_x
      write(MsgOut,'(A,I10)') '  adjacent grids to calculate density along y = ', &
                              experiments%emfit%n_grid_cut_y
      write(MsgOut,'(A,I10)') '  adjacent grids to calculate density along z = ', &
                              experiments%emfit%n_grid_cut_z
    end if

    return

  end subroutine setup_experiments_emfit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_experiments_emfit_img
  !> @brief        setup experiments information for image
  !! @authors      RV
  !! @param[in]    exp_info    : EXPERIMENTS section control parameter
  !! @param[in]    molecule : experiments information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_experiments_emfit_img(exp_info, molecule)
    ! formal arguments
    type(s_exp_info),        intent(in)    :: exp_info
    type(s_molecule),        intent(in)    :: molecule

    experiments%emfit_img%pixel_size  = exp_info%emfit_pixel_size 
    experiments%emfit_img%shift_x     = exp_info%emfit_shift_x    
    experiments%emfit_img%shift_y     = exp_info%emfit_shift_y   
    experiments%emfit_img%period      = exp_info%emfit_period
    experiments%emfit_img%sigma       = exp_info%emfit_sigma    

    emfit_icycle = -1
    
    call Euler_angle2Matrix(exp_info%emfit_roll_angle, &
        exp_info%emfit_tilt_angle , exp_info%emfit_yaw_angle, &
        experiments%emfit_img%rot_matrix)

    experiments%emfit_img%inv_rot_matrix = transpose(experiments%emfit_img%rot_matrix)
      
    ! distance in pixel where the gaussian is truncated
    experiments%emfit_img%cutoff = 	&
          ceiling((sqrt((-(2*(experiments%emfit_img%sigma**2)))&
          * log(exp_info%emfit_tolerance)))/ experiments%emfit_img%pixel_size)

    call read_header_spi(exp_info%emfit_target,&
          experiments%emfit_img%image_size)

    call alloc_experiments(experiments, ExperimentsEmfitImg, &
            experiments%emfit_img%image_size,molecule%num_atoms,0,0,0)

    call read_data_spi(exp_info%emfit_target, experiments%emfit_img%target_img)

    emfit_target_test = exp_info%emfit_target

  end subroutine setup_experiments_emfit_img

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_experimental_restraint
  !> @brief        calculate restraint energy from experimental data
  !! @authors      TM
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    inum    : pointer for restraint function
  !! @param[in]    const   : force constants
  !! @param[in]    ref     : reference value
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[out]   eexp    : restraint energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_experimental_restraint(enefunc, coord, inum, &
                                           calc_force, force, virial, eexp, cv)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eexp
    real(wp),                intent(inout) :: cv

    if (experiments%do_emfit) then

      if (experiments%emfit_type == ExperimentsEmfit) then 

        call compute_energy_experimental_restraint_emfit &
              (enefunc, coord, inum, calc_force, force, virial, eexp, cv)

      else if (experiments%emfit_type == ExperimentsEmfitImg) then
       
        call compute_energy_experimental_restraint_emfit_img &
              (enefunc, coord, inum, calc_force, force, virial, eexp, cv)
      endif
    end if

    return

  end subroutine compute_energy_experimental_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_experimental_restraint_emfit
  !> @brief        calculate restraint energy from experimental data
  !! @authors      OM, TM
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    inum    : pointer for restraint function
  !! @param[in]    const   : force constants
  !! @param[in]    ref     : reference value
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[out]   eexp    : restraint energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_experimental_restraint_emfit(enefunc, coord, inum, &
                                     calc_force, force, virial, e_emfit, cv)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: e_emfit
    real(wp),                intent(inout) :: cv

    ! local variables
    integer           :: i, j, k, n, m, natom, group_id, ifound
    integer           :: alloc_size, i1, j1, k1
    integer           :: igx_min, igy_min, igz_min, igx_max, igy_max, igz_max
    integer           :: sum_ig, ave_ig
    integer           :: idomain, axis_length(3), iaxis
    integer           :: itmp1, itmp2, icount1
    integer           :: id, omp_get_thread_num
    integer           :: max_natom_local
    real(wp)          :: yzfactor, zfactor, max_length
    real(wp)          :: drho_dx, drho_dy, drho_dz, before_allreduce(2), after_allreduce(2)
    real(wp)          :: dcc_dx, dcc_dy, dcc_dz
    real(wp)          :: norm_sim, dot_exp_sim
    real(wp)          :: f, g, d, coeff_rho, coeff_drho
    real(wp)          :: pi, vol, corrcoeff
    real(wp)          :: dot_exp_drhodx_tmp, dot_exp_drhody_tmp, dot_exp_drhodz_tmp
    real(wp)          :: dot_sim_drhodx_tmp, dot_sim_drhody_tmp, dot_sim_drhodz_tmp
    real(wp)          :: emfit_force_norm, weight
    real(wp)          :: inv_dx, inv_dy, inv_dz
    real(wp)          :: fact1, fact2, fact3
    logical           :: do_allocate
    integer,  pointer :: nx, ny, nz
    integer,  pointer :: n_grid_cut_x, n_grid_cut_y, n_grid_cut_z
    integer,  pointer :: ig(:,:)
    integer,  pointer :: numatoms(:), atom_id(:,:)
    real(wp), pointer :: norm_exp            ! is calculated in setup_emfit just once
    real(wp), pointer :: sigma
    real(wp), pointer :: simulated_map(:,:,:), target_map(:,:,:)
    real(wp), pointer :: bound_x(:), bound_y(:), bound_z(:)


!    call timer(TimerEmfit, TimerOn)

    sigma          => experiments%emfit%sigma
    nx             => experiments%emfit%nx
    ny             => experiments%emfit%ny
    nz             => experiments%emfit%nz
    ig             => experiments%emfit%ig
    bound_x        => experiments%emfit%bound_x
    bound_y        => experiments%emfit%bound_y
    bound_z        => experiments%emfit%bound_z
    norm_exp       => experiments%emfit%norm_exp
    simulated_map  => experiments%emfit%simulated_map
    target_map     => experiments%emfit%target_map
    n_grid_cut_x   => experiments%emfit%n_grid_cut_x
    n_grid_cut_y   => experiments%emfit%n_grid_cut_y
    n_grid_cut_z   => experiments%emfit%n_grid_cut_z
    inv_dx         =  1.0_wp/experiments%emfit%dx
    inv_dy         =  1.0_wp/experiments%emfit%dy
    inv_dz         =  1.0_wp/experiments%emfit%dz

    numatoms       => enefunc%restraint_numatoms
    atom_id        => enefunc%restraint_atomlist
    group_id       =  enefunc%restraint_grouplist(1,inum)
    weight         =  enefunc%restraint_const(1,inum)

    natom = numatoms(group_id)
    pi    = acos(-1.0_wp)
    vol   = experiments%emfit%dx * experiments%emfit%dy * experiments%emfit%dz


    ! perform emfit or not
    !
    emfit_icycle = emfit_icycle + 1
    if (experiments%emfit%emfit_period /= 0) then
      if (mod(emfit_icycle,experiments%emfit%emfit_period) /= 0) then
        do m = 1, natom_local
          n = list_local(m)
          force(1:3,n) = force(1:3,n) + emfit_force(1:3,n)
        end do
        corrcoeff = corrcoeff_save
        experiments%emfit%corrcoeff = corrcoeff
        e_emfit = weight * (1.0_wp - corrcoeff)
        cv = corrcoeff
        return
      end if
    else
      return
    end if

    ! coefficient for drho/d[xyz]
    coeff_drho = - sigma**2 * pi / 6.0_wp / vol

    ! coefficient for rho
    !   - sign is missing in drho/dq equstions of BJ 2012
    coeff_rho = (sigma**2 * pi / 6.0_wp)**(3.0_wp/2.0_wp) / vol

    ! calculate which grid each atom is in
    ! and the range of the grid used for calculation
    !
    igx_min = nx
    igy_min = ny
    igz_min = nz
    igx_max = 0
    igy_max = 0
    igz_max = 0

    !$omp parallel do                            &
    !$omp private(m, n)                          &
    !$omp reduction(max:igx_max,igy_max,igz_max) &
    !$omp reduction(min:igx_min,igy_min,igz_min)
    !
    do m = 1, natom

      ! calculate the index of the grid that each atom is located
      !
      n = atom_id(m,group_id)

      ig(1,n) = int((coord(1,n) - experiments%emfit%x0)*inv_dx + 0.5_wp)
      ig(2,n) = int((coord(2,n) - experiments%emfit%y0)*inv_dy + 0.5_wp)
      ig(3,n) = int((coord(3,n) - experiments%emfit%z0)*inv_dz + 0.5_wp)

      if (ig(1,n) > igx_max) igx_max = ig(1,n)
      if (ig(2,n) > igy_max) igy_max = ig(2,n)
      if (ig(3,n) > igz_max) igz_max = ig(3,n)

      if (ig(1,n) < igx_min) igx_min = ig(1,n)
      if (ig(2,n) < igy_min) igy_min = ig(2,n)
      if (ig(3,n) < igz_min) igz_min = ig(3,n)

      ! check if the calculation (grids around the atom) doesn't go outside the box
      ! comparison like igx(n) + n_grid_cut_x >= nx doesn't work for integer
      ! overflow, when igx is max integer
      !
      if (ig(1,n) < n_grid_cut_x .or. ig(1,n) >= nx - n_grid_cut_x .or. &
          ig(2,n) < n_grid_cut_y .or. ig(2,n) >= ny - n_grid_cut_y .or. &
          ig(3,n) < n_grid_cut_z .or. ig(3,n) >= nz - n_grid_cut_z ) then

        write(MsgOut,'(A,I10)') 'Gaussian kernes is extending outside the map box for atom ', n
        write(MsgOut,'(A,3F12.3)') 'Coordinate: ', coord(1,n), coord(2,n), coord(3,n)
        write(MsgOut,*) 'grid index igx(n), igy(n), igz(n): ', ig(1,n), ig(2,n), ig(3,n)
        write(MsgOut,*) 'grid index range, xmin(n), xmax(n), ymin(n), ymax(n), zmin(n), zmax(n) ', &
          int(ig(1,n),8)-n_grid_cut_x, int(ig(1,n),8)+n_grid_cut_x, int(ig(2,n),8)-n_grid_cut_y,   &
          int(ig(2,n),8)+n_grid_cut_y, int(ig(3,n),8)-n_grid_cut_z, int(ig(3,n),8)+n_grid_cut_z

        call error_msg('Compute_Energy_Experimental_Restraint_Emfit> Gaussian kernel is extending outside the map box (see "Chapter: Trouble shooting" in the user manual)')

      end if

    end do
    !$omp end parallel do


    ! min/max grid index that needs to be calculated for the whole system
    !
    igx_min = igx_min - n_grid_cut_x
    igy_min = igy_min - n_grid_cut_y
    igz_min = igz_min - n_grid_cut_z
    igx_max = igx_max + n_grid_cut_x
    igy_max = igy_max + n_grid_cut_y
    igz_max = igz_max + n_grid_cut_z


    ! make partitions for domain decomposition
    !
    icount(1)         = natom
    domain_index(:)   = 1
    ig_max_local(1,1) = igx_max
    ig_max_local(1,2) = igy_max
    ig_max_local(1,3) = igz_max
    ig_min_local(1,1) = igx_min
    ig_min_local(1,2) = igy_min
    ig_min_local(1,3) = igz_min

    do i = 1, nproc_city - 1

      ! select one domain which constains largest number of particle
      !
      idomain = maxloc(icount(1:i),dim=1)

      ! select partitioning axis (longest axis in the domain)
      !
      ig_min_local(i+1,1:3)  = ig_min_local(idomain,1:3)
      ig_max_local(i+1,1:3)  = ig_max_local(idomain,1:3)

      axis_length(1:3) = ig_max_local(idomain,1:3) - ig_min_local(idomain,1:3)
      iaxis = maxloc(axis_length(1:3),dim=1)

      ! calculate the averaged coordinates of particles in the selected axis
      ! and make partition in the selected axis
      !
      sum_ig = 0

      !$omp parallel do   &
      !$omp private(m, n) &
      !$omp reduction(+:sum_ig)
      !
      do m = 1, natom
        n = atom_id(m,group_id)
        if (domain_index(n) == idomain) then
          sum_ig = sum_ig + ig(iaxis,n)
        end if
      end do
      !$omp end parallel do

      ave_ig = aint(sum_ig/dble(icount(idomain)))
      ig_max_local(idomain,iaxis) = ave_ig
      ig_min_local(i + 1,  iaxis) = ave_ig + 1

      ! update domain assignment information
      !
      icount1 = 0

      !$omp parallel do          &
      !$omp private (m, n)       &
      !$omp reduction(+:icount1)
      !
      do m = 1, natom
        n = atom_id(m,group_id)
        if (domain_index(n) == idomain) then
          if (ave_ig < ig(iaxis,n)) then
            domain_index(n) = i + 1
            icount1 = icount1 + 1
          end if
        end if
      end do
      !$omp end parallel do

      icount(i+1)     = icount1
      icount(idomain) = icount(idomain) - icount(i+1)

    end do

    igx_max = ig_max_local(my_city_rank+1,1)
    igy_max = ig_max_local(my_city_rank+1,2)
    igz_max = ig_max_local(my_city_rank+1,3)
    igx_min = ig_min_local(my_city_rank+1,1)
    igy_min = ig_min_local(my_city_rank+1,2)
    igz_min = ig_min_local(my_city_rank+1,3)


    ! make atom list in my domain
    !
    ifound = 0
    do m = 1, natom
      n = atom_id(m,group_id)
      if (ig(1,n) < igx_min - n_grid_cut_x) cycle
      if (ig(1,n) > igx_max + n_grid_cut_x) cycle
      if (ig(2,n) < igy_min - n_grid_cut_y) cycle
      if (ig(2,n) > igy_max + n_grid_cut_y) cycle
      if (ig(3,n) < igz_min - n_grid_cut_z) cycle
      if (ig(3,n) > igz_max + n_grid_cut_z) cycle
      ifound = ifound + 1
      list_local(ifound) = n
    end do
    natom_local = ifound

!    call timer(TimerEmfit, TimerOff)

    ! MPI allgather ifound
    !
#ifdef MPI
    call mpi_allgather(ifound,          1, mpi_integer, &
                       num_atoms_local, 1, mpi_integer, mpi_comm_city, ierror)
#endif
    max_natom_local = maxval(num_atoms_local(:))


!    call timer(TimerEmfit, TimerOn)

    ! allocate memory by ifound
    !
    if (allocated(dot_exp_drhodx)) then
      do_allocate = .false.
      if (max_natom_local > size(dot_exp_drhodx(:))) then
        do_allocate = .true.
        deallocate(dot_exp_drhodx, dot_exp_drhody, dot_exp_drhodz, &
                   dot_sim_drhodx, dot_sim_drhody, dot_sim_drhodz, &
                   ig_upper, ig_lower, erfa, expa,                 &
                   derfa_x, derfa_y, derfa_z, dexpa_x, dexpa_y, dexpa_z)
      end if
    else
      do_allocate = .true.
    end if

    if (do_allocate) then
      alloc_size = int(1.05_wp*max_natom_local)
      allocate(dot_exp_drhodx(    1:alloc_size))
      allocate(dot_exp_drhody(    1:alloc_size))
      allocate(dot_exp_drhodz(    1:alloc_size))
      allocate(dot_sim_drhodx(    1:alloc_size))
      allocate(dot_sim_drhody(    1:alloc_size))
      allocate(dot_sim_drhodz(    1:alloc_size))
      allocate(ig_upper      (1:3,1:alloc_size))
      allocate(ig_lower      (1:3,1:alloc_size))
      allocate(erfa    (-n_grid_cut_x:n_grid_cut_x+1,1:alloc_size))
      allocate(expa    (-n_grid_cut_x:n_grid_cut_x+1,1:alloc_size))
      allocate(derfa_x (-n_grid_cut_x:n_grid_cut_x+1,1:alloc_size))
      allocate(derfa_y (-n_grid_cut_y:n_grid_cut_y+1,1:alloc_size))
      allocate(derfa_z (-n_grid_cut_z:n_grid_cut_z+1,1:alloc_size))
      allocate(dexpa_x (-n_grid_cut_x:n_grid_cut_x+1,1:alloc_size))
      allocate(dexpa_y (-n_grid_cut_y:n_grid_cut_y+1,1:alloc_size))
      allocate(dexpa_z (-n_grid_cut_z:n_grid_cut_z+1,1:alloc_size))
    end if


    !$omp parallel do     &
    !$omp private (m, n)
    !
    do m = 1, ifound
      n = list_local(m)
      ig_lower(1,m) = ig(1,n) - n_grid_cut_x
      ig_lower(2,m) = ig(2,n) - n_grid_cut_y
      ig_lower(3,m) = ig(3,n) - n_grid_cut_z
      ig_upper(1,m) = ig(1,n) + n_grid_cut_x
      ig_upper(2,m) = ig(2,n) + n_grid_cut_y
      ig_upper(3,m) = ig(3,n) + n_grid_cut_z

      if (igx_max < ig_upper(1,m)) ig_upper(1,m) = igx_max
      if (igy_max < ig_upper(2,m)) ig_upper(2,m) = igy_max
      if (igz_max < ig_upper(3,m)) ig_upper(3,m) = igz_max
      if (igx_min > ig_lower(1,m)) ig_lower(1,m) = igx_min
      if (igy_min > ig_lower(2,m)) ig_lower(2,m) = igy_min
      if (igz_min > ig_lower(3,m)) ig_lower(3,m) = igz_min
    end do
    !$omp end parallel do


    ! calculate error and exp functions and differences for integral and save
    ! only for the grid within the cutoff
    !   f: coefficient for x for erf calculation
    !   g: coefficient for y for exp calculation
    !
    f = sqrt(3.0_wp / (2.0_wp * sigma**2))
    g = -3.0_wp / (2.0_wp * sigma**2)

    !$omp parallel do         &
    !$omp private(j, m, n, d)
    !
    do m = 1, ifound

      n = list_local(m)

      do j = -n_grid_cut_x, n_grid_cut_x+1
        d = bound_x(ig(1,n)+j) - coord(1,n)
        erfa(j,m) = erf(f*d)
        expa(j,m) = exp(g*d*d)
      end do
      do j = -n_grid_cut_x, n_grid_cut_x
        derfa_x(j,m) = erfa(j+1,m) - erfa(j,m)
        dexpa_x(j,m) = expa(j+1,m) - expa(j,m)
      end do

      do j = -n_grid_cut_y, n_grid_cut_y+1
        d = bound_y(ig(2,n)+j) - coord(2,n)
        erfa(j,m) = erf(f*d)
        expa(j,m) = exp(g*d*d)
      end do
      do j = -n_grid_cut_y, n_grid_cut_y
        derfa_y(j,m) = erfa(j+1,m) - erfa(j,m)
        dexpa_y(j,m) = expa(j+1,m) - expa(j,m)
      end do

      do j = -n_grid_cut_z, n_grid_cut_z+1
        d = bound_z(ig(3,n)+j) - coord(3,n)
        erfa(j,m) = erf(f*d)
        expa(j,m) = exp(g*d*d)
      end do
      do j = -n_grid_cut_z, n_grid_cut_z
        derfa_z(j,m) = erfa(j+1,m) - erfa(j,m)
        dexpa_z(j,m) = expa(j+1,m) - expa(j,m)
      end do

    end do
    !$omp end parallel do


    ! calculate simulated density map
    !
    simulated_map(igx_min:igx_max, igy_min:igy_max, igz_min:igz_max) = 0.0_wp

    !$omp parallel &
    !$omp private(id, m, n, k, j, i, zfactor, yzfactor, i1, j1, k1)
#ifdef OMP
    id = omp_get_thread_num()
#else
    id = 0
#endif
    do m = 1, ifound
      n = list_local(m)
      do k = id+ig_lower(3,m), ig_upper(3,m), nthread
        k1 = k - ig(3,n)
        zfactor = coeff_rho*derfa_z(k1,m)
        do j = ig_lower(2,m), ig_upper(2,m)
          j1 = j - ig(2,n)
          yzfactor = derfa_y(j1,m)*zfactor
          do i = ig_lower(1,m), ig_upper(1,m)
            i1 = i - ig(1,n)
            simulated_map(i,j,k) = simulated_map(i,j,k) + derfa_x(i1,m)*yzfactor
          end do
        end do
      end do
      !$omp barrier
    end do
    !$omp end parallel


    ! calculate two nominators for each atom first
    !
    dot_exp_drhodx(:) = 0.0_wp
    dot_exp_drhody(:) = 0.0_wp
    dot_exp_drhodz(:) = 0.0_wp

    dot_sim_drhodx(:) = 0.0_wp
    dot_sim_drhody(:) = 0.0_wp
    dot_sim_drhodz(:) = 0.0_wp

    !$omp parallel do schedule (dynamic)                                         &
    !$omp private(i, j, k, m, n, drho_dx, drho_dy, drho_dz, fact1, fact2, fact3, &
    !$omp         dot_exp_drhodx_tmp, dot_exp_drhody_tmp, dot_exp_drhodz_tmp,    &
    !$omp         dot_sim_drhodx_tmp, dot_sim_drhody_tmp, dot_sim_drhodz_tmp,    &
    !$omp         i1, j1, k1)
    !
    do m = 1, ifound

      n = list_local(m)

      dot_exp_drhodx_tmp = 0.0_wp
      dot_exp_drhody_tmp = 0.0_wp
      dot_exp_drhodz_tmp = 0.0_wp

      dot_sim_drhodx_tmp = 0.0_wp
      dot_sim_drhody_tmp = 0.0_wp
      dot_sim_drhodz_tmp = 0.0_wp

      do k = ig_lower(3,m), ig_upper(3,m)
        k1 = k - ig(3,n)
        do j = ig_lower(2,m), ig_upper(2,m)
          j1 = j - ig(2,n)

          fact1 = derfa_y(j1,m) * derfa_z(k1,m)
          fact2 = dexpa_y(j1,m) * derfa_z(k1,m)
          fact3 = dexpa_z(k1,m) * derfa_y(j1,m)

          do i = ig_lower(1,m), ig_upper(1,m)
            i1 = i - ig(1,n)

            drho_dx = dexpa_x(i1,m) * fact1
            drho_dy = derfa_x(i1,m) * fact2
            drho_dz = derfa_x(i1,m) * fact3

            dot_exp_drhodx_tmp = dot_exp_drhodx_tmp + target_map(i,j,k) * drho_dx
            dot_exp_drhody_tmp = dot_exp_drhody_tmp + target_map(i,j,k) * drho_dy
            dot_exp_drhodz_tmp = dot_exp_drhodz_tmp + target_map(i,j,k) * drho_dz

            dot_sim_drhodx_tmp = dot_sim_drhodx_tmp + simulated_map(i,j,k) * drho_dx
            dot_sim_drhody_tmp = dot_sim_drhody_tmp + simulated_map(i,j,k) * drho_dy
            dot_sim_drhodz_tmp = dot_sim_drhodz_tmp + simulated_map(i,j,k) * drho_dz

          end do
        end do
      end do

      dot_exp_drhodx(m) = dot_exp_drhodx_tmp * coeff_drho
      dot_exp_drhody(m) = dot_exp_drhody_tmp * coeff_drho
      dot_exp_drhodz(m) = dot_exp_drhodz_tmp * coeff_drho

      dot_sim_drhodx(m) = dot_sim_drhodx_tmp * coeff_drho
      dot_sim_drhody(m) = dot_sim_drhody_tmp * coeff_drho
      dot_sim_drhodz(m) = dot_sim_drhodz_tmp * coeff_drho

    end do
    !$omp end parallel do


    norm_sim = 0.0_wp
    dot_exp_sim = 0.0_wp

    !$omp parallel do                                    &
    !$omp private(i, j, k)                               &
    !$omp reduction(+:norm_sim) reduction(+:dot_exp_sim)
    !
    do k = igz_min, igz_max
      do j = igy_min, igy_max
        do i = igx_min, igx_max
          norm_sim = norm_sim + simulated_map(i,j,k)**2
          dot_exp_sim = dot_exp_sim + target_map(i,j,k) * simulated_map(i,j,k)
        end do
      end do
    end do
    !$omp end parallel do

!    call timer(TimerEmfit, TimerOff)

#ifdef MPI
    before_allreduce(1) = norm_sim
    before_allreduce(2) = dot_exp_sim

    call mpi_allreduce(before_allreduce, after_allreduce, 2,  &
                       mpi_wp_real,  mpi_sum,                 &
                       mpi_comm_city, ierror)

    norm_sim    = after_allreduce(1)
    dot_exp_sim = after_allreduce(2)
#endif

!    call timer(TimerEmfit, TimerOn)

    norm_sim  = sqrt(norm_sim)
    corrcoeff = dot_exp_sim / (norm_exp*norm_sim)
    experiments%emfit%corrcoeff = corrcoeff
    e_emfit = weight * (1.0_wp - corrcoeff)
    cv = corrcoeff
    corrcoeff_save = corrcoeff

    if (calc_force) then

      emfit_force(:,:) = 0.0_wp
      fact1 = 1.0_wp / norm_exp
      fact2 = 1.0_wp / norm_sim
      fact3 = dot_exp_sim/norm_sim**3

      !$omp parallel do                           &
      !$omp private(m, n, dcc_dx, dcc_dy, dcc_dz) &
      !$omp reduction(+:force)
      !
      do m = 1, ifound
        n = list_local(m)

        dcc_dx = fact1*(dot_exp_drhodx(m)*fact2 - dot_sim_drhodx(m)*fact3)
        dcc_dy = fact1*(dot_exp_drhody(m)*fact2 - dot_sim_drhody(m)*fact3)
        dcc_dz = fact1*(dot_exp_drhodz(m)*fact2 - dot_sim_drhodz(m)*fact3)

        emfit_force(1,n) = weight*dcc_dx
        emfit_force(2,n) = weight*dcc_dy
        emfit_force(3,n) = weight*dcc_dz

        force(1,n) = force(1,n) + emfit_force(1,n)
        force(2,n) = force(2,n) + emfit_force(2,n)
        force(3,n) = force(3,n) + emfit_force(3,n)
      end do
      !$omp end parallel do

    end if

!    call timer(TimerEmfit, TimerOff)

    return

  end subroutine compute_energy_experimental_restraint_emfit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_experimental_restraint_emfit_img
  !> @brief        calculate restraint energy from experimental data, the data includes two set of images
  !! @authors      Rémi Vuillemot
  !! @param[in]    enefunc : potential energy functions information: this potential engergy function is extracted from CHARMM forcefield
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    inum    : pointer for restraint function
  !! @param[in]    const   : force constants
  !! @param[in]    ref     : reference value
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[out]   eexp    : restraint energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_experimental_restraint_emfit_img(enefunc, coord, inum, &
    calc_force, force, virial, eexp, cv)
    ! formal arguments
    type(s_enefunc), target, intent(inout) 	:: enefunc
    real(wp),                intent(in)    	:: coord(:,:)
    integer,                 intent(in)    	:: inum
    logical,                 intent(in)    	:: calc_force
    real(wp),                intent(inout) 	:: force(:,:)
    real(wp),                intent(inout) 	:: virial(3,3)
    real(wp),                intent(inout) 	:: eexp
    real(wp),                intent(inout) 	:: cv
  
    integer                   :: i, j, a, n_atoms, n_pix, group_id, &
                                  n, n_atoms_group, cutoff, image_size
    real(wp) :: roll_angle, tilt_angle, yaw_angle, shift_x, shift_y, pixel_size
    real(wp) :: gaussian, sigma, norm, threshold
    real(wp) :: sum_sim2, sum_exp2, sum_simpexp, force_constant	
    real(wp) :: const1, const2		
    real(wp) :: mu(2), dpsim(2),  dcc(2)
    real(wp) :: rot_matrix(4,4),  inv_rot_matrix(4,4)
    integer,  pointer ::numatoms(:), atom_id(:,:)
    real(wp), pointer :: sim_image(:,:), exp_image(:,:)
    real(wp), pointer :: rot_coord(:,:), gaussians_saved(:,:,:), emfit_img_force(:,:)
    integer,pointer	:: pixels(:,:)
    character(MaxFilename) :: outfile 
  
    n_atoms = size(coord,2)
    atom_id        => enefunc%restraint_atomlist
    group_id       =  enefunc%restraint_grouplist(1,inum)
    numatoms       => enefunc%restraint_numatoms
    n_atoms_group = numatoms(group_id)

    image_size = experiments%emfit_img%image_size
    pixel_size = experiments%emfit_img%pixel_size
    roll_angle = experiments%emfit_img%roll_angle
    tilt_angle = experiments%emfit_img%tilt_angle
    yaw_angle  = experiments%emfit_img%yaw_angle 
    shift_x    = experiments%emfit_img%shift_x   
    shift_y    = experiments%emfit_img%shift_y   
    sigma      = experiments%emfit_img%sigma
    cutoff     = experiments%emfit_img%cutoff
    exp_image  => experiments%emfit_img%target_img 
    sim_image  => experiments%emfit_img%simulated_img  
    rot_coord  => experiments%emfit_img%rot_coord   
    pixels     => experiments%emfit_img%pixels     
    emfit_img_force  => experiments%emfit_img%emfit_img_force       
    gaussians_saved  => experiments%emfit_img%gaussians_saved


    ! ------------------------------------------------------------------------
    ! PERFORM EMFIT OR NOT
    ! ------------------------------------------------------------------------
    emfit_icycle = emfit_icycle + 1
    if (experiments%emfit_img%period /= 0) then
      if (mod(emfit_icycle,experiments%emfit_img%period) /= 0) then
        do a=1, n_atoms_group
          n=  atom_id(a,group_id)
          force(1:3,n) = force(1:3,n) + emfit_img_force(1:3,n)
        end do
        cv= corrcoeff_save
        force_constant = enefunc%restraint_const(1,inum)
        eexp = force_constant * (1.0_wp - cv)
        return
      end if
    else
      return
    end if

    ! ------------------------------------------------------------------------
    ! ROTATE AND SHIFT PDB
    ! ------------------------------------------------------------------------
    do a=1, n_atoms_group
      n=  atom_id(a,group_id)
      rot_coord(1:3,a)=matmul(experiments%emfit_img%rot_matrix(1:3,1:3),coord(1:3,n))
      rot_coord(1,a) = rot_coord(1,a) - (shift_x * pixel_size)
      rot_coord(2,a) = rot_coord(2,a) - (shift_y * pixel_size)
    end do
    
    ! ------------------------------------------------------------------------
    ! SELECT PIXELS TO INTEGRATE
    ! ------------------------------------------------------------------------
    ! number of pixels considered per atom
    n_pix = int(cutoff*2)  	

    do a = 1, n_atoms_group
      pixels(1,a) = 1+nint(rot_coord(1,a)/pixel_size + (image_size/2))-cutoff ! starting pixel
      pixels(2,a) = 1+nint(rot_coord(2,a)/pixel_size + (image_size/2))-cutoff

      ! Check if the atom is extending outside the box
      if (((pixels(1,a) <= 0) .or. (pixels(2,a) <= 0)) .or. &
      (((pixels(1,a)+n_pix) > image_size) .or. ((pixels(2,a)+n_pix) > image_size))) then
        write(MsgOut,'(A)') 'Gaussian kernel is extending outside the map box for atom '
        write(MsgOut,'(A,I10)') '   pix1 : ', pixels(1,a)
        write(MsgOut,'(A,I10)') '   pix2 : ', pixels(2,a)
        call error_msg('Compute_Energy_Experimental_Restraint_Emfit_Img> Gaussian kernel is extending outside the map box ')
      endif
    end do
  
    ! ------------------------------------------------------------------------
    ! GENERATE SIM IMAGE
    ! ------------------------------------------------------------------------
    sim_image(:,:) = 0.0_wp
  
    !$omp parallel do default(none)                                 &
    !$omp private(a,i,j , mu, norm, gaussian)                   &
    !$omp shared(pixels, n_pix,image_size, pixel_size, &
    !$omp 	rot_coord, sigma, sim_image, gaussians_saved, atom_id, group_id, n_atoms_group)
    !
  
    do a=1, n_atoms_group
      do j=pixels(2,a), pixels(2,a) + n_pix
        do i=pixels(1,a), pixels(1,a) + n_pix

          mu = (/i,j/)
          mu = (mu -1- (image_size / 2)) * pixel_size
  
          norm = (rot_coord(1,a)-mu(1))**2 + (rot_coord(2,a)-mu(2))**2
          gaussian = exp(-norm / (2 * (sigma ** 2)))
          
          sim_image(i,j) = sim_image(i,j) + gaussian
  
          ! save the gaussian to avoid computing them for the gradient
          gaussians_saved(i,j,a) = gaussian    
        end do
      end do
    end do
  
    !$omp end parallel do
  
    ! ------------------------------------------------------------------------
    ! COMPUTE CC
    ! ------------------------------------------------------------------------
    sum_sim2 = 0.0
    sum_exp2 = 0.0
    sum_simpexp = 0.0

    do j=1, image_size
      
      do i=1, image_size
        sum_sim2 = sum_sim2 + sim_image(i,j) ** 2
        sum_exp2 = sum_exp2 + exp_image(i,j) ** 2
        sum_simpexp = sum_simpexp + exp_image(i,j) * sim_image(i,j)
      end do
    end do
  
    cv= (sum_simpexp / sqrt(sum_sim2 * sum_exp2))
    force_constant = enefunc%restraint_const(1,inum)
    eexp = force_constant * (1.0_wp - cv)
    corrcoeff_save = cv
  
    ! ------------------------------------------------------------------------
    ! GRADIENT COMPUTATION
    ! ------------------------------------------------------------------------
    if (calc_force) then
  
      const1 = 1/ sqrt(sum_sim2 * sum_exp2)
      const2 = sum_simpexp / (sqrt(sum_exp2) * sum_sim2**(1.5_wp))
      emfit_img_force(:,:) = 0.0_wp
  
      !$omp parallel do default(none)                                 &
      !$omp private(a,n, i, j,mu, dpsim,dcc)     &
      !$omp shared(pixels, n_pix, image_size, pixel_size, &
      !$omp 	rot_coord, sigma, sim_image, gaussians_saved, atom_id, group_id,&
      !$omp 	exp_image, const1, const2, force_constant,&
      !$omp 	emfit_img_force,inv_rot_matrix, force, experiments)
      !
  
      do a= 1, n_atoms_group
  
        n=  atom_id(a,group_id)
        dcc = 0.0_wp
  
        do j=pixels(2,a), pixels(2,a) + n_pix
          do i=pixels(1,a), pixels(1,a) + n_pix

            mu = (/i,j/)
            mu = (mu -1- (image_size / 2)) * pixel_size
            dpsim = -(rot_coord(1:2,a) - mu) * gaussians_saved(i,j,a) / (sigma ** 2)
            dcc = dcc + ((exp_image(i,j) * dpsim * const1) - ( sim_image(i,j) * dpsim * const2))

          end do
        enddo
  
        emfit_img_force(1:2,n) = dcc * force_constant
        emfit_img_force(1:3,n) = matmul(experiments%emfit_img%inv_rot_matrix(1:3,1:3), &
                                  emfit_img_force(1:3,n))
  
        force(1:3,n) = force(1:3,n) + emfit_img_force(1:3,n)
      end do
  
      !$omp end parallel do
    end if


    ! ------------------------------------------------------------------------
    ! WRITE OUTPUT IMAGE
    ! ------------------------------------------------------------------------
    if (main_rank) then
      outfile = emfit_target_test(:index(emfit_target_test, '.', back=.true.)-1) //"_sim.spi"
      if (mod(emfit_icycle,100) == 0) then 
        call write_spi(outfile, sim_image)
      endif
    endif

    return
  end subroutine compute_energy_experimental_restraint_emfit_img
  
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine Euler_angle2Matrix(roll, tilt, yaw, A)
    real(wp),                intent(in)    :: roll
    real(wp),                intent(in)    :: tilt
    real(wp),                intent(in)    :: yaw
    real(wp),                intent(inout) :: A(3,3)

    !local variables
    real(wp) :: pi
    real(wp) :: a1, b1, g1, sa, ca, sb, cb, sg, cg, cc, cs, sc, ss

    pi = acos(-1.0_wp)
    a1=pi*roll/180.0_wp
    b1=pi*tilt/180.0_wp
    g1=pi*yaw/180.0_wp
    sa= sin(a1)
    ca= cos(a1)
    sb= sin(b1)
    cb= cos(b1)
    sg= sin(g1)
    cg= cos(g1)
    cc = cb * ca
    cs = cb * sa
    sc = sb * ca
    ss = sb * sa
    A(1, 1) =  cg * cc - sg * sa
    A(1, 2) =  cg * cs + sg * ca
    A(1, 3) = -cg * sb
    A(2, 1) = -sg * cc - cg * sa
    A(2, 2) = -sg * cs + cg * ca
    A(2, 3) =  sg * sb
    A(3, 1) =  sc
    A(3, 2) =  ss
    A(3, 3) =  cb

  end subroutine Euler_angle2Matrix

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         read_header_spi
  !> @brief		read image_size in header of Spider image 
  !!			
  !! @authors           Rémi Vuillemot
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine read_header_spi(filename, image_size)
    character(MaxFilename),         intent(in)    	 :: filename
    integer,						intent(inout)    :: image_size

    !local variable
    real(4)								:: buffer
    integer ::unit_no, status, i, LENBYT, LABREC, LABBYT

    unit_no = get_unit_no()
    open(UNIT=unit_no, FILE=filename, FORM="UNFORMATTED", access="STREAM", convert='NATIVE')

    LABBYT=256
    i=1
    do while(i<=LABBYT)
      read(unit_no,IOSTAT=status) buffer
      if (status /= 0) &
        call error_msg("Read_Spider> Error while reading spider image")
      if (i==12) then
        image_size = nint(buffer)
      else if (i==22) then
	      LABBYT = ceiling(real(buffer)/4)
      endif
      i=i+1
    end do
    call close_file(unit_no)
    return

  end subroutine read_header_spi  

   !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         read_data_spi
  !> @brief		read Spider image format
  !!			
  !! @authors           Rémi Vuillemot
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine read_data_spi(filename, target_image)
    character(MaxFilename),         intent(in)    	 :: filename
    real(wp),						intent(inout)    :: target_image(:,:)

    !local variable
    real(4)								:: buffer
    integer :: i,j, im_size,unit_no, status,LABBYT, NX, LABREC, LENBYT

    im_size = size(target_image(:,1))

    unit_no = get_unit_no()
    open(UNIT=unit_no, FILE=filename, FORM="UNFORMATTED", access="STREAM", convert='NATIVE')

    i=1
    LABBYT=256
    do while(i<=LABBYT)
      read(unit_no,IOSTAT=status) buffer
      if (status /= 0) &
        call error_msg("Read_Spider> Error while reading spider image")
      if (i==22) then
	      LABBYT = ceiling(real(buffer)/4)
      endif
      if (i==12) then
	      NX = real(buffer)
      endif
      if (i==13) then
	      LABREC = real(buffer)
      endif
      if (i==23) then
	      LENBYT = real(buffer)
      endif

      i=i+1
    end do

    do i= 1, im_size
      do j= 1, im_size
        read(unit_no,IOSTAT=status) buffer
        if (status == 0) then
          target_image(j,i) = real(buffer,wp)
        else
          call error_msg("Read_Spider> Error while reading spider image")
        endif
      end do
    end do

    call close_file(unit_no)
    return

  end subroutine read_data_spi  

     !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         write_spi
  !> @brief		write Spider image format
  !!			
  !! @authors           Rémi Vuillemot
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine write_spi(filename, target_image)
    character(MaxFilename),         intent(in)    	 :: filename
    real(wp),						intent(inout)    :: target_image(:,:)

    !local variable
    integer :: i,j, unit_no, status
    integer :: LABBYT, NX,LENBYT, LABREC
    logical :: file_exists


    ! HEADER DATA
    NX = size(target_image(:,1))
    LENBYT = NX * 4
    LABREC = 1024 / LENBYT
    IF (MOD(1024,LENBYT) .NE. 0) LABREC = LABREC + 1
    LABBYT = LABREC * LENBYT


    inquire(file=filename, exist=file_exists)
    if (file_exists) then
      unit_no = get_unit_no()
      open(UNIT=unit_no, FILE=filename, FORM="UNFORMATTED", access="STREAM", convert='NATIVE', status="replace")
  
    else
      unit_no = get_unit_no()
      open(UNIT=unit_no, FILE=filename, FORM="UNFORMATTED", access="STREAM", convert='NATIVE', status="new")
  
    endif


    ! WRITE HEADER
    ! 1 NZ
    write(unit_no) real(1,4)
    ! 2 NY
    write(unit_no) real(NX,4)  
    ! 3 IREC
    write(unit_no) real(LABREC + NX,4)
    ! 4
    write(unit_no) real(0,4)
    ! 5 IFORM
    write(unit_no) real(1,4)
    ! skip 6
    do i =1, 6
      write(unit_no) real(0,4)
    enddo
    ! 12 NX
    write(unit_no) real(NX,4)
    ! 13 LABREC
    write(unit_no) real(LABREC,4)
    ! skip 7
    do i =1, 7
      write(unit_no) real(0,4)
    enddo
    ! 21 SCALE
    write(unit_no) real(1.0,4)
    ! 22 LABBYT
    write(unit_no) real(LABBYT,4)
    ! 23 LENBYT
    write(unit_no) real(LENBYT,4)
    ! skip
    do i =1, (int(LABBYT/4) -23)
      write(unit_no) real(0,4)
    enddo

    do i= 1, NX
      do j= 1, NX
        write(unit_no) real(target_image(j,i),4)
      end do
    end do

    call close_file(unit_no)
    return

  end subroutine write_spi  

end module at_experiments_mod