!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_minimize_mod
!> @brief   energy minimization
!! @authors Jaewoon Jung (JJ), Norio Takase (NT), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_minimize_mod

  use sp_output_mod
  use sp_update_domain_mod
  use sp_dynvars_mod
  use sp_energy_mod
  use sp_constraints_mod
  use sp_communicate_mod
  use sp_constraints_str_mod
  use sp_output_str_mod
  use sp_minimize_str_mod
  use sp_dynvars_str_mod
  use sp_boundary_str_mod
  use sp_pairlist_str_mod
  use sp_enefunc_str_mod
  use sp_domain_str_mod
  use sp_constraints_str_mod
  use sp_constraints_mod
  use fileio_control_mod
  use timers_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_min_info
    integer          :: method           = MinimizeMethodSD
    integer          :: nsteps           = 100
    integer          :: eneout_period    =  10
    integer          :: crdout_period    =   0
    integer          :: rstout_period    =   0
    integer          :: nbupdate_period  =  10
    logical          :: verbose          = .false.
    real(wp)         :: force_scale_init = 0.00005_wp
    real(wp)         :: force_scale_max  = 0.0001_wp
  end type s_min_info

  ! subroutines
  public  :: show_ctrl_minimize
  public  :: read_ctrl_minimize
  public  :: setup_minimize
  public  :: run_min
  public  :: steepest_descent

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_minimize
  !> @brief        show MINIMIZE section usage
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_minimize(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('min')

        write(MsgOut,'(A)') '[MINIMIZE]'
        write(MsgOut,'(A)') 'method        = SD        # [SD]'
        write(MsgOut,'(A)') 'nsteps        = 100       # number of minimization steps'
        write(MsgOut,'(A)') '# eneout_period = 10        # energy output period'
        write(MsgOut,'(A)') '# crdout_period = 0         # coordinates output period'
        write(MsgOut,'(A)') '# rstout_period = 0         # restart output period'
        write(MsgOut,'(A)') '# nbupdate_period = 10      # nonbond update period'
        write(MsgOut,'(A)') '# verbose       = NO        # output verbosly'
        write(MsgOut,'(A)') ' '


      end select

    else

      select case (run_mode)

      case ('min')

        write(MsgOut,'(A)') '[MINIMIZE]'
        write(MsgOut,'(A)') 'method        = SD        # [SD]'
        write(MsgOut,'(A)') 'nsteps        = 100       # number of minimization steps'
        write(MsgOut,'(A)') ' '

      end select


    end if

    return

  end subroutine show_ctrl_minimize
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_minimize
  !> @brief        read control parameters in MINIMIZE section
  !! @authors      JJ
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   min_info : MINIMIZE section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_minimize(handle, min_info)

    ! parameters
    character(*),            parameter     :: Section = 'Minimize'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_min_info),        intent(inout) :: min_info


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'method',      &
                               min_info%method, MinimizeMethodTypes)
    call read_ctrlfile_integer(handle, Section, 'nsteps',     &
                               min_info%nsteps)
    call read_ctrlfile_integer(handle, Section, 'eneout_period', &
                               min_info%eneout_period)
    call read_ctrlfile_integer(handle, Section, 'crdout_period', &
                               min_info%crdout_period)
    call read_ctrlfile_integer(handle, Section, 'rstout_period', &
                               min_info%rstout_period)
    call read_ctrlfile_integer(handle, Section, 'nbupdate_period', &
                               min_info%nbupdate_period)
    call read_ctrlfile_real    (handle, Section, 'force_scale_max', &
                               min_info%force_scale_max)
    call read_ctrlfile_real    (handle, Section, 'force_scale_init', &
                               min_info%force_scale_init)
    call read_ctrlfile_logical(handle, Section, 'verbose',        &
                               min_info%verbose)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Minimize> Parameters of MIN'

      write(MsgOut,'(A20,A10,A20,I10)')                                   &
            '  method          = ', MinimizeMethodTypes(min_info%method), &
            '  nsteps          = ', min_info%nsteps
      write(MsgOut,'(A20,I10,A20,I10)')                                   &
            '  eneout_period   = ', min_info%eneout_period,               &
            '  crdout_period   = ', min_info%crdout_period
      write(MsgOut,'(A20,I10,A20,I10)')                                   &
            '  rstout_period   = ', min_info%rstout_period,               &
            '  nbupdate_period = ', min_info%nbupdate_period
      write(MsgOut,'(A20,F10.3,A20,F10.3)')                               &
            '  force_scale_init= ', min_info%force_scale_init,            &
            '  force_scale_max = ', min_info%force_scale_max
      ! verbose
      if (min_info%verbose) then
        write(MsgOut,'(A)') '  verbose         =        yes'
      else
        write(MsgOut,'(A)') '  verbose         =         no'
      end if

      write(MsgOut,'(A)') ' '
    end if


    ! error check
    !
    if (main_rank) then
      if (min_info%eneout_period > 0 .and.                        &
          mod(min_info%nsteps, min_info%eneout_period) /= 0) then
        write(MsgOut,'(A)') 'Read_Ctrl_Minimize> Error in eneout_period'
        write(MsgOut,'(A)') '  mod(nsteps, eneout_period) is not ZERO'
        call error_msg
      end if

      if (min_info%crdout_period > 0 .and.                        &
          mod(min_info%nsteps, min_info%crdout_period) /= 0) then
        write(MsgOut,'(A)') 'Read_Ctrl_Minimize> Error in crdout_period'
        write(MsgOut,'(A)') '  mod(nsteps, crdout_period) is not ZERO'
        call error_msg
      end if

      if (min_info%rstout_period > 0 .and.                        &
          mod(min_info%nsteps, min_info%rstout_period) /= 0) then
        write(MsgOut,'(A)') 'Read_Ctrl_Minimize> Error in rstout_period'
        write(MsgOut,'(A)') '  mod(nsteps, rstout_period) is not ZERO'
        call error_msg
      end if

      if (min_info%nbupdate_period > 0 .and.                      &
          mod(min_info%nsteps, min_info%nbupdate_period) /= 0) then
        write(MsgOut,'(A)') 'Read_Ctrl_Minimize> Error in nbupdate_period'
        write(MsgOut,'(A)') '  mod(nsteps, nbupdate_period) is not ZERO'
        call error_msg
      end if

    end if


    return

  end subroutine read_ctrl_minimize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_minimize
  !> @brief        setup minimize information
  !> @authors      JJ
  !! @param[in]    min_info : MINIMIZE section control parameters information
  !! @param[out]   minimize : minimize information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_minimize(min_info, minimize)

    ! formal arguments
    type(s_min_info),        intent(in)    :: min_info
    type(s_minimize),        intent(inout) :: minimize


    ! initialize variables in minimize
    !
    call init_minimize(minimize)


    ! setup minimize
    !
    minimize%method          = min_info%method
    minimize%nsteps          = min_info%nsteps
    minimize%eneout_period   = min_info%eneout_period
    minimize%crdout_period   = min_info%crdout_period
    minimize%rstout_period   = min_info%rstout_period
    minimize%nbupdate_period = min_info%nbupdate_period
    minimize%verbose          = min_info%verbose
    minimize%force_scale_init = min_info%force_scale_init
    minimize%force_scale_max  = min_info%force_scale_max

    return

  end subroutine setup_minimize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_min
  !> @brief        perform energy minimization
  !! @authors      CK
  !! @param[inout] output      : output information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions
  !! @param[inout] dynvars     : dynamic variables
  !! @param[inout] minimize    : minimize information
  !! @param[inout] pairlist    : non-bond pair list
  !! @param[inout] boundary    : boundary condition
  !! @param[inout] comm        : information of communication
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_min(output, domain, enefunc, dynvars, minimize, &
                     pairlist, boundary, constraints, comm)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_domain),   target, intent(inout) :: domain
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_minimize),         intent(inout) :: minimize
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_comm),             intent(inout) :: comm

    ! Open output files
    !
    call open_output(output)

    ! MIN main loop
    !
    select case (minimize%method)

    case (MinimizeMethodSD)

      call steepest_descent (output, domain, enefunc, dynvars, minimize, &
                             pairlist, boundary, constraints, comm)

    end select

    ! close output files
    !
    call close_output(output)


    return

  end subroutine run_min


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    steepest_descent
  !> @brief        steepest descent integrator
  !> @authors      JJ
  !! @param[inout] output   : output information
  !! @param[inout] domain   : domain information
  !! @param[inout] enefunc  : potential energy functions information
  !! @param[inout] dynvars  : dynamic variables information
  !! @param[inout] minimize : minimize information
  !! @param[inout] pairlist : non-bond pair list information
  !! @param[inout] boundary : boundary conditions information
  !! @param[inout] comm     : information of communication
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine steepest_descent(output, domain, enefunc, dynvars, minimize, &
                              pairlist, boundary, constraints, comm)

    ! formal arguments
    type(s_output),          intent(inout) :: output
    type(s_domain),  target, intent(inout) :: domain
    type(s_enefunc),         intent(inout) :: enefunc
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_minimize),        intent(inout) :: minimize
    type(s_pairlist),        intent(inout) :: pairlist
    type(s_boundary),        intent(inout) :: boundary
    type(s_constraints),     intent(inout) :: constraints
    type(s_comm),            intent(inout) :: comm

    ! local variables
    real(dp)                 :: energy, energy_ref, delta_energy
    real(dp)                 :: delta_r, delta_rini, delta_rmax
    real(dp)                 :: rmsg, ene_total
    integer                  :: nsteps, ncell, nb, i, j, jx, k, l

    real(dp),        pointer :: coord(:,:,:), coord_ref(:,:,:)
    real(dp),        pointer :: force(:,:,:)
    real(wp),        pointer :: coord_pbc(:,:,:), force_omp(:,:,:,:)
    real(dp),        pointer :: force_long(:,:,:)
    real(wp),        pointer :: force_pbc(:,:,:,:)
    real(dp),        pointer :: virial_cell(:,:), virial(:,:)
    integer,         pointer :: natom(:), nsolute(:), nwater(:) 
    integer,         pointer :: water_list(:,:,:)


    natom       => domain%num_atom
    nsolute     => domain%num_solute
    nwater      => domain%num_water
    nwater      => domain%num_water 
    water_list  => domain%water_list
    coord       => domain%coord
    coord_ref   => domain%coord_ref
    coord_pbc   => domain%translated
    force       => domain%force
    force_long  => domain%force_long
    force_omp   => domain%force_omp
    force_pbc   => domain%force_pbc
    virial_cell => domain%virial_cellpair
    virial      => dynvars%virial

    ncell       =  domain%num_cell_local
    nb          =  domain%num_cell_boundary
    nsteps      =  minimize%nsteps

    delta_rini = real(minimize%force_scale_init,dp)
    delta_rmax = real(minimize%force_scale_max,dp)
    delta_r    = delta_rini


    ! Compute energy of the initial structure
    !
    if (enefunc%table%tip4) call decide_dummy(domain, constraints, coord)
    call compute_energy(domain, enefunc, pairlist, boundary, coord, &
                        .true., .true., .true., .true.,      &
                        enefunc%nonb_limiter,                &
                        dynvars%energy, coord_pbc, force,    &
                        force_long, force_omp,               &
                        force_pbc, virial_cell,              &
                        dynvars%virial, dynvars%virial_long, &
                        dynvars%virial_extern)

    call communicate_force(domain, comm, force)

    if (enefunc%table%tip4) &
      call water_force_redistribution(constraints, domain, force, virial)

    dynvars%total_pene = dynvars%energy%bond          &
                       + dynvars%energy%angle         &
                       + dynvars%energy%urey_bradley  &
                       + dynvars%energy%dihedral      &
                       + dynvars%energy%improper      &
                       + dynvars%energy%cmap          &
                       + dynvars%energy%electrostatic &
                       + dynvars%energy%van_der_waals &
                       + dynvars%energy%restraint_position
    dynvars%energy%total = dynvars%total_pene

    rmsg = 0.0_dp
    do j = 1, ncell
      do jx = 1, nsolute(j)
        rmsg = rmsg + force(1,jx,j)*force(1,jx,j) &
                    + force(2,jx,j)*force(2,jx,j) &
                    + force(3,jx,j)*force(3,jx,j)
      end do
      do l = 1, nwater(j)
        do k = 1, 3
          jx = water_list(k,l,j)
          rmsg = rmsg + force(1,jx,j)*force(1,jx,j) &
                      + force(2,jx,j)*force(2,jx,j) &
                      + force(3,jx,j)*force(3,jx,j)
        end do
      end do
    end do
    rmsg = rmsg / real(3*domain%num_atom_all,dp)

#ifdef MPI
    call mpi_allreduce(mpi_in_place, rmsg, 1, mpi_real8, mpi_sum, &
                       mpi_comm_city, ierror)
#endif
    rmsg = sqrt(rmsg)

    do j = 1, ncell
      do jx = 1, natom(j)
        coord_ref(1:3,jx,j) = coord(1:3,jx,j)
      end do
    end do

    do j = 1, ncell
      do jx = 1, nsolute(j)
        coord(1,jx,j) = coord(1,jx,j) + delta_r*force(1,jx,j)/rmsg
        coord(2,jx,j) = coord(2,jx,j) + delta_r*force(2,jx,j)/rmsg
        coord(3,jx,j) = coord(3,jx,j) + delta_r*force(3,jx,j)/rmsg
      end do
      do l = 1, nwater(j)
        do k = 1, 3
          jx = water_list(k,l,j)
          coord(1,jx,j) = coord(1,jx,j) + delta_r*force(1,jx,j)/rmsg
          coord(2,jx,j) = coord(2,jx,j) + delta_r*force(2,jx,j)/rmsg
          coord(3,jx,j) = coord(3,jx,j) + delta_r*force(3,jx,j)/rmsg
        end do
      end do
    end do

    if (constraints%fast_water) &
      call compute_settle_min(coord_ref, domain, constraints, coord)

    if (enefunc%table%tip4) call decide_dummy(domain, constraints, coord)

    dynvars%rms_gradient = rmsg

    call output_dynvars(output, enefunc, dynvars)

    ! Main loop of the Steepest descent method
    !
    do i = 1, nsteps

      dynvars%step = i

      ! save old energy
      !
      if (my_country_rank == 0) energy_ref = dynvars%total_pene

      ! Compute energy
      !
      call communicate_coor(domain, comm)
      call compute_energy(domain, enefunc, pairlist, boundary, coord, &
                          .true., .true., .true., .true.,      &
                           enefunc%nonb_limiter,               &
                          dynvars%energy, coord_pbc,           &
                          force, force_long, force_omp,        &
                          force_pbc, virial_cell,              &
                          dynvars%virial, dynvars%virial_long, &
                          dynvars%virial_extern)

      call communicate_force(domain, comm, force)
      
      if (enefunc%table%tip4) &
        call water_force_redistribution(constraints, domain, force, virial)

      ! Broad cast total energy
      !
      dynvars%total_pene = dynvars%energy%bond          &
                         + dynvars%energy%angle         &
                         + dynvars%energy%urey_bradley  &
                         + dynvars%energy%dihedral      &
                         + dynvars%energy%improper      &
                         + dynvars%energy%cmap          &
                         + dynvars%energy%electrostatic &
                         + dynvars%energy%van_der_waals &
                         + dynvars%energy%restraint_position

      dynvars%energy%total = dynvars%total_pene

      if (my_country_rank == 0) delta_energy = dynvars%energy%total - energy_ref
#ifdef MPI
      call mpi_bcast(delta_energy, 1, mpi_real8, 0, mpi_comm_country, ierror)
#endif
      if (delta_energy > 0.0_dp) then
        delta_r = 0.5_dp*delta_r
      else
        delta_r = min(delta_rmax, 1.2_dp*delta_r)
      end if

      rmsg = 0.0_dp
      do j = 1, ncell
        do jx = 1, nsolute(j)
          rmsg = rmsg + force(1,jx,j)*force(1,jx,j) &
                      + force(2,jx,j)*force(2,jx,j) &
                      + force(3,jx,j)*force(3,jx,j)
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            rmsg = rmsg + force(1,jx,j)*force(1,jx,j) &
                        + force(2,jx,j)*force(2,jx,j) &
                        + force(3,jx,j)*force(3,jx,j)
          end do
        end do
      end do

      rmsg = rmsg / real(3*domain%num_atom_all,dp)

#ifdef MPI
      call mpi_allreduce(mpi_in_place, rmsg, 1, mpi_real8, mpi_sum, &
                         mpi_comm_city, ierror)
#endif

      rmsg = sqrt(rmsg)

      do j = 1, ncell
        do jx = 1, nsolute(j)
          coord(1,jx,j) = coord(1,jx,j) + delta_r*force(1,jx,j)/rmsg
          coord(2,jx,j) = coord(2,jx,j) + delta_r*force(2,jx,j)/rmsg
          coord(3,jx,j) = coord(3,jx,j) + delta_r*force(3,jx,j)/rmsg
        end do
        do l = 1, nwater(j)
          do k = 1, 3
            jx = water_list(k,l,j)
            coord(1,jx,j) = coord(1,jx,j) + delta_r*force(1,jx,j)/rmsg
            coord(2,jx,j) = coord(2,jx,j) + delta_r*force(2,jx,j)/rmsg
            coord(3,jx,j) = coord(3,jx,j) + delta_r*force(3,jx,j)/rmsg
          end do
        end do
      end do
      do j = 1, ncell
        do jx = 1, natom(j)
          coord_ref(1:3,jx,j) = coord(1:3,jx,j)
        end do
      end do
      if (constraints%fast_water) &
        call compute_settle_min(coord_ref, domain, constraints, coord)

      if (enefunc%table%tip4) call decide_dummy(domain, constraints, coord)

      dynvars%rms_gradient = rmsg

      ! output trajectory & restart
      !
      call output_min(output, enefunc, minimize, boundary, &
                                   dynvars, delta_r, domain)

      ! update nonbond pairlist
      !
      if (minimize%nbupdate_period > 0) then

        call domain_interaction_update_min(i, minimize, domain, enefunc, &
                                           pairlist, boundary, comm)

      end if

      ! output parallel I/O restart
      !
      call output_prst_min(output, enefunc, minimize, boundary, &
                                   dynvars, domain)
    end do


    return

  end subroutine steepest_descent

end module sp_minimize_mod
