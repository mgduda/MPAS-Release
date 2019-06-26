! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module atm_core_interface

   use mpas_attlist
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_io_units

   contains


   !***********************************************************************
   !
   !  routine atm_setup_core
   !
   !> \brief   Atmosphere core setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a core_type for the atm core.
   !
   !-----------------------------------------------------------------------
   subroutine atm_setup_core(core)

      use mpas_derived_types, only : core_type
      use atm_core, only : atm_core_init, atm_core_run, atm_core_finalize

      implicit none

      type (core_type), pointer :: core

      core % core_init => atm_core_init
      core % core_run => atm_core_run
      core % core_finalize => atm_core_finalize
      core % define_packages => atm_define_packages
      core % setup_packages => atm_setup_packages
      core % setup_decompositions => atm_setup_decompositions
      core % setup_clock => atm_setup_clock
      core % setup_log => atm_setup_log
      core % get_mesh_stream => atm_get_mesh_stream
      core % setup_immutable_streams => atm_setup_immutable_streams
      core % setup_derived_dimensions => atm_setup_derived_dimensions
      core % setup_decomposed_dimensions => atm_setup_decomposed_dimensions
      core % setup_block => atm_setup_block
      core % setup_namelist => atm_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'

#include "inc/core_variables.inc"

   end subroutine atm_setup_core


   !***********************************************************************
   !
   !  routine atm_setup_domain
   !
   !> \brief   Atmosphere domain setup routine
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    18 March 2015
   !> \details 
   !>  This routine is intended to setup the necessary variables within 
   !>  a domain_type for the init atm core.
   !
   !-----------------------------------------------------------------------
   subroutine atm_setup_domain(domain)

      use mpas_derived_types, only : domain_type

      implicit none

      type (domain_type), pointer :: domain

#include "inc/domain_variables.inc"

   end subroutine atm_setup_domain


   !***********************************************************************
   !
   !  function atm_setup_packages
   !
   !> \brief   Package setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  This routine is responsible for setting up packages for the
   !>  atmosphere core. It may use ay logic based on configuration options
   !>  to set packages variables to either .true. or .false. Model fields are
   !>  not allocated until after this routine has been called.
   !
   !-----------------------------------------------------------------------
   function atm_setup_packages(configs, packages, iocontext) result(ierr)

      use mpas_dmpar
      use mpas_derived_types, only : mpas_pool_type, mpas_io_context_type
      use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_package

#ifdef DO_PHYSICS
      use mpas_atmphys_control
      use mpas_atmphys_packages
#endif

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      type (mpas_pool_type), intent(inout) :: packages
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer :: ierr

      logical, pointer :: iauActive
      character(len=StrKIND), pointer :: config_iau_option
      integer, pointer :: config_num_passive_tracers

      logical, pointer :: tracer1Active
      logical, pointer :: tracer2Active
      logical, pointer :: tracer3Active
      logical, pointer :: tracer4Active
      logical, pointer :: tracer5Active
      logical, pointer :: tracer6Active
      logical, pointer :: tracer7Active
      logical, pointer :: tracer8Active
      logical, pointer :: tracer9Active
      logical, pointer :: tracer10Active

      ierr = 0

      nullify(config_iau_option)
      call mpas_pool_get_config(configs, 'config_IAU_option', config_iau_option)

      nullify(iauActive)
      call mpas_pool_get_package(packages, 'iauActive', iauActive)

      if (trim(config_iau_option) /= 'off') then
         iauActive = .true.
      else
         iauActive = .false.
      end if

      !
      ! Set up packages for passive tracers
      !
      nullify(config_num_passive_tracers)
      call mpas_pool_get_config(configs, 'config_num_passive_tracers', config_num_passive_tracers)

      nullify(tracer1Active)
      nullify(tracer2Active)
      nullify(tracer3Active)
      nullify(tracer4Active)
      nullify(tracer5Active)
      nullify(tracer6Active)
      nullify(tracer7Active)
      nullify(tracer8Active)
      nullify(tracer9Active)
      nullify(tracer10Active)
      call mpas_pool_get_package(packages, 'tracer1Active',  tracer1Active)
      call mpas_pool_get_package(packages, 'tracer2Active',  tracer2Active)
      call mpas_pool_get_package(packages, 'tracer3Active',  tracer3Active)
      call mpas_pool_get_package(packages, 'tracer4Active',  tracer4Active)
      call mpas_pool_get_package(packages, 'tracer5Active',  tracer5Active)
      call mpas_pool_get_package(packages, 'tracer6Active',  tracer6Active)
      call mpas_pool_get_package(packages, 'tracer7Active',  tracer7Active)
      call mpas_pool_get_package(packages, 'tracer8Active',  tracer8Active)
      call mpas_pool_get_package(packages, 'tracer9Active',  tracer9Active)
      call mpas_pool_get_package(packages, 'tracer10Active', tracer10Active)

      tracer1Active  = .false.
      tracer2Active  = .false.
      tracer3Active  = .false.
      tracer4Active  = .false.
      tracer5Active  = .false.
      tracer6Active  = .false.
      tracer7Active  = .false.
      tracer8Active  = .false.
      tracer9Active  = .false.
      tracer10Active = .false.

      if (config_num_passive_tracers >= 1) then
         call mpas_log_write('----- Setting up packages for passive tracers -----')
         call mpas_log_write('Activating tracer 1')
         tracer1Active = .true.
      if (config_num_passive_tracers >= 2) then
         call mpas_log_write('Activating tracer 2')
         tracer2Active = .true.
      if (config_num_passive_tracers >= 3) then
         call mpas_log_write('Activating tracer 3')
         tracer3Active = .true.
      if (config_num_passive_tracers >= 4) then
         call mpas_log_write('Activating tracer 4')
         tracer4Active = .true.
      if (config_num_passive_tracers >= 5) then
         call mpas_log_write('Activating tracer 5')
         tracer5Active = .true.
      if (config_num_passive_tracers >= 6) then
         call mpas_log_write('Activating tracer 6')
         tracer6Active = .true.
      if (config_num_passive_tracers >= 7) then
         call mpas_log_write('Activating tracer 7')
         tracer7Active = .true.
      if (config_num_passive_tracers >= 8) then
         call mpas_log_write('Activating tracer 8')
         tracer8Active = .true.
      if (config_num_passive_tracers >= 9) then
         call mpas_log_write('Activating tracer 9')
         tracer9Active = .true.
      if (config_num_passive_tracers >= 10) then
         call mpas_log_write('Activating tracer 10')
         tracer10Active = .true.
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if
      end if


#ifdef DO_PHYSICS
      !check that all the physics options are correctly defined and that at
      !least one physics parameterization is called (using the logical moist_physics):
      call physics_namelist_check(configs)

      ierr = atmphys_setup_packages(configs,packages,iocontext)
      if(ierr /= 0) then
         call mpas_log_write('Package setup failed for atmphys in core_atmosphere', messageType=MPAS_LOG_ERR)
      endif
#endif

   end function atm_setup_packages


   !***********************************************************************
   !
   !  function atm_setup_clock
   !
   !> \brief   Simulation clock setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details 
   !>  The purpose of this routine is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This routine is called from the superstructure after the framework 
   !>  has been initialized but before any fields have been allocated and 
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   function atm_setup_clock(core_clock, configs) result(ierr)

      use mpas_derived_types, only : MPAS_Clock_type, mpas_pool_type
      use atm_core, only : atm_simulation_clock_init

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      ierr = 0

      call atm_simulation_clock_init(core_clock, configs, ierr)

   end function atm_setup_clock


   !***********************************************************************
   !
   !  function atm_setup_log
   !
   !> \brief   Log setup routine
   !> \author  Matt Hoffman
   !> \date    14 February 2017
   !> \details
   !>  The purpose of this routine is to set up the logging manager
   !>  and allow the core to specify details of the configuration.
   !
   !-----------------------------------------------------------------------
   function atm_setup_log(logInfo, domain) result(iErr)!{{{

      use mpas_derived_types, only : mpas_log_type, domain_type
      use mpas_log, only : mpas_log_init, mpas_log_open

      implicit none

      type (mpas_log_type), intent(inout), pointer :: logInfo  !< logging information object to set up
      type (domain_type), intent(in), pointer :: domain  !< domain object to provide info for setting up log manager
      integer :: iErr

      ! Local variables
      integer :: local_err

      iErr = 0

      ! Initialize log manager
      call mpas_log_init(logInfo, domain, err=local_err)
      iErr = ior(iErr, local_err)

      ! Set core specific options here
      ! (At present, there are not any.  There could eventually be choices about the file naming conventions
      !  or other settings controlling behavior.)

      ! After core has had a chance to modify log defaults, open the output log
      call mpas_log_open(err=local_err)
      iErr = ior(iErr, local_err)

      call mpas_log_write('')
#ifdef SINGLE_PRECISION
      call mpas_log_write('Using default single-precision reals')
#else
      call mpas_log_write('Using default double-precision reals')
#endif
      call mpas_log_write('')

   end function atm_setup_log!}}}


   !***********************************************************************
   !
   !  function atm_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details 
   !>  This routine returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping 
   !>  routine. At the time this routine is called, only namelist options 
   !>  are available.
   !
   !-----------------------------------------------------------------------
   function atm_get_mesh_stream(configs, stream) result(ierr)

      use mpas_kind_types, only : StrKIND
      use mpas_derived_types, only : mpas_pool_type
      use mpas_pool_routines, only : mpas_pool_get_config

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      logical, pointer :: config_do_restart

      ierr = 0

      call mpas_pool_get_config(configs, 'config_do_restart', config_do_restart)

      if (.not. associated(config_do_restart)) then
         call mpas_log_write('config_do_restart was not found when defining mesh stream.', messageType=MPAS_LOG_ERR)
         ierr = 1
      else if (config_do_restart) then
         write(stream,'(a)') 'restart'
      else
         write(stream,'(a)') 'input'
      end if

   end function atm_get_mesh_stream


   !***********************************************************************
   !
   !  function atm_setup_decompositions
   !
   !> \brief   Decomposition setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    11 March 2015
   !> \details 
   !>  This function is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function atm_setup_decompositions(decompList) result(ierr)

      use mpas_derived_types, only : mpas_decomp_list, mpas_decomp_function, MPAS_DECOMP_NOERR
      use mpas_decomp, only : mpas_decomp_create_decomp_list, mpas_decomp_register_method, &
                              mpas_uniform_decomp

      implicit none

      type (mpas_decomp_list), pointer :: decompList
      integer :: ierr

      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => mpas_uniform_decomp

      call mpas_decomp_register_method(decompList, 'uniform', decompFunc, ierr)

      if ( ierr == MPAS_DECOMP_NOERR ) then
         ierr = 0
      end if

   end function atm_setup_decompositions


   !***********************************************************************
   !
   !  function atm_setup_block
   !
   !> \brief   Block setup function
   !> \author  Doug Jacobsen, Michael Duda
   !> \date    03/18/2015
   !> \details 
   !>  This function is a wrapper function to properly setup a block to be
   !>  an atmosphere core block.
   !
   !-----------------------------------------------------------------------
   function atm_setup_block(block) result(ierr)

      use mpas_derived_types, only : block_type

      implicit none

      type (block_type), pointer :: block
      integer :: ierr

      ierr = 0

      call atm_generate_structs(block, block % structs, block % dimensions, block % packages)

   end function atm_setup_block


#include "inc/setup_immutable_streams.inc"

#include "inc/block_dimension_routines.inc"

#include "inc/define_packages.inc"

#include "inc/structs_and_variables.inc"

#include "inc/namelist_call.inc"

#include "inc/namelist_defines.inc"

end module atm_core_interface
