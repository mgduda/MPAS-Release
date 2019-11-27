module atm_core

   type core_type
       character(len=64) :: modelName
       character(len=64) :: coreName
       character(len=64) :: modelVersion
       character(len=64) :: executableName
       character(len=64) :: git_version
   end type core_type


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
   subroutine atm_setup_core()

      implicit none

      type(core_type) :: core

      write(0,*) 'Called atm_setup_core'

#include "core_variables.inc"

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
   subroutine atm_setup_domain()

      implicit none

      write(0,*) 'Called atm_setup_domain'

   end subroutine atm_setup_domain

end module atm_core
