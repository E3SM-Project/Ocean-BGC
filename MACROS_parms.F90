MODULE MACROS_parms

  !-----------------------------------------------------------------------------
  !   Backfill by imitating analogous module with name ecosys_parms
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !   This module manages the parameter variables for the module MACROS_mod.
  !   Most of the variables are not parameters in the Fortran sense. In the
  !   the Fortran sense, they are vanilla module variables.
  !
  !   This modules handles initializing the variables to default values and
  !   reading them from the namelist MACROS_parms. The values used are echoed
  !   to stdout for record keeping purposes.
  !
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !   variables/subroutines/function used from other modules
  !   The following are used extensively in this set of routines, so are used at
  !   the module level. The use statements for variables that are only needed
  !   locally are located at the module subprogram level.
  !-----------------------------------------------------------------------------

! !USES:

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !   public/private declarations
  !   all module variables are public and should have their values preserved
  !-----------------------------------------------------------------------------

  PUBLIC
  SAVE

  !-----------------------------------------------------------------------------
  !   kinds definitions
  !-----------------------------------------------------------------------------
   integer, parameter, public ::               &
      MACROS_char           = 256                    ,&
      MACROS_log            = kind(.true.)           ,&
      MACROS_i4             = selected_int_kind(6)   ,&
      MACROS_i8             = selected_int_kind(13)  ,&
      MACROS_r4             = selected_real_kind(6)  ,&
      MACROS_r8             = selected_real_kind(13)

  !-----------------------------------------------------------------------------
  !   floating point constants used across ecosystem module
  !-----------------------------------------------------------------------------

  real (MACROS_r8), parameter, private :: &
       c0 =   0.0_MACROS_r8,         &
       c1 =   1.0_MACROS_r8,         &
       c2 =   2.0_MACROS_r8,         &
       c10 = 10.0_MACROS_r8,         &
       p5  =  0.5_MACROS_r8

  REAL(KIND=MACROS_r8), PARAMETER :: &
       spd = 86400.0_MACROS_r8,        & ! number of seconds in a day
       dps = c1 / spd,          & ! number of days in a second
       yps = c1 / (365.0_MACROS_r8*spd)  ! number of years in a second

  type, public :: MACROS_indices_type
   integer (MACROS_i4) :: &
      prot_ind,         & ! MACROS index
      poly_ind,        & ! MACROS index
      lip_ind,         & ! MACROS index
      zooC_ind,        & ! MACROS index
      spC_ind,         & ! MACROS index
      diatC_ind,       & ! MACROS index
      diazC_ind,       & ! MACROS index
      phaeoC_ind         ! MACROS index

   character (MACROS_char), allocatable, dimension(:) ::  &
      short_name,      & ! short name of variable
      long_name,       & ! long name of variable
      units              ! units of variable
  end type MACROS_indices_type

  type, public :: MACROS_input_type
   real (MACROS_r8), allocatable, dimension(:,:,:) :: &
      MACROS_tracers
   real (MACROS_r8), allocatable, dimension(:,:) :: &
      cell_thickness
   integer (MACROS_i4), allocatable, dimension(:) ::  &
      number_of_active_levels
  end type MACROS_input_type

  !-------------------------------------------------------------------------
  ! Here SE removes piston velocity related, since the organics ride bubbles.
  ! SB will handle the surfactant and spray flux chemistry initially.
  ! Also flux and temperature are sacrificed since initial mechanism minimalist.
  !-------------------------------------------------------------------------

  type, public :: MACROS_output_type
   real (MACROS_r8), allocatable, dimension(:,:,:) :: &
      MACROS_tendencies
  end type MACROS_output_type

  !-----------------------------------------------------------------------------
  !   Sea air gas transfer originally diagnosed in this location for DMS
  !   but the macromolecules will scavenged on bubble surfaces
  !   in a routine likely handled by Burrows and company at PNNL, at least initially.
  !-----------------------------------------------------------------------------

  type, public :: MACROS_diagnostics_type
    real (MACROS_r8), allocatable, dimension(:,:) :: &
      diag_PROT_S_TOTAL,         &
      diag_POLY_S_TOTAL,         &
      diag_LIP_S_TOTAL,         &
      diag_PROT_R_TOTAL,         &
      diag_POLY_R_TOTAL,         &
      diag_LIP_R_TOTAL         
  end type MACROS_diagnostics_type

  !-----------------------------------------------------------------------------
  !   Floating point constants used across MACROS module.
  !   mm encourages se to sacrifice as much as possible this go around.
  !   Hence all light, temperature and prokaryotic details disappear. 
  !   Even Redfield goes the way of the great white buffalo
  !-----------------------------------------------------------------------------

  real (MACROS_r8) ::    &
    f_prot      ,  &
    f_poly      ,  &
    f_lip       ,  &
    k_C_p_base  ,  &
    zooC_avg    ,  &
    mort        ,  &
    k_prot_bac  ,  &
    k_poly_bac  ,  &
    k_lip_bac   ,  &   
    inject_scale 

  real (MACROS_r8), parameter :: &
    epsC      = 1.00e-8_MACROS_r8   ! small C concentration (mmol C/m^3)

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  SUBROUTINE MACROS_parms_init

    !---------------------------------------------------------------------------
    !   default parameter values
    !---------------------------------------------------------------------------

    f_prot      = 0.6_MACROS_r8
    f_poly      = 0.2_MACROS_r8
    f_lip       = 0.2_MACROS_r8
    k_C_p_base  = dps*0.1_MACROS_r8 
    zooC_avg    = 0.3_MACROS_r8
    mort        = 0.0_MACROS_r8
    k_prot_bac  = dps*0.1_MACROS_r8
    k_poly_bac  = dps*0.01_MACROS_r8
    k_lip_bac   = dps*1.0_MACROS_r8
    inject_scale = 1.0_MACROS_r8

    !---------------------------------------------------------------------------

  END SUBROUTINE MACROS_parms_init

  !*****************************************************************************

END MODULE MACROS_parms
