MODULE DMS_parms

  !-----------------------------------------------------------------------------
  !   Backfill by imitating analogous module with name ecosys_parms
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !   This module manages the parameter variables for the module DMS_mod.
  !   Most of the variables are not parameters in the Fortran sense. In the
  !   the Fortran sense, they are vanilla module variables.
  !
  !   This modules handles initializing the variables to default values and
  !   reading them from the namelist DMS_parms. The values used are echoed
  !   to stdout for record keeping purposes.
  !
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !   variables/subroutines/function used from other modules
  !   The following are used extensively in this tracegas, so are used at
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
      DMS_char           = 256                    ,&
      DMS_log            = kind(.true.)           ,&
      DMS_i4             = selected_int_kind(6)   ,&
      DMS_i8             = selected_int_kind(13)  ,&
      DMS_r4             = selected_real_kind(6)  ,&
      DMS_r8             = selected_real_kind(13)

  !-----------------------------------------------------------------------------
  !   floating point constants used across ecosystem module
  !-----------------------------------------------------------------------------

  real (DMS_r8), parameter, private :: &
       c0 =   0.0_DMS_r8,         &
       c1 =   1.0_DMS_r8,         &
       c2 =   2.0_DMS_r8,         &
       c10 = 10.0_DMS_r8,         &
       p5  =  0.5_DMS_r8

  REAL(KIND=DMS_r8), PARAMETER :: &
       spd = 86400.0_DMS_r8,        & ! number of seconds in a day
       dps = c1 / spd,          & ! number of days in a second
       yps = c1 / (365.0_DMS_r8*spd)  ! number of years in a second

  type, public :: DMS_indices_type
   integer (DMS_i4) :: &
      dms_ind,         & ! DMS index
      dmsp_ind,        & ! DMSP index
      no3_ind,         & ! DMSP index
      doc_ind,         & ! DMSP index
      zooC_ind,        & ! DMSP index
      spC_ind,         & ! DMSP index
      spCaCO3_ind,     & ! DMSP index
      diatC_ind,       & ! DMSP index
      diazC_ind,       & ! DMSP index
      phaeoC_ind,      & ! DMSP index
      spChl_ind,       & ! DMSP index
      diatChl_ind,     & ! DMSP index
      diazChl_ind,     & ! DMSP index
      phaeoChl_ind       ! DMSP index

   character (DMS_char), allocatable, dimension(:) ::  &
      short_name,      & ! short name of variable
      long_name,       & ! long name of variable
      units              ! units of variable
  end type DMS_indices_type

  type, public :: DMS_input_type
   real (DMS_r8), allocatable, dimension(:,:,:) :: &
      DMS_tracers
   real (DMS_r8), allocatable, dimension(:,:) :: &
      cell_thickness
   integer (DMS_i4), allocatable, dimension(:) ::  &
      number_of_active_levels
  end type DMS_input_type

  type, public :: DMS_forcing_type
   real (DMS_r8), allocatable, dimension(:) :: &
      ShortWaveFlux_surface,     &
      surfacePressure,           &
      iceFraction,               &
      windSpeedSquared10m,       &
      SST,                       &
      SSS
   real (DMS_r8), allocatable, dimension(:,:) :: &
      netFlux
   logical (DMS_log) :: &
      lcalc_DMS_gas_flux
  end type DMS_forcing_type

  type, public :: DMS_output_type
   real (DMS_r8), allocatable, dimension(:,:,:) :: &
      DMS_tendencies
  end type DMS_output_type

  type, public :: DMS_flux_diagnostics_type
    real (DMS_r8), allocatable, dimension(:) :: &
      diag_DMS_IFRAC,        &
      diag_DMS_XKW,          &
      diag_DMS_ATM_PRESS,    &
      diag_DMS_PV,           &
      diag_DMS_SCHMIDT,      &
      diag_DMS_SAT,          &
      diag_DMS_SURF,         &
      diag_DMS_WS
  end type DMS_flux_diagnostics_type

  type, public :: DMS_diagnostics_type
    real (DMS_r8), allocatable, dimension(:,:) :: &
      diag_DMS_S_DMSP,          &
      diag_DMS_S_TOTAL,         &
      diag_DMS_R_B,             &
      diag_DMS_R_PHOT,          &
      diag_DMS_R_BKGND,         &
      diag_DMS_R_TOTAL,         &
      diag_DMSP_S_PHAEO,        &
      diag_DMSP_S_NONPHAEO,     &
      diag_DMSP_S_ZOO,          &
      diag_DMSP_S_TOTAL,        &
      diag_DMSP_R_B,            &
      diag_DMSP_R_BKGND,        &
      diag_DMSP_R_TOTAL,        &
      diag_Cyano_frac,          &
      diag_Cocco_frac,          &
      diag_Eukar_frac,          &
      diag_diatS,               &
      diag_diatN,               &
      diag_phytoN,              &
      diag_coccoS,              &
      diag_cyanoS,              &
      diag_eukarS,              &
      diag_diazS,               &
      diag_phaeoS,              &
      diag_zooS,                &
      diag_zooCC,               &
      diag_RSNzoo
  end type DMS_diagnostics_type

  !-----------------------------------------------------------------------------
  !   floating point constants used across DMS module
  !-----------------------------------------------------------------------------

  real (DMS_r8) ::    &
    k_S_p_base  ,  &
    zooC_avg    ,  &
    mort        ,  &
    k_conv      ,  &
    k_S_z       ,  &
    B_preexp    ,  &
    B_exp       ,  &
    k_S_B       ,  &
    k_bkgnd     ,  &
    j_dms_perI  ,  &
    inject_scale ,  &
    T_cryo_hi      ,  &
    T_cryo_lo      ,  &
    T_lo           ,  &
    T_hi           ,  &
    Min_cyano_frac ,  &
    Max_cyano_frac ,  &
    Min_yld        ,  &
    Max_yld        ,  &
    G_phaeo_S      ,  &
    Sp_ref         ,  &
    Stress_mult    ,  &
    R             ,  &
    Rs2n_diat     ,  &
    Rs2n_phaeo    ,  &
    Rs2n_cocco    ,  &
    Rs2n_cyano    ,  &
    Rs2n_eukar    ,  &
    Rs2n_diaz

  real (DMS_r8), public ::    &
    f_qsw_par_DMS = 0.45_DMS_r8   ! PAR fraction

  real (DMS_r8), parameter :: &
    epsC      = 1.00e-8_DMS_r8   ! small C concentration (mmol C/m^3)

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  SUBROUTINE DMS_parms_init

    !---------------------------------------------------------------------------
    !   default parameter values
    !---------------------------------------------------------------------------

    k_S_p_base  = 0.1_DMS_r8 * dps
    zooC_avg    = 0.3_DMS_r8
    mort        = 0.0_DMS_r8
    k_conv      = 1.0_DMS_r8 * dps
    k_S_z       = 0.1_DMS_r8 * dps
    B_preexp    = 0.1_DMS_r8
    B_exp       = 0.5_DMS_r8
    k_S_B       = 30.0_DMS_r8 * dps
    k_bkgnd     = 0.01_DMS_r8 * dps
    j_dms_perI  = 0.005_DMS_r8 * dps
    inject_scale = 1.00_DMS_r8
    T_cryo_hi      =  1.0_DMS_r8
    T_cryo_lo      =  -1.0_DMS_r8
    T_lo           =  15.0_DMS_r8
    T_hi           =  20.0_DMS_r8
    Min_cyano_frac =  0.0_DMS_r8
    Max_cyano_frac =  0.5_DMS_r8
    Min_yld        =  0.2_DMS_r8
    Max_yld        =  0.7_DMS_r8
    G_phaeo_S      =  0.4_DMS_r8
    Sp_ref         =  0.1_DMS_r8
    Stress_mult    = 10.0_DMS_r8
    R             = 0.137_DMS_r8
    Rs2n_diat     = 0.01_DMS_r8
    Rs2n_phaeo    = 0.3_DMS_r8
    Rs2n_cocco    = 0.1_DMS_r8
    Rs2n_cyano    = 0.0_DMS_r8
    Rs2n_eukar    = 0.1_DMS_r8
    Rs2n_diaz     = 0.0_DMS_r8

    !---------------------------------------------------------------------------

  END SUBROUTINE DMS_parms_init

  !*****************************************************************************

END MODULE DMS_parms
