MODULE BGC_parms

  !-----------------------------------------------------------------------------

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
      BGC_char           = 256                    ,&
      BGC_log            = kind(.true.)           ,&
      BGC_i4             = selected_int_kind(6)   ,&
      BGC_i8             = selected_int_kind(13)  ,&
      BGC_r4             = selected_real_kind(6)  ,&
      BGC_r8             = selected_real_kind(13)

  !-----------------------------------------------------------------------------
  !   floating point constants used across ecosystem module
  !-----------------------------------------------------------------------------

  real (BGC_r8), parameter, private :: &
       c0 =   0.0_BGC_r8,         &
       c1 =   1.0_BGC_r8,         &
       c2 =   2.0_BGC_r8,         &
       c10 = 10.0_BGC_r8,         &
       p5  =  0.5_BGC_r8

  REAL(KIND=BGC_r8), PARAMETER :: &
       spd = 86400.0_BGC_r8,        & ! number of seconds in a day
       dps = c1 / spd,          & ! number of days in a second
       yps = c1 / (365.0_BGC_r8*spd)  ! number of years in a second

  INTEGER (KIND=BGC_i4), PARAMETER :: &
     autotroph_cnt   = 4

  real (BGC_r8), public :: T0_Kelvin_BGC

  !-----------------------------------------------------------------------------
  !   functional group
  !-----------------------------------------------------------------------------

  TYPE, PUBLIC :: autotroph_type
     CHARACTER(BGC_char) :: sname, lname
     LOGICAL(KIND=BGC_log) :: &
        Nfixer,                             & ! flag set to true if this autotroph fixes N2
        imp_calcifier,                      & ! flag set to true if this autotroph implicitly handles calcification
        exp_calcifier                         ! flag set to true if this autotroph explicitly handles calcification
     INTEGER (KIND=BGC_i4) :: &
        grazee_ind,                         & ! which grazee category does autotroph belong to
        temp_function,                      & ! functional form of temperature parameterization 
        Chl_ind, C_ind, Fe_ind,             & ! tracer indices for Chl, C, Fe content
        Si_ind, CaCO3_ind                     ! tracer indices for Si, CaCO3 content
     REAL(KIND=BGC_r8) :: &
        kFe, kPO4, kDOP, kNO3, kNH4, kSiO3, & ! nutrient uptake half-sat constants
        Qp,                                 & ! P/C ratio
        gQfe_0, gQfe_min,                   & ! initial and minimum fe/C ratio
        alphaPI,                            & ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
        PCref,                              & ! max C-spec. grth rate at tref (1/sec)
        thetaN_max,                         & ! max thetaN (Chl/N) (mg Chl/mmol N)
        loss_thres, loss_thres2,            & ! conc. where losses go to zero
        temp_thres,temp_thresS,temp_thresN, & ! Temp. where concentration threshold and photosynth. rate drops
        temp_optN, temp_optS,               & ! Temp. function optimal value
        mort, mort2,                        & ! linear and quadratic mortality rates (1/sec), (1/sec/((mmol C/m3))
        agg_rate_max, agg_rate_min,         & ! max and min agg. rate (1/d)
        z_umax_0,                           & ! max zoo growth rate at tref (1/sec)
        z_grz,                              & ! grazing coef. (mmol C/m^3)
        graze_zoo, graze_poc, graze_doc,    & ! routing of grazed term, remainder goes to dic
        loss_poc,                           & ! routing of loss term
        f_zoo_detr                            ! fraction of zoo losses to detrital
  END TYPE

  type, public :: BGC_indices_type
   integer (BGC_i4) :: &
      po4_ind,          & ! dissolved inorganic phosphate
      no3_ind,          & ! dissolved inorganic nitrate
      sio3_ind,         & ! dissolved inorganic silicate
      nh4_ind,          & ! dissolved ammonia
      fe_ind,           & ! dissolved inorganic iron
      o2_ind,           & ! dissolved oxygen
      dic_ind,          & ! dissolved inorganic carbon
      dic_alt_co2_ind,  & ! dissolved inorganic carbon with alternative CO2
      alk_ind,          & ! alkalinity
      doc_ind,          & ! dissolved organic carbon
      don_ind,          & ! dissolved organic nitrogen
      dofe_ind,         & ! dissolved organic iron
      dop_ind,          & ! dissolved organic phosphorus
      dopr_ind,         & ! refractory DOP
      donr_ind,         & ! refractory DON
      zooC_ind,         & ! zooplankton carbon
      spC_ind,          & ! small phyto carbon
      spChl_ind,        & ! small phyto chlorophyll
      spFe_ind,         & ! small phyto iron
      spCaCO3_ind,      & ! small phyto calcium carbonate
      diatC_ind,        & ! diatom carbon
      diatChl_ind,      & ! diatom chlorophyll
      diatFe_ind,       & ! diatom iron
      diatSi_ind,       & ! diatom silicon
      phaeoC_ind,       & ! Phaeocystis carbon
      phaeoChl_ind,     & ! Phaeocystis chlorophyll
      phaeoFe_ind,      & ! Phaeocystis iron
      diazC_ind,        & ! diazotroph carbon
      diazChl_ind,      & ! diazotroph chlorophyll
      diazFe_ind          ! diazotroph iron

   integer (BGC_i4) :: &
      sp_ind,          & ! autotroph index for small phyto
      diat_ind,        & ! autotroph index for diatoms
      diaz_ind,        & ! autotroph index for diazotrophs
      phaeo_ind          ! autotroph index for Phaeocystis

   character (BGC_char), allocatable, dimension(:) ::  &
      short_name,      & ! short name of variable
      long_name,       & ! long name of variable
      units              ! units of variable

  end type BGC_indices_type

  type, public :: BGC_input_type
   real (BGC_r8), allocatable, dimension(:,:,:) :: &
      BGC_tracers
   real (BGC_r8), allocatable, dimension(:,:) :: &
      PotentialTemperature, Salinity,            &
      cell_center_depth, cell_thickness, cell_bottom_depth
   real (BGC_r8), allocatable, dimension(:) :: &
      cell_latitude
   integer (BGC_i4), allocatable, dimension(:) ::  &
      number_of_active_levels
  end type BGC_input_type

  type, public :: BGC_forcing_type
   real (BGC_r8), allocatable, dimension(:,:) :: &
      FESEDFLUX, NUTR_RESTORE_RTAU,              &
      NO3_CLIM, PO4_CLIM, SiO3_CLIM
   real (BGC_r8), allocatable, dimension(:) :: &
      dust_FLUX_IN, &
      ShortWaveFlux_surface,     &
      surfacePressure,           &
      iceFraction,               &
      windSpeedSquared10m,       &
      atmCO2,                    &
      atmCO2_ALT_CO2,            &
      surface_pH,                &
      surface_pH_alt_co2,        &
      surfaceDepth,              &
      SST,                       &
      SSS
   real (BGC_r8), allocatable, dimension(:,:) :: &
      depositionFlux, &
      riverFlux, &
      gasFlux, &
      seaIceFlux, &
      netFlux
   logical (BGC_log) :: &
      lcalc_O2_gas_flux,   &
      lcalc_CO2_gas_flux
  end type BGC_forcing_type

  type, public :: BGC_output_type
   real (BGC_r8), allocatable, dimension(:,:,:) :: &
      BGC_tendencies
   real (BGC_r8), allocatable, dimension(:,:) :: &
      PH_PREV_3D, PH_PREV_ALT_CO2_3D
  end type BGC_output_type

  type, public :: BGC_flux_diagnostics_type
   real (BGC_r8), allocatable, dimension(:) :: &
      pistonVel_O2,        &
      SCHMIDT_O2,          &
      O2SAT,               &
      xkw,                 &
      co2star,             &
      dco2star,            &
      pco2surf,            &
      dpco2,               &
      pistonVel_CO2,       &
      SCHMIDT_CO2,         &
      co2star_alt_co2,     &
      dco2star_alt_co2,    &
      pco2surf_alt_co2,    &
      dpco2_alt_co2
  end type BGC_flux_diagnostics_type

  type, public :: BGC_diagnostics_type
! 3D stuff
    real (BGC_r8), allocatable, dimension(:,:) :: &
      diag_tot_Nfix,     &! diag array for total N fixation
      diag_O2_PRODUCTION,&! diag array for o2 production
      diag_O2_CONSUMPTION,&! diag array for o2 consumption
      diag_AOU,          &! diag array for AOU
      diag_PO4_RESTORE,  &! diag array for po4 restoring
      diag_NO3_RESTORE,  &! diag array for no3 restoring
      diag_SiO3_RESTORE, &! diag array for sio3 restoring
      diag_PAR_avg,      &! diag array for available radiation avg over mixed layer
      diag_POC_FLUX_IN,  &! diag array for poc flux into cell
      diag_POC_PROD,     &! diag array for poc production
      diag_POC_REMIN,    &! diag array for poc remineralization
      diag_POC_ACCUM,    &! diag array for poc accumulation
      diag_CaCO3_FLUX_IN,&! diag array for caco3 flux into cell
      diag_CaCO3_PROD,   &! diag array for caco3 production
      diag_CaCO3_REMIN,  &! diag array for caco3 remineralization
      diag_SiO2_FLUX_IN, &! diag array for sio2 flux into cell
      diag_SiO2_PROD,    &! diag array for sio2 production
      diag_SiO2_REMIN,   &! diag array for sio2 remineralization
      diag_dust_FLUX_IN, &! diag array for dust flux into cell
      diag_dust_REMIN,   &! diag array for dust remineralization
      diag_P_iron_FLUX_IN, &! diag array for p_iron flux into cell
      diag_P_iron_PROD,    &! diag array for p_iron production
      diag_P_iron_REMIN,   &! diag array for p_iron remineralization
      diag_auto_graze_TOT, &! diag array for total autotroph grazing
      diag_zoo_loss,       &! diag array for zooplankton loss
      diag_photoC_TOT,     &! diag array for total C fixation
      diag_photoC_NO3_TOT, &! diag array for total C fixation from NO3
      diag_DOC_prod,       &! diag array for doc production
      diag_DOC_remin,      &! diag array for doc remineralization
      diag_DON_prod,       &! diag array for don production
      diag_DON_remin,      &! diag array for don remineralization
      diag_DOFe_prod,      &! diag array for dofe production
      diag_DOFe_remin,     &! diag array for dofe remineralization
      diag_DOP_prod,       &! diag array for dop production
      diag_DOP_remin,      &! diag array for dop remineralization
      diag_Fe_scavenge,    &! diag array for iron scavenging
      diag_Fe_scavenge_rate,   &! diag array for iron scavenging rate
      diag_NITRIF,         &! diag array for nitrification
      diag_DENITRIF,       &! diag array for denitrification
      diag_DONr_remin,     &! diag array for DONrefractory remin
      diag_DOPr_remin       ! diag array for DOPrefractory remin

! more 3D stuff
    real (BGC_r8), allocatable, dimension(:,:) :: &
      diag_CO3,            &! diag array for 3D carbonate ion
      diag_HCO3,           &! diag array for 3D bicarbonate ion
      diag_H2CO3,          &! diag array for 3D carbonic acid
      diag_pH_3D,          &! diag array for 3D pH
      diag_CO3_ALT_CO2,    &! diag array for 3D carbonate ion, alternative CO2
      diag_HCO3_ALT_CO2,   &! diag array for 3D bicarbonate ion, alternative CO2
      diag_H2CO3_ALT_CO2,  &! diag array for 3D carbonic acid, alternative CO2
      diag_pH_3D_ALT_CO2,  &! diag array for 3D pH, alternative CO2
      diag_co3_sat_calc,   &! diag array for co3 concentration at calcite saturation
      diag_co3_sat_arag,   &! diag array for co3 concentration at aragonite saturation
      diag_calcToSed,      &! diag array for calcite flux sedimentary burial
      diag_pocToSed,       &! diag array for poc burial flux to sediments
      diag_ponToSed,       &! diag array for pon burial flux to sediments
      diag_popToSed,       &! diag array for pop burial flux to sediments
      diag_bsiToSed,       &! diag array for bsi burial flux to sediments
      diag_dustToSed,      &! diag array for dust burial flux to sediments
      diag_pfeToSed,       &! diag array for pFe burial flux to sediments
      diag_SedDenitrif,    &! diag array for sedimentary denitrification
      diag_OtherRemin,     &! diag array for non-oxic, non-denitr sed remin
      diag_tot_CaCO3_form   ! diag array for CaCO3 formation

! 3D stuff for each autotroph
! dim (nx,ny,nz,autotroph_cnt)
    real (BGC_r8), allocatable, dimension(:,:,:) :: &
      diag_N_lim,          &! diag array for N limitation
      diag_P_lim,          &! diag array for P limitation
      diag_Fe_lim,         &! diag array for Fe limitation
      diag_SiO3_lim,       &! diag array for SiO3 limitation
      diag_light_lim,      &! diag array for light limitation
      diag_photoC,         &! diag array for C fixation
      diag_photoC_NO3,     &! diag array for C fixation from NO3
      diag_photoFe,        &! diag array for Fe uptake
      diag_photoNO3,       &! diag array for NO3 uptake
      diag_photoNH4,       &! diag array for NH4 uptake
      diag_DOP_uptake,     &! diag array for DOP uptake
      diag_PO4_uptake,     &! diag array for PO4 uptake
      diag_auto_graze,     &! diag array for autotroph grazing
      diag_auto_loss,      &! diag array for autotroph loss
      diag_auto_agg,       &! diag array for autotroph aggregate
      diag_bSi_form,       &! diag array for Si uptake
      diag_CaCO3_form,     &! diag array for CaCO3 formation
      diag_Nfix            ! diag array for N fixation

! 2D stuff for each autotroph
! dim (nx,ny,autotroph_cnt)
    real (BGC_r8), allocatable, dimension(:,:) :: &
      diag_photoC_zint,    &! diag array for C fixation vertical integral
      diag_photoC_NO3_zint,&! diag array for C fixation from NO3 vertical integral
      diag_CaCO3_form_zint  ! diag array for CaCO3 formation vertical integral

! 2D vertical integral of photoC
! dim (nx,ny)
    real (BGC_r8), allocatable, dimension(:) :: &
      diag_photoC_TOT_zint,      &!
      diag_photoC_NO3_TOT_zint

! 2D vertical integral of conservative subterms of source sink term for nutrients
! dim (nx,ny)
    real (BGC_r8), allocatable, dimension(:) :: &
      diag_Jint_Ctot,      &! diag array for Ctot
      diag_Jint_100m_Ctot, &! diag array for Ctot, 0-100m
      diag_Jint_Ntot,      &! diag array for Ntot
      diag_Jint_100m_Ntot, &! diag array for Ntot, 0-100m
      diag_Jint_Ptot,      &! diag array for Ptot
      diag_Jint_100m_Ptot, &! diag array for Ptot, 0-100m
      diag_Jint_Sitot,     &! diag array for Sitot
      diag_Jint_100m_Sitot  ! diag array for Sitot, 0-100m

! 2D vertical integral of total chlorophyll over top 100m
! dim (nx,ny)
    real (BGC_r8), allocatable, dimension(:) :: &
      diag_Chl_TOT_zint_100m

! 2D stuff
! dim (nx,ny)
    real (BGC_r8), allocatable, dimension(:) :: &
      diag_tot_CaCO3_form_zint,&! diag array for CaCO3 formation vertical integral
      diag_tot_bSi_form,       &! diag array for Si uptake
      diag_zsatcalc,           &! diag array for calcite saturation depth
      diag_zsatarag,           &! diag array for aragonite saturation depth
      diag_O2_ZMIN,            &! diag array for O2 minimum
      diag_O2_ZMIN_DEPTH        ! diag array for depth of O2 minimum
  end type BGC_diagnostics_type

  !-----------------------------------------------------------------------------
  !   Redfield Ratios, dissolved & particulate
  !-----------------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: &
       parm_Red_D_C_P  = 117.0_BGC_r8,                 & ! carbon:phosphorus
       parm_Red_D_N_P  =  16.0_BGC_r8,                 & ! nitrogen:phosphorus
       parm_Red_D_O2_P = 170.0_BGC_r8,                 & ! oxygen:phosphorus
       parm_Remin_D_O2_P = 138.0_BGC_r8,               & ! oxygen:phosphorus
       parm_Red_P_C_P  = parm_Red_D_C_P,                 & ! carbon:phosphorus
       parm_Red_D_C_N  = parm_Red_D_C_P/parm_Red_D_N_P,  & ! carbon:nitrogen
       parm_Red_P_C_N  = parm_Red_D_C_N,                 & ! carbon:nitrogen
       parm_Red_D_C_O2 = parm_Red_D_C_P/parm_Red_D_O2_P, & ! carbon:oxygen
       parm_Remin_D_C_O2 = parm_Red_D_C_P/parm_Remin_D_O2_P, & ! carbon:oxygen
       parm_Red_P_C_O2 = parm_Red_D_C_O2,                & ! carbon:oxygen
       parm_Red_Fe_C   = 3.0e-6_BGC_r8,                & ! iron:carbon
       parm_Red_D_C_O2_diaz = parm_Red_D_C_P/150.0_BGC_r8! carbon:oxygen
                                                           ! for diazotrophs

  !----------------------------------------------------------------------------
  !   ecosystem parameters accessible via namelist input
  !----------------------------------------------------------------------------

  REAL(KIND=BGC_r8) :: &
       parm_Fe_bioavail,      & ! fraction of Fe flux that is bioavailable
       parm_o2_min,           & ! min O2 needed for prod & consump. (nmol/cm^3)
       parm_o2_min_delta,     & ! width of min O2 range (nmol/cm^3)
       parm_kappa_nitrif,     & ! nitrification inverse time constant (1/sec)
       parm_nitrif_par_lim,   & ! PAR limit for nitrif. (W/m^2)
       parm_z_mort_0,         & ! zoo linear mort rate (1/sec)
       parm_z_mort2_0,        & ! zoo quad mort rate (1/sec/((mmol C/m3))
       parm_labile_ratio,     & ! fraction of loss to DOC that routed directly to DIC (non-dimensional)
       parm_POMbury,          & ! scale factor for burial of POC, PON, and POP
       parm_BSIbury,          & ! scale factor burial of bSi
       parm_fe_scavenge_rate0,& ! base scavenging rate
       parm_f_prod_sp_CaCO3,  & !fraction of sp prod. as CaCO3 prod.
       parm_POC_diss,         & ! base POC diss len scale
       parm_SiO2_diss,        & ! base SiO2 diss len scale
       parm_CaCO3_diss          ! base CaCO3 diss len scale

  REAL(KIND=BGC_r8), DIMENSION(4) :: &
       parm_scalelen_z,       & ! depths of prescribed scalelen values
       parm_scalelen_vals       ! prescribed scalelen values

  !---------------------------------------------------------------------
  !     Misc. Rate constants
  !---------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: &
       fe_scavenge_thres1 = 0.8e-3_BGC_r8,  & !upper thres. for Fe scavenging
       dust_fescav_scale  = 1.0e9,      & !dust scavenging scale factor
       fe_max_scale2      = 1200.0_BGC_r8     !unitless scaling coeff.

  !---------------------------------------------------------------------
  !     Compute iron remineralization and flux out.
  !     dust remin gDust = 0.035 gFe      mol Fe     1e9 nmolFe
  !                        --------- *  ---------- * ----------
  !			    gDust       55.847 gFe     molFe
  !
  !     dust_to_Fe          conversion - dust to iron (nmol Fe/g Dust) 
  !---------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: &
       dust_to_Fe=0.035_BGC_r8/55.847_BGC_r8*1.0e9_BGC_r8
 
  !----------------------------------------------------------------------------
  !     Partitioning of phytoplankton growth, grazing and losses
  !
  !     All f_* variables are fractions and are non-dimensional
  !----------------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: &
      caco3_poc_min    = 0.4_BGC_r8,  & !minimum proportionality between 
                                          !   QCaCO3 and grazing losses to POC 
                                          !   (mmol C/mmol CaCO3)
      spc_poc_fac      = 0.11_BGC_r8, & !small phyto grazing factor (1/mmolC)
      f_graze_sp_poc_lim = 0.3_BGC_r8, & 
      f_photosp_CaCO3  = 0.4_BGC_r8,  & !proportionality between small phyto 
                                          !    production and CaCO3 production
      f_graze_CaCO3_remin = 0.33_BGC_r8, & !fraction of spCaCO3 grazing 
                                             !          which is remin
      f_graze_si_remin    = 0.35_BGC_r8      !fraction of diatom Si grazing 
                                             !          which is remin

  !----------------------------------------------------------------------------
  !     fixed ratios
  !----------------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: &
       r_Nfix_photo=1.25_BGC_r8         ! N fix relative to C fix (non-dim)

  !-----------------------------------------------------------------------
  !     SET FIXED RATIOS for N/C, P/C, SiO3/C, Fe/C
  !     assumes C/N/P of 117/16/1 based on Anderson and Sarmiento, 1994
  !     for diazotrophs a N/P of 45 is assumed based on Letelier & Karl, 1998
  !-----------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: &
      Q             = 0.137_BGC_r8,  & !N/C ratio (mmol/mmol) of phyto & zoo
      Qp_zoo_pom    = 0.00855_BGC_r8,& !P/C ratio (mmol/mmol) zoo & pom
      Qfe_zoo       = 3.0e-6_BGC_r8, & !zooplankton fe/C ratio
      gQsi_0        = 0.137_BGC_r8,  & !initial Si/C ratio
      gQsi_max      = 0.685_BGC_r8,  & !max Si/C ratio
      gQsi_min      = 0.0457_BGC_r8, & !min Si/C ratio
      QCaCO3_max    = 0.4_BGC_r8,    & !max QCaCO3
      ! carbon:nitrogen ratio for denitrification
      denitrif_C_N  = parm_Red_D_C_P/136.0_BGC_r8

  !----------------------------------------------------------------------------
  !     loss term threshold parameters, chl:c ratios
  !----------------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: &
      thres_z1          = 100.0e2_BGC_r8, & !threshold = C_loss_thres for z shallower than this (cm)
      thres_z2          = 150.0e2_BGC_r8, & !threshold = 0 for z deeper than this (cm)
      loss_thres_zoo    = 0.005_BGC_r8,  & !zoo conc. where losses go to zero
      CaCO3_temp_thres1 = 6.0_BGC_r8,   & !upper temp threshold for CaCO3 prod
      CaCO3_temp_thres2 = -2.0_BGC_r8,  & !lower temp threshold
      CaCO3_sp_thres    = 4.0_BGC_r8      ! bloom condition thres (mmolC/m3)

  !---------------------------------------------------------------------
  !     temperature functions
  !---------------------------------------------------------------------

  INTEGER (KIND=BGC_i4), PARAMETER ::   &
         tfnc_q10				  = 1,		 & ! classic Q10, no decline
         tfnc_quasi_mmrt 		  = 2  		   ! consider Topt, mimic macromolecular rate theory
  !---------------------------------------------------------------------
  !     fraction of incoming shortwave assumed to be PAR
  !---------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: &
       f_qsw_par = 0.45_BGC_r8   ! PAR fraction

  !---------------------------------------------------------------------
  !     Temperature parameters
  !---------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: &
       Tref = 30.0_BGC_r8, & ! reference temperature (C)
       Q_10 = 1.5_BGC_r8     ! factor for temperature dependence (non-dim)

  !---------------------------------------------------------------------
  !  DOM parameters for refractory components and DOP uptake
  !---------------------------------------------------------------------

  REAL(KIND=BGC_r8), PARAMETER :: &
       DOC_reminR  = (c1/250.0_BGC_r8) * dps,          & ! rate for semi-labile DOC 1/250days
       DON_reminR  = (c1/160.0_BGC_r8) * dps,          & ! rate for semi-labile DON 1/160days
       DOFe_reminR = (c1/160.0_BGC_r8) * dps,          & ! rate for semi-labile DOFe 1/160days
       DOP_reminR  = (c1/160.0_BGC_r8) * dps,          & ! rate for semi-labile DOP 1/160days  
       DONr_reminR = (c1/(365.0_BGC_r8*2.5_BGC_r8)) * dps, & ! timescale for refrac DON 1/2.5yrs
       DOPr_reminR = (c1/(365.0_BGC_r8*2.5_BGC_r8)) * dps, & ! timescale for refrac DOP 1/2.5yrs
       DONrefract = 0.08_BGC_r8,                       & ! fraction of DON to refractory pool
       DOPrefract = 0.03_BGC_r8                          ! fraction of DOP to refractory pool

   real (BGC_r8), parameter :: &
      epsC      = 1.00e-8, & ! small C concentration (mmol C/m^3)
      epsTinv   = 3.17e-8, & ! small inverse time scale (1/year) (1/sec)
      epsnondim = 1.00e-6    ! small non-dimensional number (non-dim)

   real (BGC_r8), parameter :: &
      cks  = 9.,  & ! constant used in Fe quota modification
      cksi = 5.     ! constant used in Si quota modification

   real (BGC_r8), parameter, public :: &
      xkw_coeff = 8.6e-9_BGC_r8     ! a = 0.31 cm/hr s^2/m^2 in (s/cm)

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  subroutine BGC_parms_init(BGC_indices, autotrophs)

! !INPUT/OUTPUT PARAMETERS:

  type(autotroph_type), dimension(autotroph_cnt), intent(inout) :: autotrophs

  type(BGC_indices_type), intent(inout) :: BGC_indices

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    integer (BGC_i4) :: auto_ind

    !---------------------------------------------------------------------------
    !  set specific autotroph indices
    !---------------------------------------------------------------------------

    BGC_indices%sp_ind   = 1
    BGC_indices%diat_ind = 2
    BGC_indices%diaz_ind = 3
    BGC_indices%phaeo_ind = 4

    !---------------------------------------------------------------------------
    !   default namelist settings
    !---------------------------------------------------------------------------

    parm_Fe_bioavail       = 1.0_BGC_r8
    parm_o2_min            = 4.0_BGC_r8
    parm_o2_min_delta      = 2.0_BGC_r8
    parm_kappa_nitrif      = 0.06_BGC_r8 * dps  ! (= 1/( days))
    parm_nitrif_par_lim    = 1.0_BGC_r8
    parm_z_mort_0          = 0.1_BGC_r8 * dps
    parm_z_mort2_0         = 0.4_BGC_r8 * dps
    parm_labile_ratio      = 0.85_BGC_r8
    parm_POMbury           = 1.4_BGC_r8         ! x1 default
    parm_BSIbury           = 0.65_BGC_r8        ! x1 default
    parm_fe_scavenge_rate0 = 3.0_BGC_r8         ! x1 default
    parm_f_prod_sp_CaCO3   = 0.055_BGC_r8       ! x1 default
    parm_POC_diss          = 88.0e2_BGC_r8
    parm_SiO2_diss         = 250.0e2_BGC_r8
    parm_CaCO3_diss        = 150.0e2_BGC_r8

    parm_scalelen_z    = (/ 130.0e2_BGC_r8, 290.0e2_BGC_r8, 670.0e2_BGC_r8, 1700.0e2_BGC_r8 /) ! x1 default
    parm_scalelen_vals = (/     1.0_BGC_r8,     3.0_BGC_r8,     5.0_BGC_r8,      9.0_BGC_r8 /) ! x1 default

    auto_ind = BGC_indices%sp_ind
    autotrophs(auto_ind)%sname         = 'sp'
    autotrophs(auto_ind)%lname         = 'Small Phyto'
    autotrophs(auto_ind)%Nfixer        = .false.
    autotrophs(auto_ind)%imp_calcifier = .true.
    autotrophs(auto_ind)%exp_calcifier = .false.
    autotrophs(auto_ind)%grazee_ind    = auto_ind
    autotrophs(auto_ind)%kFe           = 0.04e-3_BGC_r8
    autotrophs(auto_ind)%kPO4          = 0.01_BGC_r8
    autotrophs(auto_ind)%kDOP          = 0.26_BGC_r8
    autotrophs(auto_ind)%kNO3          = 0.1_BGC_r8
    autotrophs(auto_ind)%kNH4          = 0.01_BGC_r8
    autotrophs(auto_ind)%kSiO3         = 0.0_BGC_r8
    autotrophs(auto_ind)%Qp            = 0.00855_BGC_r8
    autotrophs(auto_ind)%gQfe_0        = 20.0e-6_BGC_r8
    autotrophs(auto_ind)%gQfe_min      = 3.0e-6_BGC_r8
    autotrophs(auto_ind)%alphaPI       = 0.6_BGC_r8 * dps
    autotrophs(auto_ind)%PCref         = 5.5_BGC_r8 * dps
    autotrophs(auto_ind)%thetaN_max    = 2.5_BGC_r8
    autotrophs(auto_ind)%loss_thres    = 0.04_BGC_r8
    autotrophs(auto_ind)%loss_thres2   = 0.0_BGC_r8
    autotrophs(auto_ind)%temp_thres    = -20.0_BGC_r8
    autotrophs(auto_ind)%temp_thresN   = -20.0_BGC_r8
    autotrophs(auto_ind)%temp_thresS   = -20.0_BGC_r8
    autotrophs(auto_ind)%temp_function = tfnc_q10
    autotrophs(auto_ind)%temp_optN     = 50.0_BGC_r8
    autotrophs(auto_ind)%temp_optS     = 50.0_BGC_r8
    autotrophs(auto_ind)%mort          = 0.12_BGC_r8 * dps
    autotrophs(auto_ind)%mort2         = 0.001_BGC_r8 * dps
    autotrophs(auto_ind)%agg_rate_max  = 0.9_BGC_r8
    autotrophs(auto_ind)%agg_rate_min  = 0.01_BGC_r8
    autotrophs(auto_ind)%z_umax_0      = 3.3_BGC_r8 * dps ! x1 default
    autotrophs(auto_ind)%z_grz         = 1.05_BGC_r8              
    autotrophs(auto_ind)%graze_zoo     = 0.3_BGC_r8
    autotrophs(auto_ind)%graze_poc     = 0.0_BGC_r8
    autotrophs(auto_ind)%graze_doc     = 0.15_BGC_r8
    autotrophs(auto_ind)%loss_poc      = 0.0_BGC_r8
    autotrophs(auto_ind)%f_zoo_detr    = 0.15_BGC_r8

    auto_ind = BGC_indices%diat_ind
    autotrophs(auto_ind)%sname         = 'diat'
    autotrophs(auto_ind)%lname         = 'Diatom'
    autotrophs(auto_ind)%Nfixer        = .false.
    autotrophs(auto_ind)%imp_calcifier = .false.
    autotrophs(auto_ind)%exp_calcifier = .false.
    autotrophs(auto_ind)%grazee_ind    = auto_ind
    autotrophs(auto_ind)%kFe           = 0.06e-3_BGC_r8
    autotrophs(auto_ind)%kPO4          = 0.05_BGC_r8
    autotrophs(auto_ind)%kDOP          = 0.9_BGC_r8
    autotrophs(auto_ind)%kNO3          = 0.5_BGC_r8
    autotrophs(auto_ind)%kNH4          = 0.05_BGC_r8
    autotrophs(auto_ind)%kSiO3         = 0.8_BGC_r8
    autotrophs(auto_ind)%Qp            = 0.00855_BGC_r8
    autotrophs(auto_ind)%gQfe_0        = 20.0e-6_BGC_r8
    autotrophs(auto_ind)%gQfe_min      = 3.0e-6_BGC_r8
    autotrophs(auto_ind)%alphaPI       = 0.465_BGC_r8 * dps
    autotrophs(auto_ind)%PCref         = 5.5_BGC_r8 * dps
    autotrophs(auto_ind)%thetaN_max    = 4.0_BGC_r8
    autotrophs(auto_ind)%loss_thres    = 0.04_BGC_r8
    autotrophs(auto_ind)%loss_thres2   = 0.0_BGC_r8
    autotrophs(auto_ind)%temp_thres    = -20.0_BGC_r8
    autotrophs(auto_ind)%temp_thresN   = 35.0_BGC_r8
    autotrophs(auto_ind)%temp_thresS   = 10.0_BGC_r8
    autotrophs(auto_ind)%temp_function = tfnc_q10
    autotrophs(auto_ind)%temp_optN     = 16.3_BGC_r8
    autotrophs(auto_ind)%temp_optS     = 5.0_BGC_r8
    autotrophs(auto_ind)%mort          = 0.12_BGC_r8 * dps
    autotrophs(auto_ind)%mort2         = 0.001_BGC_r8 * dps
    autotrophs(auto_ind)%agg_rate_max  = 0.9_BGC_r8
    autotrophs(auto_ind)%agg_rate_min  = 0.02_BGC_r8
    autotrophs(auto_ind)%z_umax_0      = 3.23_BGC_r8 * dps ! x1 default
    autotrophs(auto_ind)%z_grz         = 1.0_BGC_r8              
    autotrophs(auto_ind)%graze_zoo     = 0.3_BGC_r8
    autotrophs(auto_ind)%graze_poc     = 0.42_BGC_r8
    autotrophs(auto_ind)%graze_doc     = 0.15_BGC_r8
    autotrophs(auto_ind)%loss_poc      = 0.0_BGC_r8
    autotrophs(auto_ind)%f_zoo_detr    = 0.2_BGC_r8

    auto_ind = BGC_indices%diaz_ind
    autotrophs(auto_ind)%sname         = 'diaz'
    autotrophs(auto_ind)%lname         = 'Diazotroph'
    autotrophs(auto_ind)%Nfixer        = .true.
    autotrophs(auto_ind)%imp_calcifier = .false.
    autotrophs(auto_ind)%exp_calcifier = .false.
    autotrophs(auto_ind)%grazee_ind    = auto_ind
    autotrophs(auto_ind)%kFe           = 0.04e-3_BGC_r8
    autotrophs(auto_ind)%kPO4          = 0.02_BGC_r8
    autotrophs(auto_ind)%kDOP          = 0.09_BGC_r8
    autotrophs(auto_ind)%kNO3          = 1.0_BGC_r8
    autotrophs(auto_ind)%kNH4          = 0.15_BGC_r8
    autotrophs(auto_ind)%kSiO3         = 0.0_BGC_r8
    autotrophs(auto_ind)%Qp            = 0.002735_BGC_r8
    autotrophs(auto_ind)%gQfe_0        = 60.0e-6_BGC_r8
    autotrophs(auto_ind)%gQfe_min      = 12.0e-6_BGC_r8
    autotrophs(auto_ind)%alphaPI       = 0.4_BGC_r8 * dps
    autotrophs(auto_ind)%PCref         = 0.7_BGC_r8 * dps
    autotrophs(auto_ind)%thetaN_max    = 2.5_BGC_r8
    autotrophs(auto_ind)%loss_thres    = 0.022_BGC_r8
    autotrophs(auto_ind)%loss_thres2   = 0.001_BGC_r8
    autotrophs(auto_ind)%temp_thres    = 14.0_BGC_r8
    autotrophs(auto_ind)%temp_thresN   = -20.0_BGC_r8
    autotrophs(auto_ind)%temp_thresS   = -20.0_BGC_r8
    autotrophs(auto_ind)%temp_function = tfnc_q10
    autotrophs(auto_ind)%temp_optN     = 50.0_BGC_r8
    autotrophs(auto_ind)%temp_optS     = 50.0_BGC_r8
    autotrophs(auto_ind)%mort          = 0.15_BGC_r8 * dps
    autotrophs(auto_ind)%mort2         = 0.0_BGC_r8
    autotrophs(auto_ind)%agg_rate_max  = 0.0_BGC_r8
    autotrophs(auto_ind)%agg_rate_min  = 0.0_BGC_r8
    autotrophs(auto_ind)%z_umax_0      = 0.6_BGC_r8 * dps
    autotrophs(auto_ind)%z_grz         = 1.2_BGC_r8              
    autotrophs(auto_ind)%graze_zoo     = 0.3_BGC_r8
    autotrophs(auto_ind)%graze_poc     = 0.05_BGC_r8
    autotrophs(auto_ind)%graze_doc     = 0.15_BGC_r8
    autotrophs(auto_ind)%loss_poc      = 0.0_BGC_r8
    autotrophs(auto_ind)%f_zoo_detr    = 0.15_BGC_r8

    auto_ind = BGC_indices%phaeo_ind
    autotrophs(auto_ind)%sname         = 'phaeo'
    autotrophs(auto_ind)%lname         = 'Phaeocystis'
    autotrophs(auto_ind)%Nfixer        = .false.
    autotrophs(auto_ind)%imp_calcifier = .false.
    autotrophs(auto_ind)%exp_calcifier = .false.
    autotrophs(auto_ind)%grazee_ind    = BGC_indices%diat_ind
    autotrophs(auto_ind)%kFe           = 0.075e-3_BGC_r8
    autotrophs(auto_ind)%kPO4          = 0.05_BGC_r8
    autotrophs(auto_ind)%kDOP          = 0.9_BGC_r8
    autotrophs(auto_ind)%kNO3          = 0.7_BGC_r8
    autotrophs(auto_ind)%kNH4          = 0.05_BGC_r8
    autotrophs(auto_ind)%kSiO3         = 0.0_BGC_r8
    autotrophs(auto_ind)%Qp            = 0.00855_BGC_r8
    autotrophs(auto_ind)%gQfe_0        = 20.0e-6_BGC_r8
    autotrophs(auto_ind)%gQfe_min      = 3.0e-6_BGC_r8
    autotrophs(auto_ind)%alphaPI       = 0.77_BGC_r8 * dps
    autotrophs(auto_ind)%PCref         = 5.5_BGC_r8 * dps
    autotrophs(auto_ind)%thetaN_max    = 2.5_BGC_r8
    autotrophs(auto_ind)%loss_thres    = 0.04_BGC_r8
    autotrophs(auto_ind)%loss_thres2   = 0.0_BGC_r8
    autotrophs(auto_ind)%temp_thres    = -20.0_BGC_r8
    autotrophs(auto_ind)%temp_thresN   = 35.0_BGC_r8
    autotrophs(auto_ind)%temp_thresS   = 10.0_BGC_r8
    autotrophs(auto_ind)%temp_function = tfnc_quasi_mmrt 
    autotrophs(auto_ind)%temp_optN     = 16.3_BGC_r8
    autotrophs(auto_ind)%temp_optS     = 5.0_BGC_r8
    autotrophs(auto_ind)%mort          = 0.12_BGC_r8 * dps
    autotrophs(auto_ind)%mort2         = 0.001_BGC_r8 * dps
    autotrophs(auto_ind)%agg_rate_max  = 0.9_BGC_r8
    autotrophs(auto_ind)%agg_rate_min  = 0.02_BGC_r8
    autotrophs(auto_ind)%z_umax_0      = 3.23_BGC_r8 * dps ! x1 default
    autotrophs(auto_ind)%z_grz         = 1.0_BGC_r8
    autotrophs(auto_ind)%graze_zoo     = 0.3_BGC_r8
    autotrophs(auto_ind)%graze_poc     = 0.42_BGC_r8
    autotrophs(auto_ind)%graze_doc     = 0.15_BGC_r8
    autotrophs(auto_ind)%loss_poc      = 0.0_BGC_r8
    autotrophs(auto_ind)%f_zoo_detr    = 0.2_BGC_r8

  END SUBROUTINE BGC_parms_init

  !*****************************************************************************

END MODULE BGC_parms
