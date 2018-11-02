!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module BGC_mod

!BOP
! !MODULE: BGC_mod
!
! !DESCRIPTION:
!
!  Multispecies ecosystem based on Doney et al. 1996, Moore et al., 2002
!  Based on POP Global NCAR Nitrogen Ecosystem Model
!  version 0.0 (June 15th, 1998) from S.C. Doney.
!  Based on Doney et al., 1996 model.
!  Climate and Global Dynamics, NCAR
!  (doney@whoi.edu)
!
!  Version 1.0
!  Multispecies, multiple limiting nutrient version of ecosystem
!  based on mixed layer model of Moore et al.(2002).  Implemented here with
!  fixed elemental ratios and including only the diatoms and small
!  phytoplankton, with a parameterization of calcification,
!  by Keith Lindsay and Keith Moore, Fall 2001 - Spring 2002.
!  Calcification parameterization based on Moore et al. 2002.
!
!  Version 2.0, January 2003
!    Adds diazotrophs as a phytoplankton group, (based on Moore et al., 2002a)
!    Allows for variable fe/C for all phytoplankton groups
!     Allows for variable si/C for the diatoms
!     Adds explicit tracers for DON, DOP, DOFe
!     variable remin length scale for detrital soft POM and bSi f(temperature)
!     Extensive modifications to iron scavenging parameterization
!     Addition of a sedimentary dissolved iron source,
!        (implemented in ballast code as excess remin in bottom cell)
!        coded by J.K. Moore, (jkmoore@uci.edu)
!
!   Version 2.01. March 2003
!     corrected O2 bug
!     corrected grazing parameter z_grz bug at depth
!     dust dissolution at depth releases iron,
!     increased length scale for dust diss., increased hard fraction dust
!     no deep ocean reduction in scavenging rates,
!     increase bSi OC/ballast ratio 0.3 -> 0.35,
!     corrected bug in diazotroph photoadaptation, and diat and sp adapatation
!
!   Version 2.02.
!     corrected bug in Fe_scavenge (units for dust), May 2003
!     changed C/N/P ratios to 117/16/1 (Anderson & Sarmiento, 1994)
!
!   Version 2.03., July 2003
!     Remin of DOM no longer temperature dependent,
!     new iron scavenging parameterization added,
!     some dissolution of hard fraction of ballast materials added
!
!   Version 2.1, September 2003
!     modfied iron scavenging and dust dissolution at depth
!
!   Version 2.11, March 2004
!     fixed bug in iron scavenging code, replace dust and POC flux_in w/ flux_out
!
!   Version 2.12, April 2004 - Final version for GBC paper revision,
!     (Questions/comments, Keith Moore - jkmoore@uci.edu
!
!   References
!   Doney, S.C., Glover, D.M., Najjar, R.G., 1996. A new coupled, one-dimensional
!   biological-physical model for the upper ocean: applications to the JGOFS
!   Bermuda Time-Series Study (BATS) site. Deep-Sea Res. II, 43: 591-624.
!
!   Moore, JK, Doney, SC, Kleypas, JA, Glover, DM, Fung, IY, 2002. An intermediate
!   complexity marine ecosystem model for the global domain. Deep-Sea Res. II, 49:
!   403-462.
!
!   Moore, JK, Doney, SC, Glover, DM, Fung, IY, 2002. Iron cycling and nutrient
!   limitation patterns in surface waters of the world ocean. Deep-Sea Res. II,
!   49: 463-507.

! !REVISION HISTORY:

!  SVN:$Id:  $

!-----------------------------------------------------------------------
!  variables/subroutines/function used from other modules
!  The following are used extensively in this ecosys, so are used at
!  the module level. The use statements for variables that are only needed
!  locally are located at the module subprogram level.
!-----------------------------------------------------------------------

! !USES:

   use BGC_parms
   use co2calc

! !INPUT PARAMETERS:
!-----------------------------------------------------------------------
!  include ecosystem parameters
!  all variables from this modules have a parm_ prefix
!-----------------------------------------------------------------------

   implicit none
   save
   private

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   public :: &
      BGC_tracer_cnt,     &
      BGC_init,           &
      BGC_SurfaceFluxes,  &
      BGC_SourceSink

!-----------------------------------------------------------------------
!  module variables
!-----------------------------------------------------------------------

! increase tracer_cnt for phaeo (swang)
   integer (BGC_i4), parameter :: &
      BGC_tracer_cnt = 30

   integer (BGC_r8), parameter, private :: &
       c0 =   0.0_BGC_r8,         &
       c1 =   1.0_BGC_r8,         &
       c2 =   2.0_BGC_r8,         &
       c10 = 10.0_BGC_r8,         &
       p5  =  0.5_BGC_r8

!-----------------------------------------------------------------------
!  restoring climatologies for nutrients
!-----------------------------------------------------------------------

   logical (BGC_log) :: &
      lrest_po4,  & ! restoring on po4 ?
      lrest_no3,  & ! restoring on no3 ?
      lrest_sio3    ! restoring on sio3 ?

!  real (BGC_r8), dimension(ecosys_tracer_cnt) :: &
!     surf_avg                  ! average surface tracer values

   logical (BGC_log) :: &
      ecosys_qsw_distrb_const

!-----------------------------------------------------------------------

   real (BGC_r8), parameter :: &
      phlo_surf_init = 7.0_BGC_r8, & ! low bound for surface ph for no prev soln
      phhi_surf_init = 9.0_BGC_r8, & ! high bound for surface ph for no prev soln
      phlo_3d_init = 6.0_BGC_r8,   & ! low bound for subsurface ph for no prev soln
      phhi_3d_init = 9.0_BGC_r8,   & ! high bound for subsurface ph for no prev soln
      del_ph = 0.20_BGC_r8           ! delta-ph for prev soln

!-----------------------------------------------------------------------
!  derived type for implicit handling of sinking particulate matter
!-----------------------------------------------------------------------

   type sinking_particle
      real (BGC_r8) :: &
         diss,        & ! dissolution length for soft subclass
         gamma,       & ! fraction of production -> hard subclass
         mass,        & ! mass of 1e9 base units in g
         rho            ! QA mass ratio of POC to this particle class

      real (BGC_r8) :: &
         sflux_in,    & ! incoming flux of soft subclass (base units/cm^2/sec)
         hflux_in,    & ! incoming flux of hard subclass (base units/cm^2/sec)
         prod,        & ! production term (base units/cm^3/sec)
         sflux_out,   & ! outgoing flux of soft subclass (base units/cm^2/sec)
         hflux_out,   & ! outgoing flux of hard subclass (base units/cm^2/sec)
         sed_loss,    & ! loss to sediments (base units/cm^s/sec)
         remin          ! remineralization term (base units/cm^3/sec)
   end type sinking_particle

!-----------------------------------------------------------------------


!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: BGC_init
! !INTERFACE:

 subroutine BGC_init(BGC_indices, autotrophs)

! !DESCRIPTION:
!  Initialize ecosys tracer module. This involves setting metadata, reading
!  the module namelist, setting initial conditions, setting up forcing,
!  and defining additional tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

! !INPUT/OUTPUT PARAMETERS:

  type(autotroph_type), dimension(autotroph_cnt), intent(inout) :: autotrophs

  type(BGC_indices_type), intent(inout) :: BGC_indices

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (BGC_i4) :: &
      auto_ind,        &
      Chl_ind,         &
      C_ind,           &
      Fe_ind,          &
      Si_ind,          &
      CaCO3_ind

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   BGC_indices%short_name(BGC_indices%po4_ind)='PO4'
   BGC_indices%long_name(BGC_indices%po4_ind)='Dissolved Inorganic Phosphate'

   BGC_indices%short_name(BGC_indices%no3_ind)='NO3'
   BGC_indices%long_name(BGC_indices%no3_ind)='Dissolved Inorganic Nitrate'

   BGC_indices%short_name(BGC_indices%sio3_ind)='SiO3'
   BGC_indices%long_name(BGC_indices%sio3_ind)='Dissolved Inorganic Silicate'

   BGC_indices%short_name(BGC_indices%nh4_ind)='NH4'
   BGC_indices%long_name(BGC_indices%nh4_ind)='Dissolved Ammonia'

   BGC_indices%short_name(BGC_indices%fe_ind)='Fe'
   BGC_indices%long_name(BGC_indices%fe_ind)='Dissolved Inorganic Iron'

   BGC_indices%short_name(BGC_indices%o2_ind)='O2'
   BGC_indices%long_name(BGC_indices%o2_ind)='Dissolved Oxygen'

   BGC_indices%short_name(BGC_indices%dic_ind)='DIC'
   BGC_indices%long_name(BGC_indices%dic_ind)='Dissolved Inorganic Carbon'

   BGC_indices%short_name(BGC_indices%dic_alt_co2_ind)='DIC_ALT_CO2'
   BGC_indices%long_name(BGC_indices%dic_alt_co2_ind)='Dissolved Inorganic Carbon, Alternative CO2'

   BGC_indices%short_name(BGC_indices%alk_ind)='ALK'
   BGC_indices%long_name(BGC_indices%alk_ind)='Alkalinity'

   BGC_indices%short_name(BGC_indices%doc_ind)='DOC'
   BGC_indices%long_name(BGC_indices%doc_ind)='Dissolved Organic Carbon'

   BGC_indices%short_name(BGC_indices%don_ind)='DON'
   BGC_indices%long_name(BGC_indices%don_ind)='Dissolved Organic Nitrogen'

   BGC_indices%short_name(BGC_indices%dofe_ind)='DOFe'
   BGC_indices%long_name(BGC_indices%dofe_ind)='Dissolved Organic Iron'

   BGC_indices%short_name(BGC_indices%dop_ind)='DOP'
   BGC_indices%long_name(BGC_indices%dop_ind)='Dissolved Organic Phosphorus'

   BGC_indices%short_name(BGC_indices%dopr_ind)='DOPr'
   BGC_indices%long_name(BGC_indices%dopr_ind)='Refractory DOP'

   BGC_indices%short_name(BGC_indices%donr_ind)='DONr'
   BGC_indices%long_name(BGC_indices%donr_ind)='Refractory DON'

   BGC_indices%short_name(BGC_indices%zooC_ind)='zooC'
   BGC_indices%long_name(BGC_indices%zooC_ind)='Zooplankton Carbon'

! now do autotrophs

   do auto_ind = 1, autotroph_cnt

      if (auto_ind == BGC_indices%sp_ind) then
        Chl_ind = BGC_indices%spChl_ind
        C_ind   = BGC_indices%spC_ind
        Fe_ind  = BGC_indices%spFe_ind
      elseif (auto_ind == BGC_indices%diat_ind) then
        Chl_ind = BGC_indices%diatChl_ind
        C_ind   = BGC_indices%diatC_ind
        Fe_ind  = BGC_indices%diatFe_ind
      elseif (auto_ind == BGC_indices%diaz_ind) then
        Chl_ind = BGC_indices%diazChl_ind
        C_ind   = BGC_indices%diazC_ind
        Fe_ind  = BGC_indices%diazFe_ind
      elseif (auto_ind == BGC_indices%phaeo_ind) then
        Chl_ind = BGC_indices%phaeoChl_ind
        C_ind   = BGC_indices%phaeoC_ind
        Fe_ind  = BGC_indices%phaeoFe_ind
      endif

      BGC_indices%short_name(Chl_ind) = trim(autotrophs(auto_ind)%sname) // 'Chl'
      BGC_indices%long_name(Chl_ind)  = trim(autotrophs(auto_ind)%lname) // ' Chlorophyll'
      autotrophs(auto_ind)%Chl_ind = Chl_ind

      BGC_indices%short_name(C_ind) = trim(autotrophs(auto_ind)%sname) // 'C'
      BGC_indices%long_name(C_ind)  = trim(autotrophs(auto_ind)%lname) // ' Carbon'
      autotrophs(auto_ind)%C_ind = C_ind

      BGC_indices%short_name(Fe_ind) = trim(autotrophs(auto_ind)%sname) // 'Fe'
      BGC_indices%long_name(Fe_ind)  = trim(autotrophs(auto_ind)%lname) // ' Iron'
      autotrophs(auto_ind)%Fe_ind = Fe_ind

      if (autotrophs(auto_ind)%kSiO3 > c0) then
         Si_ind = BGC_indices%diatSi_ind
         BGC_indices%short_name(Si_ind) = trim(autotrophs(auto_ind)%sname) // 'Si'
         BGC_indices%long_name(Si_ind)  = trim(autotrophs(auto_ind)%lname) // ' Silicon'
         autotrophs(auto_ind)%Si_ind = Si_ind
      else
         autotrophs(auto_ind)%Si_ind = 0
      endif

      if (autotrophs(auto_ind)%imp_calcifier .or. &
          autotrophs(auto_ind)%exp_calcifier) then
         CaCO3_ind = BGC_indices%spCaCO3_ind
         BGC_indices%short_name(CaCO3_ind) = trim(autotrophs(auto_ind)%sname) // 'CaCO3'
         BGC_indices%long_name(CaCO3_ind)  = trim(autotrophs(auto_ind)%lname) // ' CaCO3'
         autotrophs(auto_ind)%CaCO3_ind = CaCO3_ind
      else
         autotrophs(auto_ind)%CaCO3_ind = 0
      endif
   end do

   BGC_indices%units(:)                       = 'mmol/m^3'
   BGC_indices%units(BGC_indices%alk_ind)     = 'meq/m^3'
   BGC_indices%units(BGC_indices%spChl_ind)   = 'mg/m^3'
   BGC_indices%units(BGC_indices%diatChl_ind) = 'mg/m^3'
   BGC_indices%units(BGC_indices%diazChl_ind) = 'mg/m^3'
   BGC_indices%units(BGC_indices%phaeoChl_ind) = 'mg/m^3'

!-----------------------------------------------------------------------
!EOC

 end subroutine BGC_init

!***********************************************************************
!BOP
! !IROUTINE: BGC_SourceSink
! !INTERFACE:

 subroutine BGC_SourceSink(autotrophs, BGC_indices, BGC_input, BGC_forcing,   &
                           BGC_output, BGC_diagnostic_fields, numLevelsMax,   &
                           numColumnsMax, numColumns, alt_co2_use_eco)

! !DESCRIPTION:
!  Compute time derivatives for ecosystem state variables
!
! !REVISION HISTORY:
!  same as module

   implicit none

! !INPUT PARAMETERS:

  type(autotroph_type), dimension(autotroph_cnt), intent(in) :: autotrophs
  type(BGC_indices_type),     intent(in ) :: BGC_indices
  type(BGC_input_type),       intent(in ) :: BGC_input
  type(BGC_forcing_type),     intent(in ) :: BGC_forcing

  integer (BGC_i4) :: numLevelsMax, numColumnsMax, numColumns
  logical (BGC_log) :: alt_co2_use_eco

! !OUTPUT PARAMETERS:

  type(BGC_output_type),      intent(inout) :: BGC_output
  type(BGC_diagnostics_type), intent(inout) :: BGC_diagnostic_fields

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   type(sinking_particle) :: &
      POC,            & ! base units = nmol C
      P_CaCO3,        & ! base units = nmol CaCO3
      P_SiO2,         & ! base units = nmol SiO2
      dust,           & ! base units = g
      P_iron            ! base units = nmol Fe

   real (BGC_r8), parameter :: &
      mpercm = 0.01_BGC_r8   !  meters/cm

   real (BGC_r8) :: &
      QA_dust_def,    & ! incoming deficit in the QA(dust) POC flux
      dust_flux_in_loc,&! local copy of incoming surface dust flux
      SED_DENITRIF,   & ! sedimentary denitrification (nmol N/cm^3/sec)
      OTHER_REMIN,    & ! organic C remin not due oxic or denitrif (nmolC/cm^3/sec)
      ZSATCALC,       & ! Calcite Saturation Depth
      ZSATARAG,       & ! Aragonite Saturation Depth
      CO3_CALC_ANOM_km1,&! CO3 concentration above calcite saturation at k-1
      CO3_ARAG_ANOM_km1 ! CO3 concentration above aragonite saturation at k-1

!  real (BGC_r8), dimension(km,numColumns) :: &
   real (BGC_r8), allocatable, dimension(:,:) :: &
      DIC_loc,        & ! local copy of model DIC
      DIC_ALT_CO2_loc,& ! local copy of model DIC_ALT_CO2
      ALK_loc,        & ! local copy of model ALK
      PO4_loc,        & ! local copy of model PO4
      NO3_loc,        & ! local copy of model NO3
      SiO3_loc,       & ! local copy of model SiO3
      NH4_loc,        & ! local copy of model NH4
      Fe_loc,         & ! local copy of model Fe
      O2_loc,         & ! local copy of model O2
      DOC_loc,        & ! local copy of model DOC
      zooC_loc,       & ! local copy of model zooC
      DON_loc,        & ! local copy of model DON
      DOFe_loc,       & ! local copy of model DOFe
      DOP_loc,        & ! local copy of model DOP
      DOPr_loc,       & ! local copy of model DOPr
      DONr_loc          ! local copy of model DONr

!  real (BGC_r8), dimension(km,numColumns,autotroph_cnt) :: &
   real (BGC_r8), allocatable, dimension(:,:,:) :: &
      autotrophChl_loc, & ! local copy of model autotroph Chl
      autotrophC_loc,   & ! local copy of model autotroph C
      autotrophFe_loc,  & ! local copy of model autotroph Fe
      autotrophSi_loc,  & ! local copy of model autotroph Si
      autotrophCaCO3_loc  ! local copy of model autotroph CaCO3

   logical (BGC_log) :: zero_mask

   real (BGC_r8) :: &
      work1,work2,work3,work4,work5,tmpTopt,tmpTmax ! temporaries

   real (BGC_r8) :: &
      f_loss_thres,   &! fraction of grazing loss reduction at depth
      ztop,           & ! depth of top of cell
      PAR_out,        & ! photosynthetically available radiation (W/m^2)
      PAR_in,         & ! photosynthetically available radiation (W/m^2)
      KPARdz,         & ! PAR adsorption coefficient (non-dim)
      PAR_avg,        & ! average PAR over mixed layer depth (W/m^2)
      DOC_prod,       & ! production of DOC (mmol C/m^3/sec)
      DOC_remin,      & ! remineralization of DOC (mmol C/m^3/sec)
      DON_remin,      & ! portion of DON remineralized
      DOFe_remin,     & ! portion of DOFe remineralized
      DOP_remin,      & ! portion of DOP remineralized
      NITRIF,         & ! nitrification (NH4 -> NO3) (mmol N/m^3/sec)
      DENITRIF,       & ! WC nitrification (NO3 -> N2) (mmol N/m^3/sec)
      RESTORE           ! restoring terms for nutrients (mmol ./m^3/sec)

   real (BGC_r8) :: &
      z_umax,         & ! max. zoo growth rate at local T (1/sec)
      C_loss_thres      ! bio-C threshold at which losses go to zero (mmol C/m^3)

   real (BGC_r8) :: &
      Tfunc,          & ! temp response function GD98 (non-dim)
      f_nut,          & ! nut limitation factor, modifies C fixation (non-dim)
      PCmax,          & ! max value of PCphoto at temperature TEMP (1/sec)
      light_lim,      & ! light limitation factor
      PCphoto,        & ! C-specific rate of photosynth. (1/sec)
      pChl              ! Chl synth. regulation term (mg Chl/mmol N)

   real (BGC_r8) :: & ! max of 39 continuation lines
      f_zoo_detr,     & ! frac of zoo losses into large detrital pool (non-dim)
      Fe_scavenge_rate,&! annual scavenging rate of iron as % of ambient
      Fe_scavenge,    & ! loss of dissolved iron, scavenging (mmol Fe/m^3/sec)
      Zprime,         & ! used to limit zoo mort at low biomass (mmol C/m^3)
      zoo_loss,       & ! mortality & higher trophic grazing on zooplankton (mmol C/m^3/sec)
      zoo_loss_doc,   & ! zoo_loss routed to doc (mmol C/m^3/sec)
      zoo_loss_dic      ! zoo_loss routed to dic (mmol C/m^3/sec)

   real (BGC_r8) :: &
      VNC,            & ! C-specific N uptake rate (mmol N/mmol C/sec)
      VPO4,           & ! C-specific PO4 uptake (non-dim)
      VDOP,           & ! C-specific DOP uptake rate (non-dim)
      VPtot,          & ! total P uptake rate (non-dim)
      VFe,            & ! C-specific Fe uptake (non-dim)
      VSiO3             ! C-specific SiO3 uptake (non-dim)

   real (BGC_r8), dimension(autotroph_cnt) :: &
      thetaC,         & ! local Chl/C ratio (mg Chl/mmol C)
      QCaCO3,         & ! CaCO3/C ratio (mmol CaCO3/mmol C)
      VNO3,           & ! NO3 uptake rate (non-dim)
      VNH4,           & ! NH4 uptake rate (non-dim)
      VNtot,          & ! total N uptake rate (non-dim)
      NO3_V,          & ! nitrate uptake (mmol NO3/m^3/sec)
      NH4_V,          & ! ammonium uptake (mmol NH4/m^3/sec)
      PO4_V,          & ! PO4 uptake (mmol PO4/m^3/sec)
      DOP_V,          & ! DOP uptake (mmol DOP/m^3/sec)
      Qfe,            & ! init fe/C ratio (mmolFe/mmolC)
      gQfe,           & ! fe/C for growth
      Qsi,            & ! initial Si/C ratio (mmol Si/mmol C)
      gQsi,           & ! diatom Si/C ratio for growth (new biomass)
      Pprime,         & ! used to limit autotroph mort at low biomass (mmol C/m^3)
      auto_graze,     & ! autotroph grazing rate (mmol C/m^3/sec)
      auto_graze_zoo, & ! auto_graze routed to zoo (mmol C/m^3/sec)
      auto_graze_poc, & ! auto_graze routed to poc (mmol C/m^3/sec)
      auto_graze_doc, & ! auto_graze routed to doc (mmol C/m^3/sec)
      auto_graze_dic, & ! auto_graze routed to dic (mmol C/m^3/sec)
      auto_loss,      & ! autotroph non-grazing mort (mmol C/m^3/sec)
      auto_loss_poc,  & ! auto_loss routed to poc (mmol C/m^3/sec)
      auto_loss_doc,  & ! auto_loss routed to doc (mmol C/m^3/sec)
      auto_loss_dic,  & ! auto_loss routed to dic (mmol C/m^3/sec)
      auto_agg,       & ! autotroph aggregation (mmol C/m^3/sec)
      photoC,         & ! C-fixation (mmol C/m^3/sec)
      photoFe,        & ! iron uptake
      photoSi,        & ! silicon uptake (mmol Si/m^3/sec)
      CaCO3_PROD,     & ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
      photoacc,       & ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
      Nfix,           & ! total Nitrogen fixation (mmol N/m^3/sec)
      Nexcrete          ! fixed N excretion

   real (BGC_r8) :: &
      remaining_P       ! used in routing P from autotrophs w/ Qp different from Qp_zoo_pom

   real (BGC_r8), dimension(autotroph_cnt) :: &
      remaining_P_dop,& ! remaining_P from mort routed to DOP pool
      remaining_P_dip   ! remaining_P from mort routed to remin

   real (BGC_r8) :: &
      DON_prod,       & ! production of dissolved organic N
      DOFe_prod,      & ! produciton of dissolved organic Fe
      DOP_prod,       & ! production of dissolved organic P
      O2_PRODUCTION,  & ! O2 production
      O2_CONSUMPTION, & ! O2 consumption
      DONr_remin,     & ! portion of refractory DON remineralized
      DOPr_remin        ! portion of refractory DOP remineralized

   real (BGC_r8) :: &
      partial_thickness_100m,  &! 
      CO3,            &! carbonate ion
      HCO3,           &! bicarbonate ion
      H2CO3,          &! carbonic acid
      CO3_ALT_CO2,    &! carbonate ion, alternative CO2
      HCO3_ALT_CO2,   &! bicarbonate ion, alternative CO2
      H2CO3_ALT_CO2,  &! carbonic acid, alternative CO2
      OMEGA_CALC,     &! solubility ratio for aragonite
      OMEGA_ARAG       ! solubility ratio for calcite

   integer (BGC_i4) :: &
      k,              & ! vertical level index
      n,              & ! tracer index
      auto_ind,       & ! autotroph functional group index
      auto_ind2,      & ! autotroph functional group index
      kmax,           & ! maximum number of vertical levels in column
      column            ! index for looping over columns

!-----------------------------------------------------------------------
!  local copies of non-autotroph indices
!-----------------------------------------------------------------------

   integer (BGC_i4) ::    &
      po4_ind,            &
      no3_ind,            &
      sio3_ind,           &
      nh4_ind,            &
      fe_ind,             &
      o2_ind,             &
      dic_ind,            &
      dic_alt_co2_ind,    &
      alk_ind,            &
      doc_ind,            &
      don_ind,            &
      dofe_ind,           &
      dop_ind,            &
      dopr_ind,           &
      donr_ind,           &
      zooC_ind

   logical (BGC_log) :: &
      lcalc_co2_terms, &! are any alt_co2 terms being time averaged
      lalt_co2_terms    ! are any alt_co2 terms being time averaged

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  initialize  all tendencies to zero
!-----------------------------------------------------------------------

   BGC_output%BGC_tendencies = c0

!-----------------------------------------------------------------------
!  allocate local copies of tracers
!-----------------------------------------------------------------------

   allocate(DIC_loc(numLevelsMax,numColumns))
   allocate(DIC_ALT_CO2_loc(numLevelsMax,numColumns))
   allocate(ALK_loc(numLevelsMax,numColumns))
   allocate(PO4_loc(numLevelsMax,numColumns))
   allocate(NO3_loc(numLevelsMax,numColumns))
   allocate(SiO3_loc(numLevelsMax,numColumns))
   allocate(NH4_loc(numLevelsMax,numColumns))
   allocate(Fe_loc(numLevelsMax,numColumns))
   allocate(O2_loc(numLevelsMax,numColumns))
   allocate(DOC_loc(numLevelsMax,numColumns))
   allocate(zooC_loc(numLevelsMax,numColumns))
   allocate(DON_loc(numLevelsMax,numColumns))
   allocate(DOFe_loc(numLevelsMax,numColumns))
   allocate(DOP_loc(numLevelsMax,numColumns))
   allocate(DOPr_loc(numLevelsMax,numColumns))
   allocate(DONr_loc(numLevelsMax,numColumns))

   allocate(autotrophChl_loc(numLevelsMax,numColumns,autotroph_cnt))
   allocate(autotrophC_loc(numLevelsMax,numColumns,autotroph_cnt))
   allocate(autotrophFe_loc(numLevelsMax,numColumns,autotroph_cnt))
   allocate(autotrophSi_loc(numLevelsMax,numColumns,autotroph_cnt))
   allocate(autotrophCaCO3_loc(numLevelsMax,numColumns,autotroph_cnt))

!-----------------------------------------------------------------------
!  assign non-autotroph indices.  this is not necessary but results in fewer
!    differences between original and new code.
!-----------------------------------------------------------------------

   po4_ind         = BGC_indices%po4_ind
   no3_ind         = BGC_indices%no3_ind
   sio3_ind        = BGC_indices%sio3_ind
   nh4_ind         = BGC_indices%nh4_ind
   fe_ind          = BGC_indices%fe_ind
   o2_ind          = BGC_indices%o2_ind
   dic_ind         = BGC_indices%dic_ind
   dic_alt_co2_ind = BGC_indices%dic_alt_co2_ind
   alk_ind         = BGC_indices%alk_ind
   doc_ind         = BGC_indices%doc_ind
   don_ind         = BGC_indices%don_ind
   dofe_ind        = BGC_indices%dofe_ind
   dop_ind         = BGC_indices%dop_ind
   dopr_ind        = BGC_indices%dopr_ind
   donr_ind        = BGC_indices%donr_ind
   zooC_ind        = BGC_indices%zooC_ind

!-----------------------------------------------------------------------
!  initialize vertical integrals and sums over all autotrophs
!-----------------------------------------------------------------------

   BGC_diagnostic_fields%diag_tot_CaCO3_form = c0
   BGC_diagnostic_fields%diag_tot_bSi_form = c0
   BGC_diagnostic_fields%diag_CaCO3_form_zint = c0
   BGC_diagnostic_fields%diag_tot_CaCO3_form_zint = c0
   BGC_diagnostic_fields%diag_photoC_zint = c0
   BGC_diagnostic_fields%diag_photoC_TOT_zint = c0
   BGC_diagnostic_fields%diag_tot_Nfix = c0
   BGC_diagnostic_fields%diag_photoC_NO3_zint = c0
   BGC_diagnostic_fields%diag_photoC_NO3_TOT = c0
   BGC_diagnostic_fields%diag_photoC_NO3_TOT_zint = c0
   BGC_diagnostic_fields%diag_Chl_TOT_zint_100m = c0
   BGC_diagnostic_fields%diag_Jint_Ctot = c0
   BGC_diagnostic_fields%diag_Jint_100m_Ctot = c0
   BGC_diagnostic_fields%diag_Jint_Ntot = c0
   BGC_diagnostic_fields%diag_Jint_100m_Ntot = c0
   BGC_diagnostic_fields%diag_Jint_Ptot = c0
   BGC_diagnostic_fields%diag_Jint_100m_Ptot = c0
   BGC_diagnostic_fields%diag_Jint_Sitot = c0
   BGC_diagnostic_fields%diag_Jint_100m_Sitot = c0

!-----------------------------------------------------------------------
!  initialize ALL diags so that land below ocean points have zeros
!-----------------------------------------------------------------------

   BGC_diagnostic_fields%diag_CO3 = c0
   BGC_diagnostic_fields%diag_HCO3 = c0
   BGC_diagnostic_fields%diag_H2CO3 = c0
   BGC_diagnostic_fields%diag_pH_3D = c0
   BGC_diagnostic_fields%diag_CO3_ALT_CO2 = c0
   BGC_diagnostic_fields%diag_HCO3_ALT_CO2 = c0
   BGC_diagnostic_fields%diag_H2CO3_ALT_CO2 = c0
   BGC_diagnostic_fields%diag_pH_3D_ALT_CO2 = c0
   BGC_diagnostic_fields%diag_co3_sat_calc = c0
   BGC_diagnostic_fields%diag_co3_sat_arag = c0
   BGC_diagnostic_fields%diag_NO3_RESTORE = c0
   BGC_diagnostic_fields%diag_NITRIF = c0
   BGC_diagnostic_fields%diag_DENITRIF = c0
   BGC_diagnostic_fields%diag_SiO3_RESTORE = c0
   BGC_diagnostic_fields%diag_PO4_RESTORE = c0
   BGC_diagnostic_fields%diag_O2_PRODUCTION = c0
   BGC_diagnostic_fields%diag_O2_CONSUMPTION = c0
   BGC_diagnostic_fields%diag_AOU = c0
   BGC_diagnostic_fields%diag_PAR_avg = c0
   BGC_diagnostic_fields%diag_zoo_loss = c0
   BGC_diagnostic_fields%diag_auto_graze_TOT = c0
   BGC_diagnostic_fields%diag_photoC_TOT = c0
   BGC_diagnostic_fields%diag_DOC_prod = c0
   BGC_diagnostic_fields%diag_DOC_remin = c0
   BGC_diagnostic_fields%diag_DON_prod = c0
   BGC_diagnostic_fields%diag_DON_remin = c0
   BGC_diagnostic_fields%diag_DOP_prod = c0
   BGC_diagnostic_fields%diag_DOP_remin = c0
   BGC_diagnostic_fields%diag_DOFe_prod = c0
   BGC_diagnostic_fields%diag_DOFe_remin = c0
   BGC_diagnostic_fields%diag_Fe_scavenge = c0
   BGC_diagnostic_fields%diag_Fe_scavenge_rate = c0
   BGC_diagnostic_fields%diag_POC_FLUX_IN = c0
   BGC_diagnostic_fields%diag_POC_PROD = c0
   BGC_diagnostic_fields%diag_POC_REMIN = c0
   BGC_diagnostic_fields%diag_CaCO3_FLUX_IN = c0
   BGC_diagnostic_fields%diag_CaCO3_PROD = c0
   BGC_diagnostic_fields%diag_CaCO3_REMIN = c0
   BGC_diagnostic_fields%diag_SiO2_FLUX_IN = c0
   BGC_diagnostic_fields%diag_SiO2_PROD = c0
   BGC_diagnostic_fields%diag_SiO2_REMIN = c0
   BGC_diagnostic_fields%diag_dust_FLUX_IN = c0
   BGC_diagnostic_fields%diag_dust_REMIN = c0
   BGC_diagnostic_fields%diag_P_iron_FLUX_IN = c0
   BGC_diagnostic_fields%diag_P_iron_PROD = c0
   BGC_diagnostic_fields%diag_P_iron_REMIN = c0
   BGC_diagnostic_fields%diag_calcToSed = c0
   BGC_diagnostic_fields%diag_bsiToSed = c0
   BGC_diagnostic_fields%diag_pocToSed = c0
   BGC_diagnostic_fields%diag_SedDenitrif = c0
   BGC_diagnostic_fields%diag_OtherRemin = c0
   BGC_diagnostic_fields%diag_ponToSed = c0
   BGC_diagnostic_fields%diag_popToSed = c0
   BGC_diagnostic_fields%diag_dustToSed = c0
   BGC_diagnostic_fields%diag_pfeToSed = c0

   BGC_diagnostic_fields%diag_zsatcalc = c0
   BGC_diagnostic_fields%diag_zsatarag = c0
   BGC_diagnostic_fields%diag_O2_ZMIN = c0
   BGC_diagnostic_fields%diag_O2_ZMIN_DEPTH = c0

   BGC_diagnostic_fields%diag_N_lim = c0
   BGC_diagnostic_fields%diag_Fe_lim = c0
   BGC_diagnostic_fields%diag_P_lim = c0
   BGC_diagnostic_fields%diag_SiO3_lim = c0
   BGC_diagnostic_fields%diag_light_lim = c0
   BGC_diagnostic_fields%diag_photoNO3 = c0
   BGC_diagnostic_fields%diag_photoNH4 = c0
   BGC_diagnostic_fields%diag_PO4_uptake = c0
   BGC_diagnostic_fields%diag_DOP_uptake = c0
   BGC_diagnostic_fields%diag_photoFe = c0
   BGC_diagnostic_fields%diag_bSi_form = c0
   BGC_diagnostic_fields%diag_CaCO3_form = c0
   BGC_diagnostic_fields%diag_Nfix = c0
   BGC_diagnostic_fields%diag_auto_graze = c0
   BGC_diagnostic_fields%diag_auto_loss = c0
   BGC_diagnostic_fields%diag_auto_agg = c0
   BGC_diagnostic_fields%diag_photoC = c0
   BGC_diagnostic_fields%diag_photoC_NO3 = c0

!-----------------------------------------------------------------------
!  loop over levels
!-----------------------------------------------------------------------

   setup_loop: do column = 1, numColumns

   kmax = BGC_input%number_of_active_levels(column)
   if (kmax < 1) cycle setup_loop

   do k = 1, kmax

!-----------------------------------------------------------------------
!  create local copies of model tracers
!  treat negative values as zero
!  apply mask to local copies
!-----------------------------------------------------------------------

!maltrud intel fails if i use c0 for max() instead of 0.0
   DIC_loc(k,column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,dic_ind))
   DIC_ALT_CO2_loc(k,column) = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,dic_alt_co2_ind))
   ALK_loc(k,column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,alk_ind))
   PO4_loc(k,column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,po4_ind))
   NO3_loc(k,column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,no3_ind))
   SiO3_loc(k,column)     = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,sio3_ind))
   NH4_loc(k,column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,nh4_ind))
   Fe_loc(k,column)       = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,fe_ind))
   O2_loc(k,column)       = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,o2_ind))
   DOC_loc(k,column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,doc_ind))
   zooC_loc(k,column)     = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,zooC_ind))
   DON_loc(k,column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,don_ind))
   DOFe_loc(k,column)     = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,dofe_ind))
   DOP_loc(k,column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,dop_ind))
   DOPr_loc(k,column)     = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,dopr_ind))
   DONr_loc(k,column)     = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,donr_ind))

   do auto_ind = 1, autotroph_cnt

      n = autotrophs(auto_ind)%Chl_ind
      autotrophChl_loc(k,column,auto_ind) = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,n))

      n = autotrophs(auto_ind)%C_ind
      autotrophC_loc(k,column,auto_ind) = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,n))

      n = autotrophs(auto_ind)%Fe_ind
      autotrophFe_loc(k,column,auto_ind) = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,n))

      n = autotrophs(auto_ind)%Si_ind
      if (n > 0) then
         autotrophSi_loc(k,column,auto_ind) = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,n))
      endif

      n = autotrophs(auto_ind)%CaCO3_ind
      if (n > 0) then
         autotrophCaCO3_loc(k,column,auto_ind) = max(0.0_BGC_r8, BGC_input%BGC_tracers(k,column,n))
      endif

   end do

   end do  !  end of setup k loop

   enddo  setup_loop  !  end of setup column loop

!-----------------------------------------------------------------------
!  HERE IS WHERE MODEL AGNOSTIC SHOULD BEGIN
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  loop over columns
!-----------------------------------------------------------------------

   column_loop: do column = 1, numColumns

   kmax = BGC_input%number_of_active_levels(column)
   if (kmax < 1) cycle column_loop
   
!-----------------------------------------------------------------------
!  various k==1 initializations
!-----------------------------------------------------------------------

   dust_flux_in_loc = max (0.0_BGC_r8, BGC_forcing%dust_FLUX_IN(column))

   call init_particulate_terms(POC, P_CaCO3, P_SiO2, dust, P_iron, &
                               QA_dust_def, dust_flux_in_loc)

   PAR_out = max (0.0_BGC_r8, BGC_forcing%ShortWaveFlux_surface(column))
   PAR_out = PAR_out*f_qsw_par

!-----------------------------------------------------------------------
!  loop over levels
!-----------------------------------------------------------------------

   do k = 1, kmax

!-----------------------------------------------------------------------
!  If any phyto box are zero, set others to zeros.
!-----------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt
      zero_mask = autotrophChl_loc(k,column,auto_ind) == c0 .or. &
                       autotrophC_loc(k,column,auto_ind) == c0 .or. &
                       autotrophFe_loc(k,column,auto_ind) == c0
      if (autotrophs(auto_ind)%Si_ind > 0) &
         zero_mask = zero_mask .or. autotrophSi_loc(k,column,auto_ind) == c0
      if (zero_mask) then
         autotrophChl_loc(k,column,auto_ind) = c0
         autotrophC_loc(k,column,auto_ind) = c0
         autotrophFe_loc(k,column,auto_ind) = c0
      end if
!maltrud maybe change these to single if .and. statements
      if (autotrophs(auto_ind)%Si_ind > 0) then
         if (zero_mask) autotrophSi_loc(k,column,auto_ind) = c0
      endif
      if (autotrophs(auto_ind)%CaCO3_ind > 0) then
         if (zero_mask) autotrophCaCO3_loc(k,column,auto_ind) = c0
      endif
   end do

!-----------------------------------------------------------------------
!  set local variables, with incoming ratios
!-----------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt
      thetaC(auto_ind) = autotrophChl_loc(k,column,auto_ind) / (autotrophC_loc(k,column,auto_ind) + epsC)
      Qfe(auto_ind) = autotrophFe_loc(k,column,auto_ind) / (autotrophC_loc(k,column,auto_ind) + epsC)
      if (autotrophs(auto_ind)%Si_ind > 0) then
         Qsi(auto_ind) = min(autotrophSi_loc(k,column,auto_ind) / (autotrophC_loc(k,column,auto_ind) + epsC), gQsi_max)
      endif
   end do

!-----------------------------------------------------------------------
!  DETERMINE NEW ELEMENTAL RATIOS FOR GROWTH (NEW BIOMASS)
!  Modify these initial ratios under low ambient iron conditions
!  Modify the initial si/C ratio under low ambient Si conditions
!-----------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt
      gQfe(auto_ind) = autotrophs(auto_ind)%gQfe_0
      if (Fe_loc(k,column) < cks * autotrophs(auto_ind)%kFe) then
         gQfe(auto_ind) = &
            max((gQfe(auto_ind) * Fe_loc(k,column) / (cks * autotrophs(auto_ind)%kFe)), &
                autotrophs(auto_ind)%gQfe_min)
      end if

      if (autotrophs(auto_ind)%Si_ind > 0) then
         gQsi(auto_ind) = gQsi_0
         if ((Fe_loc(k,column) < cksi * autotrophs(auto_ind)%kFe) .and. (Fe_loc(k,column) > c0) .and. &
                (SiO3_loc(k,column) > (cksi * autotrophs(auto_ind)%kSiO3))) then
            gQsi(auto_ind) = min((gQsi(auto_ind) * cksi * autotrophs(auto_ind)%kFe / Fe_loc(k,column)), gQsi_max)
         end if

         if (Fe_loc(k,column) == c0) then
            gQsi(auto_ind) = gQsi_max
         end if

         if (SiO3_loc(k,column) < (cksi * autotrophs(auto_ind)%kSiO3)) then
            gQsi(auto_ind) = max((gQsi(auto_ind) * SiO3_loc(k,column) / (cksi * autotrophs(auto_ind)%kSiO3)), &
                                          gQsi_min)
         end if
      endif

!-----------------------------------------------------------------------
!  QCaCO3 is the percentage of sp organic matter which is associated
!  with coccolithophores
!-----------------------------------------------------------------------

      if (autotrophs(auto_ind)%CaCO3_ind > 0) then
         QCaCO3(auto_ind) = autotrophCaCO3_loc(k,column,auto_ind) / (autotrophC_loc(k,column,auto_ind) + epsC)
         if (QCaCO3(auto_ind) > QCaCO3_max) QCaCO3(auto_ind) = QCaCO3_max
      endif
   end do

!-----------------------------------------------------------------------
!  compute PAR related quantities
!  Morel, Maritorena, JGR, Vol 106, No. C4, pp 7163--7180, 2001
!
!  0.45   fraction of incoming SW -> PAR (non-dim)
!-----------------------------------------------------------------------

   PAR_in = PAR_out
!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG7: ',k,PAR_in,PAR_out(i,j),PAR_surface(i,j,bid),sum(PAR_surface)
!endif


   work1 = max(sum(autotrophChl_loc(k,column,:), dim=1), 0.02_BGC_r8)
   if (work1 < 0.13224_BGC_r8) then
      KPARdz = 0.000919_BGC_r8*(work1**0.3536_BGC_r8)
   else
      KPARdz = 0.001131_BGC_r8*(work1**0.4562_BGC_r8)
   end if

   KPARdz = KPARdz * BGC_input%cell_thickness(k,column)

   PAR_out = PAR_in * exp(-KPARdz)
   PAR_avg = PAR_in * (c1 - exp(-KPARdz)) / KPARdz

!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG6: ',k,PAR_in,PAR_out(i,j),work1,KPARdz,PAR_avg,autotrophChl_loc(k,column,:)
!endif

!-----------------------------------------------------------------------
!  compute terms of carbonate chemistry
!-----------------------------------------------------------------------

   lcalc_co2_terms = .true.
   lalt_co2_terms = .true.
!  lcalc_co2_terms = .false.
!  lalt_co2_terms = .false.

   work3 = c0
   work4 = c0
   if (lcalc_co2_terms) then
      if (BGC_output%PH_PREV_3D(k,column) /= c0) then
         work1 = BGC_output%PH_PREV_3D(k,column) - del_ph
         work2 = BGC_output%PH_PREV_3D(k,column) + del_ph
      else
         work1 = phlo_3d_init
         work2 = phhi_3d_init
      end if
      work5 = BGC_input%cell_center_depth(k,column)*0.01_BGC_r8
!     call comp_CO3terms( k, BGC_input%cell_center_depth(k,column), .true., &
      call comp_CO3terms( k, work5, .true., &
                         BGC_input%PotentialTemperature(k,column), BGC_input%Salinity(k,column), DIC_loc(k,column),  &
                         ALK_loc(k,column), PO4_loc(k,column), SiO3_loc(k,column),           &
                         work1, work2, work3, H2CO3, HCO3, CO3)
      BGC_output%PH_PREV_3D(k,column) = work3
   else
      H2CO3 = c0
      HCO3  = c0
      CO3   = c0
      BGC_output%PH_PREV_3D(k,column) = 8.0
   endif

   if (lalt_co2_terms) then
      if (BGC_output%PH_PREV_ALT_CO2_3D(k,column) /= c0) then
         work1 = BGC_output%PH_PREV_ALT_CO2_3D(k,column) - del_ph
         work2 = BGC_output%PH_PREV_ALT_CO2_3D(k,column) + del_ph
      else
         work1 = phlo_3d_init
         work2 = phhi_3d_init
      end if
      work5 = BGC_input%cell_center_depth(k,column)*0.01_BGC_r8
!     call comp_CO3terms( k, BGC_input%cell_center_depth(k,column), .true., &
      call comp_CO3terms( k, work5, .true., &
                         BGC_input%PotentialTemperature(k,column), BGC_input%Salinity(k,column), DIC_loc(k,column),  &
                         ALK_loc(k,column), PO4_loc(k,column), SiO3_loc(k,column),           &
                         work1, work2, work4, H2CO3_ALT_CO2, HCO3_ALT_CO2, CO3_ALT_CO2)
      BGC_output%PH_PREV_ALT_CO2_3D(k,column) = work4
   else
      H2CO3_ALT_CO2 = c0
      HCO3_ALT_CO2  = c0
      CO3_ALT_CO2   = c0
      BGC_output%PH_PREV_ALT_CO2_3D(k,column) = 8.0
   endif

   BGC_diagnostic_fields%diag_CO3(k,column) = CO3
   BGC_diagnostic_fields%diag_HCO3(k,column) = HCO3
   BGC_diagnostic_fields%diag_H2CO3(k,column) = H2CO3
   BGC_diagnostic_fields%diag_pH_3D(k,column) = work3
   BGC_diagnostic_fields%diag_CO3_ALT_CO2(k,column) = CO3_ALT_CO2
   BGC_diagnostic_fields%diag_HCO3_ALT_CO2(k,column) = HCO3_ALT_CO2
   BGC_diagnostic_fields%diag_H2CO3_ALT_CO2(k,column) = H2CO3_ALT_CO2
   BGC_diagnostic_fields%diag_pH_3D_ALT_CO2(k,column) = work4

   work5 = BGC_input%cell_center_depth(k,column)*0.01_BGC_r8
!  call comp_co3_sat_vals(k, BGC_input%cell_center_depth(k,column),   &
   call comp_co3_sat_vals(k, work5,   &
      BGC_input%PotentialTemperature(k,column), BGC_input%Salinity(k,column), work1, work2) 
                            
   BGC_diagnostic_fields%diag_co3_sat_calc(k,column) = work1
   BGC_diagnostic_fields%diag_co3_sat_arag(k,column) = work2

   if (k == 1) then
      ! set to -1, i.e. depth not found yet,
      ! if mask == .true. and surface supersaturated to -1
      ZSATCALC = merge(-c1, c0, CO3 > work1)
      ZSATARAG = merge(-c1, c0, CO3 > work2)
   else
      work4 = BGC_input%cell_center_depth(k-1,column) + (BGC_input%cell_center_depth(k,column) - BGC_input%cell_center_depth(k-1,column))
      if (ZSATCALC == -c1 .and. CO3 <= work1) then
         ZSATCALC = work4 * &
            CO3_CALC_ANOM_km1 / (CO3_CALC_ANOM_km1 - (CO3 - work1))
      endif
      if (ZSATARAG == -c1 .and. CO3 <= work2) then
         ZSATARAG = work4 * &
            CO3_ARAG_ANOM_km1 / (CO3_ARAG_ANOM_km1 - (CO3 - work2))
      endif
      if (ZSATCALC == -c1 .and. k == kmax) then
         ZSATCALC = BGC_input%cell_bottom_depth(k,column)
      endif
      if (ZSATARAG == -c1 .and. k == kmax) then
         ZSATARAG = BGC_input%cell_bottom_depth(k,column)
      endif
   endif

   CO3_CALC_ANOM_km1 = CO3 - work1
   CO3_ARAG_ANOM_km1 = CO3 - work2

   if (k == kmax) then
      BGC_diagnostic_fields%diag_zsatcalc(column) = ZSATCALC
      BGC_diagnostic_fields%diag_zsatarag(column) = ZSATARAG
   endif

!-----------------------------------------------------------------------
!  Tref = 30.0 reference temperature (deg. C)
!
!  Using q10 formulation with Q10 value of 2.0 (Doney et al., 1996).
!  growth, mort and grazing rates scaled by Tfunc where they are computed
!-----------------------------------------------------------------------

   Tfunc = Q_10**(((BGC_input%PotentialTemperature(k,column) + T0_Kelvin_BGC) - (Tref + T0_Kelvin_BGC)) / c10)

!-----------------------------------------------------------------------
!  calculate the loss threshold interpolation factor
!-----------------------------------------------------------------------

   if (BGC_input%cell_center_depth(k,column) > thres_z1) then
      if (BGC_input%cell_center_depth(k,column) < thres_z2) then
         f_loss_thres = (thres_z2 - BGC_input%cell_center_depth(k,column))/(thres_z2 - thres_z1)
      else
         f_loss_thres = c0
      endif
   else
      f_loss_thres = c1
   endif

!-----------------------------------------------------------------------
!  Compute Pprime for all autotrophs, used for loss terms
!  temp_thres for phaeo is the upper limit for growth (swang)
!-----------------------------------------------------------------------

!   do auto_ind = 1, autotroph_cnt
!      C_loss_thres = f_loss_thres * autotrophs(auto_ind)%loss_thres
!      if (BGC_input%PotentialTemperature(k,column) < autotrophs(auto_ind)%temp_thres)   &
!         C_loss_thres = f_loss_thres * autotrophs(auto_ind)%loss_thres2
!      if (autotrophs(auto_ind)%temp_thres2 > c0 .and. BGC_input%PotentialTemperature(k,column) > autotrophs(auto_ind)%temp_thres2) then
!         C_loss_thres = f_loss_thres * autotrophs(auto_ind)%loss_thres2
!      endif
!      Pprime(auto_ind) = max(autotrophC_loc(k,column,auto_ind) - C_loss_thres, 0.0_BGC_r8)
!   end do

   do auto_ind = 1, autotroph_cnt
      C_loss_thres = f_loss_thres * autotrophs(auto_ind)%loss_thres
   select case (autotrophs(auto_ind)%temp_function)

      case (tfnc_q10)

	    if (BGC_input%PotentialTemperature(k,column) < autotrophs(auto_ind)%temp_thres)   &
	        C_loss_thres = f_loss_thres * autotrophs(auto_ind)%loss_thres2

      case (tfnc_quasi_mmrt)

            if (BGC_input%cell_latitude(column) >= 0.0_BGC_r8 ) then
            	tmpTmax = autotrophs(auto_ind)%temp_thresN
            else
            	tmpTmax = autotrophs(auto_ind)%temp_thresS
            end if

           if (BGC_input%PotentialTemperature(k,column) > tmpTmax) &
               C_loss_thres = f_loss_thres * autotrophs(auto_ind)%loss_thres2

   end select
      Pprime(auto_ind) = max(autotrophC_loc(k,column,auto_ind) - C_loss_thres, 0.0_BGC_r8)
   end do

!-----------------------------------------------------------------------
!  Get relative nutrient uptake rates for autotrophs,
!  min. relative uptake rate modifies C fixation in the manner
!  that the min. cell quota does in GD98.
!-----------------------------------------------------------------------
!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG2: ',k,PO4_loc(k,column),DOP_loc(k,column)
!endif


   do auto_ind = 1, autotroph_cnt
      VNO3(auto_ind) = (NO3_loc(k,column) / autotrophs(auto_ind)%kNO3) / &
         (c1 + (NO3_loc(k,column) / autotrophs(auto_ind)%kNO3) + (NH4_loc(k,column) / autotrophs(auto_ind)%kNH4))
      VNH4(auto_ind) = (NH4_loc(k,column) / autotrophs(auto_ind)%kNH4) / &
         (c1 + (NO3_loc(k,column) / autotrophs(auto_ind)%kNO3) + (NH4_loc(k,column) / autotrophs(auto_ind)%kNH4))
      VNtot(auto_ind) = VNO3(auto_ind) + VNH4(auto_ind)
      if (autotrophs(auto_ind)%Nfixer) VNtot(auto_ind) = c1
      BGC_diagnostic_fields%diag_N_lim(k,column,auto_ind) = VNtot(auto_ind)

      VFe = Fe_loc(k,column) / (Fe_loc(k,column) + autotrophs(auto_ind)%kFe)
      BGC_diagnostic_fields%diag_Fe_lim(k,column,auto_ind) = VFe

      f_nut = min(VNtot(auto_ind), VFe)

      VPO4 = (PO4_loc(k,column) / autotrophs(auto_ind)%kPO4) / &
         (c1 + (PO4_loc(k,column) / autotrophs(auto_ind)%kPO4) + (DOP_loc(k,column) / autotrophs(auto_ind)%kDOP))
      VDOP = (DOP_loc(k,column) / autotrophs(auto_ind)%kDOP) / &
         (c1 + (PO4_loc(k,column) / autotrophs(auto_ind)%kPO4) + (DOP_loc(k,column) / autotrophs(auto_ind)%kDOP))
      VPtot = VPO4 + VDOP
!!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG3: ',k,auto_ind,autotrophs(auto_ind)%kPO4,autotrophs(auto_ind)%kDOP,VPO4,VDOP,VPtot
!endif

      BGC_diagnostic_fields%diag_P_lim(k,column,auto_ind) = VPtot

      f_nut = min(f_nut, VPtot)

      if (autotrophs(auto_ind)%kSiO3 > c0) then
         VSiO3 = SiO3_loc(k,column) / (SiO3_loc(k,column) + autotrophs(auto_ind)%kSiO3)
         BGC_diagnostic_fields%diag_SiO3_lim(k,column,auto_ind) = VSiO3
         f_nut = min(f_nut, VSiO3)
      endif

!-----------------------------------------------------------------------
!     get photosynth. rate, phyto C biomass change, photoadapt
!         modify growth curve for phaeo (swang)
!-----------------------------------------------------------------------

      PCmax = autotrophs(auto_ind)%PCref * f_nut * Tfunc
      if (BGC_input%PotentialTemperature(k,column) < autotrophs(auto_ind)%temp_thres) PCmax = c0

!      if (autotrophs(auto_ind)%temp_thres2 > c0) then
!         PCmax = PCmax * min(1.0_BGC_r8,((autotrophs(auto_ind)%temp_thres2 - BGC_input%PotentialTemperature(k,column)) / &
!            (autotrophs(auto_ind)%temp_thres2 - 16.3_BGC_r8)))
!         PCmax = max(PCmax, 0.0_BGC_r8)
!      endif

    select case (autotrophs(auto_ind)%temp_function)

      case (tfnc_q10)

            PCmax = PCmax

      case (tfnc_quasi_mmrt)

            if (BGC_input%cell_latitude(column) >= 0.0_BGC_r8 ) then
            	tmpTopt = autotrophs(auto_ind)%temp_optN
            	tmpTmax = autotrophs(auto_ind)%temp_thresN
            else
            	tmpTopt = autotrophs(auto_ind)%temp_optS
            	tmpTmax = autotrophs(auto_ind)%temp_thresS
            end if
            PCmax = PCmax * min(1.0_BGC_r8,((tmpTmax - BGC_input%PotentialTemperature(k,column)) / &
            (tmpTmax - tmpTopt)))
            if (BGC_input%PotentialTemperature(k,column) > tmpTmax) PCmax = c0
    end select
    
      light_lim = (c1 - exp((-c1 * autotrophs(auto_ind)%alphaPI * thetaC(auto_ind) * PAR_avg) / &
                            (PCmax + epsTinv)))
      PCphoto = PCmax * light_lim

      BGC_diagnostic_fields%diag_light_lim(k,column,auto_ind) = light_lim

      photoC(auto_ind) = PCphoto * autotrophC_loc(k,column,auto_ind)

!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG5: ',k,auto_ind,f_nut,Tfunc,autotrophs(auto_ind)%PCref,PCmax,PotentialTemperature(k,column),light_lim,PCphoto,PAR_avg,thetaC(auto_ind)
!endif

!-----------------------------------------------------------------------
!  Get nutrient uptakes by small phyto based on calculated C fixation
!  total N uptake VNC is used in photoadaption
!-----------------------------------------------------------------------

      if (VNtot(auto_ind) > c0) then
         NO3_V(auto_ind) = (VNO3(auto_ind) / VNtot(auto_ind)) * photoC(auto_ind) * Q
         NH4_V(auto_ind) = (VNH4(auto_ind) / VNtot(auto_ind)) * photoC(auto_ind) * Q
         VNC = PCphoto * Q
      else
         NO3_V(auto_ind) = c0
         NH4_V(auto_ind) = c0
         VNC = c0
      end if
      BGC_diagnostic_fields%diag_photoNO3(k,column,auto_ind) = NO3_V(auto_ind)
      BGC_diagnostic_fields%diag_photoNH4(k,column,auto_ind) = NH4_V(auto_ind)

      if (VPtot > c0) then
         PO4_V(auto_ind) = (VPO4 / VPtot) * photoC(auto_ind) * autotrophs(auto_ind)%Qp
         DOP_V(auto_ind) = (VDOP / VPtot) * photoC(auto_ind) * autotrophs(auto_ind)%Qp
      else
         PO4_V(auto_ind) = c0
         DOP_V(auto_ind) = c0
      end if
!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG4: ',k,auto_ind,VPO4,VDOP,VPtot,PO4_V(auto_ind),DOP_V(auto_ind),photoC(auto_ind),autotrophs(auto_ind)%Qp
!endif

      BGC_diagnostic_fields%diag_PO4_uptake(k,column,auto_ind) = PO4_V(auto_ind)
      BGC_diagnostic_fields%diag_DOP_uptake(k,column,auto_ind) = DOP_V(auto_ind)

      photoFe(auto_ind) = photoC(auto_ind) * gQfe(auto_ind)
      BGC_diagnostic_fields%diag_photoFe(k,column,auto_ind) = photoFe(auto_ind)

!-----------------------------------------------------------------------
!  Get nutrient uptake by diatoms based on C fixation
!-----------------------------------------------------------------------

      if (autotrophs(auto_ind)%Si_ind > 0) then
         photoSi(auto_ind) = photoC(auto_ind) * gQsi(auto_ind)
         BGC_diagnostic_fields%diag_bSi_form(k,column,auto_ind) = photoSi(auto_ind)
         BGC_diagnostic_fields%diag_tot_bSi_form(column) =   &
             BGC_diagnostic_fields%diag_tot_bSi_form(column) + photoSi(auto_ind)
      endif

!-----------------------------------------------------------------------
!  calculate pChl, (used in photoadapt., GD98)
!  2.3   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
!  GD 98 Chl. synth. term
!-----------------------------------------------------------------------

      work1 = autotrophs(auto_ind)%alphaPI * thetaC(auto_ind) * PAR_avg
      if (work1 > c0) then
         pChl = autotrophs(auto_ind)%thetaN_max * PCphoto / work1
         photoacc(auto_ind) = (pChl * VNC / thetaC(auto_ind)) * autotrophChl_loc(k,column,auto_ind)
      else
         photoacc(auto_ind) = c0
      end if

!-----------------------------------------------------------------------
!  CaCO3 Production, parameterized as function of small phyto production
!  decrease CaCO3 as function of nutrient limitation decrease CaCO3 prod
!  at low temperatures increase CaCO3 prod under bloom conditions
!  maximum calcification rate is 40% of primary production
!-----------------------------------------------------------------------

      if (autotrophs(auto_ind)%imp_calcifier) then
         CaCO3_PROD(auto_ind) = parm_f_prod_sp_CaCO3 * photoC(auto_ind)
         CaCO3_PROD(auto_ind) = CaCO3_PROD(auto_ind) * f_nut

         if (BGC_input%PotentialTemperature(k,column) < CaCO3_temp_thres1)  &
            CaCO3_PROD(auto_ind) = CaCO3_PROD(auto_ind) *  &
                         max((BGC_input%PotentialTemperature(k,column)-CaCO3_temp_thres2), 0.0_BGC_r8) / &
                         (CaCO3_temp_thres1-CaCO3_temp_thres2)

         if (autotrophC_loc(k,column,auto_ind) > CaCO3_sp_thres)  &
            CaCO3_PROD(auto_ind) = min((CaCO3_PROD(auto_ind) *  &
                         autotrophC_loc(k,column,auto_ind) / CaCO3_sp_thres), &
                         (f_photosp_CaCO3 * photoC(auto_ind)))

         BGC_diagnostic_fields%diag_CaCO3_form(k,column,auto_ind) = CaCO3_PROD(auto_ind)
         BGC_diagnostic_fields%diag_tot_CaCO3_form(k,column) =   &
             BGC_diagnostic_fields%diag_tot_CaCO3_form(k,column) + CaCO3_PROD(auto_ind)

         work1 = BGC_input%cell_thickness(k,column) * CaCO3_PROD(auto_ind)
         BGC_diagnostic_fields%diag_CaCO3_form_zint(column,auto_ind) =  &
            BGC_diagnostic_fields%diag_CaCO3_form_zint(column,auto_ind) + work1
         BGC_diagnostic_fields%diag_tot_CaCO3_form_zint(column) =   &
            BGC_diagnostic_fields%diag_tot_CaCO3_form_zint(column) + work1
      endif

!-----------------------------------------------------------------------
!  get autotroph loss (in C units)
!  autotroph agg loss
!-----------------------------------------------------------------------

      auto_loss(auto_ind) = autotrophs(auto_ind)%mort * Pprime(auto_ind) * Tfunc

      auto_agg(auto_ind) = min((autotrophs(auto_ind)%agg_rate_max * dps) * Pprime(auto_ind),  &
                                autotrophs(auto_ind)%mort2 * Pprime(auto_ind) * Pprime(auto_ind))
      auto_agg(auto_ind) = max((autotrophs(auto_ind)%agg_rate_min * dps) * Pprime(auto_ind),  &
                                auto_agg(auto_ind))

!-----------------------------------------------------------------------
!  get grazing rate (graze_sp) on autotroph (in C units)
!  compute sum of carbon in the grazee class including auto_ind
!-----------------------------------------------------------------------

      work1 = c0
      do auto_ind2 = 1, autotroph_cnt
         if (autotrophs(auto_ind2)%grazee_ind == autotrophs(auto_ind)%grazee_ind) &
            work1 = work1 + Pprime(auto_ind2)
      end do

      z_umax = autotrophs(auto_ind)%z_umax_0 * Tfunc

! decrease grazing pressure on diat, when phaeo growth decreases with
! temperature (swang)      
      if (auto_ind == BGC_indices%diat_ind) then
         if ((BGC_input%cell_latitude(column) >= 0.0_BGC_r8) .and. &
               (BGC_input%PotentialTemperature(k,column) > autotrophs(auto_ind)%temp_optN)) then
               z_umax = z_umax * max((autotrophs(auto_ind)%temp_thresN - BGC_input%PotentialTemperature(k,column)) / &
                            (autotrophs(auto_ind)%temp_thresN - autotrophs(auto_ind)%temp_optN), 0.95_BGC_r8)
         elseif ((BGC_input%cell_latitude(column) <= 0.0_BGC_r8) .and. &
                   (BGC_input%PotentialTemperature(k,column) > autotrophs(auto_ind)%temp_optS)) then
               z_umax = z_umax * max((autotrophs(auto_ind)%temp_thresS - BGC_input%PotentialTemperature(k,column)) / &
                            (autotrophs(auto_ind)%temp_thresS - autotrophs(auto_ind)%temp_optS), 0.95_BGC_r8)
         endif 
      endif

      if (work1 > c0) then
         auto_graze(auto_ind) = (Pprime(auto_ind) / work1) * &
            z_umax * zooC_loc(k,column) * (work1 / (work1 + autotrophs(auto_ind)%z_grz))
      else
         auto_graze(auto_ind) = c0
      endif

!-----------------------------------------------------------------------
!  Get N fixation by diazotrophs based on C fixation,
!  Diazotrophs fix more than they need then 20% is excreted
!-----------------------------------------------------------------------

      if (autotrophs(auto_ind)%Nfixer) then
         work1 = photoC(auto_ind) * Q
         Nfix(auto_ind)     = (work1 * r_Nfix_photo) - NO3_V(auto_ind) - NH4_V(auto_ind)
         Nexcrete(auto_ind) = Nfix(auto_ind) + NO3_V(auto_ind) + NH4_V(auto_ind) - work1
         BGC_diagnostic_fields%diag_Nfix(k,column,auto_ind) = Nfix(auto_ind)
         BGC_diagnostic_fields%diag_tot_Nfix(k,column)  =   &
            BGC_diagnostic_fields%diag_tot_Nfix(k,column) + Nfix(auto_ind)
      endif

!-----------------------------------------------------------------------
!  CALCULATE GRAZING AND OTHER MORT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  routing of grazing and loss terms
!  all aggregation goes to POC
!  currently assumes that 33% of grazed caco3 is remineralized
!  if autotrophs(sp_ind)%graze_zoo ever changes, coefficients on routing grazed sp must change!
!  min.%C routed to POC from grazing for ballast requirements = 0.4 * Qcaco3
!  min.%C routed from sp_loss = 0.59 * QCaCO3, or P_CaCO3%rho
!  NOTE: if autotrophs(diat_ind)%graze_zoo is changed, coeff.s for poc,doc and dic must change!
!-----------------------------------------------------------------------

      auto_graze_zoo(auto_ind) = autotrophs(auto_ind)%graze_zoo * auto_graze(auto_ind)
      if (autotrophs(auto_ind)%imp_calcifier) then
         auto_graze_poc(auto_ind) = auto_graze(auto_ind) * max((caco3_poc_min * QCaCO3(auto_ind)),  &
                                       min(spc_poc_fac * max(1.0_BGC_r8,Pprime(auto_ind)),              &
                                       f_graze_sp_poc_lim))
      else
         auto_graze_poc(auto_ind) = autotrophs(auto_ind)%graze_poc * auto_graze(auto_ind)
      endif
      auto_graze_doc(auto_ind) = autotrophs(auto_ind)%graze_doc * auto_graze(auto_ind)
      auto_graze_dic(auto_ind) = auto_graze(auto_ind) &
         - (auto_graze_zoo(auto_ind) + auto_graze_poc(auto_ind) + auto_graze_doc(auto_ind))

      if (autotrophs(auto_ind)%imp_calcifier) then
         auto_loss_poc(auto_ind) = QCaCO3(auto_ind) * auto_loss(auto_ind)
      else
         auto_loss_poc(auto_ind) = autotrophs(auto_ind)%loss_poc * auto_loss(auto_ind)
      endif
      auto_loss_doc(auto_ind) = (c1 - parm_labile_ratio) * (auto_loss(auto_ind) - auto_loss_poc(auto_ind))
      auto_loss_dic(auto_ind) = parm_labile_ratio * (auto_loss(auto_ind) - auto_loss_poc(auto_ind))

!-----------------------------------------------------------------------
! P from some autotrophs w/ Qp different from Qp_zoo_pom must be routed differently than other
! elements to ensure that sinking detritus and zooplankton pools get their fixed P/C ratios.
! The remaining P is split evenly between DOP and PO4.
!-----------------------------------------------------------------------

      if (autotrophs(auto_ind)%Qp /= Qp_zoo_pom) then
         remaining_P = ((auto_graze(auto_ind) + auto_loss(auto_ind) + auto_agg(auto_ind)) * autotrophs(auto_ind)%Qp) &
                       - ((auto_graze_zoo(auto_ind)) * Qp_zoo_pom) &
                       - ((auto_graze_poc(auto_ind) + auto_loss_poc(auto_ind) + auto_agg(auto_ind)) * Qp_zoo_pom)
         remaining_P_dop(auto_ind) = (c1 - parm_labile_ratio) * remaining_P
         remaining_P_dip(auto_ind) = parm_labile_ratio * remaining_P
      endif

   end do

!-----------------------------------------------------------------------
!  get fractional factor for routing of zoo losses, based on food supply
!  more material is routed to large detrital pool when diatoms eaten
!-----------------------------------------------------------------------

   work1 = c0
   work2 = c0
   do auto_ind = 1, autotroph_cnt
      work1 = work1 + autotrophs(auto_ind)%f_zoo_detr * (auto_graze(auto_ind) + epsC * epsTinv)
      work2 = work2 + (auto_graze(auto_ind) + epsC * epsTinv)
   end do
   f_zoo_detr = work1 / work2

!-----------------------------------------------------------------------
!  0.01 small zoo threshold C concentration (mmol C/m^3)
!  zoo losses, scaled by Tfunc
!-----------------------------------------------------------------------

   C_loss_thres = f_loss_thres * loss_thres_zoo

   Zprime = max(zooC_loc(k,column) - C_loss_thres, 0.0_BGC_r8)

   zoo_loss = (parm_z_mort2_0 * Zprime**1.5_BGC_r8 + parm_z_mort_0 * Zprime) * Tfunc

   zoo_loss_doc = (c1 - parm_labile_ratio) * (c1 - f_zoo_detr) * zoo_loss
   zoo_loss_dic = parm_labile_ratio * (c1 - f_zoo_detr) * zoo_loss

!-----------------------------------------------------------------------
!  compute terms for DOM
!-----------------------------------------------------------------------

   DOC_prod = zoo_loss_doc + sum(auto_loss_doc) + sum(auto_graze_doc)
   DON_prod = Q * DOC_prod
   DOP_prod = Qp_zoo_pom * zoo_loss_doc
   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%Qp == Qp_zoo_pom) then
         DOP_prod = DOP_prod + autotrophs(auto_ind)%Qp * (auto_loss_doc(auto_ind) + auto_graze_doc(auto_ind))
      else
         DOP_prod = DOP_prod + remaining_P_dop(auto_ind)
      endif
   end do
   DOFe_prod = Qfe_zoo * zoo_loss_doc
   do auto_ind = 1, autotroph_cnt
      DOFe_prod = DOFe_prod + Qfe(auto_ind) * (auto_loss_doc(auto_ind) + auto_graze_doc(auto_ind))
   end do

   DOC_remin  = DOC_loc(k,column)  * DOC_reminR
   DON_remin  = DON_loc(k,column)  * DON_reminR
   DOFe_remin = DOFe_loc(k,column) * DOFe_reminR
   DOP_remin  = DOP_loc(k,column)  * DOP_reminR

!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG8: ',k,DOP_remin,DOP_loc(k,column),PAR_avg
!endif

!-----------------------------------------------------------------------
!  Refractory remin rate due to photochemistry
!  below euphotic zone remin rate sharply decrease
!-----------------------------------------------------------------------

   if (PAR_avg > 1.0_BGC_r8) then
      DONr_remin = DONr_loc(k,column) * DONr_reminR
      DOPr_remin = DOPr_loc(k,column) * DOPr_reminR
   else
      DONr_remin = DONr_loc(k,column) * (c1/(365.0_BGC_r8*670.0_BGC_r8)) * dps  ! 1/670 yrs
      DOPr_remin = DOPr_loc(k,column) * (c1/(365.0_BGC_r8*460.0_BGC_r8)) * dps  ! 1/460 yrs
      DOC_remin = DOC_remin * 0.0685_BGC_r8
      DON_remin = DON_remin * 0.1_BGC_r8
      DOFe_remin = DOFe_remin * 0.05_BGC_r8
      DOP_remin = DOP_remin * 0.05_BGC_r8
   end if

!-----------------------------------------------------------------------
!  large detritus C
!-----------------------------------------------------------------------

   POC%prod = f_zoo_detr * zoo_loss + sum(auto_graze_poc) &
                       + sum(auto_agg) + sum(auto_loss_poc)
!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG12: ',k,POC%prod, zoo_loss,auto_graze_poc,auto_agg,auto_loss_poc
!endif


!-----------------------------------------------------------------------
!  large detrital CaCO3
!  33% of CaCO3 is remin when phyto are grazed
!-----------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%CaCO3_ind > 0) then
         P_CaCO3%prod = ((c1 - f_graze_CaCO3_REMIN) * auto_graze(auto_ind) + &
                                  auto_loss(auto_ind) + auto_agg(auto_ind)) * QCaCO3(auto_ind)
      endif
   end do

!-----------------------------------------------------------------------
!  large detritus SiO2
!  grazed diatom SiO2, 60% is remineralized
!-----------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%Si_ind > 0) then
         P_SiO2%prod = Qsi(auto_ind) &
            * ((c1 - f_graze_si_remin) * auto_graze(auto_ind) + auto_agg(auto_ind) &
               + autotrophs(auto_ind)%loss_poc * auto_loss(auto_ind))
      endif
   end do

   dust%prod = c0

!-----------------------------------------------------------------------
!  Compute iron scavenging :
!  1) compute in terms of loss per year per unit iron (%/year/fe)
!  2) scale by sinking POMx10 + Dust + bSi + CaCO3 flux
!  3) increase scavenging at higher iron (>0.6nM)
!  4) convert to net loss per second
!-----------------------------------------------------------------------

   Fe_scavenge_rate = parm_Fe_scavenge_rate0

   Fe_scavenge_rate = Fe_scavenge_rate * &
      ((POC%sflux_out + POC%hflux_out) * 120.1_BGC_r8 + &
       (P_CaCO3%sflux_out + P_CaCO3%hflux_out) * P_CaCO3%mass + &
       (P_SiO2%sflux_out + P_SiO2%hflux_out) * P_SiO2%mass + &
       (dust%sflux_out + dust%hflux_out) * dust_fescav_scale)

   if (Fe_loc(k,column) > Fe_scavenge_thres1) &
      Fe_scavenge_rate = Fe_scavenge_rate + &
                         (Fe_loc(k,column) - Fe_scavenge_thres1) * fe_max_scale2

   Fe_scavenge = yps * Fe_loc(k,column) * Fe_scavenge_rate

   P_iron%prod = (zoo_loss * f_zoo_detr * Qfe_zoo) + Fe_scavenge

   do auto_ind = 1, autotroph_cnt
      P_iron%prod = P_iron%prod &
         + Qfe(auto_ind) * (auto_agg(auto_ind) + auto_graze_poc(auto_ind) + auto_loss_poc(auto_ind))
   end do

   call compute_particulate_terms(column, k, kmax, POC, P_CaCO3, P_SiO2, dust, P_iron, &
                                  QA_dust_def, BGC_input%PotentialTemperature(k,column), &
                                  O2_loc(k,column), NO3_loc(k,column), &
                                  SED_DENITRIF, OTHER_REMIN, BGC_input%cell_thickness(k,column), &
                                  BGC_input%cell_bottom_depth(k,column), BGC_forcing%FESEDFLUX(k,column),  &
                                  BGC_diagnostic_fields)
                                  

!-----------------------------------------------------------------------
!  nitrate & ammonium
!  nitrification in low light
!  use exponential decay of PAR across model level to compute taper factor
!------------------------------------------------------------------

   if (lrest_no3) then
      RESTORE = BGC_forcing%NUTR_RESTORE_RTAU(k,column) * &
                   (BGC_forcing%NO3_CLIM(k,column) - NO3_loc(k,column))
   else
      RESTORE = c0
   endif

   BGC_diagnostic_fields%diag_NO3_RESTORE(k,column) = RESTORE

   if (PAR_out < parm_nitrif_par_lim) then
      NITRIF = parm_kappa_nitrif * NH4_loc(k,column)
      if (PAR_in > parm_nitrif_par_lim) then
         NITRIF = NITRIF * log(PAR_out / parm_nitrif_par_lim) / (-KPARdz)
      end if
   else
      NITRIF = c0
   end if

   BGC_diagnostic_fields%diag_NITRIF(k,column) = NITRIF

!-----------------------------------------------------------------------
!  Compute denitrification under low O2 conditions
!-----------------------------------------------------------------------

   work1 = ((parm_o2_min + parm_o2_min_delta) - O2_loc(k,column)) / parm_o2_min_delta
   work1 = min(max(work1,0.0_BGC_r8),1.0_BGC_r8)

   work1 = merge(0.0_BGC_r8, work1, NO3_loc(k,column) == c0)

   DENITRIF = work1 * ((DOC_remin + POC%remin - OTHER_REMIN) / denitrif_C_N &
                       - SED_DENITRIF)

   BGC_diagnostic_fields%diag_DENITRIF(k,column) = DENITRIF

!-----------------------------------------------------------------------
!  nitrate & ammonium
!-----------------------------------------------------------------------

   BGC_output%BGC_tendencies(k,column,no3_ind) = RESTORE + NITRIF - DENITRIF - SED_DENITRIF - sum(NO3_V)

   BGC_output%BGC_tendencies(k,column,nh4_ind) = -sum(NH4_V) - NITRIF + DON_remin + DONr_remin &
      + Q * (zoo_loss_dic + sum(auto_loss_dic) + sum(auto_graze_dic) &
             + POC%remin * (c1 - DONrefract))

   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%Nfixer) &
         BGC_output%BGC_tendencies(k,column,nh4_ind) = BGC_output%BGC_tendencies(k,column,nh4_ind) + Nexcrete(auto_ind)
   end do

!-----------------------------------------------------------------------
!  dissolved iron
!-----------------------------------------------------------------------

     BGC_output%BGC_tendencies(k,column,fe_ind) = P_iron%remin &
       + (Qfe_zoo * zoo_loss_dic) + DOFe_remin - sum(photoFe) - Fe_scavenge

     do auto_ind = 1, autotroph_cnt
        BGC_output%BGC_tendencies(k,column,fe_ind) = BGC_output%BGC_tendencies(k,column,fe_ind) &
           + (Qfe(auto_ind) * (auto_loss_dic(auto_ind) + auto_graze_dic(auto_ind))) &
           + auto_graze_zoo(auto_ind) * (Qfe(auto_ind)-Qfe_zoo)
     end do

!-----------------------------------------------------------------------
!  dissolved SiO3
!-----------------------------------------------------------------------

   if (lrest_sio3) then
      RESTORE = BGC_forcing%NUTR_RESTORE_RTAU(k,column) * &
                   (BGC_forcing%SiO3_CLIM(k,column) - SiO3_loc(k,column))
   else
      RESTORE = c0
   endif

   BGC_diagnostic_fields%diag_SiO3_RESTORE(k,column) = RESTORE

   BGC_output%BGC_tendencies(k,column,sio3_ind) = RESTORE + P_SiO2%remin

   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%Si_ind > 0) then
         BGC_output%BGC_tendencies(k,column,sio3_ind) = BGC_output%BGC_tendencies(k,column,sio3_ind) - photoSi(auto_ind) &
            + Qsi(auto_ind) * (f_graze_si_remin * auto_graze(auto_ind) &
                                   + (c1 - autotrophs(auto_ind)%loss_poc) * auto_loss(auto_ind))
      endif
   end do

!-----------------------------------------------------------------------
!  phosphate
!-----------------------------------------------------------------------

   if (lrest_po4) then
      RESTORE = BGC_forcing%NUTR_RESTORE_RTAU(k,column) * &
                   (BGC_forcing%PO4_CLIM(k,column) - PO4_loc(k,column))
   else
      RESTORE = c0
   endif

   BGC_diagnostic_fields%diag_PO4_RESTORE(k,column) = RESTORE

   BGC_output%BGC_tendencies(k,column,po4_ind) = RESTORE + DOP_remin + DOPr_remin - sum(PO4_V) &
      + Qp_zoo_pom * ((c1 - DOPrefract) * POC%remin + zoo_loss_dic)

!maltrud debug
work4 = BGC_output%BGC_tendencies(k,column,po4_ind)
!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG9: ',k,POC%remin,zoo_loss_dic,zoo_loss,Zprime
!endif


   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%Qp == Qp_zoo_pom) then
         BGC_output%BGC_tendencies(k,column,po4_ind) = BGC_output%BGC_tendencies(k,column,po4_ind) &
            + autotrophs(auto_ind)%Qp * (auto_loss_dic(auto_ind) + auto_graze_dic(auto_ind))
      else
         BGC_output%BGC_tendencies(k,column,po4_ind) = BGC_output%BGC_tendencies(k,column,po4_ind) + remaining_P_dip(auto_ind)
      endif
   end do

!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG1: ',k,DOP_remin,DOPr_remin,PO4_V,sum(PO4_V),remaining_P_dip,work4,BGC_tendency(k,column,po4_ind)
!endif

!-----------------------------------------------------------------------
!  autotroph Carbon
!  autotroph Chlorophyll
!  autotroph Fe
!  autotroph Si
!  autotroph CaCO3
!-----------------------------------------------------------------------

   do auto_ind = 1, autotroph_cnt
      work1 = auto_graze(auto_ind) + auto_loss(auto_ind) + auto_agg(auto_ind)

      n = autotrophs(auto_ind)%C_ind
      BGC_output%BGC_tendencies(k,column,n) = photoC(auto_ind) - work1

      n = autotrophs(auto_ind)%Chl_ind
      BGC_output%BGC_tendencies(k,column,n) = photoacc(auto_ind) - thetaC(auto_ind) * work1

      n = autotrophs(auto_ind)%Fe_ind
      BGC_output%BGC_tendencies(k,column,n) =  photoFe(auto_ind) - Qfe(auto_ind) * work1

      n = autotrophs(auto_ind)%Si_ind
      if (n > 0) then
         BGC_output%BGC_tendencies(k,column,n) =  photoSi(auto_ind) - Qsi(auto_ind) * work1
      endif

      n = autotrophs(auto_ind)%CaCO3_ind
      if (n > 0) then
         BGC_output%BGC_tendencies(k,column,n) = CaCO3_PROD(auto_ind) - QCaCO3(auto_ind) * work1
      endif
   end do

!-----------------------------------------------------------------------
!  zoo Carbon
!-----------------------------------------------------------------------

   BGC_output%BGC_tendencies(k,column,zooC_ind) = sum(auto_graze_zoo) - zoo_loss

!-----------------------------------------------------------------------
!  dissolved organic Matter
!  from sinking remin small fraction to refractory pool
!-----------------------------------------------------------------------

   BGC_output%BGC_tendencies(k,column,doc_ind) = DOC_prod - DOC_remin

   BGC_output%BGC_tendencies(k,column,don_ind) = (DON_prod * (c1 - DONrefract)) - DON_remin

   BGC_output%BGC_tendencies(k,column,donr_ind) = (DON_prod * DONrefract) - DONr_remin &
      + (POC%remin * DONrefract * Q)

   BGC_output%BGC_tendencies(k,column,dop_ind) = (DOP_prod * (c1 - DOPrefract)) - DOP_remin &
      - sum(DOP_V)

   BGC_output%BGC_tendencies(k,column,dopr_ind) = (DOP_prod * DOPrefract) - DOPr_remin &
      + (POC%remin * DOPrefract * Qp_zoo_pom)

   BGC_output%BGC_tendencies(k,column,dofe_ind) = DOFe_prod - DOFe_remin

!-----------------------------------------------------------------------
!   dissolved inorganic Carbon
!-----------------------------------------------------------------------

   BGC_output%BGC_tendencies(k,column,dic_ind) = &
      sum(auto_loss_dic) + sum(auto_graze_dic) - sum(photoC) &
      + DOC_remin + POC%remin + zoo_loss_dic &
      + P_CaCO3%remin

   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%CaCO3_ind > 0) &
         BGC_output%BGC_tendencies(k,column,dic_ind) = BGC_output%BGC_tendencies(k,column,dic_ind) &
            + f_graze_CaCO3_REMIN * auto_graze(auto_ind) * QCaCO3(auto_ind) &
            - CaCO3_PROD(auto_ind)
   end do

   if (alt_co2_use_eco) then  ! set alt dic tendency to zero if we want it to be abiotic
     BGC_output%BGC_tendencies(k,column,dic_alt_co2_ind) = BGC_output%BGC_tendencies(k,column,dic_ind)
   else
     BGC_output%BGC_tendencies(k,column,dic_alt_co2_ind) = 0.0
   endif

!-----------------------------------------------------------------------
!  alkalinity
!-----------------------------------------------------------------------

   BGC_output%BGC_tendencies(k,column,alk_ind) = -BGC_output%BGC_tendencies(k,column,no3_ind) + &
      BGC_output%BGC_tendencies(k,column,nh4_ind) + c2 * P_CaCO3%remin

   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%CaCO3_ind > 0) &
         BGC_output%BGC_tendencies(k,column,alk_ind) = BGC_output%BGC_tendencies(k,column,alk_ind) &
            + c2 * (f_graze_CaCO3_REMIN * auto_graze(auto_ind) * QCaCO3(auto_ind) &
                    - CaCO3_PROD(auto_ind))
   end do

!-----------------------------------------------------------------------
!  oxygen
!-----------------------------------------------------------------------

   O2_PRODUCTION = c0

   do auto_ind = 1, autotroph_cnt
      if (.not. autotrophs(auto_ind)%Nfixer) then
         if (photoC(auto_ind) > c0) then
            O2_PRODUCTION = O2_PRODUCTION + photoC(auto_ind) * &
               ((NO3_V(auto_ind) / (NO3_V(auto_ind) + NH4_V(auto_ind))) / parm_Red_D_C_O2 + &
                (NH4_V(auto_ind) / (NO3_V(auto_ind) + NH4_V(auto_ind))) / parm_Remin_D_C_O2)
         end if
      else
         if (photoC(auto_ind) > c0) then
            O2_PRODUCTION = O2_PRODUCTION + photoC(auto_ind) * &
               ((NO3_V(auto_ind) / (NO3_V(auto_ind) + NH4_V(auto_ind) + Nfix(auto_ind))) / parm_Red_D_C_O2 + &
                (NH4_V(auto_ind) / (NO3_V(auto_ind) + NH4_V(auto_ind) + Nfix(auto_ind))) / parm_Remin_D_C_O2 + &
                (Nfix(auto_ind) / (NO3_V(auto_ind) + NH4_V(auto_ind) + Nfix(auto_ind))) / parm_Red_D_C_O2_diaz)
         end if
      endif
   end do

   work1 = (O2_loc(k,column) - parm_o2_min) / parm_o2_min_delta
   work1 = min(max(work1,0.0_BGC_r8),1.0_BGC_r8)
   O2_CONSUMPTION = work1 * &
      ((POC%remin + DOC_remin - (SED_DENITRIF*denitrif_C_N) - OTHER_REMIN + zoo_loss_dic &
        + sum(auto_loss_dic) + sum(auto_graze_dic)) / parm_Remin_D_C_O2 + (c2 * NITRIF))

   BGC_output%BGC_tendencies(k,column,o2_ind) = O2_PRODUCTION - O2_CONSUMPTION

!-----------------------------------------------------------------------
!  various tavg/history variables
!-----------------------------------------------------------------------

   BGC_diagnostic_fields%diag_O2_PRODUCTION(k,column) = O2_PRODUCTION

   BGC_diagnostic_fields%diag_O2_CONSUMPTION(k,column) = O2_CONSUMPTION

   work1 = O2SAT_singleValue(BGC_input%PotentialTemperature(k,column), BGC_input%Salinity(k,column))
   work1 = work1 - O2_loc(k,column)
   BGC_diagnostic_fields%diag_AOU(k,column) = work1

   BGC_diagnostic_fields%diag_PAR_avg(k,column) = PAR_avg

   BGC_diagnostic_fields%diag_zoo_loss(k,column) = zoo_loss

   work1 = sum(auto_graze)
   BGC_diagnostic_fields%diag_auto_graze_TOT(k,column) = work1

   do auto_ind = 1, autotroph_cnt
      BGC_diagnostic_fields%diag_auto_graze(k,column,auto_ind) = auto_graze(auto_ind)
      BGC_diagnostic_fields%diag_auto_loss(k,column,auto_ind) = auto_loss(auto_ind)
      BGC_diagnostic_fields%diag_auto_agg(k,column,auto_ind) = auto_agg(auto_ind)
      BGC_diagnostic_fields%diag_photoC(k,column,auto_ind) = photoC(auto_ind)
      work1 = BGC_input%cell_thickness(k,column) * photoC(auto_ind)
      BGC_diagnostic_fields%diag_photoC_zint(column,auto_ind) =  &
          BGC_diagnostic_fields%diag_photoC_zint(column,auto_ind) + work1
   end do

   work1 = sum(photoC)
   BGC_diagnostic_fields%diag_photoC_TOT(k,column) = work1
   work1 = work1 * BGC_input%cell_thickness(k,column)
   BGC_diagnostic_fields%diag_photoC_TOT_zint(column) =  &
      BGC_diagnostic_fields%diag_photoC_TOT_zint(column) + work1

   do auto_ind = 1, autotroph_cnt
      if (VNtot(auto_ind) > c0) then
         work1 = (VNO3(auto_ind) / VNtot(auto_ind)) * photoC(auto_ind)
      else
         work1 = c0
      end if

      BGC_diagnostic_fields%diag_photoC_NO3(k,column,auto_ind) = work1

      work1 = work1 * BGC_input%cell_thickness(k,column)
      BGC_diagnostic_fields%diag_photoC_NO3_zint(column,auto_ind) = &
          BGC_diagnostic_fields%diag_photoC_NO3_zint(column,auto_ind) + work1

      BGC_diagnostic_fields%diag_photoC_NO3_TOT(k,column) = &
          BGC_diagnostic_fields%diag_photoC_NO3_TOT(k,column) +  &
          BGC_diagnostic_fields%diag_photoC_NO3(k,column,auto_ind)

      BGC_diagnostic_fields%diag_photoC_NO3_TOT_zint(column) = &
          BGC_diagnostic_fields%diag_photoC_NO3_TOT_zint(column) +  &
          BGC_diagnostic_fields%diag_photoC_NO3_zint(column,auto_ind)

   end do

   BGC_diagnostic_fields%diag_DOC_prod(k,column) = DOC_prod

   BGC_diagnostic_fields%diag_DOC_remin(k,column) = DOC_remin

   BGC_diagnostic_fields%diag_DON_prod(k,column) = DON_prod

   BGC_diagnostic_fields%diag_DON_remin(k,column) = DON_remin

   BGC_diagnostic_fields%diag_DOP_prod(k,column) = DOP_prod

   BGC_diagnostic_fields%diag_DOP_remin(k,column) = DOP_remin

   BGC_diagnostic_fields%diag_DOFe_prod(k,column) = DOFe_prod

   BGC_diagnostic_fields%diag_DOFe_remin(k,column) = DOFe_remin

   BGC_diagnostic_fields%diag_Fe_scavenge(k,column) = Fe_scavenge

   BGC_diagnostic_fields%diag_Fe_scavenge_rate(k,column) = Fe_scavenge_rate

   ztop = c0
   if (k > 1) ztop = BGC_input%cell_bottom_depth(k-1,column)
   work2 = min(100.0e2_BGC_r8 - ztop, BGC_input%cell_thickness(k,column))
   partial_thickness_100m = merge(work2, 0.0_BGC_r8, work2 > c0)

   work1 = BGC_output%BGC_tendencies(k,column,dic_ind) + BGC_output%BGC_tendencies(k,column,doc_ind) &
           + BGC_output%BGC_tendencies(k,column,zooC_ind) &
           + sum(BGC_output%BGC_tendencies(k,column,autotrophs(:)%C_ind), dim=1)
   do auto_ind = 1, autotroph_cnt
      n = autotrophs(auto_ind)%CaCO3_ind
      if (n > 0) then
         work1 = work1 + BGC_output%BGC_tendencies(k,column,n)
      endif
   end do

   BGC_diagnostic_fields%diag_Jint_Ctot(column) = BGC_diagnostic_fields%diag_Jint_Ctot(column) + &
      work1*BGC_input%cell_thickness(k,column) + POC%sed_loss + P_CaCO3%sed_loss

   BGC_diagnostic_fields%diag_Jint_100m_Ctot(column) =  &
      BGC_diagnostic_fields%diag_Jint_100m_Ctot(column) + work1*partial_thickness_100m + &
      merge(POC%sed_loss + P_CaCO3%sed_loss, 0.0_BGC_r8, BGC_input%cell_bottom_depth(k,column) <= 100.0e2_BGC_r8)

   work1 = BGC_output%BGC_tendencies(k,column,no3_ind) + BGC_output%BGC_tendencies(k,column,nh4_ind) &
           + BGC_output%BGC_tendencies(k,column,don_ind) + BGC_output%BGC_tendencies(k,column,donr_ind) &
           + Q * BGC_output%BGC_tendencies(k,column,zooC_ind) &
           + Q * sum(BGC_output%BGC_tendencies(k,column,autotrophs(:)%C_ind), dim=1)
   ! add back column and sediment denitrification
   work1 = work1 + DENITRIF + SED_DENITRIF
   ! subtract out N fixation
   do auto_ind = 1, autotroph_cnt
      if (autotrophs(auto_ind)%Nfixer) work1 = work1 - Nfix(auto_ind)
   end do

   BGC_diagnostic_fields%diag_Jint_Ntot(column) = BGC_diagnostic_fields%diag_Jint_Ntot(column) + &
      work1*BGC_input%cell_thickness(k,column) + POC%sed_loss * Q

   BGC_diagnostic_fields%diag_Jint_100m_Ntot(column) =  &
      BGC_diagnostic_fields%diag_Jint_100m_Ntot(column) + work1*partial_thickness_100m + &
      merge(POC%sed_loss * Q, 0.0_BGC_r8, BGC_input%cell_bottom_depth(k,column) <= 100.0e2_BGC_r8)

   work1 = BGC_output%BGC_tendencies(k,column,po4_ind) + BGC_output%BGC_tendencies(k,column,dop_ind) &
           + BGC_output%BGC_tendencies(k,column,dopr_ind) &
           + Qp_zoo_pom * BGC_output%BGC_tendencies(k,column,zooC_ind)
   do auto_ind = 1, autotroph_cnt
      n = autotrophs(auto_ind)%C_ind
      work1 = work1 + autotrophs(auto_ind)%Qp * BGC_output%BGC_tendencies(k,column,n)
   end do

   BGC_diagnostic_fields%diag_Jint_Ptot(column) = BGC_diagnostic_fields%diag_Jint_Ptot(column) + &
      work1*BGC_input%cell_thickness(k,column) + POC%sed_loss * Qp_zoo_pom

   BGC_diagnostic_fields%diag_Jint_100m_Ptot(column) =  &
      BGC_diagnostic_fields%diag_Jint_100m_Ptot(column) + work1*partial_thickness_100m + &
      merge(POC%sed_loss * Qp_zoo_pom, 0.0_BGC_r8, BGC_input%cell_bottom_depth(k,column) <= 100.0e2_BGC_r8)

   work1 = BGC_output%BGC_tendencies(k,column,sio3_ind)
   do auto_ind = 1, autotroph_cnt
      n = autotrophs(auto_ind)%Si_ind
      if (n > 0) then
         work1 = work1 + BGC_output%BGC_tendencies(k,column,n)
      endif
   end do

   BGC_diagnostic_fields%diag_Jint_Sitot(column) = BGC_diagnostic_fields%diag_Jint_Sitot(column) + &
      work1*BGC_input%cell_thickness(k,column) + P_SiO2%sed_loss

   BGC_diagnostic_fields%diag_Jint_100m_Sitot(column) =  &
      BGC_diagnostic_fields%diag_Jint_100m_Sitot(column) + work1*partial_thickness_100m + &
      merge(P_SiO2%sed_loss, 0.0_BGC_r8, BGC_input%cell_bottom_depth(k,column) <= 100.0e2_BGC_r8)

   do auto_ind = 1, autotroph_cnt
      work1 = autotrophChl_loc(k,column,auto_ind) * partial_thickness_100m
      BGC_diagnostic_fields%diag_Chl_TOT_zint_100m(column) = &
          BGC_diagnostic_fields%diag_Chl_TOT_zint_100m(column) +  &
          autotrophChl_loc(k,column,auto_ind) * partial_thickness_100m
   end do

   enddo ! k loop

! find O2 minimum
   ! work1 = O2 at this level
   ! work2 = vertical min of O2
   ! work3 = depth of min

   k = 1
   work1 = O2_loc(k,column)
   work2 = work1
   work3 = BGC_input%cell_center_depth(k,column)

   do k = 2,kmax
      work1 = O2_loc(k,column)
      if (work1 < work2) then
         work2 = work1
         work3 = BGC_input%cell_center_depth(k,column)
      endif
   end do

   BGC_diagnostic_fields%diag_O2_ZMIN(column) = work2
   BGC_diagnostic_fields%diag_O2_ZMIN_DEPTH(column) = work3

   enddo column_loop ! i loop

   deallocate(DIC_loc)
   deallocate(DIC_ALT_CO2_loc)
   deallocate(ALK_loc)
   deallocate(PO4_loc)
   deallocate(NO3_loc)
   deallocate(SiO3_loc)
   deallocate(NH4_loc)
   deallocate(Fe_loc)
   deallocate(O2_loc)
   deallocate(DOC_loc)
   deallocate(zooC_loc)
   deallocate(DON_loc)
   deallocate(DOFe_loc)
   deallocate(DOP_loc)
   deallocate(DOPr_loc)
   deallocate(DONr_loc)

   deallocate(autotrophChl_loc)
   deallocate(autotrophC_loc)
   deallocate(autotrophFe_loc)
   deallocate(autotrophSi_loc)
   deallocate(autotrophCaCO3_loc)

!-----------------------------------------------------------------------
!EOC

 end subroutine BGC_SourceSink

!***********************************************************************
!***********************************************************************
!BOP
! !IROUTINE: init_particulate_terms
! !INTERFACE:

 subroutine init_particulate_terms(POC, P_CaCO3, P_SiO2, dust, P_iron, &
       QA_dust_def, NET_DUST_IN)

! !DESCRIPTION:
!  Set incoming fluxes (put into outgoing flux for first level usage).
!  Set dissolution length, production fraction and mass terms.
!
!  The first 6 arguments are intent(inout) in
!  order to preserve contents on other blocks.

! !INPUT/OUTPUT PARAMETERS:

   type(sinking_particle), intent(inout) :: &
      POC,          & ! base units = nmol C
      P_CaCO3,      & ! base units = nmol CaCO3
      P_SiO2,       & ! base units = nmol SiO2
      dust,         & ! base units = g
      P_iron          ! base units = nmol Fe

   real (BGC_r8), intent(inout) :: &
      QA_dust_def     ! incoming deficit in the QA(dust) POC flux

! !INPUT PARAMETERS:

   real (BGC_r8), intent(in) :: &
      NET_DUST_IN     ! dust flux

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  parameters, from Armstrong et al. 2000
!
!  July 2002, length scale for excess POC and bSI modified by temperature
!  Value given here is at Tref of 30 deg. C, JKM
!-----------------------------------------------------------------------

    POC%diss      = parm_POC_diss   ! diss. length (cm), modified by TEMP
    POC%gamma     = c0              ! not used
    POC%mass      = 12.01_BGC_r8        ! molecular weight of POC
    POC%rho       = c0              ! not used

    P_CaCO3%diss  = parm_CaCO3_diss ! diss. length (cm)
    P_CaCO3%gamma = 0.30_BGC_r8         ! prod frac -> hard subclass
    P_CaCO3%mass  = 100.09_BGC_r8       ! molecular weight of CaCO
    P_CaCO3%rho   = 0.05_BGC_r8 * P_CaCO3%mass / POC%mass ! QA mass ratio for CaCO3

    P_SiO2%diss   = parm_SiO2_diss  ! diss. length (cm), modified by TEMP
    P_SiO2%gamma  = 0.030_BGC_r8        ! prod frac -> hard subclass
    P_SiO2%mass   = 60.08_BGC_r8        ! molecular weight of SiO2
    P_SiO2%rho    = 0.05_BGC_r8 * P_SiO2%mass / POC%mass ! QA mass ratio for SiO2

    dust%diss     = 20000.0_BGC_r8      ! diss. length (cm)
    dust%gamma    = 0.97_BGC_r8         ! prod frac -> hard subclass
    dust%mass     = 1.0e9_BGC_r8        ! base units are already grams
    dust%rho      = 0.05_BGC_r8 * dust%mass / POC%mass ! QA mass ratio for dust

    P_iron%diss   = 60000.0_BGC_r8      ! diss. length (cm) - not used
    P_iron%gamma  = c0              ! prod frac -> hard subclass - not used
    P_iron%mass   = c0              ! not used
    P_iron%rho    = c0              ! not used

!-----------------------------------------------------------------------
!  Set incoming fluxes
!-----------------------------------------------------------------------

    P_CaCO3%sflux_out = c0
    P_CaCO3%hflux_out = c0

    P_SiO2%sflux_out = c0
    P_SiO2%hflux_out = c0

    if (NET_DUST_IN /= c0) then
       dust%sflux_out = (c1 - dust%gamma) * NET_DUST_IN
       dust%hflux_out = dust%gamma * NET_DUST_IN
    else
       dust%sflux_out = c0
       dust%hflux_out = c0
    endif

    P_iron%sflux_out = c0
    P_iron%hflux_out = c0

!-----------------------------------------------------------------------
!  Hard POC is QA flux and soft POC is excess POC.
!-----------------------------------------------------------------------

    POC%sflux_out = c0
    POC%hflux_out = c0

!-----------------------------------------------------------------------
!  Compute initial QA(dust) POC flux deficit.
!-----------------------------------------------------------------------

    QA_dust_def = dust%rho * &
       (dust%sflux_out + dust%hflux_out)

!-----------------------------------------------------------------------
!EOC

 end subroutine init_particulate_terms

!***********************************************************************
!BOP
! !IROUTINE: compute_particulate_terms
! !INTERFACE:

 subroutine compute_particulate_terms(column, k, kmax, POC, P_CaCO3, P_SiO2, dust, P_iron, &
       QA_dust_def, TEMP, O2_loc, NO3_loc, SED_DENITRIF, OTHER_REMIN, &
       cell_thickness, cell_bottom_depth, FESEDFLUX_loc, BGC_diagnostic_fields)

! !DESCRIPTION:
!  Compute outgoing fluxes and remineralization terms. Assumes that
!  production terms have been set. Incoming fluxes are assumed to be the
!  outgoing fluxes from the previous level.
!
!  It is assumed that there is no production of dust.
!
!  Instantaneous remineralization in the bottom cell is implemented by
!  setting the outgoing flux to zero.
!
!  For POC, the hard subclass is the POC flux qualitatively associated
!  with the ballast flux. The soft subclass is the excess POC flux.
!
!  Remineralization for the non-iron particulate pools is computing
!  by first computing the outgoing flux and then computing the
!  remineralization from conservation, i.e.
!     flux_in - flux_out + prod * dz - remin * dz == 0.
!
!  For iron, remineralization is first computed from POC remineralization
!  and then flux_out is computed from conservation. If the resulting
!  flux_out is negative or should be zero because of the sea floor, the
!  remineralization is adjusted.
!  Note: all the sinking iron is in the P_iron%sflux pool, hflux Fe not
!        explicitly tracked, it is assumed that total iron remin is
!        proportional to total POC remin.
!
!  Based upon Armstrong et al. 2000
!
!  July 2002, added temperature effect on remin length scale of
!  excess POC (all soft POM& Iron) and on SiO2.
!  new variable passed into ballast, Tfunc, main Temperature function
!  computed in ecosystem routine.  scaling factor for dissolution
!  of excess POC, Fe, and Bsi now varies with location (f(temperature)).
!
!  Added diffusive iron flux from sediments at depths < 1100m,
!  based on Johnson et al., 1999, value of 5 umolFe/m2/day,
!      this value too high, using 2 umolFe/m2/day here
!
!  Allow hard fraction of ballast to remin with long length scale 40,000m
!     thus ~ 10% of hard ballast remins over 4000m water column.
!
!  Sinking dust flux is decreased by assumed instant solubility/dissolution
!     at ocean surface from the parm_Fe_bioavail.
!
!  Modified to allow different Q10 factors for soft POM and bSI remin,
!  water TEMP is now passed in instead of Tfunc (1/2005, JKM)

! !USES:

! !INPUT PARAMETERS:

   integer (BGC_i4), intent(in) :: column, k, kmax ! vertical model level

   real (BGC_r8), intent(in) :: &
      TEMP,         & ! temperature for scaling functions bsi%diss
      O2_loc,       & ! dissolved oxygen used to modify POC%diss, Sed fluxes
      NO3_loc         ! dissolved nitrate used to modify sed fluxes

   real (BGC_r8), intent(in) :: &
      cell_thickness,         & !
      cell_bottom_depth         !

   real (BGC_r8), intent(in) :: &
      FESEDFLUX_loc          ! sediment Fe flux at level k

! !INPUT/OUTPUT PARAMETERS:

   type(sinking_particle), intent(inout) :: &
      POC,          & ! base units = nmol C
      P_CaCO3,      & ! base units = nmol CaCO3
      P_SiO2,       & ! base units = nmol SiO2
      dust,         & ! base units = g
      P_iron          ! base units = nmol Fe

   real (BGC_r8), intent(inout) :: &
      QA_dust_def,  & ! incoming deficit in the QA(dust) POC flux
      SED_DENITRIF, & ! sedimentary denitrification (umolN/cm^2/s)
      OTHER_REMIN     ! sedimentary remin not due to oxic or denitrification

   type(BGC_diagnostics_type), intent(inout) :: BGC_diagnostic_fields

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (BGC_r8) :: poc_diss, & ! diss. length used (cm)
                sio2_diss,& ! diss. length varies spatially with O2
                caco3_diss,&
                dust_diss

   real (BGC_r8) :: &
      work,               & ! temporary for summed quantities to be averaged
      TfuncS,             & ! temperature scaling from soft POM remin
      scalelength,        & ! used to scale dissolution length scales
      DECAY_Hard,         & ! scaling factor for dissolution of Hard Ballast
      DECAY_HardDust        ! scaling factor for dissolution of Hard dust

   real (BGC_r8) :: &
      decay_POC_E,        & ! scaling factor for dissolution of excess POC
      decay_SiO2,         & ! scaling factor for dissolution of SiO2
      decay_CaCO3,        & ! scaling factor for dissolution of CaCO3
      decay_dust,         & ! scaling factor for dissolution of dust
      POC_PROD_avail,     & ! POC production available for excess POC flux
      new_QA_dust_def,    & ! outgoing deficit in the QA(dust) POC flux
      flux, flux_alt,     & ! temp variables used to update sinking flux
      dz_loc, dzr_loc       ! dz, dzr at a particular i,j location

   integer (BGC_i4) :: &
      n             ! loop indices

   logical (BGC_log) :: &
      poc_error             ! POC error flag

   real (BGC_r8), parameter :: &
      mpercm = 0.01_BGC_r8   !  meters/cm

!-----------------------------------------------------------------------
!  incoming fluxes are outgoing fluxes from previous level
!-----------------------------------------------------------------------

   P_CaCO3%sflux_in = P_CaCO3%sflux_out
   P_CaCO3%hflux_in = P_CaCO3%hflux_out

   P_SiO2%sflux_in = P_SiO2%sflux_out
   P_SiO2%hflux_in = P_SiO2%hflux_out

   dust%sflux_in = dust%sflux_out
   dust%hflux_in = dust%hflux_out

   POC%sflux_in = POC%sflux_out
   POC%hflux_in = POC%hflux_out

   P_iron%sflux_in = P_iron%sflux_out
   P_iron%hflux_in = P_iron%hflux_out

!-----------------------------------------------------------------------
!  initialize loss to sediments = 0 and local copy of percent sed
!-----------------------------------------------------------------------

   P_iron%sed_loss = c0
   POC%sed_loss = c0
   P_CaCO3%sed_loss = c0
   P_SiO2%sed_loss = c0
   dust%sed_loss = c0
   SED_DENITRIF=c0
   OTHER_REMIN=c0

!-----------------------------------------------------------------------
!  compute scalelength and decay factors
!-----------------------------------------------------------------------

   if (cell_bottom_depth < parm_scalelen_z(1)) then
      scalelength = parm_scalelen_vals(1)
   else if (cell_bottom_depth >= parm_scalelen_z(size(parm_scalelen_z))) then
      scalelength = parm_scalelen_vals(size(parm_scalelen_z))
   else
      do n = 2, size(parm_scalelen_z)
         if (cell_bottom_depth < parm_scalelen_z(n)) then
            scalelength = parm_scalelen_vals(n-1) &
               + (parm_scalelen_vals(n) - parm_scalelen_vals(n-1)) &
                 * (cell_bottom_depth - parm_scalelen_z(n-1))/(parm_scalelen_z(n) - parm_scalelen_z(n-1))
            exit
         endif
      end do
   endif

   DECAY_Hard     = exp(-cell_thickness / 4.0e6_BGC_r8)
   DECAY_HardDust = exp(-cell_thickness / 1.2e7_BGC_r8)

!----------------------------------------------------------------------
!   Tref = 30.0 reference temperature (deg. C)
!-----------------------------------------------------------------------

   TfuncS = 1.5_BGC_r8**(((TEMP + T0_Kelvin_BGC) - (Tref + T0_Kelvin_BGC)) / c10)

   poc_error = .false.

   dz_loc = cell_thickness
   dzr_loc = c1 / dz_loc

   poc_diss = POC%diss
   sio2_diss = P_SiO2%diss
   caco3_diss = P_CaCO3%diss
   dust_diss = dust%diss

!-----------------------------------------------------------------------
!  increase POC diss length scale where O2 concentrations are low
!-----------------------------------------------------------------------

   if ((O2_loc >= 5.0_BGC_r8) .and. (O2_loc < 40.0_BGC_r8)) then
      poc_diss = POC%diss*(c1+(3.3_BGC_r8-c1)*(40.0_BGC_r8 - O2_loc)/35.0_BGC_r8)
   else if (O2_loc < 5.0_BGC_r8) then
      poc_diss = POC%diss * 3.3_BGC_r8
   endif

!-----------------------------------------------------------------------
!  apply scalelength factor to length scales
!-----------------------------------------------------------------------

   poc_diss = scalelength * poc_diss
   sio2_diss = scalelength * sio2_diss
   caco3_diss = scalelength * caco3_diss
   dust_diss = scalelength * dust_diss

!-----------------------------------------------------------------------
!  apply temperature dependence to sio2_diss length scale
!-----------------------------------------------------------------------

   sio2_diss = sio2_diss / TfuncS

!-----------------------------------------------------------------------
!  decay_POC_E and decay_SiO2 set locally, modified by O2
!-----------------------------------------------------------------------

   decay_POC_E = exp(-dz_loc / poc_diss)
   decay_SiO2  = exp(-dz_loc / sio2_diss)
   decay_CaCO3 = exp(-dz_loc / caco3_diss)
   decay_dust  = exp(-dz_loc / dust_diss)

!-----------------------------------------------------------------------
!  Set outgoing fluxes for non-iron pools.
!  The outoing fluxes for ballast materials are from the
!  solution of the coresponding continuous ODE across the model
!  level. The ODE has a constant source term and linear decay.
!  It is assumed that there is no sub-surface dust production.
!-----------------------------------------------------------------------

   P_CaCO3%sflux_out = P_CaCO3%sflux_in * decay_caco3 + &
      P_CaCO3%prod * ((c1 - P_CaCO3%gamma) * (c1 - decay_caco3) &
         * caco3_diss)

   P_CaCO3%hflux_out = P_CaCO3%hflux_in * DECAY_Hard + &
      P_CaCO3%prod * (P_CaCO3%gamma * dz_loc)

   P_SiO2%sflux_out = P_SiO2%sflux_in * decay_SiO2 + &
      P_SiO2%prod * ((c1 - P_SiO2%gamma) * (c1 - decay_SiO2) &
         * sio2_diss)

   P_SiO2%hflux_out = P_SiO2%hflux_in * DECAY_Hard + &
      P_SiO2%prod * (P_SiO2%gamma * dz_loc)

   dust%sflux_out = dust%sflux_in * decay_dust

   dust%hflux_out = dust%hflux_in * DECAY_HardDust

!-----------------------------------------------------------------------
!  Compute how much POC_PROD is available for deficit reduction
!  and excess POC flux after subtracting off fraction of non-dust
!  ballast production from net POC_PROD.
!-----------------------------------------------------------------------

   POC_PROD_avail = POC%prod - &
      P_CaCO3%rho * P_CaCO3%prod - &
      P_SiO2%rho * P_SiO2%prod

!-----------------------------------------------------------------------
!  Check for POC production bounds violations
!-----------------------------------------------------------------------

   if (POC_PROD_avail < c0) then
      poc_error = .true.
   endif

!-----------------------------------------------------------------------
!  Compute 1st approximation to new QA_dust_def, the QA_dust
!  deficit leaving the cell. Ignore POC_PROD_avail at this stage.
!-----------------------------------------------------------------------

   if (QA_dust_def > 0) then
      new_QA_dust_def = QA_dust_def * &
         (dust%sflux_out + dust%hflux_out) / &
         (dust%sflux_in + dust%hflux_in)
   else
      new_QA_dust_def = c0
   endif

!-----------------------------------------------------------------------
!  Use POC_PROD_avail to reduce new_QA_dust_def.
!-----------------------------------------------------------------------

   if (new_QA_dust_def > c0) then
      new_QA_dust_def = new_QA_dust_def - POC_PROD_avail * dz_loc
      if (new_QA_dust_def < c0) then
         POC_PROD_avail = -new_QA_dust_def * dzr_loc
         new_QA_dust_def = c0
      else
         POC_PROD_avail = c0
      endif
   endif

   QA_dust_def = new_QA_dust_def

!-----------------------------------------------------------------------
!  Compute outgoing POC fluxes. QA POC flux is computing using
!  ballast fluxes and new_QA_dust_def. If no QA POC flux came in
!  and no production occured, then no QA POC flux goes out. This
!  shortcut is present to avoid roundoff cancellation errors from
!  the dust%rho * dust_flux_out - QA_dust_def computation.
!  Any POC_PROD_avail still remaining goes into excess POC flux.
!-----------------------------------------------------------------------

   if (POC%hflux_in == c0 .and. POC%prod == c0) then
      POC%hflux_out = c0
   else
      POC%hflux_out = P_CaCO3%rho * &
         (P_CaCO3%sflux_out + P_CaCO3%hflux_out) + &
         P_SiO2%rho * &
         (P_SiO2%sflux_out + P_SiO2%hflux_out) + &
         dust%rho * &
         (dust%sflux_out + dust%hflux_out) - &
         new_QA_dust_def
      POC%hflux_out = max(POC%hflux_out, 0.0_BGC_r8)
   endif

   POC%sflux_out = POC%sflux_in * decay_POC_E + &
      POC_PROD_avail *((c1 - decay_POC_E) * &
      poc_diss)

!-----------------------------------------------------------------------
!  Compute remineralization terms. It is assumed that there is no
!  sub-surface dust production.
!-----------------------------------------------------------------------

   P_CaCO3%remin = P_CaCO3%prod + &
      ((P_CaCO3%sflux_in - P_CaCO3%sflux_out) + &
      (P_CaCO3%hflux_in - P_CaCO3%hflux_out)) * dzr_loc

   P_SiO2%remin = P_SiO2%prod + &
      ((P_SiO2%sflux_in - P_SiO2%sflux_out) + &
      (P_SiO2%hflux_in - P_SiO2%hflux_out)) * dzr_loc

   POC%remin = POC%prod + &
      ((POC%sflux_in - POC%sflux_out) + &
      (POC%hflux_in - POC%hflux_out)) * dzr_loc
!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG10: ',k,POC%remin, POC%prod,POC%sflux_in,POC%sflux_out,POC%hflux_in,POC%hflux_out
!endif

   dust%remin = &
      ((dust%sflux_in - dust%sflux_out) + &
      (dust%hflux_in - dust%hflux_out)) * dzr_loc

!-----------------------------------------------------------------------
!  Compute iron remineralization and flux out.
!-----------------------------------------------------------------------

   if (POC%sflux_in + POC%hflux_in == c0) then
      P_iron%remin = (POC%remin * parm_Red_Fe_C)
   else
      P_iron%remin = (POC%remin * &
         (P_iron%sflux_in + P_iron%hflux_in) / &
         (POC%sflux_in + POC%hflux_in))
   endif
   P_iron%remin = P_iron%remin +                &
                  (P_iron%sflux_in * 1.5e-5_BGC_r8)

   P_iron%sflux_out = P_iron%sflux_in + dz_loc * &
      ((c1 - P_iron%gamma) * P_iron%prod - P_iron%remin)

   if (P_iron%sflux_out < c0) then
      P_iron%sflux_out = c0
      P_iron%remin = P_iron%sflux_in * dzr_loc + &
         (c1 - P_iron%gamma) * P_iron%prod
   endif

!-----------------------------------------------------------------------
!  Compute iron release from dust remin/dissolution
!
!  dust remin gDust = 0.035 / 55.847 * 1.0e9 = 626712.0 nmolFe
!                      gFe     molFe     nmolFe
!  Also add in Fe source from sediments if applicable to this cell.
!-----------------------------------------------------------------------


   P_iron%remin = P_iron%remin &
      + dust%remin * dust_to_Fe &
      + (FESEDFLUX_loc * dzr_loc)

   P_iron%hflux_out = P_iron%hflux_in

!-----------------------------------------------------------------------
!  Bottom Sediments Cell?
!  If so compute sedimentary burial and denitrification N losses.
!  Using empirical relations from Bohlen et al., 2012 (doi:10.1029/2011GB004198) for Sed Denitrification
!  OTHER_REMIN estimates organic matter remineralized in the sediments
!      by the processes other than oxic remin and denitrification (SO4 and CO2,
!      etc..)
!      based on Soetaert et al., 1996, varies between 10% and 50%
!      0.4_r8 is a coefficient with units mmolC/cm2/yr sinking flux,
!      OTHER_REMIN is 50% above this high flux value,
!      In special case where bottom O2 has been depleted to < 1.0 uM,
!               all sedimentary remin is due to DENITRIFICATION + OTHER_REMIN
!  POC burial from Dunne et al. 2007 (doi:10.1029/2006GB002907), maximum of 80% burial efficiency imposed
!  Bsi preservation in sediments based on
!     Ragueneau et al. 2000 (doi:10.1016/S0921-8181(00)00052-7)
!  Calcite is preserved in sediments above the lysocline, dissolves below.
!       Here a constant depth is used for lysocline.
!-----------------------------------------------------------------------

 if (k == kmax) then

   flux = POC%sflux_out+POC%hflux_out

   if (flux > c0) then
      flux_alt = flux*mpercm*spd ! convert to mmol/m^2/day

      POC%sed_loss = flux * min(0.8_BGC_r8, parm_POMbury &
         * (0.013_BGC_r8 + 0.53_BGC_r8 * flux_alt*flux_alt / (7.0_BGC_r8 + flux_alt)**2))

      SED_DENITRIF = dzr_loc * flux &
         * (0.06_BGC_r8 + 0.19_BGC_r8 * 0.99_BGC_r8**(O2_loc-NO3_loc))

!maltrud prevent denitrification if NO# concentration is too low
      if (NO3_loc < 5.0_BGC_r8) SED_DENITRIF = 0.

      flux_alt = flux*1.0e-6_BGC_r8*spd*365.0_BGC_r8 ! convert to mmol/cm^2/year
      OTHER_REMIN = dzr_loc &
         * min(min(0.1_BGC_r8 + flux_alt,0.5_BGC_r8) * (flux - POC%sed_loss), &
               (flux - POC%sed_loss - (SED_DENITRIF*dz_loc*denitrif_C_N)))

!----------------------------------------------------------------------------------
!              if bottom water O2 is depleted, assume all remin is denitrif + other
!----------------------------------------------------------------------------------

      if (O2_loc < c1) then
         OTHER_REMIN = dzr_loc * &
         (flux - POC%sed_loss - (SED_DENITRIF*dz_loc*denitrif_C_N))
      endif

   endif

   flux = P_SiO2%sflux_out+P_SiO2%hflux_out
   flux_alt = flux*mpercm*spd ! convert to mmol/m^2/day
   ! first compute burial efficiency, then compute loss to sediments
   if (flux_alt > c2) then
      P_SiO2%sed_loss = 0.2_BGC_r8
   else
      P_SiO2%sed_loss = 0.04_BGC_r8
   endif
   P_SiO2%sed_loss = flux * parm_BSIbury * P_SiO2%sed_loss

   if (cell_bottom_depth < 3300.0e2_BGC_r8) then
      flux = P_CaCO3%sflux_out + P_CaCO3%hflux_out
      P_CaCO3%sed_loss = flux
   endif

!----------------------------------------------------------------------------------
!  Update sinking fluxes and remin fluxes, accounting for sediments.
!  flux used to hold sinking fluxes before update.
!----------------------------------------------------------------------------------

   flux = P_CaCO3%sflux_out + P_CaCO3%hflux_out
   if (flux > c0) then
      P_CaCO3%remin = P_CaCO3%remin &
         + ((flux - P_CaCO3%sed_loss) * dzr_loc)
   endif

   flux = P_SiO2%sflux_out + P_SiO2%hflux_out
   if (flux > c0) then
      P_SiO2%remin = P_SiO2%remin &
         + ((flux - P_SiO2%sed_loss) * dzr_loc)
   endif

   flux = POC%sflux_out + POC%hflux_out
   if (flux > c0) then
      POC%remin = POC%remin &
         + ((flux - POC%sed_loss) * dzr_loc)
   endif
!maltrud debug
!if(i==20.and.j==20.and.my_task==0)then
!write(*,*)'DEBUG11: ',k,POC%remin, POC%sed_loss,POC%sflux_out,POC%hflux_out,flux
!endif


!-----------------------------------------------------------------------
!   Remove all Piron and dust that hits bottom, sedimentary Fe source
!        accounted for by FESEDFLUX_loc elsewhere.
!-----------------------------------------------------------------------

   flux = (P_iron%sflux_out + P_iron%hflux_out)
   if (flux > c0) then
      P_iron%sed_loss = flux
   endif

   dust%sed_loss = dust%sflux_out + dust%hflux_out

!-----------------------------------------------------------------------
!   Set all outgoing fluxes to 0.0
!-----------------------------------------------------------------------

!maltrud not sure we need this if-test since it is inside a k==kmax already
!           if (k == kmax) then
      P_CaCO3%sflux_out = c0
      P_CaCO3%hflux_out = c0

      P_SiO2%sflux_out = c0
      P_SiO2%hflux_out = c0

      dust%sflux_out = c0
      dust%hflux_out = c0

      POC%sflux_out = c0
      POC%hflux_out = c0

      P_iron%sflux_out = c0
      P_iron%hflux_out = c0
!           endif

 endif  ! top level k == kmax

!-----------------------------------------------------------------------
!  Set tavg variables.
!-----------------------------------------------------------------------

   work = POC%sflux_in + POC%hflux_in
   BGC_diagnostic_fields%diag_POC_FLUX_IN(k,column) = work

   BGC_diagnostic_fields%diag_POC_PROD(k,column) = POC%prod

   BGC_diagnostic_fields%diag_POC_REMIN(k,column) = POC%remin

   work = P_CaCO3%sflux_in + P_CaCO3%hflux_in
   BGC_diagnostic_fields%diag_CaCO3_FLUX_IN(k,column) = work

   BGC_diagnostic_fields%diag_CaCO3_PROD(k,column) = P_CaCO3%prod

   BGC_diagnostic_fields%diag_CaCO3_REMIN(k,column) = P_CaCO3%remin

   work = P_SiO2%sflux_in + P_SiO2%hflux_in
   BGC_diagnostic_fields%diag_SiO2_FLUX_IN(k,column) = work

   BGC_diagnostic_fields%diag_SiO2_PROD(k,column) = P_SiO2%prod

   BGC_diagnostic_fields%diag_SiO2_REMIN(k,column) = P_SiO2%remin

   work = dust%sflux_in + dust%hflux_in
   BGC_diagnostic_fields%diag_dust_FLUX_IN(k,column) = work

   BGC_diagnostic_fields%diag_dust_REMIN(k,column) = dust%remin

   work = P_iron%sflux_in + P_iron%hflux_in
   BGC_diagnostic_fields%diag_P_iron_FLUX_IN(k,column) = work

   BGC_diagnostic_fields%diag_P_iron_PROD(k,column) = P_iron%prod

   BGC_diagnostic_fields%diag_P_iron_REMIN(k,column) = P_iron%remin

! ***********************************************************************
! - Accumulte losses of BGC tracers to sediments
! ***********************************************************************

   BGC_diagnostic_fields%diag_calcToSed(k,column) = P_CaCO3%sed_loss

   BGC_diagnostic_fields%diag_bsiToSed(k,column) = P_SiO2%sed_loss

   BGC_diagnostic_fields%diag_pocToSed(k,column) = POC%sed_loss

   work = SED_DENITRIF * cell_thickness
   BGC_diagnostic_fields%diag_SedDenitrif(k,column) = work

   work = OTHER_REMIN * cell_thickness
   BGC_diagnostic_fields%diag_OtherRemin(k,column) = work

   work = (POC%sed_loss * Q)
   BGC_diagnostic_fields%diag_ponToSed(k,column) = work

   work = (POC%sed_loss * Qp_zoo_pom)
   BGC_diagnostic_fields%diag_popToSed(k,column) = work

   BGC_diagnostic_fields%diag_dustToSed(k,column) = dust%sed_loss

   BGC_diagnostic_fields%diag_pfeToSed(k,column) = P_iron%sed_loss

!-----------------------------------------------------------------------
!EOC

 end subroutine compute_particulate_terms

!***********************************************************************
!BOP
! !IROUTINE: BGC_SurfaceFluxes
! !INTERFACE:

 subroutine BGC_SurfaceFluxes(BGC_indices, BGC_input, BGC_forcing,   &
                              BGC_flux_diagnostic_fields,   &
                              numColumnsMax, numColumns)

! !DESCRIPTION:
!  Compute surface fluxes for ecosystem state variables
!
! !REVISION HISTORY:
!  same as module

   implicit none

! !INPUT PARAMETERS:

  type(BGC_indices_type), intent(in )   :: BGC_indices
  type(BGC_input_type),   intent(in )   :: BGC_input
  type(BGC_forcing_type), intent(inout) :: BGC_forcing

  integer (BGC_i4) :: numColumnsMax, numColumns

! !OUTPUT PARAMETERS:

  type(BGC_flux_diagnostics_type), intent(inout) :: BGC_flux_diagnostic_fields

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (BGC_r8) :: &
      xkw,              &
      xkw_ice,          &
      SCHMIDT_O2,       &
      O2SAT_1atm,       &
      pistonVel_O2,     &
      O2SAT ,           &
      depth,            &
      SCHMIDT_CO2,      &
      pistonVel_CO2,    &
      phlo,             &
      phhi,             &
      ph_new,           &
      co2star,          &
      dco2star,         &
      pco2surf,         &
      dpco2

   real (BGC_r8), allocatable, dimension(:) :: &
      DIC_loc,        & ! local copy of model DIC
      DIC_ALT_CO2_loc,& ! local copy of model DIC_ALT_CO2
      ALK_loc,        & ! local copy of model ALK
      PO4_loc,        & ! local copy of model PO4
      NO3_loc,        & ! local copy of model NO3
      SiO3_loc,       & ! local copy of model SiO3
      O2_loc            ! local copy of model O2

   logical (BGC_log) :: zero_mask
   logical (BGC_log), parameter :: locmip_k1_k2_bug_fix = .true.

   integer (BGC_i4) :: &
      n,              & ! tracer index
      column            ! index for looping over columns

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  allocate local copies of tracers
!-----------------------------------------------------------------------

   allocate(DIC_loc(numColumns))
   allocate(DIC_ALT_CO2_loc(numColumns))
   allocate(ALK_loc(numColumns))
   allocate(PO4_loc(numColumns))
   allocate(NO3_loc(numColumns))
   allocate(SiO3_loc(numColumns))
   allocate(O2_loc(numColumns))

!-----------------------------------------------------------------------
!  zero out diagnostic arrays
!-----------------------------------------------------------------------

    BGC_flux_diagnostic_fields%pistonVel_O2 = c0
    BGC_flux_diagnostic_fields%pistonVel_CO2 = c0
    BGC_flux_diagnostic_fields%SCHMIDT_O2 = c0
    BGC_flux_diagnostic_fields%SCHMIDT_CO2 = c0
    BGC_flux_diagnostic_fields%O2SAT = c0
    BGC_flux_diagnostic_fields%xkw = c0
    BGC_flux_diagnostic_fields%co2star = c0
    BGC_flux_diagnostic_fields%dco2star = c0
    BGC_flux_diagnostic_fields%pco2surf = c0
    BGC_flux_diagnostic_fields%dpco2 = c0
    BGC_flux_diagnostic_fields%co2star_alt_co2 = c0
    BGC_flux_diagnostic_fields%dco2star_alt_co2 = c0
    BGC_flux_diagnostic_fields%pco2surf_alt_co2 = c0
    BGC_flux_diagnostic_fields%dpco2_alt_co2 = c0

!-----------------------------------------------------------------------
!  loop over columns
!-----------------------------------------------------------------------

   column_loop: do column = 1, numColumns

!-----------------------------------------------------------------------
!  create local copies of model tracers
!  treat negative values as zero
!-----------------------------------------------------------------------

!maltrud intel fails if i use c0 for max() instead of 0.0
      DIC_loc(column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(1,column,BGC_indices%dic_ind))
      DIC_ALT_CO2_loc(column) = max(0.0_BGC_r8, BGC_input%BGC_tracers(1,column,BGC_indices%dic_alt_co2_ind))
      ALK_loc(column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(1,column,BGC_indices%alk_ind))
      PO4_loc(column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(1,column,BGC_indices%po4_ind))
      NO3_loc(column)      = max(0.0_BGC_r8, BGC_input%BGC_tracers(1,column,BGC_indices%no3_ind))
      SiO3_loc(column)     = max(0.0_BGC_r8, BGC_input%BGC_tracers(1,column,BGC_indices%sio3_ind))
      O2_loc(column)       = max(0.0_BGC_r8, BGC_input%BGC_tracers(1,column,BGC_indices%o2_ind))

!-----------------------------------------------------------------------
!  convert total iron to bioavailable
!-----------------------------------------------------------------------

      BGC_forcing%depositionFlux(column, BGC_indices%fe_ind) =  &
          BGC_forcing%depositionFlux(column, BGC_indices%fe_ind) * parm_Fe_bioavail
      BGC_forcing%riverFlux(column, BGC_indices%fe_ind) =  &
          BGC_forcing%riverFlux(column, BGC_indices%fe_ind) * parm_Fe_bioavail
      BGC_forcing%gasFlux(column, BGC_indices%fe_ind) =  &
          BGC_forcing%gasFlux(column, BGC_indices%fe_ind) * parm_Fe_bioavail
      BGC_forcing%seaIceFlux(column, BGC_indices%fe_ind) =  &
          BGC_forcing%seaIceFlux(column, BGC_indices%fe_ind) * parm_Fe_bioavail

      if (BGC_forcing%iceFraction(column) < 0.0_BGC_r8) BGC_forcing%iceFraction(column) = 0.0_BGC_r8
      if (BGC_forcing%iceFraction(column) > 1.0_BGC_r8) BGC_forcing%iceFraction(column) = 1.0_BGC_r8

      xkw = xkw_coeff * BGC_forcing%windSpeedSquared10m(column)
      xkw_ice = (1.0_BGC_r8 - BGC_forcing%iceFraction(column)) * xkw

!-----------------------------------------------------------------------
!  compute O2 flux
!-----------------------------------------------------------------------

      if (BGC_forcing%lcalc_O2_gas_flux) then

         SCHMIDT_O2 = SCHMIDT_O2_singleValue(BGC_forcing%SST(column))
         O2SAT_1atm = O2SAT_singleValue(BGC_forcing%SST(column), BGC_forcing%SSS(column))

         pistonVel_O2 = xkw_ice * SQRT(660.0_BGC_r8 / SCHMIDT_O2)
         O2SAT = BGC_forcing%surfacePressure(column) * O2SAT_1atm
         BGC_forcing%gasFlux(column, BGC_indices%o2_ind) = pistonVel_O2 * (O2SAT - O2_loc(column))

         BGC_flux_diagnostic_fields%pistonVel_O2(column) = pistonVel_O2
         BGC_flux_diagnostic_fields%SCHMIDT_O2(column) = SCHMIDT_O2
         BGC_flux_diagnostic_fields%O2SAT(column) = O2SAT
         BGC_flux_diagnostic_fields%xkw(column) = xkw_ice
      endif

!-----------------------------------------------------------------------
!  compute CO2 flux
!-----------------------------------------------------------------------

      if (BGC_forcing%lcalc_CO2_gas_flux) then

         SCHMIDT_CO2 = SCHMIDT_CO2_singleValue(BGC_forcing%SST(column))

         pistonVel_CO2 = xkw_ice * SQRT(660.0_BGC_r8 / SCHMIDT_CO2)

         if (BGC_forcing%surface_pH(column) /= c0) then
            phlo = BGC_forcing%surface_pH(column) - del_ph
            phhi = BGC_forcing%surface_pH(column) + del_ph
         else
            phlo = phlo_surf_init
            phhi = phhi_surf_init
         end if

         depth = BGC_forcing%surfaceDepth(column)
         call co2calc_1point(depth, locmip_k1_k2_bug_fix, .true., &
                             BGC_forcing%SST(column), BGC_forcing%SSS(column), &
                             DIC_loc(column), ALK_loc(column), PO4_loc(column), SiO3_loc(column), &
                             phlo, phhi, ph_new, BGC_forcing%atmCO2(column), &
                             BGC_forcing%surfacePressure(column), co2star, &
                             dco2star, pco2surf, dpco2)

         BGC_forcing%surface_pH(column) = ph_new

         BGC_forcing%gasFlux(column, BGC_indices%dic_ind) = pistonVel_CO2 * dco2star

         BGC_flux_diagnostic_fields%co2star(column)  = co2star
         BGC_flux_diagnostic_fields%dco2star(column) = dco2star
         BGC_flux_diagnostic_fields%pco2surf(column) = pco2surf
         BGC_flux_diagnostic_fields%dpco2(column)    = dpco2
         BGC_flux_diagnostic_fields%pistonVel_CO2(column) = pistonVel_CO2
         BGC_flux_diagnostic_fields%SCHMIDT_CO2(column) = SCHMIDT_CO2

         if (BGC_forcing%surface_pH_alt_co2(column) /= c0) then
            phlo = BGC_forcing%surface_pH_alt_co2(column) - del_ph
            phhi = BGC_forcing%surface_pH_alt_co2(column) + del_ph
         else
            phlo = phlo_surf_init
            phhi = phhi_surf_init
         end if

         call co2calc_1point(depth, locmip_k1_k2_bug_fix, .true., &
                             BGC_forcing%SST(column), BGC_forcing%SSS(column), &
                             DIC_ALT_CO2_loc(column), ALK_loc(column), PO4_loc(column), SiO3_loc(column), &
                             phlo, phhi, ph_new, BGC_forcing%atmCO2_ALT_CO2(column), &
                             BGC_forcing%surfacePressure(column), co2star, &
                             dco2star, pco2surf, dpco2)

         BGC_forcing%surface_pH_alt_co2(column) = ph_new

         BGC_forcing%gasFlux(column, BGC_indices%dic_alt_co2_ind) = pistonVel_CO2 * dco2star

         BGC_flux_diagnostic_fields%co2star_alt_co2(column)  = co2star
         BGC_flux_diagnostic_fields%dco2star_alt_co2(column) = dco2star
         BGC_flux_diagnostic_fields%pco2surf_alt_co2(column) = pco2surf
         BGC_flux_diagnostic_fields%dpco2_alt_co2(column)    = dpco2

      endif

!-----------------------------------------------------------------------
!  Sum up components for each tracer to get net flux
!-----------------------------------------------------------------------

      do n = 1, BGC_tracer_cnt
         BGC_forcing%netFlux(column, n)           =   &
            BGC_forcing%depositionFlux(column, n) +   &
            BGC_forcing%gasFlux(column, n)        +   &
            BGC_forcing%riverFlux(column, n)      +   &
            BGC_forcing%seaIceFlux(column, n)
      enddo

!-----------------------------------------------------------------------
!  Apply NO & NH fluxes to alkalinity
!-----------------------------------------------------------------------

      BGC_forcing%netFlux(column, BGC_indices%alk_ind) = BGC_forcing%netFlux(column, BGC_indices%alk_ind) +  &
          BGC_forcing%netFlux(column, BGC_indices%nh4_ind) - BGC_forcing%netFlux(column, BGC_indices%no3_ind)

   enddo column_loop

   deallocate(DIC_loc)
   deallocate(DIC_ALT_CO2_loc)
   deallocate(ALK_loc)
   deallocate(PO4_loc)
   deallocate(NO3_loc)
   deallocate(SiO3_loc)
   deallocate(O2_loc)

!-----------------------------------------------------------------------
!EOC

 end subroutine BGC_SurfaceFluxes

!***********************************************************************
!*****************************************************************************
!BOP
! !IROUTINE: SCHMIDT_O2_singleValue
! !INTERFACE:

 function SCHMIDT_O2_singleValue(SST)

! !DESCRIPTION:
!  Compute Schmidt number of O2 in seawater as function of SST
!  where LAND_MASK is true. Give zero where LAND_MASK is false.
!
!  ref : Keeling et al, Global Biogeochem. Cycles, Vol. 12,
!        No. 1, pp. 141-163, March 1998
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (BGC_r8), intent(in) :: SST


! !OUTPUT PARAMETERS:

   real (BGC_r8) :: SCHMIDT_O2_singleValue

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (BGC_r8), parameter :: &
      a = 1638.0_BGC_r8, &
      b = 81.83_BGC_r8, &
      c = 1.483_BGC_r8, &
      d = 0.008004_BGC_r8

!-----------------------------------------------------------------------

   SCHMIDT_O2_singleValue = a + SST * (-b + SST * (c + SST * (-d)))

!-----------------------------------------------------------------------
!EOC

 end function SCHMIDT_O2_singleValue

!*****************************************************************************
!BOP
! !IROUTINE: O2SAT_singleValue
! !INTERFACE:

 function O2SAT_singleValue(SST, SSS)

! !DESCRIPTION:
!
!  Computes oxygen saturation concentration at 1 atm total pressure
!  in mmol/m^3 given the temperature (t, in deg C) and the salinity (s,
!  in permil) where LAND_MASK is true. Give zero where LAND_MASK is false.
!
!  FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
!  THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
!
!  *** NOTE: THE "A_3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
!  *** IT SHOULD NOT BE THERE.                                ***
!
!  O2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
!  0 permil <= S <= 42 permil
!  CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil,
!  O2SAT = 282.015 mmol/m^3
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (BGC_r8), intent(in) :: &
      SST, & ! sea surface temperature (C)
      SSS    ! sea surface salinity (psu)

! !OUTPUT PARAMETERS:

    real (BGC_r8) :: O2SAT_singleValue

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (BGC_r8) :: TS

!-----------------------------------------------------------------------
!  coefficients in expansion
!-----------------------------------------------------------------------

   real (BGC_r8), parameter :: &
      a_0 = 2.00907_BGC_r8, &
      a_1 = 3.22014_BGC_r8, &
      a_2 = 4.05010_BGC_r8, &
      a_3 = 4.94457_BGC_r8, &
      a_4 = -2.56847E-1_BGC_r8, &
      a_5 = 3.88767_BGC_r8, &
      b_0 = -6.24523E-3_BGC_r8, &
      b_1 = -7.37614E-3_BGC_r8, &
      b_2 = -1.03410E-2_BGC_r8, &
      b_3 = -8.17083E-3_BGC_r8, &
      c_0 = -4.88682E-7_BGC_r8

      TS = log( ((T0_Kelvin_BGC + 25.0_BGC_r8) - SST) / (T0_Kelvin_BGC + SST) )

      O2SAT_singleValue = exp(a_0+TS*(a_1+TS*(a_2+TS*(a_3+TS*(a_4+TS*a_5)))) + &
         SSS*( (b_0+TS*(b_1+TS*(b_2+TS*b_3))) + SSS*c_0 ))

!-----------------------------------------------------------------------
!  Convert from ml/l to mmol/m^3
!-----------------------------------------------------------------------

   O2SAT_singleValue = O2SAT_singleValue / 0.0223916_BGC_r8

!-----------------------------------------------------------------------
!EOC

 end function O2SAT_singleValue

!*****************************************************************************
!*****************************************************************************
!BOP
! !IROUTINE: SCHMIDT_CO2_singleValue
! !INTERFACE:

 function SCHMIDT_CO2_singleValue(SST)

! !DESCRIPTION:
!  Compute Schmidt number of CO2 in seawater as function of SST
!  where LAND_MASK is true. Give zero where LAND_MASK is false.
!
!  ref : Wanninkhof, J. Geophys. Res, Vol. 97, No. C5,
!  pp. 7373-7382, May 15, 1992
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (BGC_r8), intent(in) :: SST

! !OUTPUT PARAMETERS:

   real (BGC_r8) :: SCHMIDT_CO2_singleValue

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (BGC_r8), parameter :: &
      a = 2073.1_BGC_r8, &
      b = 125.62_BGC_r8, &
      c = 3.6276_BGC_r8, &
      d = 0.043219_BGC_r8

   SCHMIDT_CO2_singleValue = a + SST * (-b + SST * (c + SST * (-d)))

!-----------------------------------------------------------------------
!EOC

 end function SCHMIDT_CO2_singleValue

!*****************************************************************************

 end module BGC_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
