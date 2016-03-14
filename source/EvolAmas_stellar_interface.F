#if (_EVOL_STELL_<10) /* Simple stellar evol (MS-->remnant) or no stellar evolution */

c----------------------------------------------------------------------
      double precision function stellar_mass(id)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         stellar_mass = M_TN/Coef_ConvM
         return
      end if
#endif
      if (iTet_SE(id).LE.0) then
         stellar_mass = 0.0d0
      else
         stellar_mass = Met_SE(id)
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_radius(id)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c
c function :
c ^^^^^^^^^^
      double precision RelMR
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         stellar_radius = RelMR(M_TN/Coef_ConvM,-1.0d0,iType_CtrObj)
         !write(91,*) M_TN/Coef_ConvM,stellar_radius !###################
         return
      end if
#endif
      if (iTet_SE(id).LE.0) then
         stellar_radius = 0.0d0
      else
#if (_TYPE_FOR_MR_RELATION_>0)
         stellar_radius = RelMR(Met_SE(id),-1.0d0,_TYPE_FOR_MR_RELATION_)
#else
         stellar_radius = RelMR(Met_SE(id),-1.0d0,iTet_SE(id))
#endif
      end if
c
      end
c
c----------------------------------------------------------------------
      integer*1 function stellar_type(id)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         stellar_type = iType_CtrObj
         return
      end if
#endif
      stellar_type = iTet_SE(id)
c
      end
c
c----------------------------------------------------------------------
      integer*1 function stellar_generaltype(id)
c----------------------------------------------------------------------
c     returns 1 for MS
c             2 for post-MS (non compact)
c             3 for white dwarf
c             4 for neutron star
c             5 for stellar bh
c
      include 'EvolAmas_Etoiles_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
      stellar_generaltype = iTet_SE(id)
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_evoltimescale(id) ! in years
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c
c function :
c ^^^^^^^^^^
      double precision T_MS
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         if (iType_CtrObj.LE.0) then
            stellar_evoltimescale = 1.0d30
         else
            if (iType_CtrObj.EQ.TypeEtoile_MS) then
               stellar_evoltimescale = T_MS(M_TN/Coef_ConvM)
            else
               stellar_evoltimescale = 1.0d30
            end if
         end if
         return
      end if
#endif
      if (iTet_SE(id).LE.0) then
         stellar_evoltimescale = 1.0d30
      else
         if (iTet_SE(id).EQ.TypeEtoile_MS) then
            stellar_evoltimescale = T_MS(Met_SE(id))
         else
            stellar_evoltimescale = 1.0d30
         end if
      end if
c
      end
c
c----------------------------------------------------------------------
      logical function stellar_is_compact(id)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         stellar_is_compact = (
     $        iType_CtrObj.EQ.TypeEtoile_WD .OR.
     $        iType_CtrObj.EQ.TypeEtoile_NS .OR.
     $        iType_CtrObj.EQ.TypeEtoile_BH )
         return
      end if
#endif
      stellar_is_compact = (
     $     iTet_SE(id).EQ.TypeEtoile_WD .OR.
     $     iTet_SE(id).EQ.TypeEtoile_NS .OR.
     $     iTet_SE(id).EQ.TypeEtoile_BH )
c
      end
c
c----------------------------------------------------------------------
      logical function stellar_type_is_compact(iType)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer*1 iType
c......................................................................
c
      stellar_type_is_compact = (
     $     iType.EQ.TypeEtoile_WD .OR.
     $     iType.EQ.TypeEtoile_NS .OR.
     $     iType.EQ.TypeEtoile_BH )
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_dateofbirth(id)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         stellar_dateofbirth = BirthDate_CtrObj
         return
      end if
#endif
      stellar_dateofbirth = DNet_SE(id)
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_age(id)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_Temps_common.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         if (iType_CtrObj.LE.0) then
            stellar_age = -1.0d30
         else
            stellar_age = max(0.0d0,Tps_Amas/UneAnnee-BirthDate_CtrObj)
         end if
         return
      end if
#endif
      if (iTet_SE(id).LE.0) then
         stellar_age = -1.0d30
      else
         stellar_age = max(0.0d0,Tps_SE(id)/UneAnnee-DNet_SE(id)) !# change Tps_Amas to Tps_SE; Marc, Sep 9 2004
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_centralheliumfrac(id) ! Mass of central He (not including primordial) / total stellar mass
                                                              ! Required for collisional rejuvenation
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision T_MS,M_TMS,he_frac_TMS
c
c function :
c ^^^^^^^^^^
      double precision stellar_age
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         if (iType_CtrObj.LE.0) then
            stellar_centralheliumfrac = -1.0d30
            return
         end if
         call GetMSdata(M_TN/Coef_ConvM,T_MS,M_TMS,he_frac_TMS)
         stellar_centralheliumfrac = min(1.0d0,
     $        stellar_age(id)/T_MS) * he_frac_TMS

         return
      end if
#endif
      if (iTet_SE(id).LE.0) then
         stellar_centralheliumfrac = -1.0d30
         return
      end if
      call GetMSdata(Met_SE(id),T_MS,M_TMS,he_frac_TMS)
      stellar_centralheliumfrac = min(1.0d0,stellar_age(id)/T_MS) *
     $     he_frac_TMS
c
      end
c
#if (_EVOL_STELL_==0) 
c----------------------------------------------------------------------
      subroutine stellar_SetMass_ms(
     $     id,mass)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer id
      double precision mass ! Mass in Msun
c......................................................................
c
      Met_SE(id) = mass
c
      end
c
#endif
c
c----------------------------------------------------------------------
      subroutine stellar_SetByAge_ms(
     $     id,present_mass,metallicity,age_in_years)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_Temps_common.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
      double precision present_mass, metallicity, age_in_years
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision age_eff
c
c function :
c ^^^^^^^^^^
      double precision T_MS
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         if (present_mass.LE.0.0d0) then
            iType_CtrObj = -iType_CtrObj
            M_TN  = 0.0d0
         else
            iType_CtrObj = TypeEtoile_MS
            M_TN  = Coef_ConvM*present_mass ! convert from Msun into code units
            age_eff = min(1.01d0*T_MS(present_mass),age_in_years) ! coef 1.01 to ensure that star is due for evolution off the MS 
            BirthDate_CtrObj = Tps_Amas/UneAnnee - age_eff ! effective birthdate
         end if
         return
      end if
#endif

      if (present_mass.LE.0.0d0) then
         iTet_SE(id) = -abs(iTet_SE(id))
         Met_SE(id)  = 0.0d0
      else
         iTet_SE(id) = TypeEtoile_MS
         Met_SE(id)  = present_mass
         age_eff = min(1.01d0*T_MS(present_mass),age_in_years) ! coef 1.01 to ensure that star is due for evolution off the MS 
         DNet_SE(id) = Tps_Amas/UneAnnee - age_eff ! effective birthdate
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine stellar_SetByHeFrac_ms(
     $     id,present_mass,metallicity,he_frac)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_Temps_common.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
      double precision present_mass, metallicity, he_frac
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision T_MS,M_TMS,he_frac_TMS,age_in_years, age_eff
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         if (present_mass.LE.0.0d0) then
            iType_CtrObj = -iType_CtrObj
            M_TN  = 0.0d0
         else
            iType_CtrObj = TypeEtoile_MS
            M_TN  = Coef_ConvM*present_mass ! convert from Msun into code units
            call GetMSdata(present_mass,T_MS,M_TMS,he_frac_TMS)
            !!write(98,*) present_mass,T_MS,M_TMS,he_frac_TMS !#############
            age_in_years = he_frac/he_frac_TMS*T_MS
            age_eff = min(1.01d0*T_MS,age_in_years) ! coef 1.01 to ensure that star is due for evolution off the MS 
            BirthDate_CtrObj = Tps_Amas/UneAnnee - age_eff ! effective birthdate
         end if
         return
      end if
#endif
      if (present_mass.LE.0.0d0) then
         iTet_SE(id) = -abs(iTet_SE(id))
         Met_SE(id)  = 0.0d0
      else
         iTet_SE(id) = TypeEtoile_MS
         Met_SE(id)  = present_mass
         call GetMSdata(present_mass,T_MS,M_TMS,he_frac_TMS)
         age_in_years = he_frac/he_frac_TMS*T_MS
         age_eff = min(1.01d0*T_MS,age_in_years) ! coef 1.01 to ensure that star is due for evolution off the MS 
         DNet_SE(id) = Tps_Amas/UneAnnee - age_eff ! effective birthdate
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine stellar_SetMassType_compact(
     $     id,new_mass,new_type)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_Temps_common.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
      double precision new_mass
      integer*1 new_type
c
c function :
c ^^^^^^^^^^
      logical stellar_type_is_compact
c......................................................................
c
      if (.NOT.stellar_type_is_compact(new_type)) then
         write(0,*) '!!! Type ',new_type,
     $        ' doesn''t correspond to compact star !!!'
         call Meurt('stellar_SetMassType_compact')
      end if

#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         M_TN  = Coef_ConvM*new_mass ! convert from Msun into code units
         iType_CtrObj = new_type
         return
      end if
#endif
      iTet_SE(id) = new_type
      Met_SE(id)  = new_mass
c
      end
c
c----------------------------------------------------------------------
      subroutine stellar_Destroy(id)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         iType_CtrObj=0
         M_TN = 0.0d0
         return
      end if
#endif
      iTet_SE(id) = 0
      Met_SE(id)  = 0.0d0
c
      end
c
#endif /* (_EVOL_STELL_<10) */ /* Simple stellar evol (MS-->remnant) or no stellar evolution */

#if (_EVOL_STELL_==1) /* Simple stellar evol (MS-->remnant) */
c
c----------------------------------------------------------------------
      subroutine evolve_star(id,age_in_years)
c----------------------------------------------------------------------
      include 'EvolAmas_Etoiles_inc.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
      double precision age_in_years
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision M_evol
      integer*1 iT_evol
c......................................................................
c
#if (_PRESENCE_TN_>0) /* A central object is present */
      if (id.EQ.id_CtrObj) then
         if (iType_CtrObj.EQ.TypeEtoile_MS) then
            call Evoluer_Etoile
     $           (M_TN/Coef_ConvM,age_in_years,M_evol,iT_evol)
            M_TN=Coef_ConvM*M_evol
            iType_CtrObj=iT_evol
         end if
         return
      end if
#endif
      if (iTet_SE(id).EQ.TypeEtoile_MS) then
         call Evoluer_Etoile(Met_SE(id),age_in_years,M_evol,iT_evol)
         Met_SE(id)=M_evol
         iTet_SE(id)=iT_evol
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function natal_kick_star()
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_EvolStell_common.f'
c......................................................................
c
      natal_kick_star = NatalKickLastEvol
      NatalKickLastEvol = 0.0d0
c
      end
c
#endif /* (_EVOL_STELL_==1) */ /* Simple stellar evol (MS-->remnant) */


#if (_EVOL_STELL_==10) /* SeBa stellar evolution */
c----------------------------------------------------------------------
      subroutine create_star_from_zams(id,
     $     M_ZAMS,metallicity,age_in_years, iflag)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_TN_common.f'
c
c arguments :
c ^^^^^^^^^^^
      integer id                                       ! user-provided
      double precision M_ZAMS,metallicity,age_in_years ! user-provided
      integer iflag ! to signal errors or problems
c
c local variable :
c ^^^^^^^^^^^^^^^^
      double precision tMyr
      integer nstar
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      if (metallicity.NE.solar_metal) then
         write(0,*) '!!! Only solar metallicity is allowed (',
     $        solar_metal,') !!!'
         call Meurt('create_star_from_zams')
      end if

      iflag = 0                 ! not used yet
      nstar = 1                 ! Initialize a single star only
      tMyr = 1.0d-6*age_in_years
      call initstar(nstar, id, M_ZAMS)
      call evstar(id, tMyr)
c
      end
c
c----------------------------------------------------------------------
      subroutine create_stars_from_zams(Nstar,
     $     zams_masses,metallicities,ages_in_years, iflag)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_Etoiles_param.f'
      include 'Param_Amas.f'
c
c arguments :
c ^^^^^^^^^^^
      integer Nstar
      double precision zams_masses(Nstar),metallicities(Nstar),
     $     ages_in_years(Nstar)
      integer iflag ! to signal errors or problems
c
c local variable :
c ^^^^^^^^^^^^^^^^
      double precision tMyr
      integer id(iDimSE)
      integer i
c......................................................................
c
      do i=1,Nstar
         if (metallicities(i).NE.solar_metal) then
            write(0,*) '!!! Only solar metallicity is allowed (',
     $           solar_metal,') !!!'
            call Meurt('create_stars_from_zams')
         end if
      end do

      iflag = 0                 ! not used yet
      do i=1,Nstar
         id(i)=i
      end do
      call initstar(Nstar, id, zams_masses)
      do i=1,Nstar
         tMyr=1.0d-6*ages_in_years(i)
         call evstar(i, tMyr)
      end do
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_mass(id)
c----------------------------------------------------------------------
c     returns the present-day mass in M_sun
c
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      call getmass(id, stellar_mass)
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_radius(id)
c----------------------------------------------------------------------
c     returns the present-day radius in R_sun
c
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('stellar_radius')
      end if
#endif
      call getradius(id, stellar_radius)
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_zamsmass(id) ! May not be required
c----------------------------------------------------------------------
c     returns the ZAMS mass in M_sun
c
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      call Meurt('stellar_zamsmass NOT IMPLEMENTED YET')
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_metallicity(id) ! Not implemented yet by SeBa
c----------------------------------------------------------------------
c     returns the ZAMS mass in M_sun
c
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      call Meurt('stellar_metallicity NOT IMPLEMENTED YET')
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_age(id) ! REQUIRED!
c----------------------------------------------------------------------
c     returns total (effective) age (since "birth" on ZAMS) in years
c
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      call Meurt('stellar_age NOT IMPLEMENTED YET')
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_phaseage(id) ! May not be required
c----------------------------------------------------------------------
c     returns (effective) time spent in current phase in years
c
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      call Meurt('stellar_phaseage NOT IMPLEMENTED YET')
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_centralheliumfrac(id) ! REQUIRED for collisional "rejuvenation"
c----------------------------------------------------------------------
c     Mass of central He (not including primordial) / total stellar mass
c
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      call Meurt('stellar_centralhelium NOT IMPLEMENTED YET')
c
      end
c
c----------------------------------------------------------------------
      integer*1 function stellar_type(id)
c----------------------------------------------------------------------
c     See the list of types in ~/starlab/4.3/inc/star/star_support.h
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c
c local variable :
c ^^^^^^^^^^^^^^^^
      integer stype     
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      callgetstype(id, stype)
      stellar_type = stype

      end
c
c----------------------------------------------------------------------
      integer*1 function stellar_generaltype(id)
c----------------------------------------------------------------------
c     See the list of types in ~/starlab/4.3/inc/star/star_support.h
c
c     returns 1 for MS
c             2 for post-MS (non compact)
c             3 for white dwarf
c             4 for neutron star
c             5 for stellar bh
c
      implicit none
      include 'EvolAmas_Etoiles_param.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c
c local variable :
c ^^^^^^^^^^^^^^^^
      integer stype     
c......................................................................

      callgetstype(id, stype)
      if (stype.LT.0) then
         stellar_generaltype = 0
      else if (stype.LE.3) then
         stellar_generaltype = TypeEtoile_MS
      else if (stype.LE.8) then
         stellar_generaltype = TypeEtoile_RG
      else if (stype.LE.10) then
         stellar_generaltype = TypeEtoile_MS
      else if (stype.LE.11) then
         stellar_generaltype = TypeEtoile_RG
      else if (stype.LE.14) then
         stellar_generaltype = TypeEtoile_WD
      else if (stype.LE.18) then
         stellar_generaltype = TypeEtoile_NS
      else if (stype.LE.19) then
         stellar_generaltype = TypeEtoile_BH
      else
         stellar_generaltype = 0
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_evoltimescale(id) ! UPDATE REQUIRED
c----------------------------------------------------------------------
c     Total duration of the current phase, in years
c
      implicit none
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
      double precision T_MS
      integer*1 gentype, stellar_generaltype
      double precision stellar_mass
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
c
c VERY ROUGH!! UPDATE by use of SeBa routines required
c
      gentype = stellar_generaltype(id)
      if (gentype.EQ.TypeEtoile_MS) then
         stellar_evoltimescale = T_MS(stellar_mass(id))
      else if (gentype.EQ.TypeEtoile_RG) then
         stellar_evoltimescale = 0.05d0*T_MS(stellar_mass(id))
      else if (gentype.EQ.TypeEtoile_WD) then
         stellar_evoltimescale = 1.0d9
      else
         stellar_evoltimescale = 1.0d30
      end if
c
      end
c
c----------------------------------------------------------------------
      logical function stellar_is_compact(id)
c----------------------------------------------------------------------
c     returns TRUE for white dwarfs, neutron stars and stellar BHs,
c     used to set collisional time
c
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c
c local variable :
c ^^^^^^^^^^^^^^^^
      integer stype 
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      call getstype(id, stype)
      if (stype.GT.0 .AND. stype.LE.11) then
         stellar_is_compact = .FALSE.
      else
         stellar_is_compact = .TRUE.
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function stellar_dateofbirth(id) ! REQUIRED
c----------------------------------------------------------------------
c     returns effective date of birth, in years, i.e. the cluster time
c     at which a star should have been "born" on the ZAMS to have the
c     state of star id at the present time (good luck to implement this
c     one...)
c
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      call Meurt('stellar_dateofbirth NOT IMPLEMENTED YET')
      !stellar_dateofbirth = 0.0d0
c
      end
c
c----------------------------------------------------------------------
      subroutine stellar_SetByAge_ms(
     $     id,present_mass,metallicity,age_in_years) ! REQUIRED
c----------------------------------------------------------------------
c     Set the star of given id to MS evolutionary state defined by
c     its present mass and age in years (since birth)
c
c     Only solar metallicity is possible for the time being
c
      implicit none
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
      double precision present_mass, metallicity, age_in_years
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('stellar_SetByAge_ms')
      end if
#endif
      if (metallicity.NE.solar_metal) then
         write(0,*) '!!! Only solar metallicity is allowed (',
     $        solar_metal,') !!!'
         call Meurt('create_star_from_zams')
      end if

      call Meurt('stellar_SetByAge_ms NOT IMPLEMENTED YET')
c
      end
c
c----------------------------------------------------------------------
      subroutine stellar_SetByHeFrac_ms(
     $     id,present_mass,metallicity,he_frac) ! REQUIRED (for mergers)
c----------------------------------------------------------------------
c     Set the star of given id to MS evolutionary state defined my
c     its present mass and core helium content (in fraction of total mass)
c     This "Helium fraction" is defined to be 0 on the ZAMS (i.e. it
c     doesn't include "primordial" He)
c
c     Only solar metallicity is possible for the time being
c
      implicit none
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
      double precision present_mass, metallicity, he_frac
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('stellar_SetByHeFrac_ms')
      end if
#endif
      if (metallicity.NE.solar_metal) then
         write(0,*) '!!! Only solar metallicity is allowed (',
     $        solar_metal,') !!!'
         call Meurt('create_star_from_zams')
      end if

      call Meurt('stellar_SetByHeFrac_ms NOT IMPLEMENTED YET')
c
      end
c
c----------------------------------------------------------------------
      subroutine stellar_Destroy(id) ! REQUIRED
c----------------------------------------------------------------------
c     Destroy star with given id (whatever that means...)
c
      implicit none
      include 'EvolAmas_TN_common.f'
c
c argument :
c ^^^^^^^^^^
      integer id
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      call Meurt('stellar_Destroy NOT IMPLEMENTED YET')
c
      end
c
c----------------------------------------------------------------------
      subroutine evolve_star(id,age_in_years)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_TN_common.f'
c
c constants :
c ^^^^^^^^^^^
      double precision resol_phase
      parameter (resol_phase=0.01d0) ! should be enough...
c
c argument :
c ^^^^^^^^^^
      integer id
      double precision age_in_years
c
c local variable :
c ^^^^^^^^^^^^^^^^
      double precision tMyr, tMyr_udt
c
c functions :
c ^^^^^^^^^^^
      double precision stellar_evoltimescale, stellar_age
c......................................................................
c
#if (_PRESENCE_TN_==2) /* A central VMS is present */
      if (id.EQ.id_CtrObj) then
         write(0,*)
     $        '!!! Cannot deal with central object !!!'
         call Meurt('XXX')
      end if
#endif
      if ((age_in_years-stellar_age(id)) .GT.
     $     resol_phase*stellar_evoltimescale(id)) then
         call evstar(id, tMyr)
      else
         tMyr = 1.0d-6*age_in_years
         call update_time(id, tMyr_udt)
         if (tMyr.GE.tMyr_udt) then
            call evstar(id, tMyr)
         endif
      endif
c
      end
c
c----------------------------------------------------------------------
      double precision function natal_kick_star()
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_EvolStell_common.f'
c......................................................................
c
      natal_kick_star = NatalKickLastEvol
      NatalKickLastEvol = 0.0d0
c
      end
c
#endif /* (_EVOL_STELL_==10) */ /* SeBa stellar evolution */
