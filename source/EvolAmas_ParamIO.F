c======================================================================
c     Routines d'entrees/sorties des parametres de simulation
c     pour EvolAmas.F
c     code genere automatiquement par /Users/freitag/AMAS/util/EvolAmas_Param_rdb2IO.sh
c======================================================================
c
c----------------------------------------------------------------------
      subroutine Lire_Fich_Param_EvolAmas(NumFich)
c----------------------------------------------------------------------
c     Lit les donnees necessaire a l'execution de programme 
c     d'evolution d'amas. Le fichier doit etre sous le format :
c
c      1ere colonne
c     / 
c     NomAmasIn        << Nom de la variable
c       Amas_in.dat    << valeur 
c
c     NbPas             
c       100
c     ....
c
      include 'EvolAmas_param.f'
      include 'EvolAmas_Param_common.f' ! blocs common contenant les
                                        ! parametres
c
      integer NumFich
      character*80 Ligne
      character*1 car
c......................................................................
c
c---- Valeurs par defaut des parametres de simulation
c
      Type_CtrObj_def      = 'TN'
      NomSimul             = '_EvolAmas'
      Coef_Trlx_GW_Capt    = 1.0d0
      CollAccrFrac_MS_on_BH = 0.5d0
      CollAccrFrac_MS_on_NS = 0.0d0
      CollAccrFrac_MS_on_WD = 0.0d0
      CollDminDestr_MS_vs_RMN = 0.5d0
      CollMassLossFracDef  = 0.0d0
      FacMaxPG             = 1.5d0
      FacMinPG             = 0.5d0
      FacNbPasDetEvap      = 1.0d0
      FacNbPasDetTrelExtremes = 0.5d0
      FacNbPasInfo         = 1.0d0
      FacNbPasRecArbre     = 2.0d0
      FacNbPasRecCB3c      = 0.5d0
      FacNbPasRecRelax     = 10.0d0
      FacNbPasRecTpsAmas   = 1.0d0
      FacNbPasRecalcTH     = 0.5d0
      FacNbPasSauv         = 5.0d0
      FacNbPasSauvPart     = 1.0d0
      FacNbPasTestDemSauv  = 0.1d0
      FacNbPasTestDuplic   = 0.5d0
      Fac_SldAvg_Tscale    = 5.0d0
      Fac_dt_accr_M        = 0.01d0
      Fac_dt_accr_Npart    = 10.0d0
      FactR_CC             = 0.5d0
      FactTrelMax          = 1.0d3
      FracPas_AjustPassif  = 0.0d0
      Frac_Taccr           = 0.01d0
      Frac_Tadiab          = 0.01d0
      Frac_Tcoll           = 0.01d0
      Frac_Tevap           = 0.01d0
      Frac_Tevst           = 0.025d0
      Frac_Tkick           = 0.01d0
      Frac_Trelax          = 0.01d0
      Frac_Tstaccr         = 0.01d0
      Frac_dt_max          = 1.0d0
      Gamma_relax          = 0.14d0
      M_TN_ini_def         = 0.0d0
      MasseEtoileDef       = 1.0d0
      Mmax_AccrByStar      = 80.0
      Mu_mol_elec          = 1.13d0
      R_Mar_ini_def        = 1.0d30
      R_tid_disr_def_NB    = 1.0d-7
      RateAccrPerStar_MsunMyr = 1.0
      StellMet_def         = 0.02d0
      T_fin                = 1.0d5
      Tadiab_ini           = 0.001d0
      TailleAmas_en_pc     = 3.0d0
      Tmax_GW_Capt_Gyr     = 1.0d30
      Trlx_GW_def_Gyr      = 1.0d0
      effic_conv_lum       = 0.0d0
      exp_AccrByStar       = 0.0
      fact_b0_kick         = 2.0d0
      frac_accr_Coll       = 1.0d0
      frac_accr_Dechir     = 1.0d0
      frac_accr_EvolSt     = 1.0d0
      rNbEtoiles           = 1.0d6
      tAccrByStar_1Msun_Myr = 0.25
      DeltaRang_Paire      = 1
      NbSECoucheCB3c       = 1000
      NbSECouchePot        = 25
      NbSauv_Conserve      = 5
      Rang_CC              = 10
      TypeEtoileDef        = 1
      iRand_Seed           = 133111
      iVerbeux             = 1
      i_fin                = 2000000000
      i_ini                = 0
c
c---- Lecture du fichier de parametres
c
 1    continue
      read(NumFich,'(A)',END=2) Ligne
      call EnleverBlancs(Ligne)
c
c---- on saute les lignes blanches et les commentaires
c
      if (LEN(Ligne).EQ.0 .OR. Ligne(1:1).EQ.'#') then
	 goto 1
c
c---- on lit les parametres presents
c
      else if (Ligne.EQ.'')                   then
         read(NumFich,*) 

      else if (Ligne.EQ.'Type_CtrObj_def')    then
         read(NumFich,'(A)') Type_CtrObj_def
         call EnleverBlancs(Type_CtrObj_def)

      else if (Ligne.EQ.'NomSimul')           then
         read(NumFich,'(A)') NomSimul
         call EnleverBlancs(NomSimul)

      else if (Ligne.EQ.'Coef_Trlx_GW_Capt')  then
         read(NumFich,*) Coef_Trlx_GW_Capt

      else if (Ligne.EQ.'CollAccrFrac_MS_on_BH') then
         read(NumFich,*) CollAccrFrac_MS_on_BH

      else if (Ligne.EQ.'CollAccrFrac_MS_on_NS') then
         read(NumFich,*) CollAccrFrac_MS_on_NS

      else if (Ligne.EQ.'CollAccrFrac_MS_on_WD') then
         read(NumFich,*) CollAccrFrac_MS_on_WD

      else if (Ligne.EQ.'CollDminDestr_MS_vs_RMN') then
         read(NumFich,*) CollDminDestr_MS_vs_RMN

      else if (Ligne.EQ.'CollMassLossFracDef') then
         read(NumFich,*) CollMassLossFracDef

      else if (Ligne.EQ.'FacMaxPG')           then
         read(NumFich,*) FacMaxPG

      else if (Ligne.EQ.'FacMinPG')           then
         read(NumFich,*) FacMinPG

      else if (Ligne.EQ.'FacNbPasDetEvap')    then
         read(NumFich,*) FacNbPasDetEvap

      else if (Ligne.EQ.'FacNbPasDetTrelExtremes') then
         read(NumFich,*) FacNbPasDetTrelExtremes

      else if (Ligne.EQ.'FacNbPasInfo')       then
         read(NumFich,*) FacNbPasInfo

      else if (Ligne.EQ.'FacNbPasRecArbre')   then
         read(NumFich,*) FacNbPasRecArbre

      else if (Ligne.EQ.'FacNbPasRecCB3c')    then
         read(NumFich,*) FacNbPasRecCB3c

      else if (Ligne.EQ.'FacNbPasRecRelax')   then
         read(NumFich,*) FacNbPasRecRelax

      else if (Ligne.EQ.'FacNbPasRecTpsAmas') then
         read(NumFich,*) FacNbPasRecTpsAmas

      else if (Ligne.EQ.'FacNbPasRecalcTH')   then
         read(NumFich,*) FacNbPasRecalcTH

      else if (Ligne.EQ.'FacNbPasSauv')       then
         read(NumFich,*) FacNbPasSauv

      else if (Ligne.EQ.'FacNbPasSauvPart')   then
         read(NumFich,*) FacNbPasSauvPart

      else if (Ligne.EQ.'FacNbPasTestDemSauv') then
         read(NumFich,*) FacNbPasTestDemSauv

      else if (Ligne.EQ.'FacNbPasTestDuplic') then
         read(NumFich,*) FacNbPasTestDuplic

      else if (Ligne.EQ.'Fac_SldAvg_Tscale')  then
         read(NumFich,*) Fac_SldAvg_Tscale

      else if (Ligne.EQ.'Fac_dt_accr_M')      then
         read(NumFich,*) Fac_dt_accr_M

      else if (Ligne.EQ.'Fac_dt_accr_Npart')  then
         read(NumFich,*) Fac_dt_accr_Npart

      else if (Ligne.EQ.'FactR_CC')           then
         read(NumFich,*) FactR_CC

      else if (Ligne.EQ.'FactTrelMax')        then
         read(NumFich,*) FactTrelMax

      else if (Ligne.EQ.'FracPas_AjustPassif') then
         read(NumFich,*) FracPas_AjustPassif

      else if (Ligne.EQ.'Frac_Taccr')         then
         read(NumFich,*) Frac_Taccr

      else if (Ligne.EQ.'Frac_Tadiab')        then
         read(NumFich,*) Frac_Tadiab

      else if (Ligne.EQ.'Frac_Tcoll')         then
         read(NumFich,*) Frac_Tcoll

      else if (Ligne.EQ.'Frac_Tevap')         then
         read(NumFich,*) Frac_Tevap

      else if (Ligne.EQ.'Frac_Tevst')         then
         read(NumFich,*) Frac_Tevst

      else if (Ligne.EQ.'Frac_Tkick')         then
         read(NumFich,*) Frac_Tkick

      else if (Ligne.EQ.'Frac_Trelax')        then
         read(NumFich,*) Frac_Trelax

      else if (Ligne.EQ.'Frac_Tstaccr')       then
         read(NumFich,*) Frac_Tstaccr

      else if (Ligne.EQ.'Frac_dt_max')        then
         read(NumFich,*) Frac_dt_max

      else if (Ligne.EQ.'Gamma_relax')        then
         read(NumFich,*) Gamma_relax

      else if (Ligne.EQ.'M_TN_ini_def')       then
         read(NumFich,*) M_TN_ini_def

      else if (Ligne.EQ.'MasseEtoileDef')     then
         read(NumFich,*) MasseEtoileDef

      else if (Ligne.EQ.'Mmax_AccrByStar')    then
         read(NumFich,*) Mmax_AccrByStar

      else if (Ligne.EQ.'Mu_mol_elec')        then
         read(NumFich,*) Mu_mol_elec

      else if (Ligne.EQ.'R_Mar_ini_def')      then
         read(NumFich,*) R_Mar_ini_def

      else if (Ligne.EQ.'R_tid_disr_def_NB')  then
         read(NumFich,*) R_tid_disr_def_NB

      else if (Ligne.EQ.'RateAccrPerStar_MsunMyr') then
         read(NumFich,*) RateAccrPerStar_MsunMyr

      else if (Ligne.EQ.'StellMet_def')       then
         read(NumFich,*) StellMet_def

      else if (Ligne.EQ.'T_fin')              then
         read(NumFich,*) T_fin

      else if (Ligne.EQ.'Tadiab_ini')         then
         read(NumFich,*) Tadiab_ini

      else if (Ligne.EQ.'TailleAmas_en_pc')   then
         read(NumFich,*) TailleAmas_en_pc

      else if (Ligne.EQ.'Tmax_GW_Capt_Gyr')   then
         read(NumFich,*) Tmax_GW_Capt_Gyr

      else if (Ligne.EQ.'Trlx_GW_def_Gyr')    then
         read(NumFich,*) Trlx_GW_def_Gyr

      else if (Ligne.EQ.'effic_conv_lum')     then
         read(NumFich,*) effic_conv_lum

      else if (Ligne.EQ.'exp_AccrByStar')     then
         read(NumFich,*) exp_AccrByStar

      else if (Ligne.EQ.'fact_b0_kick')       then
         read(NumFich,*) fact_b0_kick

      else if (Ligne.EQ.'frac_accr_Coll')     then
         read(NumFich,*) frac_accr_Coll

      else if (Ligne.EQ.'frac_accr_Dechir')   then
         read(NumFich,*) frac_accr_Dechir

      else if (Ligne.EQ.'frac_accr_EvolSt')   then
         read(NumFich,*) frac_accr_EvolSt

      else if (Ligne.EQ.'rNbEtoiles')         then
         read(NumFich,*) rNbEtoiles

      else if (Ligne.EQ.'tAccrByStar_1Msun_Myr') then
         read(NumFich,*) tAccrByStar_1Msun_Myr

      else if (Ligne.EQ.'DeltaRang_Paire')    then
#if (_COMPIL_==3)
         read(NumFich,'(I)') DeltaRang_Paire
#else
         read(NumFich,*) DeltaRang_Paire
#endif

      else if (Ligne.EQ.'NbSECoucheCB3c')     then
#if (_COMPIL_==3)
         read(NumFich,'(I)') NbSECoucheCB3c
#else
         read(NumFich,*) NbSECoucheCB3c
#endif

      else if (Ligne.EQ.'NbSECouchePot')      then
#if (_COMPIL_==3)
         read(NumFich,'(I)') NbSECouchePot
#else
         read(NumFich,*) NbSECouchePot
#endif

      else if (Ligne.EQ.'NbSauv_Conserve')    then
#if (_COMPIL_==3)
         read(NumFich,'(I)') NbSauv_Conserve
#else
         read(NumFich,*) NbSauv_Conserve
#endif

      else if (Ligne.EQ.'Rang_CC')            then
#if (_COMPIL_==3)
         read(NumFich,'(I)') Rang_CC
#else
         read(NumFich,*) Rang_CC
#endif

      else if (Ligne.EQ.'TypeEtoileDef')      then
#if (_COMPIL_==3)
         read(NumFich,'(I)') TypeEtoileDef
#else
         read(NumFich,*) TypeEtoileDef
#endif

      else if (Ligne.EQ.'iRand_Seed')         then
#if (_COMPIL_==3)
         read(NumFich,'(I)') iRand_Seed
#else
         read(NumFich,*) iRand_Seed
#endif

      else if (Ligne.EQ.'iVerbeux')           then
#if (_COMPIL_==3)
         read(NumFich,'(I)') iVerbeux
#else
         read(NumFich,*) iVerbeux
#endif

      else if (Ligne.EQ.'i_fin')              then
#if (_COMPIL_==3)
         read(NumFich,'(I)') i_fin
#else
         read(NumFich,*) i_fin
#endif

      else if (Ligne.EQ.'i_ini')              then
#if (_COMPIL_==3)
         read(NumFich,'(I)') i_ini
#else
         read(NumFich,*) i_ini
#endif

         
      else
         write(0,100) Ligne 
 100     format(/,'!!! Le parametre ',A20,' n''a pas de sens !!!')
         read(NumFich,'(A)') car ! on saute une ligne
         call exit(1)
         
      end if
      goto 1

 2    continue
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_Fich_Param_EvolAmas(NumFich)
c----------------------------------------------------------------------
c     Ecrit les donnees necessaire a l'execution de programme 
c     d'evolution d'amas.
c
      include 'EvolAmas_param.f'
      include 'EvolAmas_Param_common.f' ! blocs common contenant les
                                        ! parametres
      integer NumFich
c......................................................................  
c   
      write(NumFich,1000)
 1000 format(
     $     '# --------------------------',/,
     $     '# Parametres pour EvolAmas.f',/,
     $     '# --------------------------',/,
     $     /)

      write(NumFich,1002) Type_CtrObj_def
 1002 format('Type_CtrObj_def     ',/,'  ',A    )
      write(NumFich,1003) NomSimul
 1003 format('NomSimul            ',/,'  ',A    )
      write(NumFich,1004) Coef_Trlx_GW_Capt
 1004 format('Coef_Trlx_GW_Capt   ',/,'  ',E12.4)
      write(NumFich,1005) CollAccrFrac_MS_on_BH
 1005 format('CollAccrFrac_MS_on_BH',/,'  ',E12.4)
      write(NumFich,1006) CollAccrFrac_MS_on_NS
 1006 format('CollAccrFrac_MS_on_NS',/,'  ',E12.4)
      write(NumFich,1007) CollAccrFrac_MS_on_WD
 1007 format('CollAccrFrac_MS_on_WD',/,'  ',E12.4)
      write(NumFich,1008) CollDminDestr_MS_vs_RMN
 1008 format('CollDminDestr_MS_vs_RMN',/,'  ',E12.4)
      write(NumFich,1009) CollMassLossFracDef
 1009 format('CollMassLossFracDef ',/,'  ',E12.4)
      write(NumFich,1010) FacMaxPG
 1010 format('FacMaxPG            ',/,'  ',E12.4)
      write(NumFich,1011) FacMinPG
 1011 format('FacMinPG            ',/,'  ',E12.4)
      write(NumFich,1012) FacNbPasDetEvap
 1012 format('FacNbPasDetEvap     ',/,'  ',E12.4)
      write(NumFich,1013) FacNbPasDetTrelExtremes
 1013 format('FacNbPasDetTrelExtremes',/,'  ',E12.4)
      write(NumFich,1014) FacNbPasInfo
 1014 format('FacNbPasInfo        ',/,'  ',E12.4)
      write(NumFich,1015) FacNbPasRecArbre
 1015 format('FacNbPasRecArbre    ',/,'  ',E12.4)
      write(NumFich,1016) FacNbPasRecCB3c
 1016 format('FacNbPasRecCB3c     ',/,'  ',E12.4)
      write(NumFich,1017) FacNbPasRecRelax
 1017 format('FacNbPasRecRelax    ',/,'  ',E12.4)
      write(NumFich,1018) FacNbPasRecTpsAmas
 1018 format('FacNbPasRecTpsAmas  ',/,'  ',E12.4)
      write(NumFich,1019) FacNbPasRecalcTH
 1019 format('FacNbPasRecalcTH    ',/,'  ',E12.4)
      write(NumFich,1020) FacNbPasSauv
 1020 format('FacNbPasSauv        ',/,'  ',E12.4)
      write(NumFich,1021) FacNbPasSauvPart
 1021 format('FacNbPasSauvPart    ',/,'  ',E12.4)
      write(NumFich,1022) FacNbPasTestDemSauv
 1022 format('FacNbPasTestDemSauv ',/,'  ',E12.4)
      write(NumFich,1023) FacNbPasTestDuplic
 1023 format('FacNbPasTestDuplic  ',/,'  ',E12.4)
      write(NumFich,1024) Fac_SldAvg_Tscale
 1024 format('Fac_SldAvg_Tscale   ',/,'  ',E12.4)
      write(NumFich,1025) Fac_dt_accr_M
 1025 format('Fac_dt_accr_M       ',/,'  ',E12.4)
      write(NumFich,1026) Fac_dt_accr_Npart
 1026 format('Fac_dt_accr_Npart   ',/,'  ',E12.4)
      write(NumFich,1027) FactR_CC
 1027 format('FactR_CC            ',/,'  ',E12.4)
      write(NumFich,1028) FactTrelMax
 1028 format('FactTrelMax         ',/,'  ',E12.4)
      write(NumFich,1029) FracPas_AjustPassif
 1029 format('FracPas_AjustPassif ',/,'  ',E12.4)
      write(NumFich,1030) Frac_Taccr
 1030 format('Frac_Taccr          ',/,'  ',E12.4)
      write(NumFich,1031) Frac_Tadiab
 1031 format('Frac_Tadiab         ',/,'  ',E12.4)
      write(NumFich,1032) Frac_Tcoll
 1032 format('Frac_Tcoll          ',/,'  ',E12.4)
      write(NumFich,1033) Frac_Tevap
 1033 format('Frac_Tevap          ',/,'  ',E12.4)
      write(NumFich,1034) Frac_Tevst
 1034 format('Frac_Tevst          ',/,'  ',E12.4)
      write(NumFich,1035) Frac_Tkick
 1035 format('Frac_Tkick          ',/,'  ',E12.4)
      write(NumFich,1036) Frac_Trelax
 1036 format('Frac_Trelax         ',/,'  ',E12.4)
      write(NumFich,1037) Frac_Tstaccr
 1037 format('Frac_Tstaccr        ',/,'  ',E12.4)
      write(NumFich,1038) Frac_dt_max
 1038 format('Frac_dt_max         ',/,'  ',E12.4)
      write(NumFich,1039) Gamma_relax
 1039 format('Gamma_relax         ',/,'  ',E12.4)
      write(NumFich,1040) M_TN_ini_def
 1040 format('M_TN_ini_def        ',/,'  ',E12.4)
      write(NumFich,1041) MasseEtoileDef
 1041 format('MasseEtoileDef      ',/,'  ',E12.4)
      write(NumFich,1042) Mmax_AccrByStar
 1042 format('Mmax_AccrByStar     ',/,'  ',E12.4)
      write(NumFich,1043) Mu_mol_elec
 1043 format('Mu_mol_elec         ',/,'  ',E12.4)
      write(NumFich,1044) R_Mar_ini_def
 1044 format('R_Mar_ini_def       ',/,'  ',E12.4)
      write(NumFich,1045) R_tid_disr_def_NB
 1045 format('R_tid_disr_def_NB   ',/,'  ',E12.4)
      write(NumFich,1046) RateAccrPerStar_MsunMyr
 1046 format('RateAccrPerStar_MsunMyr',/,'  ',E12.4)
      write(NumFich,1047) StellMet_def
 1047 format('StellMet_def        ',/,'  ',E12.4)
      write(NumFich,1048) T_fin
 1048 format('T_fin               ',/,'  ',E12.4)
      write(NumFich,1049) Tadiab_ini
 1049 format('Tadiab_ini          ',/,'  ',E12.4)
      write(NumFich,1050) TailleAmas_en_pc
 1050 format('TailleAmas_en_pc    ',/,'  ',E12.4)
      write(NumFich,1051) Tmax_GW_Capt_Gyr
 1051 format('Tmax_GW_Capt_Gyr    ',/,'  ',E12.4)
      write(NumFich,1052) Trlx_GW_def_Gyr
 1052 format('Trlx_GW_def_Gyr     ',/,'  ',E12.4)
      write(NumFich,1053) effic_conv_lum
 1053 format('effic_conv_lum      ',/,'  ',E12.4)
      write(NumFich,1054) exp_AccrByStar
 1054 format('exp_AccrByStar      ',/,'  ',E12.4)
      write(NumFich,1055) fact_b0_kick
 1055 format('fact_b0_kick        ',/,'  ',E12.4)
      write(NumFich,1056) frac_accr_Coll
 1056 format('frac_accr_Coll      ',/,'  ',E12.4)
      write(NumFich,1057) frac_accr_Dechir
 1057 format('frac_accr_Dechir    ',/,'  ',E12.4)
      write(NumFich,1058) frac_accr_EvolSt
 1058 format('frac_accr_EvolSt    ',/,'  ',E12.4)
      write(NumFich,1059) rNbEtoiles
 1059 format('rNbEtoiles          ',/,'  ',E12.4)
      write(NumFich,1060) tAccrByStar_1Msun_Myr
 1060 format('tAccrByStar_1Msun_Myr',/,'  ',E12.4)
      write(NumFich,1061) DeltaRang_Paire
 1061 format('DeltaRang_Paire     ',/,'  ',I10  )
      write(NumFich,1062) NbSECoucheCB3c
 1062 format('NbSECoucheCB3c      ',/,'  ',I10  )
      write(NumFich,1063) NbSECouchePot
 1063 format('NbSECouchePot       ',/,'  ',I10  )
      write(NumFich,1064) NbSauv_Conserve
 1064 format('NbSauv_Conserve     ',/,'  ',I10  )
      write(NumFich,1065) Rang_CC
 1065 format('Rang_CC             ',/,'  ',I10  )
      write(NumFich,1066) TypeEtoileDef
 1066 format('TypeEtoileDef       ',/,'  ',I10  )
      write(NumFich,1067) iRand_Seed
 1067 format('iRand_Seed          ',/,'  ',I10  )
      write(NumFich,1068) iVerbeux
 1068 format('iVerbeux            ',/,'  ',I10  )
      write(NumFich,1069) i_fin
 1069 format('i_fin               ',/,'  ',I10  )
      write(NumFich,1070) i_ini
 1070 format('i_ini               ',/,'  ',I10  )
      write(NumFich,2000)
 2000 format(
     $     /,
     $     '# --------------------------')
c
      end
c
c----------------------------------------------------------------------
c---- Valeurs par defaut des parametres de simulation
c----------------------------------------------------------------------
      block data BD_Init_Param
c----------------------------------------------------------------------
      include 'EvolAmas_param.f'
      include 'EvolAmas_Param_common.f'

      data Type_CtrObj_def      /'TN'        /
      data NomSimul             /'_EvolAmas' /
      data Coef_Trlx_GW_Capt    /1.0d0       /
      data CollAccrFrac_MS_on_BH /0.5d0       /
      data CollAccrFrac_MS_on_NS /0.0d0       /
      data CollAccrFrac_MS_on_WD /0.0d0       /
      data CollDminDestr_MS_vs_RMN /0.5d0    /
      data CollMassLossFracDef  /0.0d0       /
      data FacMaxPG             /1.5d0       /
      data FacMinPG             /0.5d0       /
      data FacNbPasDetEvap      /1.0d0       /
      data FacNbPasDetTrelExtremes /0.5d0    /
      data FacNbPasInfo         /1.0d0       /
      data FacNbPasRecArbre     /2.0d0       /
      data FacNbPasRecCB3c      /0.5d0       /
      data FacNbPasRecRelax     /10.0d0      /
      data FacNbPasRecTpsAmas   /1.0d0       /
      data FacNbPasRecalcTH     /0.5d0       /
      data FacNbPasSauv         /5.0d0       /
      data FacNbPasSauvPart     /1.0d0       /
      data FacNbPasTestDemSauv  /0.1d0       /
      data FacNbPasTestDuplic   /0.5d0       /
      data Fac_SldAvg_Tscale    /5.0d0       /
      data Fac_dt_accr_M        /0.01d0      /
      data Fac_dt_accr_Npart    /10.0d0      /
      data FactR_CC             /0.5d0       /
      data FactTrelMax          /1.0d3       /
      data FracPas_AjustPassif  /0.0d0       /
      data Frac_Taccr           /0.01d0      /
      data Frac_Tadiab          /0.01d0      /
      data Frac_Tcoll           /0.01d0      /
      data Frac_Tevap           /0.01d0      /
      data Frac_Tevst           /0.025d0     /
      data Frac_Tkick           /0.01d0      /
      data Frac_Trelax          /0.01d0      /
      data Frac_Tstaccr         /0.01d0      /
      data Frac_dt_max          /1.0d0       /
      data Gamma_relax          /0.14d0      /
      data M_TN_ini_def         /0.0d0       /
      data MasseEtoileDef       /1.0d0       /
      data Mmax_AccrByStar      /80.0        /
      data Mu_mol_elec          /1.13d0      /
      data R_Mar_ini_def        /1.0d30      /
      data R_tid_disr_def_NB    /1.0d-7      /
      data RateAccrPerStar_MsunMyr /1.0      /
      data StellMet_def         /0.02d0      /
      data T_fin                /1.0d5       /
      data Tadiab_ini           /0.001d0     /
      data TailleAmas_en_pc     /3.0d0       /
      data Tmax_GW_Capt_Gyr     /1.0d30      /
      data Trlx_GW_def_Gyr      /1.0d0       /
      data effic_conv_lum       /0.0d0       /
      data exp_AccrByStar       /0.0         /
      data fact_b0_kick         /2.0d0       /
      data frac_accr_Coll       /1.0d0       /
      data frac_accr_Dechir     /1.0d0       /
      data frac_accr_EvolSt     /1.0d0       /
      data rNbEtoiles           /1.0d6       /
      data tAccrByStar_1Msun_Myr /0.25       /
      data DeltaRang_Paire      /1           /
      data NbSECoucheCB3c       /1000        /
      data NbSECouchePot        /25          /
      data NbSauv_Conserve      /5           /
      data Rang_CC              /10          /
      data TypeEtoileDef        /1           /
      data iRand_Seed           /133111      /
      data iVerbeux             /1           /
      data i_fin                /2000000000  /
      data i_ini                /0           /
c
      end
c----------------------------------------------------------------------
