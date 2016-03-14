c======================================================================
c     constantes pour EvolAmas.F
c======================================================================
c
c---- fraction de pas avec relaxation
c
      double precision FracRelax
      parameter (FracRelax = 1.0d0)
c
c---- precision de l'integration Porb,Trelax
c
      double precision eps_dtrel
      parameter (eps_dtrel=0.1d0)
c
c---- valeur minimale du rayon d'une SE
c
      double precision R_SE_min
      parameter (R_SE_min=1.0d-6)
c
c---- derive max de l'energie par pas d'evolution (valeur relative)
c
      double precision drv_E_max
      parameter (drv_E_max=1.0d-10)
c
c---- limitation des R tres petits
c
      double precision FactSurv_Rmin,FactCorr_Rmin
      parameter (FactSurv_Rmin=0.5d0,FactCorr_Rmin=0.5d0) ! valeurs habituelles : 0.5 et 0.5
c
c---- Nom des differents types de fichier
c
      character*(*)
     $     TypeAMAS,TypeARBRE,TypeTEMPS, TypePG, TypeGRILLE, TypeRAND,
     $     TypeETOILES, TypeRAYLAG, TypeFRAYLAG, TypeANILAG, TypeLOG,
     $     TypeRELAX, TypeTN, TypeTH, TypeMAREE, TypeSEGR, TypeFSEGRLAG,
     $     TypeCOLL, TypeLC, TypeCONS, TypeCB3c, TypeSuivSE, TypeCaptGW,
     $     TypeDuplic, TypeEvap, TypeLogSauv, TypeTscale, TypeLIMEXT,
     $     TypeMSpec, TypeStEvol, TypeFSTSLAG, TypeMbinRad, TypeMbins,
     $     TypeExcColl, TypeStrEnc, TypeSubPop, TypeOFGrid, TypeNatKcks,
     $     TypeLagQuant, TypeCumPcoll
      parameter (
     $     TypeAMAS    = 'AMAS',
     $     TypeARBRE   = 'ARBRE',
     $     TypeTEMPS   = 'TEMPS',
     $     TypePG      = 'PG',
     $     TypeGRILLE  = 'GRILLE',
     $     TypeRAND    = 'RAND',
     $     TypeETOILES = 'ETOILES',
     $     TypeANILAG  = 'AniLag',
     $     TypeSEGR    = 'Segr',
     $     TypeFSEGRLAG = 'FracSegrLag',
     $     TypeCOLL    = 'Coll',
     $     TypeLC      = 'LC',
     $     TypeRAYLAG  = 'RayLag',
     $     TypeFRAYLAG = 'FracRayLag',
     $     TypeLOG     = 'Log',
     $     TypeRELAX   = 'GlobRelax',
     $     TypeTH      = 'TH',
     $     TypeTN      = 'TN'
     $     )
      parameter (
     $     TypeMAREE   = 'MAREE',
     $     TypeCONS    = 'CONS',
     $     TypeCB3c    = 'CB3c', ! chauffage par binaires 3-corps
     $     TypeSuivSE  = 'SuivSE',
     $     TypeCaptGW  = 'CaptGW',
     $     TypeDuplic  = 'Duplic',
     $     TypeEvap    = 'Evap',
     $     TypeLogSauv = 'LogSauv',
     $     TypeTscale  = 'Tscale',
     $     TypeLIMEXT  = 'LIMEXT',
     $     TypeMSpec   = 'MSpec',
     $     TypeStEvol  = 'StEvol',
     $     TypeFSTSLAG = 'FracSTSLag',
     $     TypeMbinRad = 'MbinRad',
     $     TypeMbins   = 'Mbins',
     $     TypeExcColl = 'CollExcept',
     $     TypeStrEnc  = 'StrongEncounters',
     $     TypeSubPop  = 'SubPop',
     $     TypeOFGrid  = 'OFGrid',
     $     TypeNatKcks = 'NatKicks',
     $     TypeLagQuant= 'LagQuant',
     $     TypeCumPColl= 'CumPColl'
     $     )
c
c     REGLE: N'EMPLOYER QUE DES NUM D'UNITES>=50
c

c
c---- Numero des fichiers "sequentiels" (restant ouverts)
c
      integer
     $     iFich_RayLag,iFich_Log,iFich_Relax,iFich_Segr,iFich_Evap,
     $     iFich_AniLag,iFich_Coll,iFich_SuivSE,iFich_LC,iFich_CaptGW,
     $     iFich_Duplic,iFich_LogSauv,iFich_Tscale,iFich_MSpec,
     $     iFich_StEvol,iFich_FSTSLag,iFich_MbinRad,iFich_StrEnc,
     $     iFich_ExcColl,iFich_OFGrid,iFich_NatKcks,iFich_LagQuant,
     $     iFich_CumPcoll,
     $     iFich_Prov
      parameter (
     $     iFich_RayLag =50,
     $     iFich_Log    =51,
     $     iFich_Relax  =52,
     $     iFich_Segr   =53,
     $     iFich_Evap   =54,
     $     iFich_AniLag =55,
     $     iFich_Coll   =56,
     $     iFich_SuivSE =57,
     $     iFich_LC     =58,
     $     iFich_CaptGW =59,
     $     iFich_Duplic =60,
     $     iFich_LogSauv=61,
     $     iFich_Tscale =62,
     $     iFich_MSpec  =63,
     $     iFich_StEvol =64,
     $     iFich_FSTSLag=65,
     $     iFich_MbinRad=66,
     $     iFich_StrEnc =67,
     $     iFich_ExcColl=68,
     $     iFich_OFGrid =69,
     $     iFich_NatKcks=70,
     $     iFich_LagQuant=71,
     $     iFich_CumPcoll=72,
     $     iFich_Prov   =99
     $     )
c
c---- Encore des numeros de fichier...
c
      integer                   ! Suivi des types stellaires
     $     iFich_STS1,iFich_STS2,iFich_STS3,iFich_STS4,iFich_STS5
      parameter (
     $     iFich_STS1=81,
     $     iFich_STS2=82,
     $     iFich_STS3=83,
     $     iFich_STS4=84,
     $     iFich_STS5=85
     $     )
c
c---- Noms de fichiers
c
      character*(*) NomFichCroissForcTN, NomFichGrilleSPH,
     $     NomFichListeSEPart,NomFichMRrelMS,NomFichParamEvolStell,
     $     NomFichDemSauv,
     $     NomFichDemArret,NomFichPID,NomFichINFO, NomFichListSubPop,
     $     NomFichOFGridCellPresc, NomFichNonInterac
      parameter (
     $     NomFichCroissForcTN='Param_Croiss_Forc_TN.asc',
     $     NomFichGrilleSPH   ='Grille_ResColl_SPH.asc',
     $     NomFichListeSEPart ='Liste_SE_a_suivre.asc',
     $     NomFichMRrelMS     ='MassRadiusMS.asc',
     $     NomFichParamEvolStell = 'ParamStellEvol.asc',
     $     NomFichDemSauv  ='_SAUV_DEMANDEE_',
     $     NomFichDemArret ='_STOP_',
     $     NomFichPID      ='_PID_',
     $     NomFichINFO     ='_INFO_',
     $     NomFichListSubPop = 'ListSubPop.asc',
     $     NomFichOFGridCellPresc = 'OFGridCellPresc.asc',
     $     NomFichNonInterac = 'NonInterac.asc'
     $     )
