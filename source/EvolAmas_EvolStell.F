#define _WIND_TYPE_ 3 /* 0: No winds */
                      /* 1: Winds according to LH95 */
                      /* 2: Winds according to Meynet et al. 94   */
                      /* 3: Winds according to Meynet & Meader 00 */
#define _FANCY_REMNANT_ 2 /* 2 for Belczynski data; then _WIND_TYPE_ doesn't matter */
#define _METALLICITY_ 1 /* 0: z=0.0001 1: z=0.02; 0 is default */

c----------------------------------------------------------------------
      subroutine InitEvolStell
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_EvolStell_common.f'
      include 'EvolAmas_param.f' ! for NomFichParamEvolStell
c
c local variables :
c ^^^^^^^^^^^^^^^^^ 
      logical lExist
      character*64 key
      double precision value
c
c functions :
c ^^^^^^^^^^^
      integer LongChaine
c......................................................................
c
      StringTypeEvolStell=''

#if (_FANCY_REMNANT_==0)
      StringTypeEvolStell = 'SimpleRemnants'
#elif (_FANCY_REMNANT_==1)
      StringTypeEvolStell = 'FancyRemnants'
#elif (_FANCY_REMNANT_==2)
      StringTypeEvolStell = 'KB05Remnants'
#  if (_METALLICITY_==0)
      StringTypeEvolStell =
     $     StringTypeEvolStell(1:LongChaine(StringTypeEvolStell)) //
     $     '-Z0.001'
#  elif (_METALLICITY_==1)
      StringTypeEvolStell =
     $     StringTypeEvolStell(1:LongChaine(StringTypeEvolStell)) //
     $     '-Z0.02'
#  else
%%%%% DO NOT COMPILE; INCONSISTANT CHOICE OF STELLAR EVOLUTION PARAMETERS
#  endif
#else
%%%%% DO NOT COMPILE; INCONSISTANT CHOICE OF STELLAR EVOLUTION PARAMETERS
#endif

#if (_FANCY_REMNANT_<2)
#  if (_WIND_TYPE_==0)
      StringTypeEvolStell =
     $     StringTypeEvolStell(1:LongChaine(StringTypeEvolStell)) //
     $     '-NoWinds'
#  elif (_WIND_TYPE_==1)
      StringTypeEvolStell =
     $     StringTypeEvolStell(1:LongChaine(StringTypeEvolStell)) //
     $     '-WindsLH95'
#  elif (_WIND_TYPE_==2)
      StringTypeEvolStell =
     $     StringTypeEvolStell(1:LongChaine(StringTypeEvolStell)) //
     $     '-WindsMMSSC94'
#  elif (_WIND_TYPE_==3)
      StringTypeEvolStell =
     $     StringTypeEvolStell(1:LongChaine(StringTypeEvolStell)) //
     $     '-WindsMM00'
#  else
%%%%% DO NOT COMPILE; INCONSISTANT CHOICE OF STELLAR EVOLUTION PARAMETERS
#  endif
#endif

/* Example for file 'ParamStellEvol.asc'

# Maximum mass for WD formation
Mzams_max_WD
 8.0

Mzams_max_NS
 22.0

Mzams_max_fallbackBH
 42.0

M_WD
 6.0

M_NS
 1.4

M_BH
 10.0

# Parameters for natal kicks: 2 Maxwellians of dispersion NatKick_Sig1D_1 and NatKick_Sig1D_2
# Values larger than 1e10 force "clean" ejection (without energy conservation issues)
#
# NatKick_Prob_1 is the probability of selecting the first Maxwellian.

NatKick_Sig1D_1
 265.0

NatKick_Sig1D_2
 0.0

NatKick_Prob_1
 1.0

*/
      inquire(file=NomFichParamEvolStell,exist=lExist)

      if (lExist) then
#if (_VERB_STELLEVOL_>0)
         write(0,*) '> InitEvolStell: file ',NomFichParamEvolStell,
     $        ' found'
#endif
         open(iFich_Prov,file=NomFichParamEvolStell)
         do while (.TRUE.)
            key='###'
            do while (key(1:1).EQ.'#' .OR. key(1:1).EQ.' ' .OR. key.EQ.'')
               read(iFich_Prov,'(A)',END=99) key
            end do
            read(iFich_Prov,*,END=99) value
#if (_VERB_STELLEVOL_>0)
            write(0,*) '> InitEvolStell: key=',key(1:index(key,' ')-1),
     $           ' value=',value
#endif
            if      (key.EQ.'Mzams_max_WD') then
               Mzams_max_WD = value
            else if (key.EQ.'Mzams_max_NS') then
               Mzams_max_NS = value
            else if (key.EQ.'Mzams_max_fallbackBH') then
               Mzams_max_fallbackBH = value
            else if (key.EQ.'M_max_WD'    ) then
               M_max_WD = value
            else if (key.EQ.'M_max_NS'    ) then
               M_max_NS = value
            else if (key.EQ.'Mremn_min_BH') then
               Mremn_min_BH = value
            else if (key.EQ.'M_WD'        ) then
               M_WD = value
            else if (key.EQ.'M_NS'        ) then
               M_NS = value
            else if (key.EQ.'M_BH'        ) then
               M_BH = value
            else if (key.EQ.'NatKick_Sig1D_1') then
               NatKick_Sig1D_1 = value
            else if (key.EQ.'NatKick_Sig1D_2') then
               NatKick_Sig1D_2 = value
            else if (key.EQ.'NatKick_Prob_1') then
               NatKick_Prob_1 = value
            else 
#if (_VERB_STELLEVOL_>0)
               write(0,*) '> InitEvolStell: unknown key, skipped'
#endif
               continue
            end if
         end do
 99      continue
         close(iFich_Prov)
      else
         continue
#if (_VERB_STELLEVOL_>0)
         write(0,*) '> InitEvolStell: file ',NomFichParamEvolStell,
     $        ' not found'
#endif
      end if

      lInitEvolStell = .TRUE.
c
      end
c
c----------------------------------------------------------------------
      subroutine write_as_comments_ParamEvolStell(iFich)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_EvolStell_common.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c local variables :
c ^^^^^^^^^^^^^^^^^ 
c......................................................................
c
      write(iFich,110)
 110  format('# Stellar evolution parameters',/,'#')
      write(iFich,100) 'Mzams_max_WD',Mzams_max_WD
      write(iFich,100) 'Mzams_max_NS',Mzams_max_NS
      write(iFich,100) 'M_max_WD',M_max_WD
      write(iFich,100) 'M_max_NS',M_max_NS
      write(iFich,100) 'Mzams_max_fallbackBH',Mzams_max_fallbackBH
      write(iFich,100) 'Mremn_min_BH',Mremn_min_BH
      write(iFich,100) 'M_WD',M_WD
      write(iFich,100) 'M_NS',M_NS
      write(iFich,100) 'M_BH',M_BH
      write(iFich,100) 'NatKick_Sig1D_1',NatKick_Sig1D_1
      write(iFich,100) 'NatKick_Sig1D_2',NatKick_Sig1D_2
      write(iFich,100) 'NatKick_Prob_1',NatKick_Prob_1
 100  format('#   ',A,' = ',G10.4)
c
      end
c
c----------------------------------------------------------------------
      subroutine write_ParamEvolStell(iFich)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_EvolStell_common.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c local variables :
c ^^^^^^^^^^^^^^^^^ 
c......................................................................
c
      write(iFich,110)
 110  format('# Stellar evolution parameters',/)
      write(iFich,100) 'Mzams_max_WD',Mzams_max_WD
      write(iFich,100) 'Mzams_max_NS',Mzams_max_NS
      write(iFich,100) 'M_max_WD',M_max_WD
      write(iFich,100) 'M_max_NS',M_max_NS
      write(iFich,100) 'Mzams_max_fallbackBH',Mzams_max_fallbackBH
      write(iFich,100) 'Mremn_min_BH',Mremn_min_BH
      write(iFich,100) 'M_WD',M_WD
      write(iFich,100) 'M_NS',M_NS
      write(iFich,100) 'M_BH',M_BH
      write(iFich,100) 'NatKick_Sig1D_1',NatKick_Sig1D_1
      write(iFich,100) 'NatKick_Sig1D_2',NatKick_Sig1D_2
      write(iFich,100) 'NatKick_Prob_1',NatKick_Prob_1
 100  format(A,/,'   ',E12.5,/)
c
      end
c
#if (_EVOL_STELL_==1) /* Simple stellar evol (MS-->remnant) */

c======================================================================
c     routine pour faire evoluer une etoile de la ZAMS a un age
c     "arbitraire"
c======================================================================
c----------------------------------------------------------------------
      subroutine Evoluer_Etoile(M_ZAMS,Age,M_evol,iT_evol)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_EvolStell_common.f'
      include 'EvolAmas_Etoiles_param.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision          ! ENTREE:
     $     M_ZAMS,              ! Masse sur la ZAMS (Msol)
     $     Age                  ! Age depuis la ZAMS (yrs)

      double precision          ! SORTIE:
     $     M_evol               ! Masse a l'age demande
      integer*1 iT_evol              ! Type d'etoile a l'age demande
c
c var. locales :
c ^^^^^^^^^^^^^^
      double precision C_TN
      logical lPrem
      data lPrem /.TRUE./
      save lPrem
c
c fonctions :
c ^^^^^^^^^^^
      double precision T_MS     ! Temps de seq. principale
c......................................................................
c
c     Le modele utilse est tres simple :     
c     Pas d'evolution sur la SP puis passage a l'etat de remanant
c
      if (lPrem) then
         if (.NOT.lInitEvolStell) call InitEvolStell
c$$$            ! Selon Binney + Merrifield (p.273)
c$$$         C_TN = (Msup_BH-Minf_BH)/log10(90.0d0/Mmax_form_NS)
         lPrem = .FALSE.
      end if
                                ! SEQUENCE PRINCIPALE
                                ! ===================
      if (Age.LE.0.0d0 .OR. Age.LT.T_MS(M_ZAMS)) then
         M_evol = M_ZAMS        ! Pas de perte de masse !
         iT_evol = TypeEtoile_MS
      else
                                ! REMANANT
                                ! ========
         call Remnant(M_ZAMS,M_evol,iT_evol)
c
c     Compute and remember natal kick (kept in common block)
c
#if (_NATKICKS_>0)
         call Get_NatalKick(M_zams,M_evol,iT_evol,
     $        NatalKickLastEvol)
#else
         NatalKickLastEvol = 0.0d0
#endif 
      end if
      !write(0,*) '^^^^^',M_ZAMS,M_evol,iT_evol !#################
c
      end
c

c----------------------------------------------------------------------
      subroutine Remnant(M_zams,M_remn,iT_remn)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_EvolStell_common.f'
      include 'EvolAmas_Etoiles_param.f'
c
c arguments :
c ^^^^^^^^^^^
                                ! Input:
      double precision M_zams 
                                ! Outputs:
      double precision M_remn 
      integer*1 iT_remn
c
c function :
c ^^^^^^^^^^
      double precision M_remnant
c......................................................................
c
      if (.NOT.lInitEvolStell) call InitEvolStell

      M_remn = M_remnant(M_zams)
      !write(0,*) '%%%%%',M_zams,M_remn !#################

      iT_remn = TypeEtoile_NS

#if (_FANCY_REMNANT_!=2)
      if (M_zams.LT.Mzams_max_WD .AND. M_remn.LT.M_max_WD) then
         iT_remn = TypeEtoile_WD 
      end if 

#else /* use Belczynski data */
#   if (_METALLICITY_==0) /* Z=0.0001 */
      if (M_zams.LT.6.426d0) then
         iT_remn = TypeEtoile_WD 
      end if
#   elif (_METALLICITY_==1) /* Z=0.02 */
      if (M_zams.LT.7.98d0) then
         iT_remn = TypeEtoile_WD 
      end if
#   endif
#endif

      if (M_remn.GE.M_max_NS) then
         iT_remn = TypeEtoile_BH
      end if
c
      end
c
#if (_FANCY_REMNANT_==0)
c----------------------------------------------------------------------
      double precision function M_remnant(M_zams)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_EvolStell_common.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision M_zams    
c......................................................................
c
      if (.NOT.lInitEvolStell) call InitEvolStell

      if (M_zams.LE.Mzams_max_WD) then
         M_remnant=min(0.95d0*M_zams,M_WD)
      else if (M_zams.LE.Mzams_max_NS) then
         M_remnant=min(0.95d0*M_zams,M_NS)
      else
         M_remnant=min(0.95d0*M_zams,M_BH)
      end if
c
      end
c
#endif /* (_FANCY_REMNANT_==0) */

#if (_FANCY_REMNANT_>0)
c----------------------------------------------------------------------
      double precision function M_remnant(M_zams)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      double precision M_zams
c
c var locales :
c ^^^^^^^^^^^^^
      double precision lm,lmr,lm1,lm2
      integer i1,i2
      integer DimRelMremn 
      parameter ( DimRelMremn  = 351 )
      double precision lM_remn(DimRelMremn)
      save lM_remn
c
c     This data has been generated by the macro Genere_Mremn_fancy 
c     in Stell_ingred.smac; se explanations there...
c
      double precision lMmin, dlM
      parameter (lMmin=-1.0d0, dlM = 0.01d0)

#if (_FANCY_REMNANT_<2)

      data lM_remn /
     &     -1.02228e+00, -1.01228e+00, -1.00228e+00, -9.92276e-01, 
     &     -9.82276e-01, -9.72276e-01, -9.62276e-01, -9.52276e-01, 
     &     -9.42276e-01, -9.32276e-01, -9.22276e-01, -9.12276e-01, 
     &     -9.02276e-01, -8.92276e-01, -8.82276e-01, -8.72276e-01, 
     &     -8.62276e-01, -8.52276e-01, -8.42276e-01, -8.32276e-01, 
     &     -8.22276e-01, -8.12276e-01, -8.02276e-01, -7.92276e-01, 
     &     -7.82276e-01, -7.72276e-01, -7.62276e-01, -7.52276e-01, 
     &     -7.42276e-01, -7.32276e-01, -7.22276e-01, -7.12276e-01, 
     &     -7.02276e-01, -6.92276e-01, -6.82276e-01, -6.72276e-01, 
     &     -6.62276e-01, -6.52276e-01, -6.42276e-01, -6.32276e-01, 
     &     -6.22276e-01, -6.12276e-01, -6.02276e-01, -5.92276e-01, 
     &     -5.82276e-01, -5.72276e-01, -5.62276e-01, -5.52276e-01, 
     &     -5.42276e-01, -5.32276e-01, -5.22276e-01, -5.12276e-01, 
     &     -5.02276e-01, -4.92276e-01, -4.82276e-01, -4.72276e-01, 
     &     -4.62276e-01, -4.52276e-01, -4.42276e-01, -4.32276e-01, 
     &     -4.22276e-01, -4.12276e-01, -4.02276e-01, -3.98356e-01, 
     &     -3.95336e-01, -3.92316e-01, -3.89296e-01, -3.86276e-01, 
     &     -3.83256e-01, -3.80236e-01, -3.77216e-01, -3.74196e-01, 
     &     -3.71176e-01, -3.68156e-01, -3.65136e-01, -3.62116e-01, 
     &     -3.59096e-01, -3.56076e-01, -3.53056e-01, -3.50036e-01, 
     &     -3.47016e-01, -3.43996e-01, -3.40977e-01, -3.37957e-01, 
     &     -3.34937e-01, -3.31917e-01, -3.28897e-01, -3.25877e-01, 
     &     -3.22857e-01, -3.19837e-01, -3.16817e-01, -3.13797e-01, 
     &     -3.10777e-01, -3.07757e-01, -3.04737e-01, -3.01717e-01, 
     &     -2.98697e-01, -2.95677e-01, -2.92657e-01, -2.89637e-01, 
     &     -2.86617e-01, -2.83597e-01, -2.80577e-01, -2.77557e-01, 
     &     -2.74537e-01, -2.71517e-01, -2.68497e-01, -2.65477e-01, 
     &     -2.62457e-01, -2.59437e-01, -2.56418e-01, -2.53398e-01, 
     &     -2.50378e-01, -2.47358e-01, -2.44338e-01, -2.41318e-01, 
     &     -2.38298e-01, -2.35278e-01, -2.32258e-01, -2.29238e-01, 
     &     -2.26218e-01, -2.23198e-01, -2.20178e-01, -2.17158e-01, 
     &     -2.14138e-01, -2.11118e-01, -2.08098e-01, -2.04120e-01, 
     &     -2.00118e-01, -1.96116e-01, -1.94636e-01, -1.93590e-01, 
     &     -1.92545e-01, -1.91452e-01, -1.90223e-01, -1.88970e-01, 
     &     -1.84674e-01, -1.80198e-01, -1.73271e-01, -1.66343e-01, 
     &     -1.59342e-01, -1.53144e-01, -1.49379e-01, -1.45614e-01, 
     &     -1.41849e-01, -1.37816e-01, -1.32867e-01, -1.27918e-01, 
     &     -1.22969e-01, -1.18020e-01, -1.13071e-01, -1.08122e-01, 
     &     -1.03173e-01, -9.82237e-02, -9.32746e-02, -8.82432e-02, 
     &     -8.25480e-02, -7.68527e-02, -7.11575e-02, -6.54622e-02, 
     &     -5.97669e-02, -5.40717e-02, -4.83764e-02, -4.26812e-02, 
     &     -3.69859e-02, -3.12906e-02, -2.53138e-02, -1.84779e-02, 
     &     -1.16421e-02, -4.80627e-03, 2.02956e-03, 8.86539e-03, 
     &     1.57012e-02, 2.25371e-02, 2.93729e-02, 3.62087e-02, 
     &     4.37286e-02, 5.13655e-02, 5.90023e-02, 6.66392e-02, 
     &     7.42761e-02, 8.19703e-02, 9.03388e-02, 9.87074e-02,
     &     1.07076e-01, 1.15444e-01, 1.23813e-01, 1.32182e-01,  
     &     1.40550e-01, 1.48919e-01, 1.57287e-01, 9.28516e-02, 
     &     9.54303e-02, 9.80090e-02, 1.00588e-01, 1.03166e-01, 
     &     1.05745e-01, 1.08324e-01, 1.10902e-01, 1.13481e-01, 
     &     1.16060e-01, 1.18639e-01, 1.21217e-01, 1.23796e-01, 
     &     1.26375e-01, 1.28953e-01, 1.31532e-01, 1.34111e-01, 
     &     1.36689e-01, 1.39268e-01, 1.41847e-01, 1.44425e-01, 
     &     1.47004e-01, 1.49583e-01, 1.52161e-01, 1.54740e-01, 
     &     1.57319e-01, 1.59898e-01, 1.62476e-01, 1.65055e-01, 
     &     1.67634e-01, 1.70212e-01, 1.72791e-01, 1.75370e-01, 
     &     1.77948e-01, 1.80527e-01, 1.83106e-01, 1.85684e-01, 
     &     2.00502e-01, 2.49468e-01, 2.98434e-01, 3.47400e-01, 
     &     3.96366e-01, 4.45332e-01, 4.94299e-01, 5.43265e-01, 
     &     5.92231e-01, 6.41197e-01, 6.90163e-01, 7.39129e-01,  ! end of common data
#if (_WIND_TYPE_==0) /* No winds */ 
     &     7.88095e-01, 8.08082e-01, 8.18768e-01, 8.29453e-01, 
     &     8.69889e-01, 9.09393e-01, 9.26167e-01, 9.42941e-01, 
     &     9.59715e-01, 9.76489e-01, 9.93262e-01, 1.01004e+00,
     &     1.02681e+00, 1.04358e+00, 1.06036e+00, 1.07713e+00, 
     &     1.09391e+00, 1.11068e+00, 1.12745e+00, 1.14423e+00, 
     &     1.16100e+00, 1.17777e+00, 1.19455e+00, 1.60772e+00, 
     &     1.61772e+00, 1.62772e+00, 1.63772e+00, 1.64772e+00, 
     &     1.65772e+00, 1.66772e+00, 1.67772e+00, 1.68772e+00, 
     &     1.69772e+00, 1.70772e+00, 1.71772e+00, 1.72772e+00, 
     &     1.73772e+00, 1.74772e+00, 1.75772e+00, 1.76772e+00, 
     &     1.77772e+00, 1.78772e+00, 1.79772e+00, 1.80772e+00, 
     &     1.81772e+00, 1.82772e+00, 1.83772e+00, 1.84772e+00, 
     &     1.85772e+00, 1.86772e+00, 1.87772e+00, 1.88772e+00, 
     &     1.89772e+00, 1.90772e+00, 1.91772e+00, 1.92772e+00, 
     &     1.93772e+00, 1.94772e+00, 1.95772e+00, 1.96772e+00, 
     &     1.97772e+00, 1.98772e+00, 1.99772e+00, 2.00772e+00, 
     &     2.01772e+00, 2.02772e+00, 2.03772e+00, 2.04772e+00, 
     &     2.05772e+00, 2.06772e+00, 2.07772e+00, 2.08772e+00, 
     &     2.09772e+00, 2.10772e+00, 2.11772e+00, 2.12772e+00, 
     &     2.13772e+00, 2.14772e+00, 2.15772e+00, 2.16772e+00, 
     &     2.17772e+00, 2.18772e+00, 2.19772e+00, 2.20772e+00, 
     &     2.21772e+00, 2.22772e+00, 2.23772e+00, 2.24772e+00, 
     &     2.25772e+00, 2.26772e+00, 2.27772e+00, 2.28772e+00, 
     &     2.29772e+00, 2.30772e+00, 2.31772e+00, 2.32772e+00, 
     &     2.33772e+00, 2.34772e+00, 2.35772e+00, 2.36772e+00, 
     &     2.37772e+00, 2.38772e+00, 2.39772e+00, 2.40772e+00, 
     &     2.41772e+00, 2.42772e+00, 2.43772e+00, 2.44772e+00, 
     &     2.45772e+00, 2.46772e+00, 2.47772e+00 /
#elif (_WIND_TYPE_==1) /* Winds according to LH95 */ 
     &     7.88298e-01, 8.08082e-01, 8.18768e-01, 8.29453e-01, 
     &     8.71427e-01, 9.09513e-01, 9.27597e-01, 9.45681e-01, 
     &     9.63764e-01, 9.81848e-01, 9.99932e-01, 1.01802e+00, 
     &     1.03610e+00, 1.05418e+00, 1.07227e+00, 1.05980e+00, 
     &     1.04119e+00, 1.02258e+00, 1.00397e+00, 9.85365e-01, 
     &     9.67060e-01, 9.51607e-01, 9.36154e-01, 9.20701e-01, 
     &     9.05249e-01, 8.89796e-01, 8.74278e-01, 8.55048e-01, 
     &     8.35819e-01, 8.16589e-01, 7.97360e-01, 7.78130e-01, 
     &     7.58901e-01, 7.39672e-01, 7.20442e-01, 7.01213e-01, 
     &     6.81983e-01, 6.62754e-01, 6.43524e-01, 6.24295e-01, 
     &     6.05066e-01, 5.85836e-01, 5.66607e-01, 5.47377e-01, 
     &     5.28148e-01, 5.08918e-01, 4.89689e-01, 4.70459e-01, 
     &     4.51230e-01, 4.32001e-01, 4.12771e-01, 3.93542e-01, 
     &     3.74312e-01, 3.55083e-01, 3.35853e-01, 3.16624e-01, 
     &     2.97395e-01, 2.78165e-01, 2.58936e-01, 2.39706e-01, 
     &     2.20477e-01, 2.01247e-01, 1.82018e-01, 1.62789e-01, 
     &     1.43559e-01, 1.24330e-01, 1.05100e-01, 8.58708e-02, 
     &     6.66413e-02, 4.74119e-02, 2.81824e-02, 8.95301e-03, 
     &     -1.02764e-02, -2.95059e-02, -4.87353e-02, -6.79647e-02, 
     &     -8.71942e-02, -1.06424e-01, -1.25653e-01, -1.44882e-01, 
     &     -1.64112e-01, -1.83341e-01, -2.02571e-01, -2.21800e-01, 
     &     -2.41030e-01, -2.60259e-01, -2.79489e-01, -2.98718e-01, 
     &     -3.17947e-01, -3.37177e-01, -3.56406e-01, -3.75636e-01, 
     &     -3.94865e-01, -4.14095e-01, -4.33324e-01, -4.52553e-01, 
     &     -4.71783e-01, -4.91012e-01, -5.10242e-01, -5.29471e-01, 
     &     -5.48701e-01, -5.67930e-01, -5.87160e-01, -6.06389e-01, 
     &     -6.25618e-01, -6.44848e-01, -6.64077e-01, -6.83307e-01, 
     &     -7.02536e-01, -7.21766e-01, -7.40995e-01 /
#elif (_WIND_TYPE_==2) /* Winds according to Meynet et al 94 */ 
     &     7.88095e-01, 8.08082e-01, 8.18768e-01, 8.29453e-01, 
     &     8.69889e-01, 8.95924e-01, 8.77462e-01, 8.59417e-01, 
     &     8.41806e-01, 8.24803e-01, 8.08431e-01, 7.92856e-01, 
     &     7.78111e-01, 7.64354e-01, 7.51624e-01, 7.40071e-01, 
     &     7.29743e-01, 7.20783e-01, 7.13246e-01, 7.07266e-01, 
     &     7.02905e-01, 7.00161e-01, 6.99098e-01, 6.99373e-01, 
     &     7.01040e-01, 7.03716e-01, 7.07440e-01, 7.11843e-01, 
     &     7.16951e-01, 7.22407e-01, 7.28224e-01, 7.34059e-01, 
     &     7.39913e-01, 7.45454e-01, 7.50669e-01, 7.55242e-01, 
     &     7.59147e-01, 7.62081e-01, 7.64007e-01, 7.64819e-01, 
     &     7.64472e-01, 7.63139e-01, 7.60781e-01, 7.57588e-01, 
     &     7.53526e-01, 7.48779e-01, 7.43320e-01, 7.37327e-01, 
     &     7.30779e-01, 7.23847e-01, 7.16516e-01, 7.08952e-01, 
     &     7.01145e-01, 6.93245e-01, 6.85246e-01, 6.77203e-01, 
     &     6.69113e-01, 6.60985e-01, 6.52817e-01, 6.44618e-01, 
     &     6.36387e-01, 6.28130e-01, 6.19848e-01, 6.11548e-01, 
     &     6.03229e-01, 5.94899e-01, 5.86557e-01, 5.78215e-01, 
     &     5.69873e-01, 5.61531e-01, 5.53189e-01, 5.44847e-01, 
     &     5.36505e-01, 5.28163e-01, 5.19821e-01, 5.11479e-01, 
     &     5.03137e-01, 4.94795e-01, 4.86453e-01, 4.78111e-01, 
     &     4.69769e-01, 4.61426e-01, 4.53084e-01, 4.44742e-01, 
     &     4.36400e-01, 4.28058e-01, 4.19716e-01, 4.11374e-01, 
     &     4.03032e-01, 3.94690e-01, 3.86348e-01, 3.78006e-01, 
     &     3.69664e-01, 3.61322e-01, 3.52980e-01, 3.44638e-01, 
     &     3.36296e-01, 3.27954e-01, 3.19612e-01, 3.11270e-01, 
     &     3.02928e-01, 2.94586e-01, 2.86244e-01, 2.77902e-01, 
     &     2.69559e-01, 2.61217e-01, 2.52875e-01, 2.44533e-01, 
     &     2.36191e-01, 2.27849e-01, 2.19507e-01 /
#elif (_WIND_TYPE_==3) /* Winds according to Meynet & Meader 00 */
     &     7.88095e-01, 8.08082e-01, 8.18768e-01, 8.29453e-01, 
     &     8.69889e-01, 9.09393e-01, 9.26167e-01, 9.42941e-01, 
     &     9.59715e-01, 9.76489e-01, 9.93262e-01, 1.01004e+00, 
     &     1.02681e+00, 1.04358e+00, 1.05528e+00, 1.05884e+00, 
     &     1.06232e+00, 1.06574e+00, 1.06881e+00, 1.07164e+00, 
     &     1.07387e+00, 1.07568e+00, 1.07688e+00, 1.07762e+00, 
     &     1.07788e+00, 1.07779e+00, 1.07735e+00, 1.07666e+00, 
     &     1.07576e+00, 1.07470e+00, 1.07357e+00, 1.07239e+00, 
     &     1.07127e+00, 1.07019e+00, 1.06931e+00, 1.06858e+00, 
     &     1.06818e+00, 1.06802e+00, 1.06833e+00, 1.06896e+00, 
     &     1.07007e+00, 1.07152e+00, 1.07340e+00, 1.07561e+00, 
     &     1.07823e+00, 1.08114e+00, 1.08442e+00, 1.08797e+00, 
     &     1.09186e+00, 1.09600e+00, 1.10044e+00, 1.10511e+00, 
     &     1.11005e+00, 1.11519e+00, 1.12057e+00, 1.12613e+00, 
     &     1.13189e+00, 1.13780e+00, 1.14389e+00, 1.15010e+00, 
     &     1.15645e+00, 1.16291e+00, 1.16947e+00, 1.17611e+00, 
     &     1.18283e+00, 1.18960e+00, 1.19641e+00, 1.20325e+00, 
     &     1.21009e+00, 1.21693e+00, 1.22378e+00, 1.23062e+00, 
     &     1.23746e+00, 1.24430e+00, 1.25114e+00, 1.25798e+00, 
     &     1.26483e+00, 1.27167e+00, 1.27851e+00, 1.28535e+00, 
     &     1.29219e+00, 1.29903e+00, 1.30588e+00, 1.31272e+00, 
     &     1.31956e+00, 1.32640e+00, 1.33324e+00, 1.34008e+00, 
     &     1.34692e+00, 1.35377e+00, 1.36061e+00, 1.36745e+00, 
     &     1.37429e+00, 1.38113e+00, 1.38797e+00, 1.39482e+00, 
     &     1.40166e+00, 1.40850e+00, 1.41534e+00, 1.42218e+00, 
     &     1.42902e+00, 1.43587e+00, 1.44271e+00, 1.44955e+00, 
     &     1.45639e+00, 1.46323e+00, 1.47007e+00, 1.47692e+00, 
     &     1.48376e+00, 1.49060e+00, 1.49744e+00 /
#endif /* Winds according to Meynet & Meader 00 */

#endif /* (_FANCY_REMNANT_<2) */


#if (_FANCY_REMNANT_==2)
c
c     KB05: Program remnant.c given by Chris Belczynski
c     see Belczynski et al 2002, ApJ 572, 407;  Belczynski et al 2004, ApJ 611, 1068
c
#if (_METALLICITY_==0) /* Z=0.0001 */
      data lM_remn /
     &     -1.02228e+00, -1.01228e+00, -1.00228e+00, -9.92276e-01, 
     &     -9.82276e-01, -9.72276e-01, -9.62276e-01, -9.52276e-01, 
     &     -9.42276e-01, -9.32276e-01, -9.22276e-01, -9.12276e-01, 
     &     -9.02276e-01, -8.92276e-01, -8.82276e-01, -8.72276e-01, 
     &     -8.62276e-01, -8.52276e-01, -8.42276e-01, -8.32276e-01, 
     &     -8.22276e-01, -8.12276e-01, -8.02276e-01, -7.92276e-01, 
     &     -7.82276e-01, -7.72276e-01, -7.62276e-01, -7.52276e-01, 
     &     -7.42276e-01, -7.32276e-01, -7.22276e-01, -7.12276e-01, 
     &     -7.02276e-01, -6.92276e-01, -6.82276e-01, -6.72276e-01, 
     &     -6.62276e-01, -6.52276e-01, -6.42276e-01, -6.32276e-01, 
     &     -6.22276e-01, -6.12276e-01, -6.02276e-01, -5.92276e-01, 
     &     -5.82276e-01, -5.72276e-01, -5.62276e-01, -5.52276e-01, 
     &     -5.42276e-01, -5.32276e-01, -5.22276e-01, -5.12276e-01, 
     &     -5.02276e-01, -4.92276e-01, -4.82276e-01, -4.72276e-01, 
     &     -4.62276e-01, -4.52276e-01, -4.42276e-01, -4.32276e-01, 
     &     -4.22276e-01, -4.12276e-01, -4.02276e-01, -3.92276e-01, 
     &     -3.82276e-01, -3.72276e-01, -3.62276e-01, -3.54080e-01, 
     &     -3.49945e-01, -3.45809e-01, -3.41673e-01, -3.37537e-01, 
     &     -3.33401e-01, -3.29265e-01, -3.25129e-01, -3.20993e-01, 
     &     -3.16857e-01, -3.12721e-01, -3.08585e-01, -3.04450e-01, 
     &     -3.00314e-01, -2.96178e-01, -2.92042e-01, -2.87906e-01, 
     &     -2.83770e-01, -2.79634e-01, -2.75498e-01, -2.71362e-01, 
     &     -2.67226e-01, -2.63090e-01, -2.58955e-01, -2.54819e-01, 
     &     -2.50683e-01, -2.46547e-01, -2.42411e-01, -2.38275e-01, 
     &     -2.34139e-01, -2.30003e-01, -2.25867e-01, -2.21731e-01, 
     &     -2.17595e-01, -2.13354e-01, -2.09112e-01, -2.04870e-01, 
     &     -2.00628e-01, -1.95697e-01, -1.90654e-01, -1.85611e-01, 
     &     -1.80526e-01, -1.74960e-01, -1.69393e-01, -1.63827e-01, 
     &     -1.58183e-01, -1.52489e-01, -1.46794e-01, -1.41020e-01, 
     &     -1.35118e-01, -1.29216e-01, -1.23693e-01, -1.18760e-01, 
     &     -1.13827e-01, -1.09596e-01, -1.05856e-01, -1.02116e-01, 
     &     -9.82466e-02, -9.43716e-02, -9.04268e-02, -8.64041e-02, 
     &     -8.23380e-02, -7.79670e-02, -7.35961e-02, -6.96338e-02, 
     &     -6.57185e-02, -6.12079e-02, -5.65274e-02, -5.04618e-02, 
     &     -4.39532e-02, -3.67450e-02, -2.93907e-02, -2.24349e-02, 
     &     -1.59620e-02, -1.13171e-02, -7.83601e-03, -5.50627e-03, 
     &     -3.59479e-03, -1.82753e-03, -2.57160e-04, 1.00875e-03, 
     &     2.11557e-03, 3.06643e-03, 3.51281e-03, 3.43350e-03, 
     &     2.71867e-03, 1.07712e-03, -1.69587e-03, -5.29581e-03, 
     &     -5.72727e-03, -1.50443e-05, 5.72737e-03, 1.16322e-02, 
     &     1.76137e-02, 2.36791e-02, 2.97877e-02, 3.60239e-02, 
     &     4.23366e-02, 4.87512e-02, 5.52210e-02, 6.07419e-02, 
     &     6.64399e-02, 7.30415e-02, 7.98237e-02, 8.66852e-02, 
     &     9.36236e-02, 1.00641e-01, 1.07730e-01, 1.14908e-01, 
     &     1.22157e-01, 1.29479e-01, 1.36883e-01, 1.44361e-01, 
     &     1.51913e-01, 1.32414e-01, 1.13965e-01, 1.13948e-01, 
     &     1.14020e-01, 1.15136e-01, 1.17609e-01, 1.20171e-01, 
     &     1.22832e-01, 1.25596e-01, 1.28367e-01, 1.31147e-01, 
     &     1.34023e-01, 1.36999e-01, 1.39831e-01, 1.43013e-01, 
     &     1.46051e-01, 1.49195e-01, 1.52448e-01, 1.55613e-01, 
     &     1.59014e-01, 1.62331e-01, 1.65761e-01, 1.70774e-01, 
     &     1.77481e-01, 1.84553e-01, 1.91538e-01, 1.98810e-01, 
     &     2.05764e-01, 2.13212e-01, 2.20327e-01, 2.27598e-01, 
     &     2.35036e-01, 2.42652e-01, 2.50430e-01, 2.57775e-01, 
     &     2.65268e-01, 2.73001e-01, 2.80800e-01, 2.88848e-01, 
     &     2.96744e-01, 3.04690e-01, 3.14895e-01, 3.27363e-01, 
     &     3.40204e-01, 3.60516e-01, 5.56571e-01, 7.06622e-01, 
     &     8.27579e-01, 9.30405e-01, 1.02061e+00, 1.07519e+00, 
     &     1.11689e+00, 1.15725e+00, 1.19643e+00, 1.23454e+00, 
     &     1.27168e+00, 1.30793e+00, 1.33784e+00, 1.34608e+00, 
     &     1.35424e+00, 1.36232e+00, 1.37031e+00, 1.37820e+00, 
     &     1.38599e+00, 1.39367e+00, 1.40124e+00, 1.40868e+00, 
     &     1.41598e+00, 1.42315e+00, 1.43016e+00, 1.11217e+00, 
     &     1.12327e+00, 1.13447e+00, 1.14552e+00, 1.15663e+00, 
     &     1.10021e+00, 1.05857e+00, 1.05522e+00, 1.05596e+00, 
     &     1.05821e+00, 1.06021e+00, 1.06231e+00, 1.06428e+00, 
     &     1.06714e+00, 1.06908e+00, 1.07092e+00, 1.07277e+00, 
     &     1.07498e+00, 1.07636e+00, 1.07755e+00, 1.07835e+00, 
     &     1.07876e+00, 1.07860e+00, 1.07732e+00, 1.07330e+00, 
     &     1.06136e+00, 1.04801e+00, 1.04378e+00, 1.04416e+00, 
     &     1.04686e+00, 1.05128e+00, 1.05554e+00, 1.06038e+00, 
     &     1.06554e+00, 1.06982e+00, 1.07512e+00, 1.07933e+00, 
     &     1.08423e+00, 1.08857e+00, 1.09358e+00, 1.09848e+00, 
     &     1.10288e+00, 1.10786e+00, 1.11217e+00, 1.11712e+00, 
     &     1.12141e+00, 1.12640e+00, 1.13062e+00, 1.13561e+00, 
     &     1.13980e+00, 1.14482e+00, 1.14900e+00, 1.15413e+00, 
     &     1.15827e+00, 1.16312e+00, 1.16763e+00, 1.17251e+00, 
     &     1.17668e+00, 1.18159e+00, 1.18578e+00, 1.19071e+00, 
     &     1.19493e+00, 1.19988e+00, 1.20396e+00, 1.20804e+00, 
     &     1.21302e+00, 1.21727e+00, 1.22204e+00, 1.22613e+00, 
     &     1.23022e+00, 1.23431e+00, 1.23840e+00, 1.24250e+00, 
     &     1.24659e+00, 1.25068e+00, 1.25477e+00, 1.25886e+00, 
     &     1.26295e+00, 1.26704e+00, 1.27113e+00, 1.27523e+00, 
     &     1.27932e+00, 1.28341e+00, 1.28750e+00, 1.29159e+00, 
     &     1.29568e+00, 1.29977e+00, 1.30386e+00, 1.30796e+00, 
     &     1.31205e+00, 1.31614e+00, 1.32023e+00, 1.32432e+00, 
     &     1.32841e+00, 1.33250e+00, 1.33659e+00, 1.34069e+00, 
     &     1.34478e+00, 1.34887e+00, 1.35296e+00 /

#elif (_METALLICITY_==1) /* Z=0.02 */

      data lM_remn /
     &     -1.02228e+00, -1.01228e+00, -1.00228e+00, -9.92276e-01, 
     &     -9.82276e-01, -9.72276e-01, -9.62276e-01, -9.52276e-01, 
     &     -9.42276e-01, -9.32276e-01, -9.22276e-01, -9.12276e-01, 
     &     -9.02276e-01, -8.92276e-01, -8.82276e-01, -8.72276e-01, 
     &     -8.62276e-01, -8.52276e-01, -8.42276e-01, -8.32276e-01, 
     &     -8.22276e-01, -8.12276e-01, -8.02276e-01, -7.92276e-01, 
     &     -7.82276e-01, -7.72276e-01, -7.62276e-01, -7.52276e-01, 
     &     -7.42276e-01, -7.32276e-01, -7.22276e-01, -7.12276e-01, 
     &     -7.02276e-01, -6.92276e-01, -6.82276e-01, -6.72276e-01, 
     &     -6.62276e-01, -6.52276e-01, -6.42276e-01, -6.32276e-01, 
     &     -6.22276e-01, -6.12276e-01, -6.02276e-01, -5.92276e-01, 
     &     -5.82276e-01, -5.72276e-01, -5.62276e-01, -5.52276e-01, 
     &     -5.42276e-01, -5.32276e-01, -5.22276e-01, -5.12276e-01, 
     &     -5.02276e-01, -4.92276e-01, -4.82276e-01, -4.72276e-01, 
     &     -4.62276e-01, -4.52276e-01, -4.42276e-01, -4.32276e-01, 
     &     -4.22276e-01, -4.12276e-01, -4.02276e-01, -3.92276e-01, 
     &     -3.82276e-01, -3.72276e-01, -3.62276e-01, -3.52276e-01, 
     &     -3.42276e-01, -3.32276e-01, -3.22276e-01, -3.12276e-01, 
     &     -3.10449e-01, -3.09515e-01, -3.08581e-01, -3.07646e-01, 
     &     -3.06712e-01, -3.05778e-01, -3.04844e-01, -3.03909e-01, 
     &     -3.02975e-01, -3.02041e-01, -3.01106e-01, -3.00172e-01, 
     &     -2.99238e-01, -2.98304e-01, -2.97369e-01, -2.96435e-01, 
     &     -2.95501e-01, -2.94567e-01, -2.93632e-01, -2.92698e-01, 
     &     -2.91764e-01, -2.90830e-01, -2.89895e-01, -2.88961e-01, 
     &     -2.88027e-01, -2.87093e-01, -2.86158e-01, -2.85224e-01, 
     &     -2.84290e-01, -2.83102e-01, -2.81915e-01, -2.80727e-01, 
     &     -2.79540e-01, -2.76879e-01, -2.73979e-01, -2.71080e-01, 
     &     -2.68164e-01, -2.65071e-01, -2.61977e-01, -2.58883e-01, 
     &     -2.55788e-01, -2.52693e-01, -2.49598e-01, -2.46486e-01, 
     &     -2.43347e-01, -2.40208e-01, -2.37042e-01, -2.33835e-01, 
     &     -2.30628e-01, -2.27337e-01, -2.23989e-01, -2.20640e-01, 
     &     -2.17099e-01, -2.13549e-01, -2.09951e-01, -2.06300e-01, 
     &     -2.02662e-01, -1.99108e-01, -1.95554e-01, -1.94603e-01, 
     &     -1.93951e-01, -1.92821e-01, -1.91555e-01, -1.89761e-01, 
     &     -1.87798e-01, -1.81993e-01, -1.75387e-01, -1.67065e-01, 
     &     -1.59573e-01, -1.55422e-01, -1.51240e-01, -1.47027e-01, 
     &     -1.42729e-01, -1.38321e-01, -1.33668e-01, -1.28886e-01, 
     &     -1.24050e-01, -1.19181e-01, -1.14269e-01, -1.09325e-01, 
     &     -1.04343e-01, -9.92962e-02, -9.41786e-02, -8.92743e-02, 
     &     -8.44016e-02, -7.92879e-02, -7.40299e-02, -6.86094e-02, 
     &     -6.30483e-02, -5.73498e-02, -5.15184e-02, -4.55552e-02, 
     &     -3.94612e-02, -3.32413e-02, -2.68983e-02, -2.04316e-02, 
     &     -1.38497e-02, -7.16791e-03, -3.98986e-04, 6.47754e-03, 
     &     1.34506e-02, 2.05409e-02, 2.77386e-02, 3.50628e-02, 
     &     4.24825e-02, 4.99903e-02, 5.76116e-02, 6.53307e-02, 
     &     7.31343e-02, 8.10480e-02, 8.90539e-02, 9.71608e-02, 
     &     1.05346e-01, 1.13514e-01, 1.21613e-01, 1.29796e-01, 
     &     1.38093e-01, 1.46483e-01, 1.54972e-01, 1.13962e-01, 
     &     1.14023e-01, 1.15091e-01, 1.17962e-01, 1.20731e-01, 
     &     1.23816e-01, 1.26830e-01, 1.29992e-01, 1.33278e-01, 
     &     1.36448e-01, 1.39910e-01, 1.43419e-01, 1.46936e-01, 
     &     1.50499e-01, 1.54181e-01, 1.57990e-01, 1.61640e-01, 
     &     1.65638e-01, 1.71332e-01, 1.79117e-01, 1.86579e-01, 
     &     1.94623e-01, 2.02446e-01, 2.10908e-01, 2.18891e-01, 
     &     2.27019e-01, 2.35298e-01, 2.43527e-01, 2.51876e-01, 
     &     2.60349e-01, 2.68598e-01, 2.77260e-01, 2.85726e-01, 
     &     2.94741e-01, 3.03108e-01, 3.12143e-01, 3.20897e-01, 
     &     3.29579e-01, 3.39382e-01, 3.53397e-01, 4.39621e-01, 
     &     5.07481e-01, 5.76216e-01, 6.41947e-01, 6.96423e-01, 
     &     7.25327e-01, 7.55393e-01, 7.66843e-01, 7.30654e-01, 
     &     7.41125e-01, 7.52487e-01, 7.63968e-01, 7.74043e-01, 
     &     7.86774e-01, 8.00226e-01, 8.12990e-01, 8.24481e-01, 
     &     8.37894e-01, 8.51481e-01, 8.64943e-01, 8.75762e-01, 
     &     8.88546e-01, 9.01208e-01, 9.13394e-01, 9.25377e-01, 
     &     9.35458e-01, 9.46888e-01, 9.58014e-01, 9.69072e-01, 
     &     9.79765e-01, 9.90217e-01, 1.00032e+00, 1.00350e+00, 
     &     1.00714e+00, 1.01009e+00, 1.01306e+00, 1.01590e+00, 
     &     1.01930e+00, 1.00406e+00, 9.32626e-01, 8.61699e-01, 
     &     7.98506e-01, 8.11284e-01, 8.23033e-01, 8.37153e-01, 
     &     8.53103e-01, 8.68561e-01, 8.82060e-01, 8.97013e-01, 
     &     9.11184e-01, 9.23568e-01, 9.34776e-01, 9.47402e-01, 
     &     9.59755e-01, 9.71640e-01, 9.83181e-01, 9.93368e-01, 
     &     1.00137e+00, 1.00516e+00, 1.00828e+00, 1.01096e+00, 
     &     1.01384e+00, 1.01668e+00, 1.01970e+00, 1.02239e+00, 
     &     1.02503e+00, 1.02781e+00, 1.03030e+00, 1.03241e+00, 
     &     1.03537e+00, 1.03768e+00, 1.03957e+00, 1.04236e+00, 
     &     1.04412e+00, 1.04620e+00, 1.04786e+00, 1.05045e+00, 
     &     1.05198e+00, 1.05356e+00, 1.05530e+00, 1.05736e+00, 
     &     1.05904e+00, 1.06034e+00, 1.06172e+00, 1.06319e+00, 
     &     1.06437e+00, 1.06627e+00, 1.06768e+00, 1.06904e+00, 
     &     1.07040e+00, 1.07177e+00, 1.07313e+00, 1.07449e+00, 
     &     1.07586e+00, 1.07722e+00, 1.07858e+00, 1.07994e+00, 
     &     1.08131e+00, 1.08267e+00, 1.08403e+00, 1.08540e+00, 
     &     1.08676e+00, 1.08812e+00, 1.08949e+00, 1.09085e+00, 
     &     1.09221e+00, 1.09358e+00, 1.09494e+00, 1.09630e+00, 
     &     1.09767e+00, 1.09903e+00, 1.10039e+00, 1.10175e+00, 
     &     1.10312e+00, 1.10448e+00, 1.10584e+00, 1.10721e+00, 
     &     1.10857e+00, 1.10993e+00, 1.11130e+00 /
#endif /* metallicity choice */
#endif /* (_FANCY_REMNANT_==2) */

c
c......................................................................
c
      lm = log10(M_zams)
      if (lm.LT.lMmin) then ! would only occur for T>>T_Hubble !
         M_remnant=0.95d0*M_zams
         return 
      end if

      !call locate (lM_zams,DimRelMremn,lm,i1)
      i1 = INT((lm-lMmin)/dlM) + 1
      i2 = i1+1
      if (i1.GE.DimRelMremn) then ! linear extension of the relation (for very massive stars...)
         i1=DimRelMremn-5
         i2=DimRelMremn
      end if
      lm1 = lmmin+(i1-1)*dlm 
      lm2 = lmmin+(i2-1)*dlm 
      lmr = lM_remn(i1) + (lm-lm1)*
     $        (lM_remn(i2)-lM_remn(i1))/(lm2-lm1)
      M_remnant = min(0.95d0*M_zams,10**lmr)
c
      end
c
#endif  /* (_FANCY_REMNANT_>0) */

#if (_NATKICKS_>0)
c----------------------------------------------------------------------
      subroutine Get_NatalKick(M_zams,M_remn,iT_remn,Vkick_kms)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_EvolStell_common.f'
      include 'EvolAmas_Etoiles_param.f'
c
c arguments :
c ^^^^^^^^^^^
                                ! Inputs:
      double precision M_zams 
      double precision M_remn 
      integer*1 iT_remn
                                ! Output:
      double precision Vkick_kms ! kick velocity in km/s
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision sigma_kick
c
c functions :
c ^^^^^^^^^^^
      double precision random, RandomGauss
c......................................................................
c
c---- Only neutron stars and black holes get kicks
c
      if (iT_remn.NE.TypeEtoile_NS .AND. iT_remn.NE.TypeEtoile_BH)
     $     then
         Vkick_kms =0.0d0
      else if (iT_remn.EQ.TypeEtoile_BH .AND.
     $        M_zams.LT.Mzams_max_fallbackBH) then
         Vkick_kms =0.0d0
      else
c
c---- We assume a bimodal Maxwellian distribution for the natal kicks
c
         if (random().LT.NatKick_Prob_1) then
            sigma_kick = NatKick_Sig1D_1
         else
            sigma_kick = NatKick_Sig1D_2
         end if
         Vkick_kms = sigma_kick * ! Is there a more efficient way?
     $        sqrt(RandomGauss()**2+RandomGauss()**2+RandomGauss()**2)
                                ! scale kick with 1/mass
         Vkick_kms = Vkick_kms * (M_NS/M_remn)
      end if
c
      end
c
#endif /* (_NATKICKS_>0) */
#endif /* (_EVOL_STELL_==1) */ /* Simple stellar evol (MS-->remnant) */

c----------------------------------------------------------------------
      double precision function T_MS(M_ZAMS) 
c----------------------------------------------------------------------
c     Fournit le tps passe sur la seq. princ. en annees
c
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      double precision M_ZAMS 
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision M_TMS,fHe_core
c......................................................................
c
      call GetMSdata(M_ZAMS,T_MS,M_TMS,fHe_core)
c
      end
c

c-------------------------------------------------------------------------------
      subroutine GetMSdata(M_ZAMS,T_MS,M_TMS,MfracCoreHe)
c-------------------------------------------------------------------------------
c     gives quantities describing the MS life of a star of ZAMS mass M_ZAMS
c     M_ZAMS -- initial MS mass
c     M_TMS -- terminal MS mass
c     T_MS -- MS lifetime in years
c     MfracCoreHe -- He core mass at terminal MS / terminal MS mass
c
c     Sources: Hurley et al. 00, Belczynski et al. 02
c
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      double precision M_ZAMS ! input 
      double precision T_MS,M_TMS,MfracCoreHe ! output
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer DimMSdata
      parameter ( DimMSdata  = 50 )
      double precision logM_ZAMS(DimMSdata)
      double precision logT_MS(DimMSdata)
      double precision logM_TMS(DimMSdata)
      double precision fHe_core(DimMSdata)
      save logM_ZAMS, logT_MS, logM_TMS, fHe_core

      data logM_ZAMS /
     &     -1.00000e+00, -9.38776e-01, -8.77551e-01, -8.16327e-01,
     &     -7.55102e-01, -6.93878e-01, -6.32653e-01, -5.71429e-01,
     &     -5.10204e-01, -4.48980e-01, -3.87755e-01, -3.26531e-01,
     &     -2.65306e-01, -2.04082e-01, -1.42857e-01, -8.16327e-02,
     &     -2.04082e-02, 4.08163e-02, 1.02041e-01, 1.63265e-01,
     &     2.24490e-01, 2.85714e-01, 3.46939e-01, 4.08163e-01,
     &     4.69388e-01, 5.30612e-01, 5.91837e-01, 6.53061e-01,
     &     7.14286e-01, 7.75510e-01, 8.36735e-01, 8.97959e-01,
     &     9.59184e-01, 1.02041e+00, 1.08163e+00, 1.14286e+00,
     &     1.20408e+00, 1.26531e+00, 1.32653e+00, 1.38776e+00,
     &     1.44898e+00, 1.51020e+00, 1.57143e+00, 1.63265e+00,
     &     1.69388e+00, 1.75510e+00, 1.81633e+00, 1.87755e+00,
     &     1.93878e+00, 2.00000e+00 /

#if (_METALLICITY_==0)
c
c     Z=0.0001
c
      data logT_MS /
     &     1.23377e+01, 1.22153e+01, 1.20928e+01, 1.19703e+01,
     &     1.18479e+01, 1.17253e+01, 1.16011e+01, 1.14770e+01,
     &     1.13510e+01, 1.12174e+01, 1.10795e+01, 1.09206e+01,
     &     1.07403e+01, 1.05362e+01, 1.03090e+01, 1.00740e+01,
     &     9.84619e+00, 9.63193e+00, 9.43470e+00, 9.24993e+00,
     &     9.07536e+00, 8.90889e+00, 8.74904e+00, 8.59518e+00,
     &     8.44729e+00, 8.30512e+00, 8.16859e+00, 8.03773e+00,
     &     7.91260e+00, 7.79315e+00, 7.67936e+00, 7.57117e+00,
     &     7.46853e+00, 7.37126e+00, 7.27939e+00, 7.19292e+00,
     &     7.11188e+00, 7.03634e+00, 6.96639e+00, 6.90208e+00,
     &     6.84347e+00, 6.79052e+00, 6.74315e+00, 6.70119e+00,
     &     6.66439e+00, 6.63242e+00, 6.60490e+00, 6.58141e+00,
     &     6.56151e+00, 6.54476e+00 /
      data logM_TMS /
     &     -1.00000e+00, -9.38776e-01, -8.77551e-01, -8.16327e-01,
     &     -7.55102e-01, -6.93878e-01, -6.32653e-01, -5.71429e-01,
     &     -5.10204e-01, -4.48980e-01, -3.87755e-01, -3.26531e-01,
     &     -2.65306e-01, -2.04082e-01, -1.42857e-01, -8.16327e-02,
     &     -2.04082e-02, 4.08163e-02, 1.02041e-01, 1.63265e-01,
     &     2.24490e-01, 2.85714e-01, 3.46939e-01, 4.08163e-01,
     &     4.69388e-01, 5.30612e-01, 5.91837e-01, 6.53061e-01,
     &     7.14286e-01, 7.75510e-01, 8.36698e-01, 8.97832e-01,
     &     9.58954e-01, 1.02011e+00, 1.08125e+00, 1.14236e+00,
     &     1.20345e+00, 1.26449e+00, 1.32550e+00, 1.38645e+00,
     &     1.44734e+00, 1.50814e+00, 1.56884e+00, 1.62943e+00,
     &     1.68986e+00, 1.75011e+00, 1.81014e+00, 1.86991e+00,
     &     1.92939e+00, 1.98854e+00 /
      data fHe_core /
     &     2.06170e-01, 2.09495e-01, 2.12821e-01, 2.16146e-01,
     &     2.19471e-01, 2.22696e-01, 2.24817e-01, 2.26938e-01,
     &     2.28166e-01, 2.25975e-01, 2.22813e-01, 2.14788e-01,
     &     2.04762e-01, 1.93253e-01, 1.81995e-01, 1.72996e-01,
     &     1.67349e-01, 1.64738e-01, 1.63398e-01, 1.60561e-01,
     &     1.55316e-01, 1.48402e-01, 1.43353e-01, 1.42460e-01,
     &     1.44940e-01, 1.49737e-01, 1.56023e-01, 1.63276e-01,
     &     1.71570e-01, 1.79912e-01, 1.88848e-01, 1.98316e-01,
     &     2.08304e-01, 2.18824e-01, 2.29884e-01, 2.41504e-01,
     &     2.53707e-01, 2.66521e-01, 2.79973e-01, 2.94091e-01,
     &     3.08905e-01, 3.24445e-01, 3.40738e-01, 3.57815e-01,
     &     3.75701e-01, 3.94424e-01, 4.04472e-01, 4.11402e-01,
     &     4.18415e-01, 4.25508e-01 /

#elif (_METALLICITY_==1)
c
c     Z=0.02 (solar metallicity)
c
      data logT_MS /
     &     1.25630e+01, 1.24406e+01, 1.23181e+01, 1.21957e+01,
     &     1.20732e+01, 1.19507e+01, 1.18273e+01, 1.17039e+01,
     &     1.15793e+01, 1.14500e+01, 1.13178e+01, 1.11710e+01,
     &     1.10058e+01, 1.08176e+01, 1.06005e+01, 1.03629e+01,
     &     1.01213e+01, 9.88779e+00, 9.67039e+00, 9.47533e+00,
     &     9.29047e+00, 9.11068e+00, 8.93461e+00, 8.76381e+00,
     &     8.59766e+00, 8.43603e+00, 8.27914e+00, 8.12733e+00,
     &     7.98106e+00, 7.84081e+00, 7.70706e+00, 7.58076e+00,
     &     7.46316e+00, 7.35215e+00, 7.24884e+00, 7.15363e+00,
     &     7.06658e+00, 6.98761e+00, 6.91660e+00, 6.85331e+00,
     &     6.79742e+00, 6.74852e+00, 6.70613e+00, 6.66972e+00,
     &     6.63871e+00, 6.61252e+00, 6.59057e+00, 6.57232e+00,
     &     6.55723e+00, 6.54485e+00 /
      data logM_TMS /
     &     -1.00000e+00, -9.38776e-01, -8.77551e-01, -8.16327e-01,
     &     -7.55102e-01, -6.93878e-01, -6.32653e-01, -5.71429e-01,
     &     -5.10204e-01, -4.48980e-01, -3.87755e-01, -3.26531e-01,
     &     -2.65306e-01, -2.04082e-01, -1.42857e-01, -8.16327e-02,
     &     -2.04082e-02, 4.08163e-02, 1.02041e-01, 1.63265e-01,
     &     2.24490e-01, 2.85714e-01, 3.46939e-01, 4.08163e-01,
     &     4.69388e-01, 5.30612e-01, 5.91837e-01, 6.53061e-01,
     &     7.14286e-01, 7.75510e-01, 8.36712e-01, 8.96113e-01,
     &     9.55178e-01, 1.01500e+00, 1.07458e+00, 1.13370e+00,
     &     1.19220e+00, 1.25000e+00, 1.30689e+00, 1.36267e+00,
     &     1.41706e+00, 1.46972e+00, 1.52029e+00, 1.56831e+00,
     &     1.61334e+00, 1.65493e+00, 1.69268e+00, 1.72631e+00,
     &     1.75571e+00, 1.78082e+00 /
      data fHe_core /
     &     2.06630e-01, 2.00334e-01, 1.94039e-01, 1.87743e-01,
     &     1.81447e-01, 1.75208e-01, 1.69598e-01, 1.63989e-01,
     &     1.58416e-01, 1.52986e-01, 1.47517e-01, 1.41852e-01,
     &     1.35881e-01, 1.29579e-01, 1.23054e-01, 1.17010e-01,
     &     1.12421e-01, 1.09674e-01, 1.08340e-01, 1.06912e-01,
     &     1.04719e-01, 1.01955e-01, 1.15789e-01, 1.28233e-01,
     &     1.37738e-01, 1.46150e-01, 1.54249e-01, 1.62402e-01,
     &     1.70797e-01, 1.79532e-01, 1.88664e-01, 1.97949e-01,
     &     2.07625e-01, 2.17901e-01, 2.28640e-01, 2.39819e-01,
     &     2.51414e-01, 2.63420e-01, 2.75799e-01, 2.88500e-01,
     &     3.01447e-01, 3.14538e-01, 3.27641e-01, 3.40592e-01,
     &     3.53199e-01, 3.65256e-01, 3.76557e-01, 3.86920e-01,
     &     3.96212e-01, 4.04324e-01 /
#else
      ERROR: _METALLICITY_ not defined
#endif
      double precision lM
      integer j1,j2
c......................................................................
c
      lM=log10(M_ZAMS)
      
      if (lM.LE.LogM_ZAMS(1)) then ! extrapolation lineaire
         j1=1
         j2=2
      else if (lM.GE.LogM_ZAMS(DimMSdata)) then ! extrapolation lineaire
         j1=DimMSdata-1
         j2=DimMSdata
      else                      ! interpolation lineaire
         call locate(logM_ZAMS,DimMSdata,lM,j1)
         j2=j1+1
      end if
      
      T_MS = logT_MS(j1) +
     $     (logT_MS(j2)-logT_MS(j1))*(lM-logM_ZAMS(j1))/
     $     (logM_ZAMS(j2)-logM_ZAMS(j1))
      T_MS = 10**T_MS

      M_TMS = logM_TMS(j1) +
     $     (logM_TMS(j2)-logM_TMS(j1))*(lM-logM_ZAMS(j1))/
     $     (logM_ZAMS(j2)-logM_ZAMS(j1))
      M_TMS = min(M_ZAMS,10**M_TMS)

      if (lM.LE.LogM_ZAMS(1)) then 
         MfracCoreHe = fHe_Core(1)
      else
         MfracCoreHe = fHe_Core(j1) +
     $        (fHe_Core(j2)-fHe_Core(j1))*(lM-logM_ZAMS(j1))/
     $        (logM_ZAMS(j2)-logM_ZAMS(j1))
         MfracCoreHe = min(MfracCoreHe,0.9d0)
      endif

c
      end
c
c----------------------------------------------------------------------
      subroutine Init_StEvol_ASC(iFich)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100)
 100  format(
     $     '# Stellar evolution events',/,
     $     '# 1:iPas_Evol 2:iSE 3:Age_in_yrs 4:R ',
     $     '5:iType_before 6:M_before 7:V2_before ',
     $     '8:iType_after  9:M_after 10:V2_after')
c
      end
c

c----------------------------------------------------------------------
      subroutine Init_NatKicks_ASC(iFich)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100)
 100  format(
     $     '# Additional info for natal kicks events')
      call OutputFileHeader(iFich, 'iPas_Evol iSE Vkick_kms Vkick ' // 
     $     'Vkick_rad Vkick_tg T_ini J_ini T_fin, J_fin, V_esc')
c
      end
c

c----------------------------------------------------------------------
      block data EvolStellDefValues
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_EvolStell_param.f'
      include 'EvolAmas_EvolStell_common.f'
      data
     $     Mzams_max_WD, Mzams_max_NS,
     $     Mzams_max_fallbackBH,
     $     M_max_WD, M_max_NS, Mremn_min_BH,
     $     M_WD, M_NS, M_BH,
     $     NatKick_Sig1D_1, NatKick_Sig1D_2,
     $     NatKick_Prob_1
     $     /Mzams_max_WD_def, Mzams_max_NS_def,
     $     Mzams_max_fallbackBH_def,
     $     M_max_WD_def, M_max_NS_def, Mremn_min_BH_def,
     $     M_WD_def, M_NS_def, M_BH_def,
     $     NatKick_Sig1D_1_def, NatKick_Sig1D_2_def,
     $     NatKick_Prob_1_def/
      data lInitEvolStell /.FALSE./
      end 
