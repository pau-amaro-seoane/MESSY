c======================================================================
c     routines liees au trou noir pour code d'evolution d'amas
c======================================================================

c----------------------------------------------------------------------
      double precision function Rmar_TN(M_stel,R_stel)
c----------------------------------------------------------------------
      include 'EvolAmas_TN_inc.f'
      include 'EvolAmas_Dimensions_common.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision M_stel,R_stel
c......................................................................
c
      Rmar_TN = Coef_ConvR*R_stel*(M_TN/(Coef_ConvM*M_stel))**0.333333d0
c
      end
c
c----------------------------------------------------------------------
      double precision function Rhoriz_eff_TN()
c----------------------------------------------------------------------
      include 'EvolAmas_TN_inc.f'
      include 'EvolAmas_common.f'
c......................................................................
c
c     Selon Miralda-Escude & Gould, un coef 8 serait plus approprie
c     pour l'horizon EFFECTIF
c     
      Rhoriz_eff_TN = 8.0d0*M_TN/VitLum2
c
      end
c
#if (_NO_CENTRAL_VMS_==0)
c----------------------------------------------------------------------
      double precision function R_VMS()
c----------------------------------------------------------------------
c     Radius of the central very massive star (aka VMS)
c
      include 'EvolAmas_TN_inc.f'
      include 'EvolAmas_Dimensions_common.f'
c
c function:
c ^^^^^^^^^
      double precision stellar_radius
c......................................................................
c
      R_VMS = Coef_ConvR*stellar_radius(id_CtrObj)
c
      end
c
#endif
#if (_NO_CROISS_TN_FORC_==0)
c----------------------------------------------------------------------
      double precision function M_TN_Croiss_Forc(Tps)
c----------------------------------------------------------------------
      include 'EvolAmas_TN_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision Tps
c
c var locale :
c ^^^^^^^^^^^^
      integer k_hunt
      save k_hunt
      data k_hunt /1/
c......................................................................
c
      if (iType_TN_CF.EQ.PARAM) then
         if (T_TN_CF.GT.0.0d0) then ! Sigurdsson et al. 1995, ApJ 446, 75
            if (Tps.LT.T_TN_CF) then
               M_TN_Croiss_Forc = M_TN_fin_CF * 
     $              (3.0d0-2.0d0*(Tps/T_TN_CF)) *(Tps/T_TN_CF)**2 
            else
               M_TN_Croiss_Forc = M_TN_fin_CF
c$$$c!#### extension
c$$$            M_TN_Croiss_Forc = M_TN_fin_CF*
c$$$     $           exp( log(M_TN_fin_CF_ext/M_TN_fin_CF)
c$$$     $           * (Tps-T_TN_CF)/(T_TN_CF_ext-T_TN_CF) )
c$$$c!####
            end if
         else                   ! Linear growth
            if (Tps.LT.-T_TN_CF) then
               M_TN_Croiss_Forc = -M_TN_fin_CF * (Tps/T_TN_CF)
            else
               M_TN_Croiss_Forc = M_TN_fin_CF
            end if
         end if
      else if (iType_TN_CF.EQ.EXPL) then
c$$$         call splint_h(T_expl_TN_CF,M_expl_TN_CF,Y2_expl_TN_CF,
c$$$     $        N_expl_TN_CF,k_hunt,Tps,M_TN_Croiss_Forc)
         call splint(T_expl_TN_CF,M_expl_TN_CF,Y2_expl_TN_CF,
     $        N_expl_TN_CF,Tps,M_TN_Croiss_Forc)
      else
         M_TN_Croiss_Forc = 0.0d0
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Lire_Param_TN_Croiss_Forc(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_TN_inc.f'
      include 'EvolAmas_common.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
c
c constantes :
c ^^^^^^^^^^^^
      double precision fT_TN_CF_def,fM_TN_fin_CF_def
      parameter (fT_TN_CF_def=0.2d0, fM_TN_fin_CF_def=0.1d0)
c
c var locales :
c ^^^^^^^^^^^^^
      character*80 Ligne
      character*64 Type
      double precision fT_TN_CF,fM_TN_fin_CF,x,y
c......................................................................
c
c---- Valeurs par defaut
c
      T_TN_CF = fT_TN_CF_def*UT_FP
      M_TN_fin_CF = fM_TN_fin_CF_def*UM_Ncorps
c
c---- Lecture du fichier
c
      Type='PARAM'
 1    continue
      read(iFich,'(A)',END=2) Ligne
      call EnleverBlancs(Ligne)
                                ! on saute les lignes blanches et les commentaires
      if (LEN(Ligne).EQ.0 .OR. Ligne(1:1).EQ.'#') then
	 goto 1
      end if

      if (Ligne(1:1).EQ.'@') then
         Type=Ligne(2:index(Ligne,' ')-1)
         read(iFich,'(A)') Ligne
      end if

      if (Type.EQ.'PARAM') then
         iType_TN_CF=PARAM
                                ! Modele parametrique simple
         do while (.TRUE.)
                                ! on lit les parametres presents
            if (Ligne.EQ.'T_TN_CF') then
               read(iFich,*) T_TN_CF
            else if (Ligne.EQ.'M_TN_fin_CF') then
               read(iFich,*) M_TN_fin_CF
            else if (Ligne.EQ.'fT_TN_CF') then
               read(iFich,*) fT_TN_CF
               T_TN_CF = fT_TN_CF*UT_FP
            else if (Ligne.EQ.'fM_TN_fin_CF') then
               read(iFich,*) fM_TN_fin_CF
               M_TN_fin_CF = fM_TN_fin_CF*UM_Ncorps
c$$$c!#### extension
c$$$            else if (Ligne.EQ.'T_TN_CF_ext') then
c$$$               read(iFich,*) T_TN_CF_ext
c$$$            else if (Ligne.EQ.'M_TN_fin_CF_ext') then
c$$$               read(iFich,*) M_TN_fin_CF_ext
c$$$c!####
               
            else
               write(0,100) Ligne 
 100           format(/,'!!! Le parametre ',A20,' n''a pas de sens !!!')
               call Meurt('Lire_Param_TN_Croiss_Forc')
            
            end if
            read(iFich,'(A)',END=3) Ligne
         end do

 3       continue
         write(0,110) T_TN_CF,M_TN_fin_CF
 110     format('> Parametres pour croissance forcee du TN',
     $        ' T_TN_CF,M_TN_fin_CF = ',2(E12.4,' '))
         write(0,115) T_TN_CF_ext,M_TN_fin_CF_ext
 115     format('> Parametres pour croissance forcee du TN',
     $        ' T_TN_CF_ext,M_TN_fin_CF_ext = ',2(E12.4,' '))

      else if (Type.EQ.'EXPL') then
         iType_TN_CF=EXPL
                                ! Specification explicite d'une loi T-->M_TN
         N_expl_TN_CF = 0
         do while (.TRUE.)
            read(iFich,*,END=9) x,y
            if (N_expl_TN_CF.LT.Nmax_expl_TN_CF) then
               if (N_expl_TN_CF.EQ.0 .OR.
     $              x.GT.T_expl_TN_CF(N_expl_TN_CF)) then
                  N_expl_TN_CF = N_expl_TN_CF + 1
                  T_expl_TN_CF(N_expl_TN_CF) = x
                  M_expl_TN_CF(N_expl_TN_CF) = y
               end if
            else
               write(0,*) '!!! Le nb max de points pour la relation ',
     $              'T-->M_TN est ',Nmax_expl_TN_CF,' !!!'
               call Meurt('Lire_Param_TN_Croiss_Forc')
            end if
         end do
 9       continue
         call spline(T_expl_TN_CF,M_expl_TN_CF,N_expl_TN_CF,
     $        1.0d30,1.0d30,Y2_expl_TN_CF)
         
         write(0,120) N_expl_TN_CF
 120     format('> Relation T-->M_TN lue (',I3,' elements)')

      end if

 2    continue
c     
      end
c 
#endif /* (_NO_CROISS_TN_FORC_==0) */
c----------------------------------------------------------------------
      subroutine read_TN(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_TN_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      read(iFich) M_TN
c
      end
c
c----------------------------------------------------------------------
      subroutine EcritTN_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_TN_inc.f'
      include 'fxdr.inc'
c
c argument :
c ^^^^^^^^^^
      integer iFich_XDR
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err, iType
c......................................................................
c
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichTN_XDR)
c
      call XDR_Ecrit_d (iFich_XDR,'M_TN',M_TN)
      call XDR_Ecrit_d (iFich_XDR,'M_ctr',M_ctr)
      call XDR_Ecrit_d (iFich_XDR,'BirthDate_CtrObj',BirthDate_CtrObj)
      iType = iType_CtrObj
      call XDR_Ecrit_i (iFich_XDR,'iType_CtrObj',iType)
c
      end
c
c----------------------------------------------------------------------
      subroutine LitTN_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_TN_inc.f'
      include 'fxdr.inc'
c
c argument :
c ^^^^^^^^^^
      integer iFich_XDR
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      character*20 NomType,NomVar
      character*80 EnTete
      integer idum, ixdr_err
c
c fonction :
c ^^^^^^^^^^
      integer ixdr_read_string
c......................................................................
c
      ixdr_err = ixdr_read_string(iFich_XDR,EnTete)
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichTN_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' de trou noir XDR valable !!!'
         call Meurt('LitTN_XDR')
      end if
      M_ctr=0.0d0
      M_TN=0.0d0
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.LT.0) goto 99

         if (NomVar.EQ.'M_TN') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,M_TN)
         else if (NomVar.EQ.'M_ctr') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,M_ctr)
         else if (NomVar.EQ.'BirthDate_CtrObj') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,BirthDate_CtrObj)
         else if (NomVar.EQ.'iType_CtrObj') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,idum)
            iType_CtrObj = idum
         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitTN_XDR')
         end if

      end do

 99   continue
      M_ctr=max(M_ctr,M_TN)
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_FichLC_ASC(iFich)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100)
 100     format(
     $     '# Donnees sur les etoiles dechirees ou ',
     $     'traversant l''horizon',/,
     $     '# Valeurs de |iFlagLC|: FULL=1,EMPTY=2,COLL=3,GW=11;',
     $     ' <0 si traversee directe',/,
     $     '# 1:iPas_Evol 2:iSE 3:R_SE 4:E_SE 5:J_SE 6:iTet_SE',
     $     ' 7:Met_SE 8:Theta_disp 9:Theta_dechir 10:dt_renc',
     $     ' 11:Porb_quad 12:Beta_dechir 13:iFlagLC 14:dE 15:Pot_Ctr',
     $     ' 16:Trel_renc 17:Trel_loc 18: Trel_orb 19:Peri 20:Apo'
     $        )
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_FichCaptGW_ASC(iFich)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100)
 100     format(
     $     '# Donnees complementaires sur les etoiles captees par ',
     $     'emission d''ondes grav.'/,
     $     '# 1:iPas_Evol 2:a 3:e 4:theta2 5:T_GW_sur_Porb',
     $     ' 6:Trel_sur_Porb'
     $        )
c
      end
c

c----------------------------------------------------------------------
      block data bd_init_TN
c----------------------------------------------------------------------
      include 'EvolAmas_TN_inc.f'
      include 'EvolAmas_Etoiles_param.f'
      
      data M_TN /0.0d0/
      data M_ctr /0.0d0/
      data BirthDate_CtrObj /0.0d0/
      data iType_CtrObj /TypeEtoile_BH/
c
      end
c----------------------------------------------------------------------
