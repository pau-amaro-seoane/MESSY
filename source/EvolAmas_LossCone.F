#define _TYPE_DECHIR_ 1
#define _TYPE_DETEC_LC_ 1
#define _ECRIT_INPUT_PA_ 0
#define _ECRIT_PA_ 0
c----------------------------------------------------------------------
      subroutine Dechirement_par_TN(TypeEtoile,Mini,beta,
     $     Mfin,dEorb)
c----------------------------------------------------------------------
      implicit none
      include 'Param_Amas.f'
      include 'EvolAmas_Etoiles_param.f'
c
c arguments :
c ^^^^^^^^^^^
      integer*1 TypeEtoile           ! Type d'etoile (MS,WD,NS,RG...)
      double precision Mini,beta  ! Masse stellaire initiale (en Msol), parametre de "profondeur" de la rencontre
      double precision Mfin,dEorb ! Masse stellaire subsistante, variation d'energie orbitale specifique
c......................................................................
c
c---- Pour l'instant, on traite tous les type sur le meme pied 
c     Seuls les rayon et masses sont pris en compte
c
c$$$      if (TypeEtoile.NE.TypeEtoile_MS) then
c$$$         write(0,*)'!!! On ne sait traiter que les dechirements',
c$$$     $        ' d''etoiles de sequ. princ. !!!'
c$$$         call Meurt('Dechirement_par_TN')
c$$$      end if

#if (_TYPE_DECHIR_==0 || _DECHIR_NULS_>1)
c
c---- On ne change rien !!!!!! (pour tests)
c
      Mfin = Mini
      dEorb = 0.0d0
#elif (_TYPE_DECHIR_==1)
c
c---- Destruction complete
c
      dEorb = 0.0d0
      if (beta.GT.0.8d0) then
         Mfin = 0.0d0
      else
         Mfin = Mini
      end if
#else
      write(0,*) '!!! Valeur invalide de _TYPE_DECHIR_'
      call Meurt('Dechirement_par_TN')
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_Proba_LC_Porb(Theta_lc,Theta_diff,Theta_V,Prob_LC)
c----------------------------------------------------------------------
c     Determine la probabilite de diffusion dans le loss-cone pendant
c     une periode orbitale. Methode basee sur le calcul de la surface
c     de la "lunule" de recouvrement entre le cone de diffusion (ouverture
c     Theta_diff) et le cone de perte (ouverture Theta_lc) avec une
c     separation angulaire Theta_V entre ces deux cones (direction de la
c     vitesse par-rapport a la direction du centre de l'amas)
c
      implicit none
      include 'CstesMath.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     Theta_lc,Theta_diff,Theta_V, ! parametres d'entree (voir ci-dessus)
     $     Prob_LC                      ! La probabilite de tomber dans le LC
c
c var locales :
c ^^^^^^^^^^^^^
      double precision t1,t2,t2t2t2,alpha,l,
     $     arg_sqrt
c......................................................................
c
c---- Reglons d'abord certains cas particuliers faciles
c
      if ((Theta_lc+Theta_diff).LT.Theta_V)
     $     then                 ! Pas de recouvrement
         Prob_LC = 0.0d0
      else if ((Theta_lc+Theta_V).LT.Theta_diff)
     $     then                 ! Cone de diff englobant cone de pertes 
                                ! C'est la limite "full loss-cone"
         Prob_LC = (Theta_lc/Theta_diff)**2
      else if ((Theta_diff+Theta_V).LT.Theta_lc)
     $     then                 ! Cone de pertes englobant cone de diff
                                ! C'est la limite "empty loss-cone"
         Prob_LC = 1.0d0
c
c---- Reste le cas d'un recouvrement partiel...
c
      else
         Prob_LC = 0.0d0
                                ! Demi-lunule du disque de perte
                                ! ------------------------------
         t1 = Theta_lc
         t2 = Theta_diff
         t2t2t2 = Theta_V**2+t1**2-t2**2
                                ! distance centre disque - corde d'intersection
         l = 0.5d0*t2t2t2/Theta_V ! et si Theta_V ~= 0 ?
                                ! demi-corde
#if (_SURVEILLE_SQRT_ > 0)
         arg_sqrt = t1**2-l**2
         if (arg_sqrt.LT.0.0d0) then
            write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
            call Meurt('Calc_Proba_LC_Porb')
         end if
         l = sqrt(arg_sqrt)
#else
         l = sqrt(t1**2-l**2) ! faire eventuellement un test
#endif
                                ! demi-angle d'ouverture de la lunule vue du centre du disque
         alpha = asin(l/t1) ! faire eventuellement un test
         if (t2t2t2.LT.0.0d0) then ! "grande" lunule
            alpha = rPi-alpha
         end if
         alpha = 2.0d0*alpha
         Prob_LC = Prob_LC + 0.5d0*t1**2*(alpha-sin(alpha))
                                ! Demi-lunule du disque de diffusion
                                ! ----------------------------------
         t1 = Theta_diff
         t2 = Theta_lc
         t2t2t2 = Theta_V**2+t1**2-t2**2                            
                                ! demi-angle d'ouverture de la lunule vue du centre du disque
         alpha = asin(l/t1) ! faire eventuellement un test
         if (t2t2t2.LT.0.0d0) then ! "grande" lunule
            alpha = rPi-alpha
         end if
         alpha = 2.0d0*alpha
         Prob_LC = Prob_LC + 0.5d0*t1**2*(alpha-sin(alpha))
                                ! On divise par la surface du disque de diffusion
                                ! pour obtenir la probabilite
         Prob_LC = Prob_LC/(rPi*Theta_diff**2)
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine ParcoursAlea_LC(Vpart_renc,Vrenc,Vlc,delta,n_pas,
     $     n_pas_adap,Vtg,V2,iIssue)
c----------------------------------------------------------------------
c     Version avec pas adaptatif pour plus d'efficacite
c
      implicit none
      include 'CstesMath.f'      
      include 'EvolAmas_param.f'
      include 'EvolAmas_Compteurs_common.f'
      include 'EvolAmas_GW_common.f'
      include 'EvolAmas_TN_common.f'
c
c arguments :
c ^^^^^^^^^^^
      ! ENTREE
      double precision
     $     Vpart_renc(3),       ! vitesse de la particule dans le referentiel de la rencontre
     $     Vrenc(3),            ! vitesse du ref. de la rencontre par-rap. a l'amas
     $     Vlc,                 ! vitesse tangentielle correspondant au LC (ref. de l'amas)
     $     delta                ! amplitude d'un pas elementaire
      integer n_pas             ! nb de pas (max) du parcours aleatoire
      ! SORTIE
      integer n_pas_adap        ! nb de "pas adaptes" realises
      double precision Vtg, V2  ! vitesse tangentielle et vitesse**2 dans referentiel de l'amas
                                ! lors du dernier pas (correspondant a chute dans le LC)
      integer iIssue            ! Indique le resultat du parcours :
                                !        0 ---> Etoile survit
                                !        1 ---> Etoile entre dans le LC
                                !        2 ---> Etoile captee sur orbite a forte 
                                !               emission d'ondes grav.
c
c constantes :
c ^^^^^^^^^^^^
      double precision minus
      parameter (minus=1.0d-20)
      double precision coef_secu, delta_max, fac_delta_adap
      parameter (coef_secu=0.5d0) ! definit la taille max du pas en fct de la
                                  ! proximite au LC
      parameter (delta_max=0.1d0*rPi) ! pas maximum
      parameter (fac_delta_adap=1.0d0) ! A mettre a 1.0d0 sauf pour experimentations!
      double precision cf_GW_1, cf_GW_2
      parameter (cf_GW_1=73.0d0/24.0d0,cf_GW_2=37.0d0/96.0d0)
      double precision DemiPi2
      parameter (DemiPi2=0.25d0*rPi*rPi)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision
     $     beta,                ! orientation aleatoire du pas de parcours brownien
     $     sd,cd,sb,cb,st,ct,sp,cp,stcp,stsp, ! sin et cos d'angles
     $     Vp_norme,Vp_tg,      ! norme et composant tg de la vitesse de la particule (ref renc)
     $     R2,                  ! distance quadratique restant a parcourir
     $     delta_adap,          ! pas (angulaire) adapte
     $     CSurV, delta2, Prob_dechir, Vtg_x, Vtg_y,
     $     arg_sqrt
      logical lPrem
      data lPrem /.TRUE./
      save lPrem
      logical lInit_GW
c
c fonctions :
c ^^^^^^^^^^^
      double precision random, RandomGauss
      logical lTest_GW_Capture
c......................................................................
c
c---- initialisations 
c
      lInit_GW = .TRUE.
#if (_ECRIT_INPUT_PA_>0) /* ecrire les arguments d'entree */
      if (lPrem) then
         lPrem=.FALSE.
         open(77,file='PA_INPUT.asc')
         write(77,100)
 100     format('# 1:Vp_renc_x 2:Vp_renc_y 3:Vp_renc_z ',
     $        '4:Vrenc_x 5:Vrenc_y 6:Vrenc_z 7:Vlc 8:delta 9:n_pas')
      end if
      write(77,'(8(E12.5," "),I10)') Vpart_renc,Vrenc,Vlc,delta,n_pas
#endif
#if (_ECRIT_PA_>0) /* ecrire le parcours aleatoire dans fichier */
      open(88,file='PA.asc')
      write(88,200) Vlc, delta, n_pas,
     $     delta_max, coef_secu
 200  format(
     $     '# Parcours aleatoire adaptatif',/,
     $     '# Vlc, delta, n_pas : ',2(E12.5,' '),I10,/,
     $     '# delta_max, coef_secu : ',2(E12.5,' '),/,
     $     '# 1:Vx 2:Vy 3:Vz 4:delta' )
#endif
      iIssue = 0
      delta2=delta**2
      Vp_tg = Vpart_renc(1)**2 + Vpart_renc(2)**2
      Vp_norme = sqrt(Vp_tg + Vpart_renc(3)**2)
      Vp_tg = sqrt(Vp_tg)

      ct = Vpart_renc(3)/Vp_norme
      st = Vp_tg/Vp_norme
      stcp = Vpart_renc(1)/Vp_norme
      stsp = Vpart_renc(2)/Vp_norme

      R2 = n_pas*delta**2
      CSurV = coef_secu/Vp_norme
      n_pas_adap = 0

      do while (R2.GT.0.0d0)

         n_pas_adap = n_pas_adap+1
         Vtg = sqrt( (Vp_norme*stcp+Vrenc(1))**2 +
     $        (Vp_norme*stsp+Vrenc(2))**2 )
         V2 = Vtg**2+(Vp_norme*ct+Vrenc(3))**2 ! Vitesse**2 dans ref amas
                                ! valeur de la longueur angulaire du pas suivant
         delta_adap = max(fac_delta_adap*delta, min(
     $        fac_delta_adap*delta_max,fac_delta_adap*CSurV*(Vtg-Vlc),
     $        sqrt(R2)) )

c$$$         Vlc_=sin(sqrt(theta2_GW))*sqrt(V2) !#######
c$$$         delta_adap = max(0.05d0*sqrt(theta2_GW), ! valeur de la longueur angulaire du pas suivant
c$$$     $        min(delta_max,CSurV*(Vtg-Vlc_),sqrt(R2)) )

c
c     DETECTION DE L'ENTREE DANS LE LOSS CONE
c
#if (_TYPE_DETEC_LC_==1)
         if (Vtg.LT.Vlc) then
            iIssue = 1
            goto 9 
         end if
#elif (_TYPE_DETEC_LC_==2) /* essai 20.4.2001 */
         if (Vtg.LT.Vlc) then
            Prob_dechir=(Vlc/delta)**2
            if (random().LT.Prob_dechir) then
               iIssue = 1
               goto 9 
            end if
         end if
#elif (_TYPE_DETEC_LC_==3) /* essai 20.4.2001 */
         call Calc_Proba_LC_Porb(Vlc,delta,Vtg,Prob_dechir)
         if (Prob_dechir.GT.1.0d-10) then
            if (random().LT.Prob_dechir) then 
                                ! il y a entre dans le LC
                                ! il faut determiner une valeur de la vitesse tg dans le LC (au cas ou...)
               if ((Vtg+Vlc).LT.delta) then ! LC entierement dans le disque de diffusion
                  Vtg = 0.50*Vlc ! arbitraire!!
               else if ((Vtg+delta).LT.Vlc) then ! disque de diffusion entierement dans le LC
                  Vtg = Vtg     ! Pour que ca soit clair!
               else             ! recouvrement partiel
                  Vtg = 0.5d0*(Vtg+Vlc-delta) 
               end if
               iIssue = 1
               goto 9 
            end if
         end if
#elif (_TYPE_DETEC_LC_==4) /* essai 26.4.2001 : on tire la "position" selon une gaussienne */
         if (Vtg.LT.Vlc .OR. (Vtg-Vlc.LT.3.0d0*delta)) then
            Vtg_x = Vtg + rRacDeux1*delta*RandomGauss()
            Vtg_y =       rRacDeux1*delta*RandomGauss()
            Vtg = sqrt(Vtg_x*Vtg_x+Vtg_y*Vtg_y)
            if (Vtg.LT.Vlc) then
               iIssue = 1
               goto 9 
            end if
         end if 
#else
         write(0,*) '!!! Valeur de _TYPE_DETEC_LC_ incorrecte !!!'
         call Meurt('ParcoursAlea_LC')
#endif

#if (_ECRIT_PA_>0) /* ecrire le parcours aleatoire dans fichier */
         write(88,210) Vp_norme*stcp+Vrenc(1),Vp_norme*stsp+Vrenc(2),
     $        Vp_norme*ct+Vrenc(3),delta_adap
 210     format(4(E12.5,' '))
#endif

#if (_CAPT_GW_>0) /* On regarde si l'etoile peut etre "captee" par emission d'ondes grav. */
c
c     CAPTAGE SUR UNE ORBITE A FORTE RADIATION GRAV.
c
         if (lTest_GW_Capture(lInit_GW,V2,Vtg**2)) then
#   if (_ECRIT_CAPT_GW_>0)
            write(iFich_CaptGW,300) iPas_Evol,
     $           a_GW, e_GW, theta2_GW, T_GW_sur_Porb,
     $           Trel_over_Porb_GW                 
 300        format(I10,' ',E11.4,' ',E20.13,3(' ',E11.4))
            call flush_MD(iFich_CaptGW)
#   endif
            iIssue = 2
            goto 9
         end if
#endif /* if (CAPT_GW>0) */

         R2=R2-delta_adap**2    ! diminution de la distance quad restant a parcourir
         cd = cos(delta_adap)
         sd = sin(delta_adap)
         beta = r2Pi*random()
         cb = cos(beta)
         sb = sin(beta)
         if (st.GT.minus) then
            cp = stcp/st
            sp = stsp/st
            stcp = cd*stcp-sd*cb*ct*cp-sd*sb*sp
            stsp = cd*stsp-sd*cb*ct*sp+sd*sb*cp
         else
            st = 0.0d0
            ct = ct/abs(ct)
            stcp = 0.0d0
            stsp = 0.0d0
         end if
         ct = cd*ct + sd*cb*st
#if (_SURVEILLE_SQRT_ > 0)
         arg_sqrt = 1.0d0-ct**2
         if (arg_sqrt.LT.0.0d0) then
            write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
            call Meurt('ParcoursAlea_LC')
         end if
         st = sqrt(arg_sqrt)
#else
         st = sqrt(1.0d0-ct**2)
#endif
      end do
c
c---- sortie de la boucle sans entree dans le LC
c
      Vtg = -Vtg
c
c---- sortie de la boucle avec entree dans le LC
c
 9    continue
c
      end
c
