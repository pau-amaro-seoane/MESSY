#if (_COLLISIONS_>0)

#define _TRAITEMENT_DBRD_M_HAUT_ 1
#define _TRAITEMENT_DBRD_M_BAS_  1
#define _FIT_LAMBDA_CRIT_MERG_ 1

c----------------------------------------------------------------------
      subroutine Collision_Stellaire (Type1,Type2,
     $     Mini1,Mini2,Vrel_norm,b_norm,
     $     Mfin1,Mfin2,Vrel_norm_p,Theta)
c----------------------------------------------------------------------
      include 'EvolAmas_Collisions_inc.f'
#if (_TYPE_COLL_>=10)
      include 'EvolAmas_Collisions_SPH_inc.f'
#endif
c     
c     arguments :
c     ^^^^^^^^^^^
      integer*1 Type1,Type2          ! le type des 2 etoiles

      double precision          ! Les conditions initiales
     $                          ! ------------------------
     $     Mini1,Mini2,         ! Mini1,Mini2 : Les masses en masses solaires
     $     Vrel_norm,b_norm     ! Vrel_norm : La vitesse relative a l'infini normalisee par la vitesse d'evasion
                                ! b_norm : le parametre d'impact normalise par R1+R2

c     La vitesse d'evasion est par definition sqrt(2*(Mini1+Mini2)/(Rini1+Rini2))

      double precision          ! Le resultat de la collison
     $                          ! --------------------------
     $     Mfin1,Mfin2,         ! Mfin1,Mfin2 : Les masses en masses solaires
     $     Vrel_norm_p,Theta    ! Vrel_norm_p : La vitesse relative a l'infini (in units of the initial Vesc)
     $                          ! Theta : angle de deviation de la trajectoire (en valeur absolue)

c     ATTENTION: Ces grandeurs sont calculees pour le ref du CM de la paire
c     APRES la collision. Ce ref est en general different de celui
c     du CM AVANT la collision. Il faut donc faire l'hypothese simplificatrice
c     qu'il s'agit du meme ref pour se servir des resultats de cette routine !!
c     
c     constantes :
c     ^^^^^^^^^^^^
#if (_TYPE_COLL_==10) /* Grille SPH */
      double precision C_dbrd_min, C_dbrd_max
      parameter (C_dbrd_min=1.05d0, C_dbrd_max=0.95d0) ! Changed from 1% to 5% on Dec 14 2003
#endif
c     
c     var locales :
c     ^^^^^^^^^^^^^
      integer*1 Type1_prv,Type2_prv
      double precision Mini1_prv,Mini2_prv,M_reste,dmin_norm, M_prov,
     $     lambda, dmin_tidal_disrup
      logical lEchange12
#if (_TYPE_COLL_==2) /* Formules de Davies */
      double precision
     $     btilde, q, beta_rel, Delta, Eta, b2, P1,P2,P3,P4,P5, d1,d2
#elif (_TYPE_COLL_==10) /* Grille SPH */
      integer iStat_Coll_Grille_SPH
      double precision factM1_grid, factM2_grid, rap, dMoverM,
     $     q_ini_grid, q_fin_grid, q_ini_real, q_fin_real, Mtot_fin
      
      logical l_interp, l_dbrd_M, lMassesScaled
#endif
c     
c     fonction :
c     ^^^^^^^^^^
      double precision random, RelMR, lambda_crit_merger_toy

#if (_HARDEN_SPH_COLL_TREATMENT_>0)
      double precision lambda_crit_merger
#endif
c......................................................................
c     
c---- conventionnellement, l'etoile 1 est la moins massive
c     
      if(Mini1.LT.Mini2) then
         lEchange12=.FALSE.
         Type1_prv=Type1
         Type2_prv=Type2
         Mini1_prv=Mini1
         Mini2_prv=Mini2
      else
         lEchange12=.TRUE.
         Type1_prv=Type2
         Type2_prv=Type1
         Mini1_prv=Mini2
         Mini2_prv=Mini1
      end if
#if (_TYPE_COLL_==10) /* Grille SPH */
      factM1_grid=1.0d0
      factM2_grid=1.0d0
      lMassesScaled = .FALSE.
#endif
c     
c---- Resultat par defaut : rien ne change
c     
      Mfin1 = Mini1_prv
      Mfin2 = Mini2_prv
      Vrel_norm_p = Vrel_norm
      Theta = 0.0d0
c     
c==== On passe en revue les differents types de collisions mais
c     seuls les cas MS+MS sont traites pour l'instant
c     =========================================================
c     
                                !  %%%%%%%%%%%
                                !  % MS + MS %
                                !  %%%%%%%%%%%
      if (Type1_prv.EQ.TypeEtoile_MS .AND. Type2_prv.EQ.TypeEtoile_MS) then
c     
#if (_TYPE_COLL_==0)
c     
c---- On ne change rien !!!!!! (pour tests)
c     
#elif (_TYPE_COLL_==1)
c     
c---- Modele ultra simple : coalescence en-dessous de Vrel_Lim
c     destruction totale en-dessus
c     
         if (Vrel_norm.GT.Fact_Vrel_Lim) then ! destruction totale
            Mfin1=0.0d0
            Mfin2=0.0d0
            Vrel_norm_p=0.0d0
            Theta=0.0d0
         else                   ! coalescence
            Mfin1=0.0d0
            Mfin2=Mini1_prv+Mini2_prv
            Mfin2=(1.0d0-CollMassLossFracDef)*Mfin2 ! Apply flat mass-loss rate
            Vrel_norm_p=0.0d0
            Theta=0.0d0
         end if

#elif (_TYPE_COLL_==2) /* formules de Davies pour Rauch(1999) */
                                ! rapport de masses
         q = Mini1_prv/Mini2_prv
                                ! vitesse relative / vitesse de la lumiere * V_sol/Vesc (normalisation par Vesc par-rapport a Davies)
         Beta_rel = Vrel_norm / c_en_Vesc_sol
                                ! distance mini d'approche pour rencontre hyperbolique
         btilde = 2.0d0*(b_norm*Vrel_norm)**2 /
     $        (1.0d0+sqrt(1.0d0+4.0d0*b_norm**2*Vrel_norm**4))
         b2 = btilde**2
         P1 = 5.5530d0  - 2.8607d0 *btilde - 10.5255d0*b2
         P2 = 0.57366d0 + 0.18251d0*btilde - 0.13155d0*b2
         P3 = 1.0d-6*(3.057d0 - 4.825d0*btilde + 2.616d0*b2)
         P4 = -16.138d0 - 7.26d0 *btilde
         P5 = 0.417d0   - 15.96d0*btilde
         Eta = q**0.9d0         ! Eta = (dM2/M2)/(dM1/M1)
         Delta = q**0.2d0*Beta_rel**1.25d0*exp(P1) ! Delta = (dM1+dM2)/(M1+M2)
         d1 = Delta*(Mini1_prv+Mini2_prv)/(Mini1_prv+Eta*Mini2_prv) ! d1 = dM1/M1     
         d2 = Eta*d1            ! d2 = dM2/M2
         Mfin1 = max(0.0d0,(1.0d0-d1)*Mini1_prv)
         Mfin2 = max(0.0d0,(1.0d0-d2)*Mini2_prv)
         if (Mfin1.GT.0.0d0 .AND. Mfin2.GT.0.0d0) then ! les 2 etoiles ont survecu
            Vrel_norm_p = Vrel_norm*(  1.0d0-min( 1.0d0,
     $           beta_rel**(-2.75d0)*exp(P4)+beta_rel**(-0.5d0)*exp(P5)
     $           ) )
            Theta = (rPi-2.0d0/(P2+P3/beta_rel**2))*q**(-0.3d0)
         else                   ! une des 2 etoiles est detruite : on admet qu'il s'agit de coalescence
            Vrel_norm_p = 0.0d0
            Theta = 0.0d0
         end if

         if (Vrel_norm.LT.Fact_Vrel_Lim) then ! coalescence
            Mfin2=Mfin1+Mfin2
            Mfin1=0.0d0
            Vrel_norm_p=0.0d0
            Theta=0.0d0
         end if

#elif (_TYPE_COLL_==3) /* Rencontre gravitionnelle pure */
         Vrel_norm_p=Vrel_norm
         Mfin1 = Mini1_prv
         Mfin2 = Mini2_prv
         Theta = atan2(0.5d0,b_norm*Vrel_norm**2)

#elif (_TYPE_COLL_==4) /* very simple parametrisation for merger vs fly-bys; flat mass loss rate */
                                ! lambda=d_min/(R1+R2) where d_min is the pericentre distance for point-mass encounter
         lambda = 2.0d0*(b_norm*Vrel_norm)**2 /
     $        (1.0d0+sqrt(1.0d0+4.0d0*b_norm**2*Vrel_norm**4))
         if (lambda.LT.lambda_crit_merger_toy(Mini1_prv,Mini2_prv,Vrel_norm)) then ! Merger
            Mfin1=0.0d0
            Mfin2=(1.0d0-CollMassLossFracDef)*(Mini1_prv+Mini2_prv) ! Apply flat mass-loss rate
            Vrel_norm_p=0.0d0
            Theta=0.0d0
         else                   ! Fly-by with same fractional mass loss (!)
            Mfin1=(1.0d0-CollMassLossFracDef)*Mini1_prv
            Mfin2=(1.0d0-CollMassLossFracDef)*Mini2_prv
            Vrel_norm_p=Vrel_norm
            Theta=0.0d0
         end if
         
#elif (_TYPE_COLL_==9)
c     
c---- Modele non-physique, pour tests
c     
         Vrel_norm_p=random()*Vrel_norm
         M_reste = random()*(Mini1_prv+Mini2_prv)
         Mfin1=random()*M_reste
         Mfin2=M_reste-Mfin1
         Theta=random()*rPi

#elif (_TYPE_COLL_==10)
c     
c---- collisions selon resultats de MES simultaions SPH
c     
c     Il faut d'abord verifier que l'on se trouve dans le domaine de masses
c     couvert par la grille
c     
         l_dbrd_M = .FALSE.
         l_interp = .TRUE.

         if      (Mini1_prv.GT.Msup_grille_SPH .OR.
     $        Mini2_prv.GT.Msup_grille_SPH)
     $        then              ! Debordement par le haut
            n_dbrdM_CollSPH_h = n_dbrdM_CollSPH_h+1
            l_dbrd_M = .TRUE.

#   if (_WRITE_COLL_EXCEPTIONS_>0)
            write(iFich_ExcColl,100) iPas_Evol,
     $           Mini1_prv,Mini2_prv,Vrel_norm,b_norm
 100        format(
     $           '= iPas_evol : ',I10,/,
     $           '- M_too_large',/,
     $           ' Mi1,Mi2,V,b =',4(' ',E12.5))
#   endif

#   if (_TRAITEMENT_DBRD_M_HAUT_==1)
c     
c     On traite ces cas en se rapportant a des cas de la grille avec masse(s)
c     plus petite(s). Ainsi on calcule une collision "remplacante".
c     
            M_prov = Mini2_prv
            Mini2_prv = C_dbrd_max*Msup_grille_SPH
            factM2_grid = Mini2_prv/M_prov
            M_prov = Mini1_prv
            Mini1_prv = max(C_dbrd_min*Minf_grille_SPH,
     $           factM2_grid*Mini1_prv)
            factM1_grid = Mini1_prv/M_prov
            lMassesScaled = .TRUE.

#   if (_WRITE_COLL_EXCEPTIONS_>0)
            write(iFich_ExcColl,105) Mini1_prv,Mini2_prv
 105        format(
     $           '+ M_rescaled',/,
     $           ' Mi1,Mi2 =',2(' ',E12.5))
#   endif

#   else
            l_interp = .FALSE.
#   endif /* if (_TRAITEMENT_DBRD_M_HAUT_==1) */

         else if (Mini1_prv.LT.Minf_grille_SPH .OR.
     $           Mini2_prv.LT.Minf_grille_SPH)
     $           then           ! Debordement par le bas
            n_dbrdM_CollSPH_b = n_dbrdM_CollSPH_b+1
            l_dbrd_M = .TRUE.

#   if (_WRITE_COLL_EXCEPTIONS_>0)
            write(iFich_ExcColl,110) iPas_Evol,
     $           Mini1_prv,Mini2_prv,Vrel_norm,b_norm
 110        format(
     $           '= iPas_evol : ',I10,/,
     $           '- M_too_small',/,
     $           ' Mi1,Mi2,V,b =',4(' ',E12.5))
#   endif

#   if (_TRAITEMENT_DBRD_M_BAS_==1)
c     
c     On traite ces cas en se rapportant a des cas de la grille avec masse(s)
c     plus grande(s). Ainsi on calcule une collision "remplacante".
c     
            M_prov = Mini1_prv
            Mini1_prv = C_dbrd_min*Minf_grille_SPH
            factM1_grid = Mini1_prv/M_prov
            M_prov = Mini2_prv
            Mini2_prv = min(C_dbrd_max*Msup_grille_SPH,factM1_grid*Mini2_prv)
            factM2_grid = Mini2_prv/M_prov
            lMassesScaled = .TRUE.

#   if (_WRITE_COLL_EXCEPTIONS_>0)
            write(iFich_ExcColl,105) Mini1_prv,Mini2_prv
#   endif

#   else
            l_interp = .FALSE.
#   endif /* if (_TRAITEMENT_DBRD_M_BAS_==1) */
         end if
c     
c---- Interpolation du resultat de la collision dans la grille SPH
c     
         if (l_interp) then
            call ResColl_selon_Grille_SPH(
     $           Mini1_prv,Mini2_prv,Vrel_norm,b_norm,
     $           Mfin1,Mfin2,Vrel_norm_p,Theta,iStat_Coll_Grille_SPH )

#   if (_HARDEN_SPH_COLL_TREATMENT_==3)
c     
c     Check whether it is a fly-by or a merger (at relatively low velocity)
c     
            if (Vrel_norm .LT. Nu_Max_Merger) then
                                ! d_min/(R1+R2)
               dmin_norm = 2.0d0*(b_norm*Vrel_norm)**2 /
     $              (1.0d0+sqrt(1.0d0+4.0d0*b_norm**2*Vrel_norm**4))
               if (dmin_norm .LT.
     $              lambda_crit_merger(Mini1_prv,Mini2_prv,Vrel_norm))
     $              then        ! Fitting formula says it's a merger...
                  if (Vrel_norm_p.NE.0.0d0)
     $                 then     ! ...but blind interpolation thought it was a fly-by

#   if (_WRITE_COLL_EXCEPTIONS_>0)
                     write(iFich_ExcColl,170) iPas_Evol,
     $                    Mini1_prv,Mini2_prv,Vrel_norm,b_norm,Mfin1,Mfin2
 170                 format(
     $                    '= iPas_evol : ',I10,/,
     $                    '- Not_flyb_but_mrg',/,
     $                    ' Mi1,Mi2,V,b,Mf1,Mf2 =',6(' ',E12.5))
#   endif
                     
                     Mfin2 = Mfin1 + Mfin2
                     Mfin1 = 0.0d0
                     Vrel_norm_p = 0.0d0
                     Theta = 0.0d0

#   if (_WRITE_COLL_EXCEPTIONS_>0)
                     write(iFich_ExcColl,175) Mfin1,Mfin2
 175                 format(
     $                    '+ Force_mrg',/,
     $                    ' Mf1,Mf2 =',2(' ',E12.5))
#   endif
                  end if
                  
               else             ! Fitting formula says it's a fly-by...
                  if (Vrel_norm_p.EQ.0.0d0)
     $                 then     ! ...but blind interpolation thought it was a merger.

#   if (_WRITE_COLL_EXCEPTIONS_>0)
                     write(iFich_ExcColl,180) iPas_Evol,
     $                    Mini1_prv,Mini2_prv,Vrel_norm,b_norm,Mfin1,Mfin2
 180                 format(
     $                    '= iPas_evol : ',I10,/,
     $                    '- Not_mrg_but_flyb',/,
     $                    ' Mi1,Mi2,V,b,Mf1,Mf2 =',6(' ',E12.5))
#   endif
                                ! Give small relative velocity
                     Vrel_norm_p = 0.01d0*Vrel_norm
                                ! Zero deflection (for lack of better prescription)
                     Theta = 0.0d0

#   if (_WRITE_COLL_EXCEPTIONS_>0)
                     write(iFich_ExcColl,185) Mfin1,Mfin2,Vrel_norm_p,Theta
 185                 format(
     $                    '+ Force_flyb',/,
     $                    ' Vp,Theta =',4(' ',E12.5))
#   endif
                  end if
               end if
            end if
c     
c     Force conservation of mass ratio in case of "fly-bys"
c     
            if (Vrel_norm_p.GT.0.0d0) then
               dMoverM = (Mini1_prv+Mini2_prv-Mfin1-Mfin2)/
     $              (Mini1_prv+Mini2_prv)
               Mtot_fin = (1.0d0-dMoverM)*(Mini2_prv+Mini1_prv)
               q_ini_real = Mini1_prv/Mini2_prv
               q_fin_real = q_ini_real
               Mfin2 = Mtot_fin/(1.0d0+q_fin_real)
               Mfin1 = Mtot_fin-Mfin2
            end if
                  
#   endif /* (_HARDEN_SPH_COLL_TREATMENT_==3) */

c     
c     Correct back to real initial masses when masses have been rescaled
c     
            if (lMassesScaled) then

#   if (_WRITE_COLL_EXCEPTIONS_>0)
               write(iFich_ExcColl,120) Mfin1,Mfin2
 120           format(
     $              '+ M_scaled_back',/,
     $              ' Mfg1,Mfg2 =',2(' ',E12.5))
#   endif
               dMoverM = (Mini1_prv+Mini2_prv-Mfin1-Mfin2)/(Mini1_prv+Mini2_prv) ! Fractional mass-loss 
               q_ini_grid = Mini1_prv/Mini2_prv
               if (Mfin2.GT.0.0d0) then
                  q_fin_grid = Mfin1/Mfin2
                  Mini1_prv = Mini1_prv/factM1_grid
                  Mini2_prv = Mini2_prv/factM2_grid
                  q_ini_real = Mini1_prv/Mini2_prv
                  q_fin_real = q_fin_grid * q_ini_real/q_ini_grid
                  Mtot_fin = (1.0d0-dMoverM)*(Mini2_prv+Mini1_prv)
                  Mfin2 = Mtot_fin/(1.0d0+q_fin_real)
                  Mfin1 = Mtot_fin-Mfin2
               else
                  Mfin1 = 0.0d0
                  Mfin2 = 0.0d0
               end if
               
#   if (_WRITE_COLL_EXCEPTIONS_>0)
               write(iFich_ExcColl,125) Mfin1,Mfin2
 125           format(
     $              ' Mf1,Mf2 =',2(' ',E12.5))
#   endif
            end if
            
            if (iStat_Coll_Grille_SPH.NE.i_interp) then
               
               if (iStat_Coll_Grille_SPH.EQ.i_extrap) then
                  n_extrap_CollSPH = n_extrap_CollSPH+1
                  
#   if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,130) iPas_Evol,
     $                 Mini1_prv,Mini2_prv,Vrel_norm,b_norm
 130              format(
     $                 '= iPas_evol : ',I10,/,
     $                 '- Extrap',/,
     $                 ' Mi1,Mi2,V,b =',4(' ',E12.5))
#   endif
c     
c     Que fait-on si la grille ne fournit pas de resultat pour ces
c     conditions initiales ?
c     - Pour l'instant, seul le cas de grand parametre d'impact peut etre
c     "sauve" en admettant une rencontre 2-corps ponctuels (pas d'effet
c     collisionnel). Les autres cas sont traite par une routine d'erreur.
c     
               else if (iStat_Coll_Grille_SPH.EQ.i_distant) then
                  n_distant_CollSPH = n_distant_CollSPH+1
                  
#   if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,135) iPas_Evol,
     $                 Mini1_prv,Mini2_prv,Vrel_norm,b_norm
 135              format(
     $                 '= iPas_evol : ',I10,/,
     $                 '- Distant',/,
     $                 ' Mi1,Mi2,V,b =',4(' ',E12.5))
#   endif

               else if (iStat_Coll_Grille_SPH.EQ.i_inconnu) then
                  n_inconnu_CollSPH = n_inconnu_CollSPH+1

#   if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,140) iPas_Evol,
     $                 Mini1_prv,Mini2_prv,Vrel_norm,b_norm
 140              format(
     $                 '= iPas_evol : ',I10,/,
     $                 '- Unknown',/,
     $                 ' Mi1,Mi2,V,b =',4(' ',E12.5))
#   endif
                  call Traiter_ResColl_SPH_inconnu(
     $                 Mini1_prv,Mini2_prv,Vrel_norm,b_norm)
                  Mfin1=Mini1_prv
                  Mfin2=Mini2_prv
                  Vrel_norm_p=Vrel_norm
                  Theta=-1.0d-6 ! petite valeur negative pour pouvoir identifier ces cas
               else
                  write(0,*) '!!! Value of iStat_Coll_Grille_SPH (',
     $                 iStat_Coll_Grille_SPH,') not valid !!!'
                  call Meurt('Collision_Stellaire')
               end if
            end if

         end if                 ! if (l_interp) then

#   if (_HARDEN_SPH_COLL_TREATMENT_>0)
c     
c     Preserve the "hierarchy": Mfin2 >= Mfin1
c     
               if (Mfin2.LT.Mfin1) then

#      if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,145) iPas_Evol,
     $                 Mini1_prv,Mini2_prv,Vrel_norm,b_norm,Mfin1,Mfin2
 145              format(
     $                 '= iPas_evol : ',I10,/,
     $                 '- Reverse',/,
     $                 ' Mi1,Mi2,V,b,Mf1,Mf2 =',6(' ',E12.5))
#      endif

                  M_prov = Mfin1
                  Mfin1 = Mfin2
                  Mfin2 = M_prov
               end if

#      if (_HARDEN_SPH_COLL_TREATMENT_==1)
c     
c     Do not allow that the small star grows when the large one has a decreased mass
c     
               if (Mfin2.LT.Mini2_prv .AND. Mfin1.GT.Mini1_prv) then

#         if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,150) iPas_Evol,
     $                 Mini1_prv,Mini2_prv,Vrel_norm,b_norm,Mfin1,Mfin2
 150              format(
     $                 '= iPas_evol : ',I10,/,
     $                 '- Unphys_growth (light star)',/,
     $                 ' Mi1,Mi2,V,b,Mf1,Mf2 =',6(' ',E12.5))
#         endif

                  Mfin2 = Mfin2 + Mfin1 - Mini1_prv
                  Mfin1 = Mini1_prv

#         if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,155) Mfin1,Mfin2
 155              format(
     $                 '+ Corr_Masses',/,
     $                 ' Mf1,Mf2 =',2(' ',E12.5))
               end if
#         endif
#      endif /* (_HARDEN_SPH_COLL_TREATMENT_==1) */

#      if (_HARDEN_SPH_COLL_TREATMENT_==2)
c     
c     Do not allow one star to grow when the other looses mass
c     (this nearly never happens in the SPH simulations)
c     
               if ((Mfin2.LT.Mini2_prv .AND. Mfin1.GT.Mini1_prv) .OR.
     $              (Mfin1.LT.Mini1_prv .AND. Mfin2.GT.Mini2_prv)) then

#         if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,150) iPas_Evol,
     $                 Mini1_prv,Mini2_prv,Vrel_norm,b_norm,Mfin1,Mfin2
 150              format(
     $                 '= iPas_evol : ',I10,/,
     $                 '- Unphys_growth (massive star)',/,
     $                 ' Mi1,Mi2,V,b,Mf1,Mf2 =',6(' ',E12.5))
#         endif
c     
c     One assumes that the mass ratio is conserved
c     (quite OK except for close, high-velocity fly-bys)
c     
                  Mtot_fin = Mfin1 + Mfin2
                  q_ini_real = Mini1_prv/Mini2_prv
                  Mfin2 = Mtot_fin/(1.0d0+q_ini_real)
                  Mfin1 = Mtot_fin - Mfin2

#         if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,155) Mfin1,Mfin2
 155              format(
     $                 '+ Corr_Masses',/,
     $                 ' Mf1,Mf2 =',2(' ',E12.5))
#         endif
               end if
#      endif /* (_HARDEN_SPH_COLL_TREATMENT_==2) */


#      if (_HARDEN_SPH_COLL_TREATMENT_<3)
c     
c     Check more carefully whether there is merging or fly-by
c     
               if (Mfin1.LT.Mdebris_coll_max) then

#         if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,160) iPas_Evol,
     $                 Mini1_prv,Mini2_prv,Vrel_norm,b_norm,Mfin1,Mfin2
 160              format(
     $                 '= iPas_evol : ',I10,/,
     $                 '- Fin_mass_too_small',/,
     $                 ' Mi1,Mi2,V,b,Mf1,Mf2 =',6(' ',E12.5))
#         endif

                  Mfin2 = Mfin1 + Mfin2
                  Mfin1 = 0.0d0
                  Vrel_norm_p = 0.0d0
                  Theta = 0.0d0 
                  if (Mfin2.LT.Mdebris_coll_max) then
                     Mfin2 = 0.0d0
                  end if

#         if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,165) Mfin1,Mfin2
 165              format(
     $                 '+ Corr_Masses',/,
     $                 ' Mf1,Mf2 =',2(' ',E12.5))
#         endif

               end if

               if (Vrel_norm .LT. Nu_Max_Merger) then
                                ! d_min/(R1+R2)
                  dmin_norm = 2.0d0*(b_norm*Vrel_norm)**2 /
     $                 (1.0d0+sqrt(1.0d0+4.0d0*b_norm**2*Vrel_norm**4))
                  if (dmin_norm .LT.
     $                 lambda_crit_merger(Mini1_prv,Mini2_prv,Vrel_norm))
     $                 then     ! Fitting formula says it's a merger...
                     if (Vrel_norm_p.NE.0.0d0)
     $                    then  ! ...but blind interpolation thought it was a fly-by

#         if (_WRITE_COLL_EXCEPTIONS_>0)
                        write(iFich_ExcColl,170) iPas_Evol,
     $                       Mini1_prv,Mini2_prv,Vrel_norm,b_norm,Mfin1,Mfin2
 170                    format(
     $                       '= iPas_evol : ',I10,/,
     $                       '- Not_flyb_but_mrg',/,
     $                       ' Mi1,Mi2,V,b,Mf1,Mf2 =',6(' ',E12.5))
#         endif

                        Mfin2 = Mfin1 + Mfin2
                        Mfin1 = 0.0d0
                        Vrel_norm_p = 0.0d0
                        Theta = 0.0d0

#         if (_WRITE_COLL_EXCEPTIONS_>0)
                        write(iFich_ExcColl,175) Mfin1,Mfin2
 175                    format(
     $                       '+ Force_mrg',/,
     $                       ' Mf1,Mf2 =',2(' ',E12.5))
#         endif
                     end if
                  else          ! Fitting formula says it's a fly-by...
                     if (Vrel_norm_p.EQ.0.0d0)
     $                    then  ! ...but blind interpolation thought it was a merger
                                ! We decide to keep the initial mass ratio and preserve fractional mass loss

#         if (_WRITE_COLL_EXCEPTIONS_>0)
                        write(iFich_ExcColl,180) iPas_Evol,
     $                       Mini1_prv,Mini2_prv,Vrel_norm,b_norm,Mfin1,Mfin2
 180                    format(
     $                       '= iPas_evol : ',I10,/,
     $                       '- Not_mrg_but_flyb',/,
     $                       ' Mi1,Mi2,V,b,Mf1,Mf2 =',6(' ',E12.5))
#         endif
                        dMoverM = (Mini1_prv+Mini2_prv-Mfin1-Mfin2)/(Mini1_prv+Mini2_prv)
                        Mtot_fin = (1.0d0-dMoverM)*(Mini2_prv+Mini1_prv)
                        q_ini_real = Mini1_prv/Mini2_prv
                        q_fin_real = q_ini_real
                        Mfin2 = Mtot_fin/(1.0d0+q_fin_real)
                        Mfin1 = Mtot_fin-Mfin2
                                ! Give small relative velocity
                        Vrel_norm_p = 0.01d0*Vrel_norm
                                ! Zero deflection (far lack of better prescription)
                        Theta = 0.0d0

#         if (_WRITE_COLL_EXCEPTIONS_>0)
                        write(iFich_ExcColl,185) Mfin1,Mfin2,Vrel_norm_p,Theta
 185                    format(
     $                       '+ Force_flyb',/,
     $                       ' Mf1,Mf2,Vp,Theta =',4(' ',E12.5))
#         endif
                     end if
                  end if
               end if
#      endif /* (_HARDEN_SPH_COLL_TREATMENT_<3) */
#   endif /* (_HARDEN_SPH_COLL_TREATMENT_>0) */

#else 
               write(0,*) '!!! Valeur invalide de _TYPE_COLL_'
               call Meurt('Collision_Stellaire')
#endif /* various treatments of collision */
                                ! on ne fait pas croitre les etoiles
                                ! deja tres massives (pour eviter l'emballement)
               if (Mfin2.GT.Mmerg_coll_max) then ! on ne change rien !!!

#   if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,190) iPas_Evol,
     $                 Mini1_prv,Mini2_prv,Vrel_norm,b_norm,Mfin1,Mfin2
 190              format(
     $                 '= iPas_evol : ',I10,/,
     $                 '- Fin_mass_too_high',/,
     $                 ' Mi1,Mi2,V,b,Mf1,Mf2 =',6(' ',E12.5))
#   endif
                  Mfin1 = Mini1_prv
                  Mfin2 = Mini2_prv
                  Vrel_norm_p = Vrel_norm
                  Theta = 0.0d0

#   if (_WRITE_COLL_EXCEPTIONS_>0)
                  write(iFich_ExcColl,195) Mfin1,Mfin2
 195              format(
     $                 '+ Corr_Masses',/,
     $                 ' Mf1,Mf2 =',2(' ',E12.5))
#   endif

               end if
                                !  %%%%%%%%%%%
                                !  % WD + MS %
                                !  %%%%%%%%%%%
            else if (
     $              (Type1_prv.EQ.TypeEtoile_WD .AND.
     $              Type2_prv.EQ.TypeEtoile_MS)
     $              .OR.
     $              (Type2_prv.EQ.TypeEtoile_WD .AND.
     $              Type1_prv.EQ.TypeEtoile_MS)
     $              ) then
#if ( _TRAITEMENT_COLL_MS_RMN_>0)
               dmin_norm = 2.0d0*(b_norm*Vrel_norm)**2 /
     $              (1.0d0+sqrt(1.0d0+4.0d0*b_norm**2*Vrel_norm**4))
#endif
#if ( _TRAITEMENT_COLL_MS_RMN_==1) /* Destruction complete de l'etoile de MS */
               if (dmin_norm.LT.CollDminDestr_MS_vs_RMN) then
                  Vrel_norm_p=0.0d0
                  Theta=0.0d0
                  if (Type1_prv.EQ.TypeEtoile_MS) then
                     Mfin1 = 0.0d0
                  else
                     Mfin2 = 0.0d0
                  end if
               end if
#elif ( _TRAITEMENT_COLL_MS_RMN_==2) /* Tidal disruption and partial accretion onto WD */
               if (Type1_prv.EQ.TypeEtoile_MS) then
                  call Solve_Roche
     $                 (Mini1_prv/Mini2_prv,dmin_tidal_disrup)
                  if (dmin_norm.LT.dmin_tidal_disrup) then
                     Vrel_norm_p=0.0d0
                     Theta=0.0d0
                     Mfin1 = 0.0d0
                     Mfin2 = Mini2_prv + CollAccrFrac_MS_on_WD*Mini1_prv
                  end if
               else
                  call Solve_Roche
     $                 (Mini2_prv/Mini1_prv,dmin_tidal_disrup)
                  if (dmin_norm.LT.dmin_tidal_disrup) then
                     Vrel_norm_p=0.0d0
                     Theta=0.0d0
                     Mfin2 = 0.0d0
                     Mfin1 = Mini1_prv + CollAccrFrac_MS_on_WD*Mini2_prv
                  end if
               end if
#endif
               continue

                                !  %%%%%%%%%%%
                                !  % NS + MS %
                                !  %%%%%%%%%%%
            else if (
     $              (Type1_prv.EQ.TypeEtoile_NS .AND.
     $              Type2_prv.EQ.TypeEtoile_MS)
     $              .OR.
     $              (Type2_prv.EQ.TypeEtoile_NS .AND.
     $              Type1_prv.EQ.TypeEtoile_MS)
     $              ) then
#if ( _TRAITEMENT_COLL_MS_RMN_>0)
               dmin_norm = 2.0d0*(b_norm*Vrel_norm)**2 /
     $              (1.0d0+sqrt(1.0d0+4.0d0*b_norm**2*Vrel_norm**4))
#endif
#if ( _TRAITEMENT_COLL_MS_RMN_==1) /* Destruction complete de l'etoile de MS */
               if (dmin_norm.LT.CollDminDestr_MS_vs_RMN) then
                  Vrel_norm_p=0.0d0
                  Theta=0.0d0
                  if (Type1_prv.EQ.TypeEtoile_MS) then
                     Mfin1 = 0.0d0
                  else
                     Mfin2 = 0.0d0
                  end if
               end if
#elif ( _TRAITEMENT_COLL_MS_RMN_==2) /* Tidal disruption and partial accretion onto BH */
               if (Type1_prv.EQ.TypeEtoile_MS) then
                  call Solve_Roche(Mini1_prv/Mini2_prv,dmin_tidal_disrup)
                  if (dmin_norm.LT.dmin_tidal_disrup) then
                     Vrel_norm_p=0.0d0
                     Theta=0.0d0
                     Mfin1 = 0.0d0
                     Mfin2 = Mini2_prv + CollAccrFrac_MS_on_NS*Mini1_prv
                  end if
               else
                  call Solve_Roche(Mini2_prv/Mini1_prv,dmin_tidal_disrup)
                  if (dmin_norm.LT.dmin_tidal_disrup) then
                     Vrel_norm_p=0.0d0
                     Theta=0.0d0
                     Mfin2 = 0.0d0
                     Mfin1 = Mini1_prv + CollAccrFrac_MS_on_NS*Mini2_prv
                  end if
               end if
#endif
               continue

                                !  %%%%%%%%%%%
                                !  % BH + MS %
                                !  %%%%%%%%%%%
            else if (
     $              (Type1_prv.EQ.TypeEtoile_BH .AND. Type2_prv.EQ.TypeEtoile_MS)
     $              .OR.
     $              (Type2_prv.EQ.TypeEtoile_BH .AND. Type1_prv.EQ.TypeEtoile_MS)
     $              ) then
#if ( _TRAITEMENT_COLL_MS_RMN_>0)
               dmin_norm = 2.0d0*(b_norm*Vrel_norm)**2 /
     $              (1.0d0+sqrt(1.0d0+4.0d0*b_norm**2*Vrel_norm**4))
#endif
#if ( _TRAITEMENT_COLL_MS_RMN_==1) /* Destruction complete de l'etoile de MS */
               if (dmin_norm.LT.CollDminDestr_MS_vs_RMN) then
                  Vrel_norm_p=0.0d0
                  Theta=0.0d0
                  if (Type1_prv.EQ.TypeEtoile_MS) then
                     Mfin1 = 0.0d0
                  else
                     Mfin2 = 0.0d0
                  end if
               end if
#elif ( _TRAITEMENT_COLL_MS_RMN_==2) /* Tidal disruption and partial accretion onto BH */
               if (Type1_prv.EQ.TypeEtoile_MS) then
                  call Solve_Roche(Mini1_prv/Mini2_prv,dmin_tidal_disrup)
                  if (dmin_norm.LT.dmin_tidal_disrup) then
                     Vrel_norm_p=0.0d0
                     Theta=0.0d0
                     Mfin1 = 0.0d0
                     Mfin2 = Mini2_prv + CollAccrFrac_MS_on_BH*Mini1_prv
                  end if
               else
                  call Solve_Roche(Mini2_prv/Mini1_prv,dmin_tidal_disrup)
                  if (dmin_norm.LT.dmin_tidal_disrup) then
                     Vrel_norm_p=0.0d0
                     Theta=0.0d0
                     Mfin2 = 0.0d0
                     Mfin1 = Mini1_prv + CollAccrFrac_MS_on_BH*Mini2_prv
                  end if
               end if
#endif
               continue

                                !  %%%%%%%%%%%
                                !  % WD + WD % improbable !
                                !  %%%%%%%%%%%
            else if (
     $              (Type1_prv.EQ.TypeEtoile_WD .AND. Type2_prv.EQ.TypeEtoile_WD)
     $              ) then
               continue         ! collision nulle pour l'instant

                                !  %%%%%%%%%%%
                                !  % NS + WD % improbable !
                                !  %%%%%%%%%%%
            else if (
     $              (Type1_prv.EQ.TypeEtoile_NS .AND. Type2_prv.EQ.TypeEtoile_WD)
     $              .OR.
     $              (Type2_prv.EQ.TypeEtoile_NS .AND. Type1_prv.EQ.TypeEtoile_WD)
     $              ) then
               continue         ! collision nulle pour l'instant

                                !  %%%%%%%%%%%
                                !  % BH + WD % improbable !
                                !  %%%%%%%%%%%
            else if (
     $              (Type1_prv.EQ.TypeEtoile_BH .AND. Type2_prv.EQ.TypeEtoile_WD)
     $              .OR.
     $              (Type2_prv.EQ.TypeEtoile_BH .AND. Type1_prv.EQ.TypeEtoile_WD)
     $              ) then
               continue         ! collision nulle pour l'instant

                                !  %%%%%%%%%%%
                                !  % NS + NS % tres improbable !
                                !  %%%%%%%%%%%
            else if (
     $              (Type1_prv.EQ.TypeEtoile_NS .AND. Type2_prv.EQ.TypeEtoile_NS)
     $              ) then
               continue         ! collision nulle pour l'instant

                                !  %%%%%%%%%%%
                                !  % BH + NS % tres improbable !
                                !  %%%%%%%%%%%
            else if (
     $              (Type1_prv.EQ.TypeEtoile_BH .AND. Type2_prv.EQ.TypeEtoile_NS)
     $              .OR.
     $              (Type2_prv.EQ.TypeEtoile_BH .AND. Type1_prv.EQ.TypeEtoile_NS)
     $              ) then
               continue         ! collision nulle pour l'instant

                                !  %%%%%%%%%%%
                                !  % BH + BH % excessivement improbable ! sauf capture par GW-emission
                                !  %%%%%%%%%%%
            else if (
     $              (Type1_prv.EQ.TypeEtoile_BH .AND. Type2_prv.EQ.TypeEtoile_BH)
     $              ) then
               continue         ! collision nulle pour l'instant

            else                ! type d'etoiles inconnus !!
               write(0,*) '!!! Pas de recette pour les collisions entre ces',
     $              ' types d''etoiles(',Type1_prv,',',Type2_prv,') !!!'
               call Meurt('Collision_Stellaire')
            end if

            if (lEchange12) then
               M_prov = Mfin1
               Mfin1  = Mfin2
               Mfin2  = M_prov
            end if
                                ! On supprime les etoiles de trop faible masse
            if (Mfin1.LT.Mdebris_coll_max) then
               Mfin1=0.0d0
               Vrel_norm_p=0.0d0
               Theta=0.0d0
            end if
            if (Mfin2.LT.Mdebris_coll_max) then
               Mfin2=0.0d0
               Vrel_norm_p=0.0d0
               Theta=0.0d0
            end if
c     
            end

c----------------------------------------------------------------------
      subroutine Init_FichColl_ASC(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Collisions_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100)
 100  format (
     $     '# Proprietes des collisions dans evolution d''amas',/,
     $     '# 1: iPas_evol 2: R_renc ',
     $     '3: iSE1 4: iSE2 5: Type1 6: Type2 7: Mini1 ',
     $     '8: Mini2 9: Vrel_norm 10: b_norm ',
     $     '11: Mfin1 12: Mfin2 13: Vrel_norm_p ',
     $     '14: Theta 15: Pot1 16: Pot2 ',
     $     '17: Age_ini1 18: Age_ini2 19: Age_fin1 20: Age_fin2' )
c
      end
c
#if (_TYPE_COLL_>=10)  /* collisions selon resultats de MES simultaions SPH */
c----------------------------------------------------------------------
      subroutine Lire_Grille_Collisions_SPH(iFich,Type)
c----------------------------------------------------------------------
c     Lecture du fichier ASCII contenant la grille (interpolee) de resultats
c     de simulations SPH de collisions stellaires
c
      include 'EvolAmas_Collisions_inc.f'
      include 'EvolAmas_Collisions_SPH_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich      
      character*3 Type          ! 'BIN' ou 'ASC'
c
c var locales :
c ^^^^^^^^^^^^^
      logical lStat
      integer i1,i3,i4,ig,N1,N3,N4
c
c fonction :
c ^^^^^^^^^^
      integer indice_grille4
c......................................................................
c
c---- Lecture de la grille
c
      if (Type.EQ.'BIN') then ! fichier binaire
         call Lit_Tab_Interp_Res_Coll_BIN1(iFich,lStat)
      else
         call Lit_Tab_Interp_Res_Coll_ASC1(iFich,lStat)
      end if
c
c---- Determination des limites en masse
c
      Minf_grille_SPH = -1.0d0
      Msup_grille_SPH = -1.0d0

      N1=N_interp(1)
      N3=N_interp(3)
      N4=N_interp(4)

      do i1=1,N1
         do i3=1,N3
            do i4=1,N4
               ig = indice_grille4(i1,i1,i3,i4)
               if (tab_stat_interp(ig,1).NE.i_valeur_non_specifiee) then
                  Minf_grille_SPH = 10.0d0**( Xmin_interp(1)+(i1-1)*
     $                 (Xmax_interp(1)-Xmin_interp(1))/(N1-1) )
                  goto 9
               end if
            end do
         end do
      end do
 9    continue

      do i1=N1,1,-1
         do i3=1,N3
            do i4=1,N4
               ig = indice_grille4(i1,i1,i3,i4)
               if (tab_stat_interp(ig,1).NE.i_valeur_non_specifiee) then
                  Msup_grille_SPH = 10.0d0**( Xmin_interp(1)+ (i1-1)*
     $                 (Xmax_interp(1)-Xmin_interp(1))/(N1-1) )
                  goto 19
               end if
            end do
         end do
      end do
 19   continue

#if (_VERB_COLL_>0)
      write(0,100) Minf_grille_SPH,Msup_grille_SPH
 100  format('** Masses stellaires extremes dans la grille SPH : ',
     $     F4.2,' ',F5.1)
#endif
      n_dbrdM_CollSPH_b = 0
      n_dbrdM_CollSPH_h = 0 
      n_extrap_CollSPH  = 0
      n_distant_CollSPH = 0 
      n_inconnu_CollSPH = 0
c
      end
c
c----------------------------------------------------------------------
      subroutine ResColl_selon_Grille_SPH(
     $        Mini1,Mini2,Vrel_norm,b_norm,
     $        Mfin1,Mfin2,Vrel_norm_p,Theta,iStat_Coll_Grille_SPH )
c----------------------------------------------------------------------
c     Determine le resultat d'une collision grace a la grille de
c     resultats de collisions SPH
c
      include 'EvolAmas_Collisions_inc.f'
      include 'EvolAmas_Collisions_SPH_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision          ! conditions initiales
     $     Mini1,Mini2,Vrel_norm,b_norm
      double precision          ! conditions finales
     $     Mfin1,Mfin2,Vrel_norm_p,Theta
      integer                   ! signale problemes et cas speciaux
     $     iStat_Coll_Grille_SPH
c
c constantes :
c ^^^^^^^^^^^^
      double precision
     $     lim_sup_dM,lim_inf_dM,dM_min,lim_sup_dE,lim_inf_dE,dE_min,
     $     dTheta_max,q_max,lim_inf_q
      parameter (
     $     lim_sup_dM=1.01d0, lim_inf_dM=-0.01d0, dM_min=1.0d-6,
     $     lim_sup_dE=1.05d0, lim_inf_dE=-0.05d0, dE_min=1.0d-6,
     $     dTheta_max=1.01d0*rPi,
     $     q_max=3.0d0, lim_inf_q=-0.01d0 )
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      double precision b2v2, dM,dE,dTheta,q,
     $     arg_sqrt
      real VarGrille(Nvar), ResGrille(Nfct)
      integer i,iStat(Nfct)
c......................................................................
c
c---- Initialisations (rencontre 2-corps ponctuels)
c
      Mfin1 = Mini1
      Mfin2 = Mini2
      Vrel_norm_p = Vrel_norm
      Theta = atan2(0.5d0,b_norm*Vrel_norm**2)
c
c---- Passage aux variables de la grille
c
      VarGrille(1)=log10(Mini1)
      VarGrille(2)=log10(Mini2)
      VarGrille(3)=log10(Vrel_norm)
      b2v2=(b_norm*Vrel_norm)**2
      VarGrille(4)=                    ! distance minimale entre points-masse /(R1+R2)
     $     2.0d0*b2v2/(1.0d0+sqrt(1.0d0+4.0d0*b2v2*Vrel_norm**2))
c
c---- Interpolation dans la grille
c
      call Interp_ResColl_Grille(VarGrille,ResGrille,iStat)
                                ! l'interpolation/extrapolation a-t-elle ete possible ?
      iStat_Coll_Grille_SPH = i_interp
      do i=1,Nfct
         if (iStat(i).EQ.i_valeur_non_specifiee) then
            iStat_Coll_Grille_SPH = i_inconnu
            goto 9
         else if (iStat(i).EQ.i_valeur_extrapolee) then
            iStat_Coll_Grille_SPH = i_extrap
         end if
      end do
 9    continue
      
      if (iStat_Coll_Grille_SPH.EQ.i_inconnu) then
         if (VarGrille(4).GT.dmin_limite) then
c
c     Rencontre a grande distance :
c     on admet une interaction purement gravitationnelle
c
            iStat_Coll_Grille_SPH = i_distant
         end if
      else
c
c     "Traduction" du resultat fournit par la grille 
c
c     ResGrille(1) : log10(dM/(M1+M2))
c     ResGrille(2) : dE_orb/E_orb
c     ResGrille(3) : dTheta=theta_coll-theta_grav
c                    ou theta_coll est l'angle de deflection
c                    et theta_grav la valeur pour rencontre 2-corps
c                    ponctuels                                 
c     ResGrille(4) : q=M1fin/M2fin
c
c     Mefions-nous des valeurs fournies par la grille...
c
         dM     = 10.0d0**ResGrille(1)
         dE     = ResGrille(2)
         dTheta = ResGrille(3)
         q      = ResGrille(4)

         if (dM.GT.1.0d0) then
            if (dM.LT.lim_sup_dM) then
               dM=1.0d0
            else
               write(0,*) '!!! Valeur trop grande de dM : ',dM,' !!!'
               call Erreur_ResColl_Grille_SPH(VarGrille,ResGrille,iStat)
            end if
         else if (dM.LT.dM_min) then
            if (dM.GT.lim_inf_dM) then
               dM = 0.0d0       ! Ne pinaillons pas !!
            else
               write(0,*) '!!! Valeur negative de dM : ',dM,' !!!'
               call Erreur_ResColl_Grille_SPH(VarGrille,ResGrille,iStat)
            end if
         end if

         if (dE.GT.1.0d0) then
            if (dE.LT.lim_sup_dE) then
               dE=1.0d0
            else
               write(0,*) '!!! Valeur trop grande de dE : ',dE,' !!!'
               call Erreur_ResColl_Grille_SPH(VarGrille,ResGrille,iStat)
            end if
         else if (dE.LT.dE_min) then
            if (dE.GT.lim_inf_dE) then
               dE=0.0d0
            else
               write(0,*) '!!! Valeur negative de dE : ',dE,' !!!'
               call Erreur_ResColl_Grille_SPH(VarGrille,ResGrille,iStat)
            end if
         end if

         if (abs(dTheta).GT.dTheta_max) then
            write(0,*) '!!! Valeur trop grande de dTheta : ',dTheta,
     $           ' !!!'
            call Erreur_ResColl_Grille_SPH(VarGrille,ResGrille,iStat)
         end if

         if (q.GT.q_max) then
            write(0,*) '!!! Valeur trop grande de q : ',q,' !!!'
               call Erreur_ResColl_Grille_SPH(VarGrille,ResGrille,iStat)
         else if (q.LT.0.0d0) then 
            if (q.GT.lim_inf_q) then
               q=0.0d0
            else
               write(0,*) '!!! Valeur negative de q : ',q,' !!!'
               call Erreur_ResColl_Grille_SPH(VarGrille,ResGrille,iStat)
            end if
         end if
c
c     Angle de deviation total
c
         Theta = Theta + dTheta
c         
c     Determination des masses
c
         Mfin2 = (Mini1+Mini2)*(1.0d0-dM)/(1.0d0+q)
         Mfin1 = q*Mfin2
c
c     Determination de la vitesse relative a l'infini
c
         if (Mfin1.GT.0.0d0) then
#if (_SURVEILLE_SQRT_ > 0)
            arg_sqrt = (1.0d0-dE)*
     $           Mini1*Mini2*(Mfin1+Mfin2)/
     $           (Mfin1*Mfin2*(Mini1+Mini2))
            if (arg_sqrt.LT.0.0d0) then
               write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
               call Meurt('ResColl_selon_Grille_SPH')
            end if
            Vrel_norm_p=Vrel_norm*sqrt(arg_sqrt)
#else
            Vrel_norm_p=Vrel_norm*sqrt( (1.0d0-dE)*
     $           Mini1*Mini2*(Mfin1+Mfin2)/
     $           (Mfin1*Mfin2*(Mini1+Mini2)) )
#endif
         else
            Vrel_norm_p = 0.0d0
         end if

      end if ! (iStat_Coll_Grille_SPH.EQ.i_inconnu)
c
      end
c
c----------------------------------------------------------------------
      subroutine Erreur_ResColl_Grille_SPH(VarGrille,ResGrille,iStat)
c----------------------------------------------------------------------
      include 'EvolAmas_Collisions_inc.f'
      include 'EvolAmas_Collisions_SPH_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      real VarGrille(Nvar), ResGrille(Nfct)
      integer iStat(Nfct)
c......................................................................
c
      write(0,100) VarGrille,ResGrille,iStat
 100  format(
     $     '!!! Erreur d''interpolation dans la grille !!!',/,
     $     '    Variables : ',4(E11.3,' '),/,
     $     '    Resultats : ',4(E11.3,' '),/,
     $     '    Status    : ',4(I3,' ') )
      call Meurt('Erreur_ResColl_Grille_SPH')
c
      end
c
c----------------------------------------------------------------------
      subroutine Traiter_ResColl_SPH_inconnu(
     $     Mini1,Mini2,Vrel_norm,b_norm)
c----------------------------------------------------------------------
      include 'EvolAmas_Collisions_inc.f'
      include 'EvolAmas_Collisions_SPH_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision                  ! conditions initiales de la collision
     $     Mini1,Mini2,Vrel_norm,b_norm ! pour laquelle on n'a pas de prescription
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      logical lFirst, lOK
      save lFirst
      data lFirst /.TRUE./
c......................................................................
c
#if (1==2)
      if (lFirst) then
         call open_append(iFich_CollUnknw,
     $        'MSCollisionsWithUnknownResult.asc',lOK)
         if (.NOT.lOK) then
            write(0,*) '!!! Impossible to open ',
     $           'MSCollisionsWithUnknownResult.asc !!!'
            call Meurt('Traiter_ResColl_SPH_inconnu')
         end if
         call OutputFileHeader(iFich_CollUnknw,
     $        'iPas_Evol Mini1 Mini2 Vrel_norm b_norm')
         lFirst = .FALSE.
      end if

      write(iFich_CollUnknw,110) iPas_Evol,Mini1,Mini2,Vrel_norm,b_norm
 110  format(I11,4(E11.4,' '))
#endif
      write(0,100) iPas_Evol,Mini1,Mini2,Vrel_norm,b_norm
 100  format(
     $     '!!! iPas_Evol=,',I11,
     $     ': Collision de resultat inconnu !!!',/,
     $     '    conditions initiales :',/,
     $     '    M1,M2,vrel_norm,b_norm = ',4(E11.4,' '))
      !!call Meurt('Traiter_ResColl_SPH_inconnu')
c
      end
c
# if (_FIT_LAMBDA_CRIT_MERG_>0)
c----------------------------------------------------------------------
      double precision function
     $     lambda_crit_merger(Mstar1,Mstar2,Vrel_over_Vesc)
c----------------------------------------------------------------------
c     Gives the critical value of lambda = d_min/(R1+R2) for merger 
c     as opposed to fly-by at low velocities. 
c
      implicit none
      include 'EvolAmas_Collisions_param.f'
      include 'EvolAmas_Etoiles_param.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     Mstar1, Mstar2,      ! Stellar masses in M_sun
     $     Vrel_over_Vesc       ! Relative velocity at infinity, in units of sqrt(2*G*(M1+M2)/(R1+R2))
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision logvi,   ! log10(Vrel_over_Vesc)
     $     M1i, M2i, R1tot, R2tot, R1half, R2half, Rratio, nu
#if (_FIT_LAMBDA_CRIT_MERG_==1)
      double precision c(0:8)
      data c /
     $      5.248820e-01, -1.068468e-01,  6.842266e-03,
     $     -2.026948e+00,  5.248322e-01,  6.158136e-04,
     $      1.320273e-01,  5.258783e-01, -4.891898e-03 /
#elif (_FIT_LAMBDA_CRIT_MERG_==2)
      double precision c(0:11)
      data c /
     $      5.061803e-01, -1.630387e-02,  1.129050e-02,
     $     -6.619561e-02, -2.004342e+00, -9.141949e-03,
     $      3.824192e-03,  4.474605e-01,  1.631038e-01,
     $      4.427096e-03, -5.458607e-03,  4.330138e-01 /
#elif (_FIT_LAMBDA_CRIT_MERG_==3)
      double precision c(0:8)
      data c /
     $     -2.718644e+00, -1.676938e+00,  3.118058e-02,
     $     -3.014788e-01,  5.174333e-01, -5.608443e-03,
     $      3.550490e+00,  1.047162e+00, -1.921033e-02 /
#else
      Value of _FIT_LAMBDA_CRIT_MERG_ not valid
#endif
      save c

c
c functions :
c ^^^^^^^^^^^
      double precision RelMR, T_MS, Rhalf_over_R_MS
c......................................................................
c
      if (Vrel_over_Vesc .GT. Nu_Max_Merger) then ! Fit doesn't apply to large velocities
         lambda_crit_merger = -1.0d30 ! This should make it clear...
         return
      end if

      M1i = min(Mstar1,Mstar2)
      M2i = max(Mstar1,Mstar2)
c
c---- The fit is done for the HALF-MASS radii of the stars in the definition of lambda and the escape velocity.
c     One needs to convert to total radii. For consistency, one uses the same radii 
c     as used in the SPH simulations.
c
      !R1tot = RelMR(M1i,0.5d0*T_MS(M1i),TypeEtoile_MS) ! This should be adapted to more general stellar evolution
      !R2tot = RelMR(M2i,0.5d0*T_MS(M2i),TypeEtoile_MS)
      R1tot = RelMR(M1i,-1.d0,TypeEtoile_MS) ! This should be adapted to more general stellar evolution
      R2tot = RelMR(M2i,-1.d0,TypeEtoile_MS)
      R1half = Rhalf_over_R_MS(M1i)*R1tot
      R2half = Rhalf_over_R_MS(M2i)*R2tot
      Rratio = (R1half+R2half) / (R1tot+R2tot)
      nu = Vrel_over_Vesc *sqrt(Rratio)
c
c---- Fit of SPH results by Atakan Gürkan, NU, Nov 2003 
c
#if (_FIT_LAMBDA_CRIT_MERG_==1)
      logvi = log10(nu)
      lambda_crit_merger = c(0) +
     $     c(1)*(M1i/M2i) + c(2)*(M2i) + c(3)*logvi +
     $     c(4)*logvi*(M1i/M2i) + c(5)*logvi*(M2i) + c(6)*logvi*logvi +
     $     c(7)*logvi*logvi*(M1i/M2i) + c(8)*logvi*logvi*(M2i)
#elif (_FIT_LAMBDA_CRIT_MERG_==2)
      logvi = log10(nu)
      lambda_crit_merger = c(0) +
     $     c(1)*M1i + c(2)*M2i + c(3)*(M1i/M2i) + c(4)*logvi +
     $     c(5)*logvi*M1i + c(6)*logvi*M2i + c(7)*logvi*(M1i/M2i) +
     $     c(8)*logvi*logvi + c(9)*logvi*logvi*M1i +
     $     c(10)*logvi*logvi*M2i + c(11)*logvi*logvi*(M1i/M2i)
#elif (_FIT_LAMBDA_CRIT_MERG_==3)
      lambda_crit_merger = c(0) +
     $     c(1)*(M1i/M2i) + c(2)*(M2i) + c(3)*nu +
     $     c(4)*nu*(M1i/M2i) + c(5)*nu*(M2i) +
     $     (c(6) + c(7)*(M1i/M2i) + c(8)*M2i)*nu**-0.18d0
#endif
      lambda_crit_merger = lambda_crit_merger * Rratio
c
      end
#endif /* if (_FIT_LAMBDA_CRIT_MERG_>0) */

#endif /* if (_TYPE_COLL_>=10) */

c----------------------------------------------------------------------
      double precision function
     $     lambda_crit_merger_toy(Mstar1,Mstar2,Vrel_over_Vesc)
c----------------------------------------------------------------------
c     Gives the critical value of lambda = d_min/(R1h+R2h) for merger 
c     as opposed to fly-by at low velocities. R1h and R2h are the
c     HALF-MASS radii of the stars!!
c
      implicit none
      include 'EvolAmas_Collisions_param.f'
      include 'EvolAmas_Etoiles_param.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     Mstar1, Mstar2,      ! Stellar masses in M_sun
     $     Vrel_over_Vesc       ! Relative velocity at infinity, in units of sqrt(2*G*(M1+M2)/(R1h+R2h))
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision logvi,   ! log10(Vrel_over_Vesc)
     $     M1i, M2i, R1tot, R2tot, R1half, R2half, Rratio, nu
c
c functions :
c ^^^^^^^^^^^
      double precision RelMR, T_MS, Rhalf_over_R_MS
c......................................................................
c
      M1i = min(Mstar1,Mstar2)
      M2i = max(Mstar1,Mstar2)
c
c---- The fit is done for the HALF-MASS radii of the stars in the definition of lambda and the escape velocity.
c     One needs to convert to total radii. For consistency, one uses the same radii 
c     as used in the SPH simulations.
c
      !R1tot = RelMR(M1i,0.5d0*T_MS(M1i),TypeEtoile_MS) ! This should be adapted to more general stellar evolution
      !R2tot = RelMR(M2i,0.5d0*T_MS(M2i),TypeEtoile_MS)
      R1tot = RelMR(M1i,-1.d0,TypeEtoile_MS) ! This should be adapted to more general stellar evolution
      R2tot = RelMR(M2i,-1.d0,TypeEtoile_MS)
      R1half = Rhalf_over_R_MS(M1i)*R1tot
      R2half = Rhalf_over_R_MS(M2i)*R2tot
      Rratio = (R1half+R2half) / (R1tot+R2tot)
      nu = Vrel_over_Vesc *sqrt(Rratio)
c
c---- Meta-fit to fit of SPH results by Atakan Gürkan, NU, Nov 2003 (!!!)
c
      lambda_crit_merger_toy =-1.14d0-0.5d0*(log10(nu)-3.00d0)
c
c---- The fit is done for the HALF-MASS radii of the stars in the definition of lambda.
c     One needs to convert to total radii. For consistency, one uses the same radii 
c     as used in the SPH simulations.
c
      lambda_crit_merger_toy = lambda_crit_merger_toy * Rratio
c
      end
#endif /* if (_COLLISIONS_>0) */

