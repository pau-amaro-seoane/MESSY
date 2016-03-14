c----------------------------------------------------------------------
      subroutine Init_FichEvap_ASC(iFich)
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
     $     '# Evaporation data',/,
     $     '# Type_Evap = 1 -> Evap. due to 2-body encounter',/,
     $     '#           = 2 -> Evap. due to potential change ',
     $     '(noticed by Determiner_SE_liees)',/,
     $     '#           = 3 -> Evap. due to potential change ',
     $     '(noticed at moment of pair selection)',/,
     $     '# 1:iPas_Evol 2:Type_Evap 3:iSE 4:M 5:R 6:Ecin 7:Egrav 8:J' )
c
      end
c
c----------------------------------------------------------------------
      subroutine Determiner_SE_liees
c----------------------------------------------------------------------
c     Procedure enlevant de l'amas toutes les particules qui ne sont 
c     plus liees; procede iterativement
c
      include 'EvolAmas_inc.f'
      include 'EvolAmas_Grille_param.f'
      include 'EvolAmas_Grille_common.f'
c
c constantes :
c ^^^^^^^^^^^^
      integer Nb_iter_max
      parameter (Nb_iter_max=20)
      double precision dM_rel_max
      parameter (dM_rel_max=0.25d0) !0.05d0
c
c variables :
c ^^^^^^^^^^^
      integer i_iter, i_evap, i_evap_tot, iSE, iRang_dum, iPG
      double precision M_amas_prec, dM_rel, Ri,Mi,Ti,Ji,A,B,Pot,Apo
c
c fonctions :
c ^^^^^^^^^^^
      double precision R_tronc_mar
c......................................................................
c
      M_amas_prec=M_amas
      i_evap_tot=0

      do i_iter=1,Nb_iter_max
         i_evap=0
c
c==== On prend les SE liees l'une apres l'autre et on regarde si sa distance
c     apo-centrique dans le potentiel courant depasse ou non le rayon de maree
c
         do iSE=1,NbSE
            if (R_SE(iSE).LT.R_perte) then
               Ri=R_SE(iSE)
               Mi=M_SE(iSE)
               Ti=T_SE(iSE)
               Ji=J_SE(iSE)
                                ! on retire la SE de l'arbre
               call Extraire_SE_Arbre(iSE,iRang_dum,A,B) 
               Pot = B + (A+Mi*(1.0d0-FracMasseInt))/Ri
                                ! on retire la SE de la grille
               iPG = -1
               call Retirer_de_PG(Ri,Ti,Mi,Ji,iPG)  
                                ! on retire la SE de l'amas
               B0 = B0 - Mi/Ri  ! mise a jour de B(R=0)
               M_amas = M_amas - Mi
                                ! calcul apocentre
               call CalcApo_SE_ext2(Mi,Ti+Pot,Ji,Apo)
c
c---- La SE peut-elle quitter l'amas ?
c
#if (_TRONC_MAREE_>0)
               if (Apo.LE.0.0d0 .OR. Apo.GT.R_tronc_mar()) then
#else
               if (Apo.LE.0.0d0) then
#endif
                  call PerteSE_par_Evap(iSE,Pot)
                  i_evap=i_evap+1
                                ! on determine si il y a lieu de reconstruire la grille
                  if (NSE_PG(iPG).LT.Nmin_PG) then
                     flag_Grille = -1
#if (_VERBOSITUDE_>5)
                     write(0,100)iPG,NSE_PG(iPG),Nmin_PG
 100                 format('|| Nb min de SE par cellule atteint NSE_PG(',
     $                    I4,') = ',I3,'< ',I3)
#endif    
                  end if
                  iPG=-1        ! car il ne faut pas re-inserer (au meme endroit) dans la grille
               else
                                ! on re-insere la SE dans l'amas
                  B0 = B0 + Mi/Ri ! mise a jour de B(R=0)
                  M_amas = M_amas + Mi
               end if
                                ! on re-insere la SE dans l'arbre 
               call Inserer_SE_Arbre(iSE,iRang_dum,A,B)
                                ! on re-insere la SE dans la grille
               call Mettre_dans_PG(Ri,Ti,Mi,Ji,iPG)
               
            end if
         end do ! Fin de la boucle sur les SE
         i_evap_tot=i_evap_tot+i_evap

         if (i_evap.EQ.0) goto 9 ! convergence atteinte

      end do ! Fin de la boucle d'iterations
c
c---- Pas de convergence
c
      write(0,*) '!!! Convergence du nb de SE liees pas atteinte !!!'
      call meurt('Determiner_SE_liees')

 9    continue
      dM_rel = (M_amas_prec-M_amas)/M_amas_prec
#if (_VERBOSITUDE_>=5)
      write(0,200) dM_rel, i_iter, i_evap_tot
 200  format('> Determiner_SE_liees: ',
     $     'Pertes par evaporation dM_rel =',E12.4,' (',I3,
     $     ' iterations, ',I5,' evaporations)')
#endif
c
c---- Perte de masse catastrophique ?
c
      if (dM_rel .GT. dM_rel_max) then
         write(0,*) '!!! Evaporation catastrophique !!!'  
         call meurt('Determiner_SE_liees')
      end if
c
      end
c

c----------------------------------------------------------------------
      subroutine PerteSE_par_Evap(iSE,Pot)
c----------------------------------------------------------------------
c     Perte d'une SE par evaporation: modifie certaines grandeurs de
c     l'amas et de la SE en consequence
c
      include 'EvolAmas_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iSE
      double precision Pot      ! energie grav specifique de la SE
c......................................................................
c
c---- bilan de masse et d'energie (La masse totale de l'amas est mise
c     a jour lors de la re-introduction des SE)
c
#if (_VERB_PERTES_>0)
      write(0,*) '~~~ Perte de la SE ',iSE,' par Evaporation'
#endif
      dMejec_evap = dMejec_evap + M_SE(iSE)
      dTtot_evap = dTtot_evap - M_SE(iSE)*T_SE(iSE)
      dEtot_evap = dEtot_evap - M_SE(iSE)*(T_SE(iSE)+Pot)
      if (abs(dEtot_evap).GT.100.0d0) then
         write(0,*) '!!! dEtot_evap = ',dEtot_evap,' !!!'
         call Meurt('PerteSE_par_Evap') 
      end if

#if (_ECRIT_EVAP_>0)
      write(iFich_Evap,100) iPas_Evol,2,iSE,M_SE(iSE),
     $     R_SE(iSE),T_SE(iSE),Pot,J_SE(iSE)
 100  format(I10,' ',I2,' ',I8,5(' ',E11.4))
      call flush_MD(iFich_Evap)
#endif

      NbSE_perte = NbSE_perte+1
      NbSE_evap = NbSE_evap+1
      rNbEtoiles = rNbEtoiles-Net_par_SE
c
c---- Modification de certaines proprietes de la particules pour la
c     "marquee" comme perdue
c
      R_SE(iSE) = R_evap
      Tps_SE(iSE) = Tps_perte   ! Necessaire pour le calcul correcte du tps median
                                ! des SE restantes
c
      end
c

