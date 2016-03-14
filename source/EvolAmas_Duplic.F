#if (_DUPLICATION_>0)
#define _VERB_DUPLIC_ 10
#define _TEST_ENER_DUPLIC_ 1
c======================================================================
c     Routines de duplications des SE pour maintenir une resolution 
c     de bon aloi
c======================================================================

c----------------------------------------------------------------------
      subroutine Init_FichDuplic_ASC(iFich)
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
     $     '# Donnees sur les duplications de particules',/,
     $     '# les lignes commencant par % indiquent le debut ',
     $     'd''un bloc de donnees de duplication.',/,
     $     '# Les lignes qui suivent indiquent le numero de la',
     $     ' SE et de son clone')
c
      end
c
c----------------------------------------------------------------------
      logical function lTest_Duplic()
c----------------------------------------------------------------------
c     indique si il faut dupliquer les SE
c
      include 'EvolAmas_inc.f'
c......................................................................
c
      lTest_Duplic =  ((NbSE-NbSE_perte).LT.0.5d0*NbSE_ini)
c     $     .AND. (Net_par_SE.GE.0.5d0)
      lTest_Duplic = (lTest_Duplic .AND. 2*(NbSE-NbSE_perte).LE.iDimSE)
c
      end
c
c----------------------------------------------------------------------
      subroutine Dupliquer_les_SE
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_Grille_param.f'
      include 'EvolAmas_Grille_common.f'
      include 'EvolAmas_PeriApo_inc.f'
      include 'Flag_common.f'
c
c constantes :
c ^^^^^^^^^^^^
      integer*1 iLibre, iOrig, iClone
      parameter (iLibre=0, iOrig=1, iClone=2)
c
c var locales :
c ^^^^^^^^^^^^^
      integer*1 Clone(iDimSE)
      integer NbSE_new, i_evap, i_evap_tot, iSE, iSE_clone1,
     $     iSE_clone2, iRang_dum, iPG
      double precision M_amas_prec, dM_rel,
     $     Ri,Mi,Ti,Ji,Ei,Tpsi,Meti,DNeti,
     $     A,B, Peri,Apo, Pot
      integer*1 iTeti
c
c fonctions :
c ^^^^^^^^^^^
      double precision R_tronc_mar
c......................................................................
c
      NbSE_new = 2*(NbSE-NbSE_perte)
      iDuplic=iDuplic+1

#if (_VERB_DUPLIC_>0)
      write(0,200) iDuplic,iPas_Evol, NbSE-NbSE_perte, NbSE_new
 200  format(
     $     '++++++++++++++++++++++++++++++++++++++++++++',/,
     $     '++ Duplication ',I4,' iPas_Evol= ',I10 ,' ++',/,
     $     '++        ',I7 ,' --> ',I7 ,' SE          ++',/,
     $     '++++++++++++++++++++++++++++++++++++++++++++')
#endif
#if (_TEST_ENER_DUPLIC_>0)
      call Calc_EnerGlob
      write(0,300) Etot_amas,Ecin_amas,Egrav_amas,Estell_grav_amas
 300  format('++ Duplication: ',
     $     'Etot_amas,Ecin_amas,Egrav_amas,Estell_grav_amas= ',
     $     4(E12.5,' '))
#endif
      if (NbSE_new.GT.iDimSE) then
         write(0,*)
     $        '!!! Trop de particules pour faire une duplication !!!'
         call Meurt('Dupliquer_les_SE')
      end if

#if (_ECRIT_DUPLIC_>0)
      write(iFich_Duplic,100) iPas_Evol, iDuplic
 100  format('% iPas_Evol: ',I10,' iDuplic: ',I4)
#endif
c
c==== Phase 1: clonage simple sans deplacer les particules
c     ---------------------------------------------------- 
c
#if (_VERB_DUPLIC_>5)
      write(0,*)
     $     '++ Duplication: clonage direct...'

#endif
      do iSE=1,NbSE
         if (R_SE(iSE).LT.R_perte) then
            Clone(iSE) = iOrig
         else
            Clone(iSE) = iLibre
         end if
      end do

      iSE_clone1=0
      iSE_clone2=0
      do iSE=1,NbSE

         if (Clone(iSE).EQ.iOrig) then ! SE a dupliquer
            Mi = 0.5d0*M_SE(iSE)
            Ji = J_SE(iSE)
            Ti = T_SE(iSE)
            Ri = R_SE(iSE)
            Tpsi = Tps_SE(iSE)
            Meti = stellar_mass(iSE)
            DNeti = stellar_birthdate(iSE)
            iTeti = stellar_type(iSE)
            Clone(iSE) = iLibre

            iSE_clone1=min(iSE,iSE_clone2+1)
                                ! On cherche la prochaine 'case libre'
                                ! pour y mettre le clone 1
            do while (Clone(iSE_clone1).NE.iLibre)
               iSE_clone1=iSE_clone1+1
            end do
            Clone(iSE_clone1) = iClone
            M_SE(iSE_clone1) = Mi
            J_SE(iSE_clone1) = Ji
            T_SE(iSE_clone1) = Ti
            R_SE(iSE_clone1) = Ri
            Tps_SE(iSE_clone1) = Tpsi
            XXX Met_SE(iSE_clone1) = Meti XXX Should be changed to avoid using *et_SE arrays
            XXX DNet_SE(iSE_clone1) = DNeti
            XXX iTet_SE(iSE_clone1) = iTeti

            iSE_clone2=iSE_clone1+1
                                ! On cherche la prochaine 'case libre'
                                ! pour y mettre le clone 2
            do while (Clone(iSE_clone2).NE.iLibre)
               iSE_clone2=iSE_clone2+1
            end do
            Clone(iSE_clone2) = iClone
            M_SE(iSE_clone2) = Mi
            J_SE(iSE_clone2) = Ji
            T_SE(iSE_clone2) = Ti
            R_SE(iSE_clone2) = Ri
            Tps_SE(iSE_clone2) = Tpsi
            XXX Met_SE(iSE_clone2) = Meti
            XXX DNet_SE(iSE_clone2) = DNeti
            XXX iTet_SE(iSE_clone2) = iTeti

#if (_ECRIT_DUPLIC_>0)
            write(iFich_Duplic,110) iSE, iSE_clone1, iSE_clone2
 110        format(I7,' ',I7,' ',I7)
#endif
         end if
      end do
c
c---- Array indices which do not correspond to 'real' particles
c     but are left empty should be ladelled as such (in order not
c     to be counted as, say, evaporated particles). This means that,
c     when duplications are made, everything is forgotten about 'lost'
c     particles: they are not in the simulation anymore.
c
c     This is probably some kind over overkill as NbSE is adjusted
c     and all 'real' particles have indices below this value
c
      do iSE=1,NbSE
         if (R_SE(iSE).GE.R_perte) then
            R_SE(iSE)=R_vide
            M_SE(iSE)=0.0d0
            J_SE(iSE)=0.0d0
            T_SE(iSE)=0.0d0
            Tps_SE(iSE)=0.0d0
            call stellar_Destroy(iSE)
         end if
      end do
c
c---- Mise a jour des proprietes globales (nombres)
c
      NbSE = NbSE_new 
                                ! C'est un peu genant de devoir mettre les compteurs de pertes a 0 !
      NbSE_perte = 0
      NbSE_evap = 0
      NbSE_dest_coll = 0
      NbSE_dest_dechir = 0
      NbSE_dest_disp = 0 
      NbSE_dest_captGW = 0
                                ! on recalcule les coefficients de conversion d'unites
      call Calc_ConvUnites
#if (_VERB_DUPLIC_>5)
      write(0,220) Net_par_SE
 220  format(
     $     '++ Duplication: Net_par_SE = ',E12.5)
#endif
c
c---- On refait le tri, l'arbre et la grille...
c
#if (_VERB_DUPLIC_>5)
      write(0,*)
     $     '++ Duplication: tri, reconstruction grille et arbre (1)'

#endif
      call TriAmas
      flag_Tri = 1   
      call Construire_Grille_PG(NbSECouchePot)
      flag_Grille = 1   
      call Construire_Arbre_Amas
      flag_Arbre = 1  
#if (_TEST_ENER_DUPLIC_>0)
      call Calc_EnerGlob
      write(0,300) Etot_amas,Ecin_amas,Egrav_amas,Estell_grav_amas
#endif
c
c==== Phase 2: tirage de nouvelles positions orbitales pour "casser" les paires
c     -------------------------------------------------------------------------
c
#if (_VERB_DUPLIC_>5)
      write(0,*) '++ Duplication: tirages orbitaux...'

#endif
      i_evap = 0
      do iSE=1,NbSE
         if (R_SE(iSE).LT.R_perte) then
            Ri=R_SE(iSE)
            Mi=M_SE(iSE)
            Ti=T_SE(iSE)
            Ji=J_SE(iSE)
                                ! on retire la SE de l'arbre
            call Extraire_SE_Arbre(iSE,iRang_dum,A,B) 
            Pot = B + (A+Mi*(1.0d0-FracMasseInt))/Ri
            Ei = Ti + Pot
                                ! on retire la SE de la grille
            iPG = -1
            call Retirer_de_PG(Ri,Ti,Mi,Ji,iPG)  
                                ! on retire la SE de l'amas
            B0 = B0 - Mi/Ri     ! mise a jour de B(R=0)
            M_amas = M_amas - Mi
                                ! calcul apocentre
            call CalcApo_SE_ext2(Mi,Ei,Ji,Apo)
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
 100              format('|| Nb min de SE par cellule atteint NSE_PG(',
     $                 I4,') = ',I3,'< ',I3)
#endif
               end if
               iPG=-1           ! car il ne faut pas re-inserer (au meme endroit) dans la grille
            else 
                                ! calcul du pericentre
               call CalcPeri_SE_ext2(Mi,Ei,Ji,Peri) 
                                ! Tirage aleatoire d'une position orbitale
               call Realiser_TirOrb_direct(
     $              Mi,Ei,Ji,Peri,Apo,A_peri,A_apo, Ri,Ti)
               Ti = 0.5d0*( Ti  + (Ji/Ri)**2 )
               R_SE(iSE) = Ri
               T_SE(iSE) = Ti
                                ! on re-insere la SE dans l'amas
               B0 = B0 + Mi/Ri  ! mise a jour de B(R=0)
               M_amas = M_amas + Mi
            end if
                                ! on re-insere la SE dans l'arbre 
            call Inserer_SE_Arbre(iSE,iRang_dum,A,B)
                                ! on re-insere la SE dans la grille
            call Mettre_dans_PG(Ri,Ti,Mi,Ji,iPG)
               
         end if
      end do
#if (_VERB_DUPLIC_>5)
      write(0,210) i_evap
 210  format(
     $     '++ Duplication: ',I4,' evaporations')
#endif
c
c---- On refait le tri, l'arbre et la grille par acquis de conscience!
c
#if (_VERB_DUPLIC_>5)
      write(0,*)
     $     '++ Duplication: tri, reconstruction grille et arbre (2)'

#endif
      call TriAmas
      flag_Tri = 1   
      call Construire_Grille_PG(NbSECouchePot)
      flag_Grille = 1   
      call Construire_Arbre_Amas
      flag_Arbre = 1 
#if (_TEST_ENER_DUPLIC_>0)
      call Calc_EnerGlob
      write(0,300) Etot_amas,Ecin_amas,Egrav_amas,Estell_grav_amas
#endif
#if (_VERB_DUPLIC_>5)
      write(0,*)
     $     '++ Duplication terminee'

#endif          
c     
      end
c
#endif
