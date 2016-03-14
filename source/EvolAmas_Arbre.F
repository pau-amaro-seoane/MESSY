#undef _DEBUG_PARCOURS_
#define _TYPE_EXTRACTION_ 2 /* 1 : laisse des noeuds vides 2: n'en laisse pas */
c=====================================================================
c     Routines implementant un arbre binaire pour le calcul 
c     du potentiel d'amas spheriques.
c======================================================================
c
c----------------------------------------------------------------------
      subroutine Initialiser_Arbre
c----------------------------------------------------------------------
c
      include 'EvolAmas_Arbre_inc.f'
c
      integer i
c......................................................................
c
      do i=1,NbSE
         iSE_arbre(i)     = 0
         iPere_arbre(i)   = 0
         iFilsG_arbre(i)  = 0
         iFilsD_arbre(i)  = 0
         iNbDecG_arbre(i) = 1   ! correspond en fait aux noeuds isoles (utile ?)
         dA_arbre (i)     = 0.0d0
         dB_arbre (i)     = 0.0d0
         iNoeud_SE(i) = -1      ! pour marquer que les SE sont detachees
                                ! de l arbre (facultatif)
      end do
      B0 = 0.0d0
      NbNoeuds = 0
      NbNoeuds_Libres = 0
c
      end
c
c----------------------------------------------------------------------
      subroutine Tri_Avec_Arbre(lEnlevFant)
c----------------------------------------------------------------------
c     Tri les SE (en R) dans 'iListeRay' grace a l'arbre binaire
c
c     Methode : 'PostOrder Traversal' 
c
c     Reference : Knuth, _The art of computer programming_, Vol I,
c                 Addison-Wesley, 1968, pp 317-318 (Cote 36.1 KNU 1)
c
       include 'EvolAmas_Arbre_inc.f'
c
      logical lEnlevFant        ! si .TRUE. => prendre garde aux
                                ! noeuds fantomes
c
      integer iLst,iNoeud,i,j
      integer Pop_Pile
c......................................................................
c
c %%%  T1  %%%        
      call Init_Pile
      iNoeud = iAdam
      iLst=0
c
c %%%  T2  %%%
 2    continue              
c
       if (iNoeud.NE.0) then
c
c %%%  T3  %%%      
          call Push_Pile(iNoeud)
          iNoeud = iFilsG_arbre(iNoeud)
          goto 2
c
c %%%  T4  %%%  
       else
          if (NElemPile.GT.0) then
             iNoeud = Pop_Pile()
c
c %%%  T5  %%%      
             iLst = iLst+1
             iListeRay(iLst) = iSE_arbre(iNoeud)
             iNoeud_SE(iListeRay(iLst)) = iNoeud
             iNoeud = iFilsD_arbre(iNoeud)
             goto 2
c
          end if
c
       end if
c
c---- A present, enlever les Noeuds-fantomes si demande
c
       if (lEnlevFant) then
          j = 0
          do i=1,iLst
             if (iListeRay(i).LE.NbSE) then
                                ! inserer seulement si il ne s agit pas
                                ! d un "noeud-fantome" contenant une
                                ! SE "nulle" 
                j = j+1
                iListeRay(j) = iListeRay(i)
             end if
          end do

       end if
c
       end  
c
c----------------------------------------------------------------------
      subroutine Construire_Arbre_Amas
c----------------------------------------------------------------------
c     Construction d'un arbre binaire le plus 'plein' possible
c
c     les SE doivent etre deja triees dans 'iListeRay'
c     ================================================
c
c     L'arbre est cree selon le principe recursif :
c
c        Inserer pere, 
c          inserer fils droit  (+recursion), 
c          inserer fils gauche (+recursion)
c
      include 'EvolAmas_Arbre_inc.f'
c
c---- tableaux contenant les sommes des m et m/r
c
      double precision Somme_M,Somme_MSurR,Somme_P
      common /TabSomme/
     $     Somme_M(0:iDimSE), Somme_MSurR(0:iDimSE), Somme_P(0:iDimSE)
c    
      integer i_milieu,i_inf,i_sup,iD,iG,iSE
      integer Pop_Pile
c......................................................................
c
c---- initialisation des donnees
c
      call Init_Pile
c
      call Initialiser_Arbre
c
      call Remplir_Tab_Somme    ! remplit Somme_M, Somme_MSurR
c
c---- construction de l'arbre
c     -----------------------
c
c---- on place les limites de la liste triee sur la pile
c
      iAdam = (NbSE+1)/2
      call Push_Pile(NbSE)
      call Push_Pile(1)
      call Push_Pile(iAdam)
c
c---- et on commence la boucle 'recursive'
c
 1    continue
c
      i_milieu = Pop_Pile()
      i_inf = Pop_Pile()        ! on prend les limites de la
      i_sup = Pop_Pile()        ! "tranche" sur la pile...
      iSE = iListeRay(i_milieu)
      iSE_arbre(i_milieu) = iSE ! on l insere dans l arbre
      iNoeud_SE(iSE) = i_milieu
      if (R_SE(iSE).LT.Revap) then ! SE non evaporee
         B0 = B0 + M_SE(iSE)/R_SE(iSE) ! mise a jour de B(R=0)
      end if
c
                                ! et on calcule les coeff dA et dB du potentiel
                                ! grace aux tableaux de sommes :
      dA_arbre(i_milieu)    = Somme_M(i_milieu)  - Somme_M(i_inf-1)
      dB_arbre(i_milieu)    = Somme_MSurR(i_sup) - Somme_MSurR(i_milieu-1)
                                ! ainsi que le nb de descendants a gauche du noeud
                                ! (y.c. lui-meme) :
      iNbDecG_arbre(i_milieu) = i_milieu - i_inf + 1
c
      if (i_sup.GT.i_milieu) then
                                ! si la tranche de droite n est pas nulle...
         iD = (i_sup+i_milieu+1)/2
                                ! ...on la place sur la pile
         call Push_Pile(i_sup) 
         call Push_Pile(i_milieu+1)
         call Push_Pile(iD)
                                ! ...et on fait les liens avec le noeud pere
         iFilsD_arbre(i_milieu) = iD
         iPere_arbre(iD) = i_milieu

      else                      ! Sinon, on indique l absence de fils
         iFilsD_arbre(i_milieu) = 0
      end if
c
      if (i_milieu.GT.i_inf) then
                                ! pareil pour la tranche gauche
         iG = (i_milieu-1+i_inf)/2
         call Push_Pile(i_milieu-1) 
         call Push_Pile(i_inf)
         call Push_Pile(iG)
c
         iFilsG_arbre(i_milieu) = iG
         iPere_arbre(iG) = i_milieu
      else
         iFilsG_arbre(i_milieu) = 0
      end if
c
      if (NElemPile.GT.0) goto 1
                                ! on continue tant qu il y a 
                                ! des elements sur la pile
c
      NbNoeuds = NbSE
      iPere_arbre(iAdam) = 0
c
      end
c   
c----------------------------------------------------------------------
      subroutine Remplir_Tab_Somme
c----------------------------------------------------------------------
c     Remplissage des tableaux Somme_M et Somme_MSurR pour
c     'Construire_Arbre_Amas'
c
      include 'EvolAmas_Arbre_inc.f'
c   
      double precision Somme_M,Somme_MSurR,Somme_P
      common /TabSomme/
     $     Somme_M(0:iDimSE), Somme_MSurR(0:iDimSE), Somme_P(0:iDimSE)
c
      integer i,iSE
      double precision mSE
c......................................................................
c
      Somme_M(0)     = 0.0d0
      Somme_MsurR(0) = 0.0d0
      Somme_P(0)     = 0.0d0
      do i=1,NbSE
         iSE = iListeRay(i)
         if (R_SE(iSE).LT.Revap) then ! SE liee
            mSE = M_SE(iSE)
         else
            mSE = 0.0d0
         end if
         Somme_M(i)     = Somme_M(i-1) + mSE
         Somme_MsurR(i) = Somme_MsurR(i-1) + mSE/R_SE(iSE)
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine Liberer_Noeud(iNoeud)
c----------------------------------------------------------------------
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iNoeud
c......................................................................
c
      !NbNoeuds = NbNoeuds-1
      NbNoeuds_Libres = NbNoeuds_Libres + 1
      if (NbNoeuds_Libres .GT. MaxNoeuds_Libres) then
         write(0,*) '!!! Debordement du tableau Noeuds_Libres !!!'
         call Meurt('Liberer_Noeud')
      end if
      Noeuds_Libres(NbNoeuds_Libres) = iNoeud
                                ! on vide le noeud
      iPere_arbre(iNoeud) = -1
      iFilsG_arbre(iNoeud) = 0
      iFilsD_arbre(iNoeud) = 0
      iSE_arbre(iNoeud) = 0
      dA_arbre(iNoeud) = 0.0d0
      dB_arbre(iNoeud) = 0.0d0
      iNbDecG_arbre(iNoeud) = 0
c
      end
c
c----------------------------------------------------------------------
      subroutine Demander_Noeud(iNoeud)
c----------------------------------------------------------------------
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iNoeud
c......................................................................
c
      !NbNoeuds = NbNoeuds+1
      if (NbNoeuds_Libres.GT.0) then
         iNoeud = Noeuds_Libres(NbNoeuds_Libres)
         NbNoeuds_Libres = NbNoeuds_Libres - 1
      else
         if (NbNoeuds.GE.iDimArbre) then
            write(0,*) '!!! Plus de place pour un nouveau noeud !!!'
            call Meurt('Demander_Noeud')
         end if
         NbNoeuds = NbNoeuds+1
         iNoeud = NbNoeuds
      end if
                                ! initialisation du noeud
      iPere_arbre(iNoeud) = 0
      iFilsG_arbre(iNoeud) = 0
      iFilsD_arbre(iNoeud) = 0
      iSE_arbre(iNoeud) = 0
      dA_arbre(iNoeud) = 0.0d0
      dB_arbre(iNoeud) = 0.0d0
      iNbDecG_arbre(iNoeud) = 0
c
      end
c
#if (_TYPE_EXTRACTION_==1)
c----------------------------------------------------------------------
      subroutine Extraire_SE_Arbre(iSE,iRang,A,B)
c----------------------------------------------------------------------
c     Extrait la SE iSE de l'arbre sans modifier les grandeurs de
c     l'amas (M_amas,B0, ...).
c     Determine son rang et les coefficients du potentiel, A et B
c     Le Noeud est laisse vide, ou plutot est associe a une SE "phantome"
c
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^                   ! EN ENTREE
      integer iSE
                                ! EN SORTIE
      integer iRang
      double precision A,B
c......................................................................
c
      call Extraire_SE_Arbre_(iSE,iRang,A,B)
c
c---- on fait du noeud courant un noeud vide
c
      R_SE(NbNoeuds+1) = R_SE(iSE) ! Le nb de noeuds dans l'arbre reste constant
      M_SE(NbNoeuds+1) = 0.0d0 ! Par contre, on cree une SE fantome.
                               ! Tant que l'on ne reintroduit pas la SE dans l'arbre,
                               ! il y aura Nb(SE reelles et fantomes) = NbNoeuds+1
      
      iSE_arbre(iNoeud_SE(iSE)) = NbNoeuds+1 ! on fait pointer le noeud sur la SE 'nulle'
      iNoeud_SE(iSE) = -1   ! on 'marque' la SE comme detachee de l'arbre
c
      end
c 
#elif (_TYPE_EXTRACTION_==2)
c----------------------------------------------------------------------
      subroutine Extraire_SE_Arbre(iSE,iRang,A,B)
c----------------------------------------------------------------------
c     Extrait la SE iSE de l'arbre sans modifier les grandeurs de
c     l'amas (M_amas,B0, ...).
c     Determine son rang et les coefficients du potentiel, A et B
c     Le noeud n'est pas laisse vide mais re-attribue. Cette routine
c     ne conduit donc pas a la creation de noeuds vides.
c
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^                   ! EN ENTREE
      integer iSE
                                ! EN SORTIE
      integer iRang
      double precision A,B
c
c var locales :
c ^^^^^^^^^^^^
      integer iFG,iFD, iFils, iPere, iNoeud_vide, iNoeud_Inf
      double precision dA,dB
c......................................................................
c
c---- On commence par extraire la SE en creant un noeud vide
c
      iNoeud_vide=iNoeud_SE(iSE)
      iFG = iFilsG_arbre(iNoeud_vide)
      iFD = iFilsD_arbre(iNoeud_vide)
      iPere = iPere_arbre(iNoeud_vide)
      call Extraire_SE_Arbre_(iSE,iRang,A,B)
c
c---- Si le noeud vide n'a pas de descendant, on peut le supprimer
c
      if (iFG.EQ.0 .AND. iFD.EQ.0) then
         if (iNoeud_vide.EQ.iFilsG_arbre(iPere))
     $        then              ! le noeud est fils gauche
            iFilsG_arbre(iPere) = 0
         else                   ! le noeud est fils droit
            iFilsD_arbre(iPere) = 0
         end if
c
c---- Si le noeud n'a qu'un descendant direct, on met ce dernier a sa place
c
      else if (iFG.EQ.0 .OR. iFD.EQ.0) then
         if (iFG.EQ.0) then     ! seulement fils droit
            iFils = iFD
         else                   ! seulement fils gauche
            iFils = iFG
         end if
         iPere_arbre(iFils) = iPere
         if (iPere.EQ.0) then   ! on supprime la racine !!
            iAdam = iFils
         else
            if (iNoeud_vide.EQ.iFilsG_arbre(iPere))
     $           then           ! le noeud est fils gauche
               iFilsG_arbre(iPere) = iFils
            else                ! le noeud est fils droit
               iFilsD_arbre(iPere) = iFils
            end if
         end if
c
c---- Si il y a deux sous-arbres, il faut faire "remonter" l'element
c     le plus a droite du sous-arbre gauche (on pourrait faire le 
c     contraire, mais bon !...)
c
      else
                                ! on commence par rechercher l'element en question
         iNoeud_inf = iFG
         do while (iFilsD_arbre(iNoeud_inf).GT.0)
            iNoeud_inf = iFilsD_arbre(iNoeud_inf)
         end do
                                ! on le detache de l'arbre et 
                                ! on rattache son eventuel descendant gauche a l'arbre
         iPere = iPere_arbre(iNoeud_inf)
         iFils = iFilsG_arbre(iNoeud_inf)
         if (iPere.EQ.iNoeud_vide) then
            iFilsG_arbre(iPere) = iFils
         else
            iFilsD_arbre(iPere) = iFils
         end if
         if (iFils.NE.0) iPere_arbre(iFils) = iPere
                                ! on remonte jusqu'au noeud vide en metant a jour
                                ! les champs des noeuds rencontres
         dA = M_SE(iSE_arbre(iNoeud_inf))
         dB = M_SE(iSE_arbre(iNoeud_inf))/R_SE(iSE_arbre(iNoeud_inf))
         do while (iPere.NE.iNoeud_vide)
            dB_arbre(iPere) = dB_arbre(iPere) - dB
            iPere = iPere_arbre(iPere)
         end do
         !dA_arbre(iNoeud_vide) = dA_arbre(iNoeud_vide) - dA !####
         !iNbDecG_arbre(iNoeud_vide) = iNbDecG_arbre(iNoeud_vide) - 1
                                ! On met la particule trouvee dans le noeud vide
         iSE_arbre(iNoeud_vide) = iSE_arbre(iNoeud_inf)
         iNoeud_SE(iSE_arbre(iNoeud_vide)) = iNoeud_vide
         dB_arbre(iNoeud_vide) = dB_arbre(iNoeud_vide) + dB
         
         iNoeud_vide = iNoeud_inf
c
      end if
c
c---- On libere le noeud que l'on vient de vider
c
      call Liberer_Noeud(iNoeud_vide)
c
      end
c 
#endif /* if (_TYPE_EXTRACTION_==1) */
c
c----------------------------------------------------------------------
      subroutine Extraire_SE_Arbre_(iSE,iRang,A,B)
c----------------------------------------------------------------------
c     Routine appelee par les 2 versions de Extraire_SE_Arbre
c     Suppression d'une SE dans l'arbre en laissant le noeud.
c
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^                   ! EN ENTREE
      integer iSE
                                ! EN SORTIE
      integer iRang
      double precision A,B
c
c var locales :
c ^^^^^^^^^^^^
      integer iNoeud_vis,iNoeud_pere
      double precision dA,dB
c......................................................................
c
c---- initialisations pour parcours de l'arbre
c
      iNoeud_vis = iNoeud_SE(iSE)
      iNoeud_Pere = iPere_arbre(iNoeud_vis)
      A = -M_ctr -dA_arbre(iNoeud_vis)
      B = B_LimExt -dB_arbre(iNoeud_vis) + M_SE(iSE)/R_SE(iSE)
      iRang = iNbDecG_arbre(iNoeud_vis)
      dA = -M_SE(iSE)
      dB = -M_SE(iSE)/R_SE(iSE)
c
c---- on remonte l'arbre en mettant a jour les coeff. du potentiel
c     et le rang
c
      dA_arbre(iNoeud_vis) = dA_arbre(iNoeud_vis) + dA
      dB_arbre(iNoeud_vis) = dB_arbre(iNoeud_vis) + dB
      iNbDecG_arbre(iNoeud_vis) = iNbDecG_arbre(iNoeud_vis) - 1

      do while (iNoeud_pere.GT.0)
         if (iFilsD_arbre(iNoeud_pere).EQ.iNoeud_vis)
     $        then              ! on remonte depuis la droite
            A = A - dA_arbre(iNoeud_pere)
            iRang = iRang + iNbDecG_arbre(iNoeud_pere)
            dB_arbre(iNoeud_pere) = dB_arbre(iNoeud_pere) + dB
         else                   ! on remonte depuis la gauche
            B = B - dB_arbre(iNoeud_pere)
            dA_arbre(iNoeud_pere) = dA_arbre(iNoeud_pere) + dA
            iNbDecG_arbre(iNoeud_pere) = iNbDecG_arbre(iNoeud_pere) - 1
         end if
         iNoeud_vis = iNoeud_pere
         iNoeud_pere = iPere_arbre(iNoeud_vis)
      end do
      
      iSE_arbre(iNoeud_SE(iSE)) = 0

      Nb_Parcours_Arbre = Nb_Parcours_Arbre+1
c
      end
c 
c----------------------------------------------------------------------
      subroutine Inserer_SE_Arbre(iSE,iRang,A,B)
c----------------------------------------------------------------------
c     Insere la SE iSE dans l'arbre sans modifier les grandeurs de
c     l'amas (M_amas,B0, ...).
c     Determine son rang et les coefficients du potentiel, A et B
c
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^                   ! EN ENTREE
      integer iSE
                                ! EN SORTIE
      integer iRang
      double precision A,B
c
c var locales :
c ^^^^^^^^^^^^
      double precision R_ins,dB_ins, M_ins
      integer iNoeud_cree, iNoeud_vis, iNoeud_prec, Direction
c......................................................................
c
      R_ins = R_SE(iSE)
      if (R_ins.GE.R_perte) then
         M_ins = 0.0d0          ! SE "perdue"
      else
         M_ins = M_SE(iSE)
#if (_LIMITE_EXT_>0)
         if (R_ins.GT.R_LimExt) then
            write(0,*)
     $           '!!! Cannot insert super-star beyond external limit (',
     $           R_ins,R_LimExt,') !!!'
            call Meurt('Inserer_SE_Arbre')
         end if
#endif
      end if
      call Demander_Noeud(iNoeud_cree)
c
      dB_ins = M_ins/R_ins
c
c---- initialisation
c   
      A = -M_ins
      B = B_LimExt
      iRang = 1
      iFilsD_arbre (iNoeud_cree) = 0
      iFilsG_arbre (iNoeud_cree) = 0           
      iSE_arbre    (iNoeud_cree) = iSE
      iNbDecG_arbre(iNoeud_cree) = 1
      dA_arbre     (iNoeud_cree) = M_ins
      dB_arbre     (iNoeud_cree) = dB_ins
      iNoeud_SE(iSE) = iNoeud_cree
      iNoeud_prec = 0
      Direction = 0
c
c---- Traitons specialement le cas d'un arbre auparavant sans element
c
      if (NbNoeuds.EQ.1) then
         iAdam = iNoeud_cree            
         iPere_arbre(iNoeud_cree) = 0
      else
c
c---- Recheche de l'emplacement dans l'arbre
c
         iNoeud_vis = iAdam
         
         do while (iNoeud_vis.GT.0)
            
            iNoeud_prec = iNoeud_vis
            
            if (R_SE(iSE_arbre(iNoeud_vis)).LT.R_ins) then
                                ! \
                                !  \| Descendre a droite
                                !  ~
               dB_arbre(iNoeud_vis) = dB_arbre(iNoeud_vis) + dB_ins
               A = A - dA_arbre(iNoeud_vis)
               iRang = iRang + iNbDecG_arbre(iNoeud_vis)
               iNoeud_vis = iFilsD_arbre(iNoeud_vis)
               Direction = Droite
               
            else                !   /
                                ! |/  Descendre a gauche
                                !  ~
               dA_arbre(iNoeud_vis)      = dA_arbre(iNoeud_vis)
     $              + M_ins
               iNbDecG_arbre(iNoeud_vis) = iNbDecG_arbre(iNoeud_vis)
     $              + 1
               B = B - dB_arbre(iNoeud_vis)
               iNoeud_vis = iFilsG_arbre(iNoeud_vis)
               Direction = Gauche
               
            end if
            
         end do
c
c---- Plus de descendant; il faut inserer !
c
         iPere_arbre(iNoeud_cree) = iNoeud_prec
         if (Direction.EQ.Gauche) then
            iFilsG_arbre(iNoeud_prec) = iNoeud_cree
         else
            iFilsD_arbre(iNoeud_prec) = iNoeud_cree
         end if
c
      end if
c
      Nb_Parcours_Arbre = Nb_Parcours_Arbre+1
c
      end
c      
c----------------------------------------------------------------------
      double precision function Pot_Amas_Arbre_i(iNoeud)
c----------------------------------------------------------------------
c     retourne le potentiel pour la SE du noeud iNoeud
c
      include 'EvolAmas_Arbre_inc.f'
c
      integer iNoeud
c
      integer N_Cour,N_Pere,iSE
      double precision A_noeud,B_noeud
c......................................................................
c
      N_Cour = iNoeud           ! Copie pour ne pas modifier le 
                                ! parametre d appel
      N_Pere = iPere_arbre(N_Cour)
      iSE    = iSE_arbre(N_Cour)
      A_noeud = -M_ctr - dA_arbre(N_Cour)
      if (R_SE(iSE).LT.Revap)   ! SE liee
     $     A_noeud = A_noeud + M_SE(iSE)
      B_noeud = B_LimExt - dB_arbre(N_Cour) 
c
 2    continue
      if (N_Cour.NE.iAdam) then ! On remonte la lignee jusqu a atteindre
                                ! la racine
         if (iFilsD_arbre(N_Pere).EQ.N_Cour) then
                                ! le noeud est descendant a droite
            A_noeud = A_noeud - dA_arbre(N_Pere)
         else
                                ! le noeud est descendant a gauche
            B_noeud = B_noeud - dB_arbre(N_Pere)
         end if
         N_Cour = N_Pere
         N_Pere = iPere_arbre(N_Cour)
c     
         goto 2
      end if
c
      Pot_Amas_Arbre_i = A_Noeud/R_SE(iSE) + B_Noeud
c
      Nb_Parcours_Arbre = Nb_Parcours_Arbre+1
c
      end
c   
c----------------------------------------------------------------------
      double precision function Pot_Amas_Arbre_r(R)
c----------------------------------------------------------------------
c     retourne le potentiel a l'endroit R
c
      include 'EvolAmas_Arbre_inc.f'
      include 'EvolAmas_Compteurs_common.f'
c
      double precision R
c
      integer iRang_pt
      double precision R_equiv_LimExt
      double precision A_pt,B_pt
      common /AB_PotAmasArbre/ A_pt,B_pt ,iRang_pt
                                ! variables pouvant etre utiles 
                                ! => possibilite de partage via common
      integer iNoeud
c......................................................................
c
c---- In case R is beyond the 'artificial' external limit...
c
#if (_LIMITE_EXT_>0)
      if (R.GT.R_LimExt) then   ! one assumes a 'virtual' particle of mass 
                                ! Mext_LimExt at position R_equiv_LimExt
         R_equiv_LimExt = - Mext_LimExt/B_LimExt
         if (R.GT.R_equiv_LimExt) then 
            A_pt = -M_amas -Mext_LimExt ! total mass of cluster (including non-modelled external part)
            B_pt = 0.0d0
         else
            A_pt = -M_amas
            B_pt = B_LimExt
         end if
         iRang_pt = NbSE-NbSE_perte+1
         Pot_Amas_Arbre_r = A_pt/R + B_pt
         return
      end if
#endif

      A_pt = -M_ctr
      B_pt = B_LimExt
      iNoeud = iAdam
      iRang_pt = 0
c
 1    continue
c     
      if (R_SE(iSE_arbre(iNoeud)).LT.R) then
                                ! Si le R de la SE courante est inferieur
                                ! au R recherche, on part dans
                                ! le sous-arbre de DROITE
         A_pt = A_pt - dA_arbre(iNoeud)
         iRang_pt = iRang_pt + iNbDecG_arbre(iNoeud)
                                ! Contribution des SE 'internes'
         iNoeud = iFilsD_arbre(iNoeud)
c     
      else
                                ! On part dans le sous-arbre de GAUCHE
         B_pt = B_pt - dB_arbre(iNoeud)
                                ! Contribution des SE 'externes'
         iNoeud = iFilsG_arbre(iNoeud)
      end if
c     
      if (iNoeud.NE.0) then
         goto 1                 ! Si il n y a pas de noeud-fils, on 
      end if                    ! a termine la descente de l arbre
c
      Pot_Amas_Arbre_r = A_pt/R + B_pt
c
      Nb_Parcours_Arbre = Nb_Parcours_Arbre+1
c     
      end
c
c----------------------------------------------------------------------
      subroutine Determiner_PotRang_Arbre_r(R,A,B,iRang)
c----------------------------------------------------------------------
c     Determine les proprietes du potentiel et le rang en R par parcours
c     de l'arbre
c
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision R        ! en ENTREE
      double precision A,B      ! en SORTIE
      integer iRang
c
c var locales :
c ^^^^^^^^^^^^^
      integer iNoeud
      double precision R_equiv_LimExt
c......................................................................
c
c---- In case R is beyond the 'artificial' external limit...
c
#if (_LIMITE_EXT_>0)
      if (R.GT.R_LimExt) then   ! one assumes a 'virtual' particle of mass 
                                ! Mext_LimExt at position R_equiv_LimExt
         R_equiv_LimExt = - Mext_LimExt/B_LimExt
         if (R.GT.R_equiv_LimExt) then 
            A = -M_amas -Mext_LimExt-M_ctr ! total mass of cluster (including non-modelled external part)
            B = 0.0d0
         else
            A = -M_amas-M_ctr
            B = B_LimExt
         end if
         iRang = NbSE-NbSE_perte+1
         return
      end if
#endif

      A = -M_ctr
      B = B_LimExt
      iNoeud = iAdam
      iRang = 0

      do while (iNoeud.NE.0)

         if (R_SE(iSE_arbre(iNoeud)).LT.R) then
                                ! Si le R de la SE courante est inferieur
                                ! au R recherche, on part dans
                                ! le sous-arbre de DROITE
            A = A - dA_arbre(iNoeud)
            iRang = iRang + iNbDecG_arbre(iNoeud)
                                ! Contribution des SE 'internes'
            iNoeud = iFilsD_arbre(iNoeud)     
         else
                                ! On part dans le sous-arbre de GAUCHE
            B = B - dB_arbre(iNoeud)
                                ! Contribution des SE 'externes'
            iNoeud = iFilsG_arbre(iNoeud)
         end if

      end do
c
      Nb_Parcours_Arbre = Nb_Parcours_Arbre+1
c     
      end
c
c----------------------------------------------------------------------
      double precision function Somme_MsurR_Arbre()
c----------------------------------------------------------------------
c     Se sert de l'arbre pour calculer la somme des M/R
c
      include 'EvolAmas_Arbre_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      integer iNoeud
c......................................................................
c
      Somme_MsurR_Arbre = -B_LimExt
      iNoeud = iAdam

      do while (iNoeud.NE.0)
         Somme_MsurR_Arbre = Somme_MsurR_Arbre + dB_arbre(iNoeud)
         iNoeud = iFilsG_arbre(iNoeud)
      end do
c
      Nb_Parcours_Arbre = Nb_Parcours_Arbre+1
c     
      end
c
c----------------------------------------------------------------------
      subroutine Pot_Arbre(R,U,dUdR)
c----------------------------------------------------------------------
c     fournit le potentiel et sa derivee en R
c
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision R,       ! en entree
     $     U,dUdR               ! en sortie
c
c bloc common particulier :
c ^^^^^^^^^^^^^^^^^^^^^^^^^
      integer iRang_pt
      double precision A_pt,B_pt
      common /AB_PotAmasArbre/ A_pt,B_pt,iRang_pt ! determine par Pot_Amas_Arbre_r
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
c
c fonction :
c ^^^^^^^^^^
      double precision Pot_Amas_Arbre_r
c......................................................................
c
      U = Pot_Amas_Arbre_r(R)
      dUdR = -A_pt/R**2
c
      end
c
c----------------------------------------------------------------------
      subroutine Determiner_OrbCirc_pour_E(Ener,Masse,R_circ,J_circ)
c----------------------------------------------------------------------
      include 'EvolAmas_Arbre_inc.f'
c
c argument(s) :
c ^^^^^^^^^^^^^
      double precision Ener,Masse    ! INPUT
      double precision R_circ,J_circ ! OUTPUT
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      double precision A_visit, B_visit, R_visit, A_int, B_int, R_ext,
     $     arg_sqrt
c
c fonction(s) :
c ^^^^^^^^^^^^^
      logical DescDroite, DescGauche
      double precision Somme_MSurR_Arbre
c......................................................................
c
      A_int = 1.0d0
      B_int = 0.0d0 ! useless
      R_ext = 1.0d40
c
c---- descente dans l'arbre
c
      call Commencer_Parcours_arbre(.TRUE.) ! on commence a la racine...
 
 1    continue

      iSE_PA = iSE_arbre(iNoeud_PA)
      B_visit = B_PA - dB_arbre(iNoeud_PA) + M_SE(iSE_PA)/R_SE(iSE_PA)
      A_visit = A_PA - FracMasseInt*Masse - dA_arbre(iNoeud_PA)
      R_visit = R_SE(iSE_PA)
      !write(0,*) 'xx$$ ',iNoeud_PA,iSE_PA,A_visit,B_visit,R_visit

      if (Ener.GT.0.5d0*A_visit/R_visit+B_visit) then
         A_int=A_visit
         B_int=B_visit
         if (DescDroite()) goto 1 ! Si on ne peut pas descendre
                                  ! a droite, on a termine
      else 
         R_ext=R_visit
         if (DescGauche()) goto 1
      end if
      
      if (A_int.GT.0.0d0) then  ! Rcirc < min(Rnoeud) => aucun noeud trouve
         A_int = -M_ctr -FracMasseInt*Masse
         B_int = -Somme_MSurR_Arbre()
      end if

#if (_LIMITE_EXT_==1)
      if (R_ext.GE.1.0d40) then
         B_visit = B_LimExt
         A_visit = -M_amas  
         R_visit = -(Mext_LimExt/B_LimExt)
         if (Ener.GT.0.5d0*A_visit/R_visit+B_visit) then
            A_int=A_visit
            B_int=B_visit
         else
            R_ext=R_visit
         end if
      end if
#endif

      R_circ = min(R_ext,-A_int/(2.0d0*(B_int-Ener)))
#if (_SURVEILLE_SQRT_ > 0)
      arg_sqrt = -A_int*R_circ
      if (arg_sqrt.LT.0.0d0) then
         write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
         call Meurt('Determiner_OrbCirc_pour_E')
      end if
      J_circ = sqrt(arg_sqrt)
#else
      J_circ = sqrt(-A_int*R_circ)
#endif
c
      return
c
      end
c  
c======================================================================
c     routines de parcours et de recherche dans arbre binaire
c     pour amas a la Henon
c======================================================================
c
c
c----------------------------------------------------------------------
      subroutine Commencer_Parcours_arbre(lAdam)
c----------------------------------------------------------------------
c     Permet de commencer un parcours dans l'arbre...
c     Si lAdam = .TRUE., on commence depuis la racine (Noeud iAdam)
c     Sinon, on commence depuis iNoeud_PA dont il faut 
c     que les proprietes (A_PA, B_PA, iRang_PA) soient a jour
c
      include 'EvolAmas_Arbre_inc.f'
c
c argument(s) :
c ^^^^^^^^^^^^^
      logical lAdam
c......................................................................
c
      if (lAdam) then 
         iNoeud_PA = iAdam
         A_PA = -M_ctr
         B_PA = B_LimExt
         iRang_PA = 0
      else
         B_PA = B_PA + dB_arbre(iNoeud_PA) - M_SE(iSE_PA)/R_SE(iSE_PA)
         A_PA = A_PA + dA_arbre(iNoeud_PA)
         iRang_PA = iRang_PA - iNbDecG_arbre(iNoeud_PA)
      end if
      iSE_PA = iSE_arbre(iNoeud_PA) !### Modif 15/01/02
c
      end
c
c----------------------------------------------------------------------
      subroutine Terminer_Parcours_arbre
c----------------------------------------------------------------------
c     Paufine les valeurs de A_PA, B_PA, iRang_PA
c     pour le noeud courant, une fois le parcours termine
c     Fixe les proprietes de la SE "courante" comme etant celle associee
c     au noeud courant
c
      include 'EvolAmas_Arbre_inc.f'
c......................................................................
c
      iSE_PA = iSE_arbre(iNoeud_PA)
      iRang_PA = iRang_PA + iNbDecG_arbre(iNoeud_PA)
      B_PA = B_PA - dB_arbre(iNoeud_PA) + M_SE(iSE_PA)/R_SE(iSE_PA)
      A_PA = A_PA - dA_arbre(iNoeud_PA)
                                ! les coefficient du pot. ainsi calcules
                                ! sont valables pour R_SE(iSE_PA) <= R < SE suivante
c
      Nb_Parcours_Arbre = Nb_Parcours_Arbre+1
c
      end
c
c----------------------------------------------------------------------
      logical function DescDroite()
c----------------------------------------------------------------------
c     descend dans l'arbre a droite
c
      include 'EvolAmas_Arbre_inc.f'
c......................................................................
c
      if  (iFilsD_arbre(iNoeud_PA).GT.0)
     $     then                 ! Il y a un descendant a droite ...
         iRang_PA = iRang_PA + iNbDecG_arbre(iNoeud_PA)
         A_PA = A_PA - dA_arbre(iNoeud_PA)
         iNoeud_PA = iFilsD_arbre(iNoeud_PA)
         DescDroite = .TRUE.
      else                      ! Il n y a pas de descendant a droite ...
         DescDroite = .FALSE.
      end if
c
      end
c        
c----------------------------------------------------------------------
      logical function DescGauche()
c----------------------------------------------------------------------
c     descend dans l'arbre a gauche
c
      include 'EvolAmas_Arbre_inc.f'
c......................................................................
c
      if  (iFilsG_arbre(iNoeud_PA).GT.0)
     $     then                 ! Il y a un descendant a gauche ...
         B_PA = B_PA - dB_arbre(iNoeud_PA)
         iNoeud_PA = iFilsG_arbre(iNoeud_PA)
         DescGauche = .TRUE.
      else                      ! Il n y a pas de descendant a gauche ...
         DescGauche = .FALSE.
      end if
c
      end
c       
c----------------------------------------------------------------------
      logical function Monter()
c----------------------------------------------------------------------
c     monte dans l'arbre 
c
      include 'EvolAmas_Arbre_inc.f'
c
      integer iP
c......................................................................
c    
      iP = iPere_arbre(iNoeud_PA)
      if (iP.GT.0) then        ! Il y a un pere ...
         if (iFilsG_arbre(iP).EQ.iNoeud_PA)
     $        then              ! on monte depuis la gauche
            iNoeud_PA = iP
            B_PA = B_PA + dB_arbre(iNoeud_PA)
         else                   ! on monte depuis la droite
            iNoeud_PA = iP
            A_PA = A_PA + dA_arbre(iNoeud_PA)
            iRang_PA = iRang_PA - iNbDecG_arbre(iNoeud_PA)
         end if

         Monter = .TRUE.
      else                      ! Il n y a pas de pere ...
         Monter = .FALSE.
      end if
c
      end
c         
c----------------------------------------------------------------------
      logical function DepArbre(iDir)
c----------------------------------------------------------------------
c     se deplace dans l'arbre dans la direction iDir
c
      include 'EvolAmas_Arbre_inc.f'
c
c argument(s) :
c ^^^^^^^^^^^^^
      integer iDir
c
c fonction(s) :
c ^^^^^^^^^^^^^
      logical Monter, DescDroite, DescGauche
c......................................................................
c
      if (iDir.EQ.Droite) then      ! On part a droite
         DepArbre = DescDroite()
      else if (iDir.EQ.Gauche) then ! On part a gauche
         DepArbre = DescGauche()
      else if (iDir.EQ.Haut) then   ! On monte ...
         DepArbre = Monter()
      else
         write(0,*) '!!! Direction non valide dans "DepArbre" !!!'
         stop
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Trouver_SE_Rang(iRang_demande,
     $     iSE_trouve,A_trouve,B_trouve)
c----------------------------------------------------------------------
c     Trouve la SE qui est la (iRang)eme depuis le centre de l'amas
c
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^
                                ! EN ENTREE
      integer iRang_demande
                                ! EN SORTIE
      integer iSE_trouve
      double precision A_trouve,B_trouve
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^     
      integer NouvRang
c
c fonction(s) :
c ^^^^^^^^^^^^^
      logical DescDroite, DescGauche
c......................................................................
c
c---- descente dans l'arbre
c
      call Commencer_Parcours_arbre(.TRUE.) ! on commence a la racine...
 
 1    continue
      NouvRang = iRang_PA + iNbDecG_arbre(iNoeud_PA)

      if (NouvRang.LT.iRang_demande) then
                                ! Descendre a droite
         if (DescDroite()) then
            goto 1
         else ! ANORMAL
            goto 999
         end if

      else if ( (NouvRang.GT.iRang_demande) .OR.
     $        (iSE_arbre(iNoeud_PA).GT.NbSE) ! afin de sauter les noeuds fantomes
     $        ) then         
                                ! Descendre a gauche
         if (DescGauche()) then
            goto 1
         else ! ANORMAL
            goto 999
         end if
            
      end if
c
c---- une SE de rang iRang_demande a ete trouvee
c     on determine ses proprietes
c
      call Terminer_Parcours_arbre ! on fixe les proprietes de la SE trouve
                                   ! (qui est la SE courante)
      iSE_trouve = iSE_PA
      A_trouve = A_PA
      B_trouve = B_PA

#if (_DEBUG_PARCOURS_>1) /* test de l'algorithme */
      if (iRang_PA.NE.iRang_demande) then
         write(0,150) iRang_PA,iRang_demande
 150     format('!!! iRand_trouve(=',I6.6,') <> iRang_demande (=',I6.6,') !!!')
         stop
      end if
#endif
c
      return
c
 999  continue ! On arrive ici si la SE n'a pas pu etre trouvee, ce qui est
               ! forcement ANORMAL
      write(0,200)iRang_demande,iNoeud_PA,iSE_PA,iRang_PA,
     $     iFilsD_arbre(iNoeud_PA),iFilsG_arbre(iNoeud_PA)
 200  format('!!! La SE de rang ',I6,' n''a pas pu etre trouvee !!!',/,
     $     '!!! Dernier noeud visite iNoeud,iSE,iRang,iFilsD,iFilsG = ',
     $     5(I6,' '),' !!!')
      call Meurt('ParcoursArbre')
c
      end
c
c----------------------------------------------------------------------
      subroutine Trouver_SE_R(R_demande,
     $     iSE_trouve,iRang_trouve,A_trouve,B_trouve)
c----------------------------------------------------------------------
c     Trouve la SE de rang i telle que R(i) <= R_demande < R(i+1)
c
c     Attention  : le noeud_PA peut etre un noeud VIDE mais les valeurs
c     ===========  de iRang_PA, A_PA, B_PA correspondent
c                  bien a la SE remplissant la condition ci-dessus   
c                  Pour connaitre le noeud associe a cette SE, il semble 
c                  necesssaire de reparcourir l'arbre avec Trouver_SE_rang !!
c
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^
                                ! EN ENTREE
      double precision R_demande
                                ! EN SORTIE
      integer iSE_trouve,iRang_trouve
      double precision A_trouve,B_trouve
c
c fonction(s) :
c ^^^^^^^^^^^^^
      logical DescDroite, DescGauche
      double precision Somme_MSurR_Arbre
c......................................................................
c
      iNoeud_memo=0
c
c---- descente dans l'arbre
c
      call Commencer_Parcours_arbre(.TRUE.) ! on commence a la racine...
 
 1    continue

      if (R_SE(iSE_arbre(iNoeud_PA)).LT.R_demande) then
                                ! On memorise le noeud et on
                                ! cherche a descendre a droite
         call MemoNoeud
         if (DescDroite()) goto 1 ! Si on ne peut pas descendre
                                  ! a droite, on a termine

      else 
         if (DescGauche()) goto 1
      end if
c
c---- Le noeud recherche est memorise dans Noeud_Memo;
c     on determine ses proprietes
c
      if (iNoeud_memo.GE.1) then ! l`arbre contient un noeud au moins avec
                                 ! R < R_demande
         call RappelleMemoNoeud
         call Terminer_Parcours_arbre ! on fixe les propietes de la SE trouve
                                      ! (qui est la SE courante)
      else                      ! R_demande < min(Rnoeud) => aucun noeud trouve
         iNoeud_PA = 0
         A_PA = -M_ctr
         B_PA = -Somme_MSurR_Arbre()
      end if

      iSE_trouve = iSE_PA
      iRang_trouve = iRang_PA
      A_trouve = A_PA
      B_trouve = B_PA

#if (_DEBUG_PARCOURS_>1) /* test de l'algorithme */
      if (iNoeud_PA.GE.1) then
         if (R_SE(iSE_arbre(iNoeud_PA)).GT.R_demande) then
            write(0,150) R_SE(iSE_arbre(iNoeud_PA)),R_demande
 150        format('!!! R_final(=',E10.3,') < R_cherche (=',E10.3,') !!!')
            call Meurt('ParcoursArbre')
         end if
      end if
#endif
c
      return
c
      end
c
c------------------------------------------------------------
      subroutine MemoNoeud
c------------------------------------------------------------
c
      include 'EvolAmas_Arbre_inc.f'
c............................................................
c
      iNoeud_memo = iNoeud_PA
      iRang_memo = iRang_PA
      A_memo = A_PA
      B_memo = B_PA
c
      end
c
c------------------------------------------------------------
      subroutine RappelleMemoNoeud
c------------------------------------------------------------
c
      include 'EvolAmas_Arbre_inc.f'
c............................................................
c
      iNoeud_PA = iNoeud_memo
      iRang_PA = iRang_memo
      A_PA = A_memo
      B_PA = B_memo
c
      end
c
c----------------------------------------------------------------------
      subroutine RemonterArbre
c----------------------------------------------------------------------
c     Remonte dans l'arbre depuis le noeud iNoeud_PA jusqu'a la 
c     racine. Ce faisant, on determine les proprietes de la SE courante
c     associee a iNoeud_PA
c
      include 'EvolAmas_Arbre_inc.f'
c
c var locales :
c ^^^^^^^^^^^^
      integer iNoeud_vis,iNoeud_pere
c......................................................................
c
c---- initialisation
c
      iSE_PA = iSE_arbre(iNoeud_PA)
      A_PA = -dA_arbre(iNoeud_PA)
      B_PA = -dB_arbre(iNoeud_PA) + M_SE(iSE_PA)/R_SE(iSE_PA)
      iRang_PA = iNbDecG_arbre(iNoeud_PA)

      iNoeud_vis = iNoeud_PA
      iNoeud_pere = iPere_arbre(iNoeud_vis)
c
c---- on remonte l'arbre en mettant a jour les coeff. du potentiel
c     et le rang
c
      do while (iNoeud_pere.GT.0)
         if (iFilsD_arbre(iNoeud_pere).EQ.iNoeud_vis)
     $        then              ! on remonte depuis la droite
            A_PA = A_PA - dA_arbre(iNoeud_pere)
            iRang_PA = iRang_PA + iNbDecG_arbre(iNoeud_pere)
         else                   ! on remonte depuis la gauche
            B_PA = B_PA - dB_arbre(iNoeud_pere)
         end if
         iNoeud_vis = iNoeud_pere
         iNoeud_pere = iPere_arbre(iNoeud_vis)
      end do
c
      Nb_Parcours_Arbre = Nb_Parcours_Arbre+1
c
      end
c
c----------------------------------------------------------------------
      subroutine Modif_Pot_Arbre(iNoeud,dM)
c----------------------------------------------------------------------
c     Modifie les coef dA et dB dans l'arbre de facon correspondante
c     a une variation de masse dM du Noeud iNoeud.
c     Doit aussi modifier la masse M_SE de la particule !!!
c     Les grandeurs globales (B0, M_amas) ne sont pas modifiees.
c
      include 'EvolAmas_Arbre_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iNoeud
      double precision dM
c
c var locales :
c ^^^^^^^^^^^^
      integer iSE,iNoeud_vis,iNoeud_pere
      double precision dA,dB
c......................................................................
c
c---- initialisations pour parcours de l'arbre
c
      iNoeud_vis = iNoeud
      iSE=iSE_Arbre(iNoeud)
      M_SE(iSE)=M_SE(iSE)+dM !!!! Ceci est dangereux mais indispensable !!!!
      iNoeud_Pere = iPere_arbre(iNoeud_vis)
      dA = dM
      dB = dM/R_SE(iSE)
c
c---- on remonte l'arbre en mettant a jour les coeff. du potentiel
c
      dA_arbre(iNoeud_vis) = dA_arbre(iNoeud_vis) + dA
      dB_arbre(iNoeud_vis) = dB_arbre(iNoeud_vis) + dB

      do while (iNoeud_pere.GT.0)
         if (iFilsD_arbre(iNoeud_pere).EQ.iNoeud_vis)
     $        then              ! on remonte depuis la droite
            dB_arbre(iNoeud_pere) = dB_arbre(iNoeud_pere) + dB
         else                   ! on remonte depuis la gauche
            dA_arbre(iNoeud_pere) = dA_arbre(iNoeud_pere) + dA
         end if
         iNoeud_vis = iNoeud_pere
         iNoeud_pere = iPere_arbre(iNoeud_vis)
      end do
c
      Nb_Parcours_Arbre = Nb_Parcours_Arbre+1
c
      end
c 
c======================================================================
c     pile FILO pour arbre binaire d'amas a la Henon
c======================================================================

c----------------------------------------------------------------------
      subroutine Init_Pile
c----------------------------------------------------------------------
c
      include 'EvolAmas_Arbre_inc.f'
c......................................................................
c
      NElemPile = 0
c
      end
c
c----------------------------------------------------------------------
      subroutine Push_Pile(iElem)
c----------------------------------------------------------------------
c
      include 'EvolAmas_Arbre_inc.f'
c
      integer iElem
c......................................................................
c
      if (NElemPile.EQ.iDimPile) then
         write(*,'(/,A)') ' Pile trop petite !! On stoppe.'
         stop
      else
         NElemPile=NElemPile+1
         iPile(NElemPile) = iElem
      end if
c
      end
c
c----------------------------------------------------------------------
      integer function Pop_Pile()
c----------------------------------------------------------------------
c
      include 'EvolAmas_Arbre_inc.f'
c......................................................................
c
      if (NElemPile.EQ.0) then
         Pop_Pile = 0
      else
         Pop_Pile = iPile(NElemPile)
         iPile(NElemPile) = 0
         NElemPile = NElemPile-1
      end if
c
      end
c
c----------------------------------------------------------------------
      logical function Pop_Pile2(iVal)
c----------------------------------------------------------------------
c
      include 'EvolAmas_Arbre_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iVal              ! retourne la valeur au sommet de la pile 
c......................................................................
c
      if (NElemPile.EQ.0) then
         Pop_Pile2 = .FALSE.
      else
         Pop_Pile2 = .TRUE.
         iVal = iPile(NElemPile)
         iPile(NElemPile) = 0
         NElemPile = NElemPile-1
      end if
c
      end
c
