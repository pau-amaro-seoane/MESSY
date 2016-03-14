c======================================================================
c     includes pour ParcoursArbreBin.f
c======================================================================

      implicit none

      include 'Param_Amas.f'
      include 'VarSE.f'
      include 'EvolAmas_Pertes_inc.f'
      include 'EvolAmas_Arbre_common.f'
      include 'EvolAmas_TN_common.f'
      include 'EvolAmas_LimExt_common.f'
c
c---- directions :
c
      integer Droite,Gauche,Haut
      parameter (Droite=1, Gauche=-1, Haut=2)
c
c---- bloc common des proprietes du noeud parcouru
c
      double precision
     $     A_PA,B_PA
      integer
     $     iNoeud_PA,iRang_PA,iSE_PA
      common /Noeud_PA_common/
     $     A_PA,B_PA,
     $     iNoeud_PA,iRang_PA,iSE_PA
c
c--- bloc common de memorisation d'un noeud "prometteur"
c  
      double precision
     $     A_memo,B_memo
      integer
     $     iNoeud_memo,iRang_memo
      common /Noeud_memo_common/
     $     A_memo,B_memo,
     $     iNoeud_memo,iRang_memo
c
c---- implementation pile "FILO"
c
      integer iDimPile
      parameter (iDimPile = 8*iDimArbre)
      integer
     $     NElemPile,iPile(iDimPile)
      common /Pile/
     $     NElemPile,iPile
