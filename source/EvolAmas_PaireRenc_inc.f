c======================================================================
c     proprietes de la paire de SE participant a une rencontre
c====================================================================== 
c
c---- "designation" des deux SE 
c
      integer iSE_PR(2),iRang_PR(2)
      common /Design_PR_common/ iSE_PR,iRang_PR
c
c---- Potentiel a la position des deux SE
c
      double precision A_PR(2),B_PR(2)
      common /Pot_PR_common/ A_PR,B_PR
c
c---- grandeurs modifiees par la rencontre
c
      double precision
     $     M_PR(2),R_PR(2),T_PR(2),J_PR(2),Tps_PR(2)
      common /Modif_PR_common/
     $     M_PR,R_PR,T_PR,J_PR,Tps_PR

