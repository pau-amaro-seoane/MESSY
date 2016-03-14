      double precision Beta_Dechir_crit
      parameter (Beta_Dechir_crit=0.8d0) ! 0.8

      double precision Beta_Dechir
      common /Dechir_common/ Beta_Dechir

      character*2
     $     SAFE,
     $     RLX_DISRP, COLL_DISRP, KICK_DISRP,
     $     RLX_CROSS, COLL_CROSS, KICK_CROSS,
     $     RLX_GWCPT, COLL_GWCPT, KICK_GWCPT
      parameter (
     $     SAFE='',
     $     RLX_DISRP='RD', COLL_DISRP='CD', KICK_DISRP='KD',
     $     RLX_CROSS='RC', COLL_CROSS='CC', KICK_CROSS='KC',
     $     RLX_GWCPT='RG', COLL_GWCPT='CG', KICK_GWCPT='KD'
     $     )

      integer
     $     iDISRP, iCROSS, iGWCPT,
     $     iRLX,   iCOLL,  iKICK,
     $     iSAFE,
     $     iRLX_DISRP, iCOLL_DISRP, iKICK_DISRP,
     $     iRLX_CROSS, iCOLL_CROSS, iKICK_CROSS,
     $     iRLX_GWCPT, iCOLL_GWCPT, iKICK_GWCPT
      parameter (
     $     iDISRP=1, iCROSS=2,  iGWCPT=3,
     $     iRLX=100, iCOLL=200, iKICK=300,
     $     iSAFE=0,
     $     iRLX_DISRP=101, iCOLL_DISRP=201, iKICK_DISRP=301,
     $     iRLX_CROSS=102, iCOLL_CROSS=202, iKICK_CROSS=302,
     $     iRLX_GWCPT=103, iCOLL_GWCPT=203, iKICK_GWCPT=303
     $     )
