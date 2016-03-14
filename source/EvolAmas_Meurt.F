c----------------------------------------------------------------------
      subroutine Meurt(Message)
c----------------------------------------------------------------------
c     routine appelee en cas d'erreur grave constatee
c     sauvegarde un certain nb de chose dans le but d'autopsier
c     le cadavre
c
      include 'EvolAmas_inc.f'
c
c argument(s) :
c ^^^^^^^^^^^^^
      character*(*) Message
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      character*32 NomFich
c......................................................................
c
      write(0,110)Message,iPas_Evol
 110  format(
     $     '!!! UNE ERREUR S''EST PRODUITE (',A,') !!!',/,
     $     '!!! iPas_Evol = ',I10,' !!!')

                                ! Sauvegardes...

                                ! LA PAIRE DE RELAX
      write(NomFich,200)iPas_Evol,'PaireRenc'
      open(iFich_Prov,file=NomFich)
      call Ecrire_PaireRenc_ASC(iFich_Prov)
      close(iFich_Prov)
                                ! LA SE ACTIVE
      write(NomFich,200)iPas_Evol,'SEact'
      open(iFich_Prov,file=NomFich)
      call Ecrire_SEact_ASC(iFich_Prov)
      close(iFich_Prov)
                                           
 200  format('_RIP_%',I10.10,'%',A,'.asc')
                                ! DONNEES GENERALES
      NomSimul='_RIP_'
      call Sauvegarde_EvolAmas(.TRUE.)

      iexit_stat = max(iexit_stat,1)
      call Quitter
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_PaireRenc_ASC(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
c
c argument(s) :
c ^^^^^^^^^^^^^
      integer iFich
c
c functions :
c ^^^^^^^^^^^
      double precision stellar_mass, stellar_radius
      integer*1 stellar_type
c......................................................................
c
      write(iFich,100) 
 100  format(
     $     '# Proprietes de la paire de SE en rencontre',/,
     $     '# 1: Nom_Grandeur 2: Valeur_SE1 3: Valeur_SE2')
      write(iFich,110) 'iSE_PR',iSE_PR
      write(iFich,110) 'irang_PR',iRang_PR
      write(iFich,120) 'iTet_PR',
     $     stellar_type(iSE_PR(1)), stellar_type(iSE_PR(2))
      write(iFich,130) 'M_PR', M_PR
      write(iFich,130) 'R_PR', R_PR
      write(iFich,130) 'T_PR', T_PR
      write(iFich,130) 'J_PR', J_PR
      write(iFich,130) 'Met_PR', 
     $     stellar_mass(iSE_PR(1)), stellar_mass(iSE_PR(2))
      write(iFich,130) 'Ret_PR',
     $     stellar_radius(iSE_PR(1)), stellar_radius(iSE_PR(2))
      write(iFich,130) 'Tps_PR', Tps_PR
      write(iFich,130) 'A_PR', A_PR
      write(iFich,130) 'B_PR', B_PR  
 110  format(A8,2(' ',I8))
 120  format(A8,2(' ',I2))
 130  format(A8,2(' ',E12.5))

C     J'ai utilise:
c     cat _qq |tr ',' '\012' | gawk '{print "      write(iFich,130) ''"$1"'', "$1""}' 

      end

c----------------------------------------------------------------------
      subroutine Ecrire_SEact_ASC(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
c
c argument(s) :
c ^^^^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100) 
 100  format(
     $     '# Proprietes de SE active',/,
     $     '# 1: Nom_Grandeur 2: Valeur') 
      write(iFich,120) 'Num_SE_act', Num_SE_act
      write(iFich,110) 'iSE_act', iSE_act
      write(iFich,110) 'iRang_act', iRang_act
      write(iFich,110) 'iPG_act', iPG_act   
      write(iFich,130) 'Pot_act', Pot_act
      write(iFich,130) 'Peri_act', Peri_act
      write(iFich,130) 'Apo_act', Apo_act                     

 110  format(A10,' ',I8)
 120  format(A10,' ',I2)
 130  format(A10,' ',E17.10)  
c
      end
         

