c======================================================================
c     Routines de suivi de SE particulieres lors de l'evolution d'un
c     amas.
c======================================================================
c
      subroutine rien
      end
c

#if (_ECRIT_SUIVI_SE_PART_>0)
c----------------------------------------------------------------------
      subroutine Init_SuiviSE
c----------------------------------------------------------------------
      include 'EvolAmas_SuiviSE_inc.f'
c......................................................................
c
      Nmax_SEsuiv = Dim_Liste_SEsuiv
      N_SEsuiv = 0
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_FichSuiviSE_ASC(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_SuiviSE_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100)
 100  format(
     $     '# Donnees de suivi de SE particulieres',/,
     $     '# 1:iPas_Evol 2:iSE 3:iRang 4:R_old 5:R_new 6:Peri_new',
     $     ' 7:Apo_new 8:T_old 9:T_new 10:J_old 11:J_new' )
c
      end
c
c----------------------------------------------------------------------
      subroutine Lire_Liste_SE_a_suivre(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_SuiviSE_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c var locale :
c ^^^^^^^^^^^^
      character*1 car1
c......................................................................
c
      do while (.TRUE.)
         read(iFich,'(A1)',END=99,ERR=999) car1
         if (car1.NE.'#') then
            if (N_SEsuiv.LT.Nmax_SEsuiv) then
               backspace(iFich)
               N_SEsuiv = N_SEsuiv +1
               read(iFich,*,ERR=999) Liste_SEsuiv(N_SEsuiv)
               write(0,*) '> Lire_Liste_SE_a_suivre: SE a suivre : ',
     $              Liste_SEsuiv(N_SEsuiv)
            else
               write(0,*) '!!! Lire_Liste_SE_a_suivre: liste pleine !!!'
               goto 99
            end if
         end if
      end do
 99   continue
      return
c
 999  continue ! erreur de lecture
      write(0,*) '!!! Erreur de lecture !!!'
      call Meurt('Lire_Liste_SE_a_suivre')
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_SuiviSE(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_SuiviSE_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100) iPas_evol,
     $     iSE_act,iRang_act,R_SE(iSE_act),R_PR(Num_SE_act),
     $     Peri_act,Apo_act,
     $     T_SE(iSE_act),T_PR(Num_SE_act),J_SE(iSE_act),J_PR(Num_SE_act)
 100  format(I10,' ',I6,' ',I6,8(' ',E11.4))
c
      end
c
c----------------------------------------------------------------------
      logical function SEact_A_Suivre() ! Determine si il faut suivre la 
                                        ! SE active
c----------------------------------------------------------------------
      include 'EvolAmas_SuiviSE_inc.f'
c
c constantes :
c ^^^^^^^^^^^^
      double precision ecc_lim
      parameter (ecc_lim=0.995d0)
c
c var locales :
c ^^^^^^^^^^^^^
      integer i
c
c fonctions :
c ^^^^^^^^^^^
      double precision Rmar_TN, RelMR, stellar_mass
      integer*1 stellar_generaltype

c......................................................................
c
      SEact_A_Suivre = .FALSE.
      if (N_SEsuiv.LE.0) return !#####
      if (R_PR(Num_SE_act).GE.R_perte) return
c
c---- On recherche d'abord la SE dans la liste des SE suivies
c
      do i=1,N_SEsuiv
         if (Liste_SEsuiv(i).EQ.iSE_act) then
            SEact_A_Suivre = .TRUE.
            return
         end if
      end do
      return !####
c
c---- On admet eventuellement la SE comme SE a suivre
c     si on a encore de la place dans la liste
c
      if (N_SEsuiv.LT.Nmax_SEsuiv) then
                                ! La condition pour etre suivie
                                ! Ici, on s'interesse aux SE "descendant" tres
                                ! pres du centre
         Mstar = stellar_mass(iSE_act)
         if ( Apo_act.GT.0.0d0 .AND.
     $        Peri_act.LT.2.0d0*Rmar_TN(Mstar,
     $        RelMR(Mstar,-1.0d0,stellar_generaltype(iSE_act)))
c$$$     $        ((Apo_act-Peri_act) .GE. ecc_lim*(Apo_act+Peri_act))
     $        ) then
            SEact_A_Suivre = .TRUE.
            N_SEsuiv = N_SEsuiv + 1
            Liste_SEsuiv(N_SEsuiv) = iSE_act
         end if
      end if
c
      end
c
#endif
