c======================================================================
c     routines de calcul des distances peri- et apo-centrique dans
c     un amas spherique. Avec ou sans potentiel exterieur.
c======================================================================
c
c----------------------------------------------------------------------
      subroutine CalcPeri_SE_ext(Rse,Mse,Ese,Jse,Peri)
c----------------------------------------------------------------------
c     Calcule la distance peri-centrique pour une SE
c     non inclue dans l'arbre de potentiel dont l'orbite passe par R, 
c     la masse est M, l'energie E et le moment cinetique J (valeurs specifiques)
c
      include 'PeriApo_SE2_inc.f'
c     
c arguments :
c ^^^^^^^^^^^
      double precision
     $     Rse,Mse,Ese,Jse,     ! INPUT
     $     Peri                 ! OUTPUT
c
c var locales :
c ^^^^^^^^^^^^^
      integer i,iSE_loc,Rang_parcours
      double precision A_parcours,B_parcours,Corr_B,Rmax
c
c fonction :
c ^^^^^^^^^^
      double precision Vr2_r2
c......................................................................
c
c---- initialisations
c
      iErr_PeriApo = iOK_err_PeriApo
      J22_SE = 0.5d0*Jse**2
      Ener_SE = Ese
      Masse_SE = Mse
      Ray_SE = Rse
      R_inf = -1.0d0
      R_sup = -1.0d0
      Noeud_inf = 0
      Noeud_sup = 0
      Rang_parcours = 0
      A_parcours = -M_ctr
      B_parcours = B_LimExt
      Noeud_loc = iAdam
      Corr_B = 0.0d0
c
c---- recherche d'un encadrement 
c
 1    continue
         iSE_loc = iSE_Arbre(Noeud_loc)
         R_loc = R_SE(iSE_loc)

         if (R_loc.GE.Revap) then
                                ! Essayer de descendre a gauche
            Noeud_loc = iFilsG_arbre(Noeud_loc)
            if (Noeud_loc.GT.0) goto 1
         else
            
            A_loc = A_parcours - dA_arbre(Noeud_loc)
            Rang_loc = Rang_parcours + iNbDecG_arbre(Noeud_loc)
            Corr_B = M_SE(iSE_loc)/R_loc
            B_loc = B_parcours - dB_arbre(Noeud_loc) + Corr_B
                                ! A_loc,B_loc sont les coefficients du potentiel
                                ! keplerien local U(r)=A/r+B pour r>=R_loc
            if (R_loc.GE.Rse .OR. Vr2_r2().GT.0.0d0) then
                                ! Borne superieure possible
                                ! Essayer de descendre a gauche
               Noeud_sup = Noeud_loc
               R_sup = R_loc
               A_sup = A_loc
               B_sup = B_loc
               Rang_sup = Rang_loc
               Noeud_loc = iFilsG_arbre(Noeud_loc)
               if (Noeud_loc.GT.0) then
                  B_parcours = B_loc - Corr_B
                  goto 1
               end if
            else
                                ! Borne inferieure possible
                                ! Essayer de descendre a droite
               Noeud_inf = Noeud_loc
               R_inf = R_loc
               A_inf = A_loc
               B_inf = B_loc
               Rang_inf = Rang_loc
               Noeud_loc = iFilsD_arbre(Noeud_loc)
               if (Noeud_loc.GT.0) then
                  A_parcours = A_loc
                  Rang_parcours = Rang_loc
                  goto 1
               end if
            end if

         end if

      continue
c
c---- test/correction de l'encadrement
c
      if (Noeud_inf.LE.0) then  ! pas de borne inf. trouvee

         A_inf = -M_ctr
         B_inf = B_loc-Corr_B   ! -[somme des Mi/Ri] (sans la SE)
         Rang_inf = 0
         A_loc = A_inf
         B_loc = B_inf
         R_loc = Rse
         do i=1,Max_cR
            R_loc = cR_inf*R_loc
            if (Vr2_r2().LT.0.0d0) goto 19
         end do
         write(0,*) '!!! Impossible de trouver R_inf avec Vrad2 < 0 !!!'
         if (lStop_err_PeriApo) then
            call Meurt('CalcPeri_SE_ext')
         else
            iErr_PeriApo = iRinfNotFound_err_PeriApo
            Mse_err_PeriApo = Mse
            Ese_err_PeriApo = Ese
            Jse_err_PeriApo = Jse
            Rloc_err_PeriApo =  R_loc
            Qloc_err_PeriApo =  1d30
            dQdRloc_err_PeriApo =  1d30
         end if
 19      continue
         R_inf = R_loc

      end if
      if (Noeud_sup.LE.0) then  ! pas de borne sup. trouvee
         A_sup = A_loc          ! c'est -la masse totale de l'amas (sans la SE)
         B_sup = 0.0d0
         R_sup = Rse
         Rang_sup = NbSE+1
      else if (R_sup .GT. Rse) then ! R_sup=min des R_SE > R_apo > Rse > R_peri
         A_sup = -M_ctr
         B_sup = B_loc-Corr_B
         R_sup = Rse
         Rang_sup = 0
      end if
c
c---- calcul du pericentre par resolution de Vrad**2(Rperi)=0 
c     entre Rinf et Rsup
c
      call Resoud_Rperi(Peri)
      if (iErr_PeriApo.EQ.iNegDelta_err_PeriApo) goto 999
                                ! verifions que Peri soit bien
                                ! entre R_inf et min(R-sup,Ray_SE)
      Rmax=min(Ray_SE,R_sup)
      
      if (Peri.GT.Rmax) then
         if (Peri-Rmax.LT.Tol_R*Rmax) then
            Peri = epsR*R_inf+epsR1*Rmax
         else
            write(0,*) '!!! Rmax < Peri : ',Rmax,Peri,' !!!'
            if (lStop_err_PeriApo) then
               call Meurt('CalcPeri_SE_ext')
            else
               iErr_PeriApo = iRsupGTPeri_err_PeriApo
               Mse_err_PeriApo = Mse
               Ese_err_PeriApo = Ese
               Jse_err_PeriApo = Jse
               Rloc_err_PeriApo = R_loc
               Qloc_err_PeriApo =  1d30
               dQdRloc_err_PeriApo =  1d30
            end if
         end if
      else if (Peri.LT.R_inf) then
         if (R_inf-Peri.LT.Tol_R*R_inf) then
            Peri = epsR1*R_inf+epsR*Rmax
         else
            write(0,*) '!!! R_inf > Peri : ',R_inf,Peri,' !!!'
            if (lStop_err_PeriApo) then
               call Meurt('CalcPeri_SE_ext')
            else
               iErr_PeriApo = iRinfLTPeri_err_PeriApo
               Mse_err_PeriApo = Mse
               Ese_err_PeriApo = Ese
               Jse_err_PeriApo = Jse
               Rloc_err_PeriApo =  R_loc
               Qloc_err_PeriApo =  1d30
               dQdRloc_err_PeriApo =  1d30
            end if
         end if
      end if

      A_peri = A_inf            ! proprietes du potentiel de l'amas
      B_peri = B_inf            ! au peri-centre
      iRang_peri = Rang_inf

 999  continue
      lStop_err_PeriApo = .TRUE. ! Resore default behaviour: one stops on errors; this is safer
c
      end
c
c----------------------------------------------------------------------
      subroutine CalcApo_SE_ext(Rse,Mse,Ese,Jse,Apo)
c----------------------------------------------------------------------
c     Calcule la distance apo-centrique pour une SE
c     non inclue dans l'arbre de potentiel dont l'orbite passe par R, 
c     la masse est M, l'energie E et le moment cinetique J (valeurs specifiques)
c
      include 'PeriApo_SE2_inc.f'
      include 'EvolAmas_Pertes_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     Rse,Mse,Ese,Jse,     ! INPUT
     $     Apo                 ! OUTPUT
c
c var locales :
c ^^^^^^^^^^^^^
      integer i,iSE_loc,Rang_parcours
      double precision A_parcours,B_parcours,Corr_B,Rmin, A, DeuxBE,
     $     Q, dQdR
c
c fonction :
c ^^^^^^^^^^
      double precision Vr2_r2
c......................................................................
c
      iErr_PeriApo = iOK_err_PeriApo
c
c---- SE libre ?
c
      if (Ese.GT.0.0d0) then
         Apo = -1.0d0
         goto 999
      end if
c
c---- initialisations
c
      J22_SE = 0.5d0*Jse**2
      Ener_SE = Ese
      Masse_SE = Mse
      Ray_SE = Rse

#if (_LIMITE_EXT_==1)
c
c     test if Apo < R_tilde : iff Q(R_tilde)<0 and dQdR(R_tilde)<0
c     where R_tilde is the position of the "effective external particle"
c
      R_loc = R_LimExt
      A_loc = -M_amas-M_ctr
      B_loc = B_LimExt
      call Calc_Q_dQ(Q,dQdR) 
      if (.NOT.(Q.LT.0.0d0 .AND. dQdR.LT.0.0d0)) then      
         R_loc = -Mext_LimExt/B_LimExt
         A_loc = -M_amas-Mext_LimExt-M_ctr
         B_loc = 0.0d0
         call Calc_Q_dQ(Q,dQdR)
         if (.NOT.(Q.LT.0.0d0 .AND. dQdR.LT.0.0d0)) then ! Apo > R_tilde
            A_apo = A_loc
            B_apo = B_loc
         else ! R_trunc < Apo < R_tilde
            A_apo = -M_amas-M_ctr
            B_apo = B_LimExt
         end if
         A_inf = A_apo
         B_inf = B_apo
         call Resoud_Rapo(Apo)
         if (iErr_PeriApo.EQ.iNegDelta_err_PeriApo) goto 999
         iRang_apo = NbSE-NbSE_perte-1
         goto 999
      endif 
#endif
c
c---- initialisations
c
      R_loc = -1.0d0
      R_sup = -1.0d0
      Noeud_inf = 0
      Noeud_sup = 0
      Rang_parcours = 0
      A_parcours = -M_ctr
      B_parcours = B_LimExt
      Noeud_loc = iAdam
      Corr_B = 0.0d0
c
c---- recherche d'un encadrement 
c
 1    continue
         iSE_loc = iSE_Arbre(Noeud_loc)
         R_loc = R_SE(iSE_loc)

         if (R_loc.GE.Revap) then
                                ! Essayer de descendre a gauche
            Noeud_loc = iFilsG_arbre(Noeud_loc)
            if (Noeud_loc.GT.0) goto 1
         else

            A_loc = A_parcours - dA_arbre(Noeud_loc)
            Rang_loc = Rang_parcours + iNbDecG_arbre(Noeud_loc)
            Corr_B = M_SE(iSE_loc)/R_loc
            B_loc = B_parcours - dB_arbre(Noeud_loc) + Corr_B
                                ! A_loc,B_loc sont les coefficients du potentiel
                                ! keplerien local U(r)=A/r+B pour r>=R_loc
            if (R_loc.LE.Rse .OR. Vr2_r2().GT.0.0d0) then
                                ! Borne inferieure possible
                                ! Essayer de descendre a droite
               Noeud_inf = Noeud_loc
               R_inf = R_loc
               A_inf = A_loc
               B_inf = B_loc
               Rang_inf = Rang_loc
               Noeud_loc = iFilsD_arbre(Noeud_loc)
               if (Noeud_loc.GT.0) then
                  A_parcours = A_loc
                  Rang_parcours = Rang_loc
                  goto 1
               end if
            else
                                ! Borne superieure possible
                                ! Essayer de descendre a gauche
               Noeud_sup = Noeud_loc
               R_sup = R_loc
               A_sup = A_loc
               B_sup = B_loc
               Rang_sup = Rang_loc
               Noeud_loc = iFilsG_arbre(Noeud_loc)
               if (Noeud_loc.GT.0) then
                  B_parcours = B_loc - Corr_B
                  goto 1
               end if
            end if

         end if

      continue
c
c---- test/correction de l'encadrement
c
      if (Noeud_inf.LE.0) then  ! pas de borne inf. trouvee
         A_inf = -M_ctr
         B_inf = B_loc - Corr_B
         Rang_inf = 0
         R_inf = Rse
      else if (Noeud_sup.LE.0) then ! pas de borne sup. trouvee
         A_sup = A_inf
         B_sup = B_inf
         Rang_sup = NbSE+1
         
         A_loc = A_sup
         B_loc = B_sup
         R_loc = Rse
         do i=1,Max_cR
            R_loc = cR_sup*R_loc
            if (Vr2_r2().LT.0.0d0) goto 19
         end do
         write(0,*)
     $        '!!! Impossible de trouver R_sup avec Vrad2 < 0 !!!'
         if (lStop_err_PeriApo) then
            call Meurt('CalcApo_SE_ext')
         else
            iErr_PeriApo = iRsupNotFound_err_PeriApo
            Mse_err_PeriApo = Mse
            Ese_err_PeriApo = Ese
            Jse_err_PeriApo = Jse
            Rloc_err_PeriApo = R_loc
            Qloc_err_PeriApo =  1d30
            dQdRloc_err_PeriApo = 1d30
         end if
 19      continue
         R_sup = R_loc
         
         if (R_inf.LT.Rse) then ! R_inf=max des R_SE < Rse < R_apo
            R_sup = Rse
            Rang_sup = NbSE+1
         end if
         
      end if
c
c---- calcul de l'apocentre par resolution de Vrad**2(Rapo)=0 
c     entre Rinf et Rsup
c
 29   continue
      call Resoud_Rapo(Apo)
      if (iErr_PeriApo.EQ.iNegDelta_err_PeriApo) goto 999
                                ! verifions que Apo soit bien
                                ! entre max(R_inf,Ray_SE) et R_sup
      Rmin=max(Ray_SE,R_inf)
      
      if (Apo.LT.Rmin) then
         if (Rmin-Apo.LT.Tol_R*Rmin) then
            Apo = epsR*Rmin+epsR1*R_sup
         else
            write(0,*) '!!! Rmin > Apo : ',Rmin,Apo,' !!!'
            if (lStop_err_PeriApo) then
               call Meurt('CalcApo_SE_ext')
            else
               iErr_PeriApo = iRinfGTApo_err_PeriApo
               Mse_err_PeriApo = Mse
               Ese_err_PeriApo = Ese
               Jse_err_PeriApo = Jse
               Rloc_err_PeriApo = R_loc
               Qloc_err_PeriApo =  1d30
               dQdRloc_err_PeriApo = 1d30
            end if
         end if
      else if (Apo.GT.R_sup) then
         if (Apo-R_sup.LT.Tol_R*R_sup) then
            Apo = epsR1*Rmin+epsR*R_sup
         else
            write(0,*) '!!! R_sup < Apo : ',R_sup,Apo,' !!!'
            if (lStop_err_PeriApo) then
               call Meurt('CalcApo_SE_ext')
            else
               iErr_PeriApo = iRsupLTApo_err_PeriApo
               Mse_err_PeriApo = Mse
               Ese_err_PeriApo = Ese
               Jse_err_PeriApo = Jse
               Rloc_err_PeriApo = R_loc
               Qloc_err_PeriApo = 1d30
               dQdRloc_err_PeriApo = 1d30
            end if
         end if
      end if

      A_apo = A_inf             ! proprietes du potentiel de l'amas
      B_apo = B_inf             ! a l'apo-centre
      iRang_apo = Rang_inf

 999  continue
      lStop_err_PeriApo = .TRUE. ! Resore default behaviour: one stops on errors; this is safer
c
      end
c     
c----------------------------------------------------------------------
      subroutine CalcPeri_SE_ext2(Mse,Ese,Jse,Peri)
c----------------------------------------------------------------------
c     Calcule la distance peri-centrique pour une particule
c     non inclue dans l'arbre de potentiel dont la masse est M,
c     l'energie E et le moment cinetique J (valeurs specifiques)
c     Cette version ne necessite pas de connaitre un R sur l'orbite
c     a priori. Retourne Peri=-10.0 si il n'existe pas d'orbite
c     pour ces valeurs de E et J
c     
      include 'PeriApo_SE2_inc.f'
      include 'Flag_common.f'
c     
c     arguments :
c     ^^^^^^^^^^^
      double precision
     $     Mse,Ese,Jse,         ! INPUT
     $     Peri                 ! OUTPUT
c     
c     var locales :
c     ^^^^^^^^^^^^^
      integer i,iSE_loc,Rang_parcours
      double precision A_parcours,B_parcours,Corr_B,R1,R2,tiny_dQdR,
     $     Q_loc,dQdR_loc
c     
c     bloc common part :
c     ^^^^^^^^^^^^^^^^^^
      double precision R_circ,E_circ
      common /PeriApo_circ/  R_circ,E_circ
c......................................................................
c
c---- initialisations
c     
      iErr_PeriApo = iOK_err_PeriApo
      J22_SE = 0.5d0*Jse**2
      Ener_SE = Ese
      Masse_SE = Mse
      R_inf = -1.0d0
      R_sup = -1.0d0
      Noeud_inf = 0
      Noeud_sup = 0
      Rang_parcours = 0
      A_parcours = -M_ctr
      B_parcours = B_LimExt
      Noeud_loc = iAdam
      Corr_B = 0.0d0
c     
c---- recherche d'un encadrement 
c     
 1    continue
      iSE_loc = iSE_Arbre(Noeud_loc)
      R_loc = R_SE(iSE_loc)

      if (R_loc.GE.Revap) then
                                ! Essayer de descendre a gauche
         Noeud_loc = iFilsG_arbre(Noeud_loc)
         if (Noeud_loc.GT.0) goto 1
      else

         A_loc = A_parcours - dA_arbre(Noeud_loc)
         Rang_loc = Rang_parcours + iNbDecG_arbre(Noeud_loc)
         Corr_B = M_SE(iSE_Loc)/R_loc
         B_loc = B_parcours - dB_arbre(Noeud_loc) + Corr_B
                                ! A_loc,B_loc sont les coefficients du potentiel
                                ! keplerien local U(r)=A/r+B pour r>=R_loc
         call Calc_Q_dQ(Q_loc,dQdR_loc)
         if (dQdR_loc.LE.0.0d0 .OR. Q_loc.GT.0.0d0) then
                                ! Borne superieure possible
                                ! Essayer de descendre a gauche
            Noeud_sup = Noeud_loc
            R_sup = R_loc
            A_sup = A_loc
            B_sup = B_loc
            Q_sup = Q_loc
            Rang_sup = Rang_loc
            Noeud_loc = iFilsG_arbre(Noeud_loc)
            if (Noeud_loc.GT.0) then
               B_parcours = B_loc - Corr_B
               goto 1
            end if
         else
                                ! Borne inferieure possible
                                ! Essayer de descendre a droite
            Noeud_inf = Noeud_loc
            R_inf = R_loc
            A_inf = A_loc
            B_inf = B_loc
            Q_inf = Q_loc
            Rang_inf = Rang_loc
            Noeud_loc = iFilsD_arbre(Noeud_loc)
            if (Noeud_loc.GT.0) then
               A_parcours = A_loc
               Rang_parcours = Rang_loc
               goto 1
            end if
         end if

      end if

      continue
c     
c---- test/correction de l'encadrement
c     
      if (Noeud_inf.LE.0) then  ! pas de borne inf. trouvee
         A_inf = -M_ctr
         B_inf = B_loc-Corr_B   
         Rang_inf = 0
         A_loc = A_inf
         B_loc = B_inf
         R_loc = R_sup
         do i=1,Max_cR
            R_loc = cR_inf*R_loc
            call Calc_Q_dQ(Q_loc,dQdR_loc)
            if (Q_loc.LT.0.0d0 .AND. dQdR_loc.GT.0.0d0 ) goto 19
         end do
         write(0,*) '!!! Impossible de trouver R_inf !!!'
         write(0,*) '    M_SE,E_SE,J_SE = ',Mse,Ese,Jse
         write(0,*) '    R_loc,Q_loc,dQdR_loc = ',R_loc,Q_loc,dQdR_loc
         if (lStop_err_PeriApo) then
            call Meurt('CalcPeri_SE_ext2')
         else
            iErr_PeriApo = iRinfNotFound_err_PeriApo
            Mse_err_PeriApo = Mse
            Ese_err_PeriApo = Ese
            Jse_err_PeriApo = Jse
            Rloc_err_PeriApo = R_loc
            Qloc_err_PeriApo = Q_loc
            dQdRloc_err_PeriApo = dQdR_loc
         end if
 19      continue
         R_inf = R_loc
         Q_inf = Q_loc
         
      else if (Noeud_sup.LE.0) then ! pas de borne sup. trouvee
         A_sup = A_inf
         B_sup = B_inf   
         Rang_sup = NbSE+1
         A_loc = A_sup
         B_loc = B_sup
         R_loc = R_inf
         do i=1,Max_cR
            R_loc = cR_sup*R_loc
            call Calc_Q_dQ(Q_loc,dQdR_loc)
            if (Q_loc.GT.0.0d0 .OR. dQdR_loc.LT.0.0d0) goto 29
         end do
         write(0,*) '!!! Impossible de trouver R_sup !!!'
         if (lStop_err_PeriApo) then
            call Meurt('CalcPeri_SE_ext2')
         else
            iErr_PeriApo = iRsupNotFound_err_PeriApo
            Mse_err_PeriApo = Mse
            Ese_err_PeriApo = Ese
            Jse_err_PeriApo = Jse
            Rloc_err_PeriApo = R_loc
            Qloc_err_PeriApo = Q_loc
            dQdRloc_err_PeriApo = dQdR_loc
         end if
 29      continue
         R_sup = R_loc
         Q_sup = Q_loc
      end if
c     
c---- Si Q(R_sup) < 0.0, il faut trouver une borne sup avec Q>0
c     si elle n'existe pas cela signifie qu'il n'y a pas d'orbite physique
c     
      if (Q_sup.LE.0.0d0) then
                                ! on cherche en R_sup avec Q>0 en s'approchant du rayon
                                ! circulaire Rc tq dQdR(Rc)=0 par bissection
         tiny_dQdR = eps_dQdR*max(1.0d0,min(abs((Q_inf-Ese)/R_inf),
     $        abs((Q_sup-Ese)/R_sup)))
         R1 = R_inf
         R2 = R_sup
         A_loc = A_inf
         B_loc = B_inf
         do i=1,Max_Bissec
            R_loc=0.5d0*(R1+R2)
            call Calc_Q_dQ(Q_sup,dQdR_loc)
            if (Q_sup.GT.0.0d0) then 
               goto 39
            else if (abs(dQdR_loc).LT.tiny_dQdR)
     $              then        ! on a atteint le rayon circulaire
                                ! mais Vrad2(Rcirc)<0 => orbite non physique
               R_circ = R_loc
               E_circ = Ener_SE-Q_sup 
               Peri = -10.0d0
               call write_QdQ_err_PeriApo
               goto 999
            end if
            if (dQdR_loc.GT.0.0d0) then
               R1=R_loc
            else
               R2=R_loc
            end if
         end do
         write(0,*) '!!! Nb de pas de bissection max effectues',
     $        ' a la recherche de R_sup !!!'
         call write_QdQ_err_PeriApo

         if (lStop_err_PeriApo) then
            call Meurt('CalcPeri_SE_ext2')
         else
            iErr_PeriApo = iRsupMaxBissec_err_PeriApo
            Mse_err_PeriApo = Mse
            Ese_err_PeriApo = Ese
            Jse_err_PeriApo = Jse
            Rloc_err_PeriApo = R_loc
            Qloc_err_PeriApo = Q_loc
            dQdRloc_err_PeriApo = dQdR_loc
         end if

 39      continue               ! on sort ici lorsqu'on a trouve un R_sup convenable
         R_sup = R_loc
      end if
c     
c---- calcul du pericentre par resolution de Vrad**2(Rperi)=0 
c     entre Rinf et Rsup
c     
      call Resoud_Rperi(Peri)
      if (iErr_PeriApo.EQ.iNegDelta_err_PeriApo) goto 999
                                ! verifions que Peri soit bien
                                ! entre R_inf et R_sup
      if (Peri.GT.R_sup) then
         if (Peri-R_sup.LT.Tol_R*R_sup) then
            Peri = epsR*R_inf+epsR1*R_sup
         else
            write(0,*) '!!! R_sup < Peri : ',R_sup,Peri,' !!!'
            if (lStop_err_PeriApo) then
               call Meurt('CalcPeri_SE_ext2')
            else
               iErr_PeriApo = iRsupGTPeri_err_PeriApo
               Mse_err_PeriApo = Mse
               Ese_err_PeriApo = Ese
               Jse_err_PeriApo = Jse
               Rloc_err_PeriApo = R_loc
               Qloc_err_PeriApo = Q_loc
               dQdRloc_err_PeriApo = dQdR_loc
            end if
         end if
      else if (Peri.LT.R_inf) then
         if (R_inf-Peri.LT.Tol_R*R_inf) then
            Peri = epsR1*R_inf+epsR*R_sup
         else
            write(0,*) '!!! R_inf > Peri : ',R_inf,Peri,' !!!'
            if (lStop_err_PeriApo) then
               call Meurt('CalcPeri_SE_ext2')
            else
               iErr_PeriApo = iRinfLTPeri_err_PeriApo
               Mse_err_PeriApo = Mse
               Ese_err_PeriApo = Ese
               Jse_err_PeriApo = Jse
               Rloc_err_PeriApo = R_loc
               Qloc_err_PeriApo = Q_loc
               dQdRloc_err_PeriApo = dQdR_loc
            end if
         end if
      end if

      A_peri = A_inf            ! proprietes du potentiel de l'amas
      B_peri = B_inf            ! au peri-centre
      iRang_peri = Rang_inf

 999  continue
      lStop_err_PeriApo = .TRUE. ! Resore default behaviour: one stops on errors; this is safer
c     
      end
c     
c----------------------------------------------------------------------
      subroutine CalcApo_SE_ext2(Mse,Ese,Jse,Apo)
c----------------------------------------------------------------------
c     Calcule la distance apo-centrique pour une particule
c     non inclue dans l'arbre de potentiel dont la masse est M,
c     l'energie E et le moment cinetique J (valeurs specifiques)
c     Cette version ne necessite pas de connaitre un R sur l'orbite
c     a priori. Retourne Apo=-10.0 si il n'existe pas d'orbite
c     pour ces valeurs de E et J
c     
      include 'PeriApo_SE2_inc.f'
      include 'EvolAmas_Pertes_inc.f'
      include 'Flag_common.f'
c     
c     arguments :
c     ^^^^^^^^^^^
      double precision
     $     Mse,Ese,Jse,         ! INPUT
     $     Apo                  ! OUTPUT
c     
c     var locales :
c     ^^^^^^^^^^^^^
      integer i,iSE_loc,Rang_parcours
      double precision A_parcours,B_parcours,Corr_B,R1,R2,tiny_dQdR,
     $     Q_loc,dQdR_loc,
     $     Q, dQdR
c     
c     bloc common part :
c     ^^^^^^^^^^^^^^^^^^
      double precision R_circ,E_circ
c......................................................................
c
      iErr_PeriApo = iOK_err_PeriApo
c     
c---- SE libre ?
c     
      if (Ese.GT.0.0d0) then
         Apo = -1.0d0
         goto 999
      end if
c
c---- initialisations
c
      J22_SE = 0.5d0*Jse**2
      Ener_SE = Ese
      Masse_SE = Mse

#if (_LIMITE_EXT_==1)
c
c     test if Apo < R_tilde : iff Q(R_tilde)<0 and dQdR(R_tilde)<0
c     where R_tilde is the position of the "effective external particle"
c
      R_loc = R_LimExt
      A_loc = -M_amas-M_ctr
      B_loc = B_LimExt
      call Calc_Q_dQ(Q,dQdR) 
      if (.NOT.(Q.LT.0.0d0 .AND. dQdR.LT.0.0d0)) then      
         R_loc = -Mext_LimExt/B_LimExt
         A_loc = -M_amas-Mext_LimExt-M_ctr
         B_loc = 0.0d0
         call Calc_Q_dQ(Q,dQdR)
         if (.NOT.(Q.LT.0.0d0 .AND. dQdR.LT.0.0d0)) then ! Apo > R_tilde
            A_apo = A_loc
            B_apo = B_loc
         else ! R_trunc < Apo < R_tilde
            A_apo = -M_amas-M_ctr
            B_apo = B_LimExt
         end if
         A_inf = A_apo
         B_inf = B_apo
         call Resoud_Rapo(Apo)
         if (iErr_PeriApo.EQ.iNegDelta_err_PeriApo) goto 999
         iRang_apo = NbSE-NbSE_perte-1
         goto 999
      endif 
#endif

c     
c---- initialisations
c     
      R_inf = -1.0d0
      R_sup = -1.0d0
      Noeud_inf = 0
      Noeud_sup = 0
      Rang_parcours = 0
      A_parcours = -M_ctr
      B_parcours = B_LimExt
      Noeud_loc = iAdam
      Corr_B = 0.0d0
c     
c---- recherche d'un encadrement 
c     
 1    continue

      iSE_loc = iSE_Arbre(Noeud_loc)
      R_loc = R_SE(iSE_loc)

      if (R_loc.GE.Revap) then
                                ! Essayer de descendre a gauche
         Noeud_loc = iFilsG_arbre(Noeud_loc)
         if (Noeud_loc.GT.0) goto 1
      else

         A_loc = A_parcours - dA_arbre(Noeud_loc)
         Rang_loc = Rang_parcours + iNbDecG_arbre(Noeud_loc)
         Corr_B = M_SE(iSE_loc)/R_loc
         B_loc = B_parcours - dB_arbre(Noeud_loc) + Corr_B
                                ! A_loc,B_loc sont les coefficients du potentiel
                                ! keplerien local U(r)=A/r+B pour r>=R_loc
         call Calc_Q_dQ(Q_loc,dQdR_loc)
         if (dQdR_loc.GE.0.0d0 .OR. Q_loc.GT.0.0d0) then
                                ! Borne inferieure possible
                                ! Essayer de descendre a droite
            Noeud_inf = Noeud_loc
            R_inf = R_loc
            A_inf = A_loc
            B_inf = B_loc
            Q_inf = Q_loc
            Rang_inf = Rang_loc
            Noeud_loc = iFilsD_arbre(Noeud_loc)
            if (Noeud_loc.GT.0) then
               A_parcours = A_loc
               Rang_parcours = Rang_loc
               goto 1
            end if
         else
                                ! Borne superieure possible
                                ! Essayer de descendre a gauche
            Noeud_sup = Noeud_loc
            R_sup = R_loc
            A_sup = A_loc
            B_sup = B_loc
            Q_sup = Q_loc
            Rang_sup = Rang_loc
            Noeud_loc = iFilsG_arbre(Noeud_loc)
            if (Noeud_loc.GT.0) then
               B_parcours = B_loc - Corr_B
               goto 1
            end if
         end if

      end if

      continue
c     
c---- test/correction de l'encadrement
c     
      if (Noeud_inf.LE.0) then  ! pas de borne inf. trouvee
         A_inf = -M_ctr
         B_inf = B_loc-Corr_B   
         Rang_inf = 0
         A_loc = A_inf
         B_loc = B_inf
         R_loc = R_sup
         do i=1,Max_cR
            R_loc = cR_inf*R_loc
            call Calc_Q_dQ(Q_loc,dQdR_loc)
            if (Q_loc.GT.0.0d0 .OR. dQdR_loc.GT.0.0d0 ) goto 19
         end do
         write(0,*) '!!! Impossible de trouver R_inf !!!'
         if (lStop_err_PeriApo) then
            call Meurt('CalcApo_SE_ext2')
         else
            iErr_PeriApo = iRinfNotFound_err_PeriApo
            Mse_err_PeriApo = Mse
            Ese_err_PeriApo = Ese
            Jse_err_PeriApo = Jse
            Rloc_err_PeriApo = R_loc
            Qloc_err_PeriApo = Q_loc
            dQdRloc_err_PeriApo = dQdR_loc
         end if
 19      continue
         R_inf = R_loc
         Q_inf = Q_loc
         
      else if (Noeud_sup.LE.0) then ! pas de borne sup. trouvee
         A_sup = A_inf
         B_sup = B_inf   
         Rang_sup = NbSE+1
         A_loc = A_sup
         B_loc = B_sup
         R_loc = R_inf
         do i=1,Max_cR
            R_loc = cR_sup*R_loc
            call Calc_Q_dQ(Q_loc,dQdR_loc)
            if (Q_loc.LT.0.0d0 .AND. dQdR_loc.LT.0.0d0) goto 29
         end do
         write(0,*) '!!! Impossible de trouver R_sup !!!'
         if (lStop_err_PeriApo) then
            call Meurt('CalcApo_SE_ext2')
         else
            iErr_PeriApo = iRsupNotFound_err_PeriApo
            Mse_err_PeriApo = Mse
            Ese_err_PeriApo = Ese
            Jse_err_PeriApo = Jse
            Rloc_err_PeriApo = R_loc
            Qloc_err_PeriApo = Q_loc
            dQdRloc_err_PeriApo = dQdR_loc
         end if
 29      continue
         R_sup = R_loc
         Q_sup = Q_loc
      end if
c     
c---- Si Q(R_inf) < 0.0, il faut trouver une borne inf avec Q>0
c     si elle n'existe pas cela signifie qu'il n'y a pas d'orbite physique
c     
      if (Q_inf.LE.0.0d0) then
                                ! on cherche en R_inf avec Q>0 en s'approchant du rayon
                                ! circulaire Rc tq dQdR(Rc)=0 par bissection
         tiny_dQdR = eps_dQdR*max(1.0d0,min(abs((Q_inf-Ese)/R_inf),
     $        abs((Q_sup-Ese)/R_sup)))

         R1 = R_inf
         R2 = R_sup
         A_loc = A_inf
         B_loc = B_inf
         do i=1,Max_Bissec
            R_loc=0.5d0*(R1+R2)
            call Calc_Q_dQ(Q_sup,dQdR_loc)
            if (Q_sup.GT.0.0d0) then 
               goto 39
            else if (abs(dQdR_loc).LT.tiny_dQdR)
     $              then        ! on a atteint le rayon circulaire
                                ! mais Vrad2(Rcirc)<0 => orbite non physique
               R_circ = R_loc
               E_circ = Ener_SE-Q_sup 
               !write(0,*) '>>>> ',R_loc,dQdR_loc,tiny_dQdR
               Apo = -10.0d0
               call write_QdQ_err_PeriApo
               return
            end if
            if (dQdR_loc.GT.0.0d0) then
               R1=R_loc
            else
               R2=R_loc
            end if
         end do
         write(0,*) '!!! Nb de pas de bissection max effectues',
     $        ' a la recherche de R_sup bbb !!!'
         call write_QdQ_err_PeriApo
         if (lStop_err_PeriApo) then
            call Meurt('CalcApo_SE_ext2')
         else
            iErr_PeriApo = iRsupMaxBissec_err_PeriApo
            Mse_err_PeriApo = Mse
            Ese_err_PeriApo = Ese
            Jse_err_PeriApo = Jse
            Rloc_err_PeriApo = R_loc
            Qloc_err_PeriApo = Q_loc
            dQdRloc_err_PeriApo = dQdR_loc
         end if

 39      continue               ! on sort ici lorsqu'on a trouve un R_sup convenable
         R_inf = R_loc
      end if
c     
c---- calcul de l'apocentre par resolution de Vrad**2(Rapo)=0 
c     entre Rinf et Rsup
c     
      call Resoud_Rapo(Apo)
      if (iErr_PeriApo.EQ.iNegDelta_err_PeriApo) goto 999
                                ! verifions que Apo soit bien
                                ! entre R_inf et R_sup
      if (Apo.LT.R_inf) then
         if (R_inf-Apo.LT.Tol_R*R_inf) then
            Apo = epsR*R_inf+epsR1*R_sup
         else
            write(0,*) '!!! R_inf > Apo : ',R_inf,Apo,' !!!'
            if (lStop_err_PeriApo) then
               call write_QdQ_err_PeriApo
               call Meurt('CalcApo_SE_ext2')
            else
               iErr_PeriApo = iRinfGTApo_err_PeriApo
               Mse_err_PeriApo = Mse
               Ese_err_PeriApo = Ese
               Jse_err_PeriApo = Jse
               Rloc_err_PeriApo = R_loc
               Qloc_err_PeriApo = Q_loc
               dQdRloc_err_PeriApo = dQdR_loc
            end if
         end if
      else if (Apo.GT.R_sup) then
         if (Apo-R_sup.LT.Tol_R*R_sup) then
            Apo = epsR1*R_inf+epsR*R_sup
         else
            write(0,*) '!!! R_sup < Apo : ',R_sup,Apo,' !!!'
            if (lStop_err_PeriApo) then
               call write_QdQ_err_PeriApo
               call Meurt('CalcApo_SE_ext2')
            else
               iErr_PeriApo = iRsupLTApo_err_PeriApo
               Mse_err_PeriApo = Mse
               Ese_err_PeriApo = Ese
               Jse_err_PeriApo = Jse
               Rloc_err_PeriApo = R_loc
               Qloc_err_PeriApo = Q_loc
               dQdRloc_err_PeriApo = dQdR_loc
            end if
         end if
      end if

      A_apo = A_inf             ! proprietes du potentiel de l'amas
      B_apo = B_inf             ! au apo-centre
      iRang_apo = Rang_inf

 999  continue
      lStop_err_PeriApo = .TRUE. ! Resore default behaviour: one stops on errors; this is safer
c     
      end
c
c----------------------------------------------------------------------
      double precision function Vr2_r2()
c----------------------------------------------------------------------
      include 'PeriApo_SE2_inc.f'
c......................................................................
c
      Vr2_r2 = -J22_SE -(A_loc - FracMasseInt*Masse_SE +
     $     (B_loc-Ener_SE)*R_loc)*R_loc
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_Q_dQ(Q,dQdR)
c----------------------------------------------------------------------
c     calcule Q=1/2*Vrad**2 et sa derivee radiale au noeud local
c
      include 'PeriApo_SE2_inc.f'
      include 'Flag_common.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision Q,dQdR
c
c var locales :
c ^^^^^^^^^^^^^
      double precision UnSurR,ASurR,J22R2
c......................................................................
c
      UnSurR = 1.0d0/R_loc
      J22R2  = J22_SE*UnSurR**2
      ASurR  = (A_loc-FracMasseInt*Masse_SE)*UnSurR
      Q = Ener_SE-B_loc-ASurR-J22R2
      dQdR = (ASurR+2.0d0*J22R2)*UnSurR
c
      end
c
c----------------------------------------------------------------------
      subroutine Resoud_Rperi(Peri)
c----------------------------------------------------------------------
c     il s'agit de resoudre une equation polynomiale de degre 2 :
c     1/2*(J/r)**2 + A/r + 2*(B_E) = 0
c
      include 'PeriApo_SE2_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision Peri
c
c var locales :
c ^^^^^^^^^^^^^
      double precision A,delta
c......................................................................
c
      A = A_inf -FracMasseInt*Masse_SE
      delta = A*A-4.0d0*J22_SE*(B_inf-Ener_SE)
      
      if (delta.LE.0.0d0) 
     $     then                 ! peut-etre une orbite circulaire...
         if (delta.GE.Tol_circul*A*A) then
           Peri = -2.0d0*J22_SE/A
         else                   ! Il y a une bulle ...
            write(0,*)'!!!  delta,A,B,J2,E = ',
     $           delta,A,B_inf,J22_SE,Ener_SE,' !!!'
            if (lStop_err_PeriApo) then
               call Meurt('Resoud_Rperi')
            else
               iErr_PeriApo = iNegDelta_err_PeriApo
               Ese_err_PeriApo = Ener_SE
               Jse_err_PeriApo = sqrt(2.0d0*J22_SE)
               A_err_PeriApo = A
               B_err_PeriApo = B_inf
            end if
         end if
      
      else
        Peri = 2.0d0*J22_SE / (-A+sqrt(delta))
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Resoud_Rapo(Apo)
c----------------------------------------------------------------------
c     il s'agit de resoudre une equation polynomiale de degre 2 :
c     1/2*(J/r)**2 + A/r + 2*(B_E) = 0
c
      include 'PeriApo_SE2_inc.f'
      include 'EvolAmas_Compteurs_common.f' !###########
c
c argument :
c ^^^^^^^^^^
      double precision Apo
c
c var locales :
c ^^^^^^^^^^^^^
      double precision A,delta,DeuxBE
c......................................................................
c
      A = A_inf -FracMasseInt*Masse_SE
      DeuxBE = 2.0d0*(B_inf-Ener_SE)
      delta = A*A-2.0d0*J22_SE*DeuxBE
c      if (iPas_Evol.GT.39268710) then
c         write(0,*)'!!!  iPas_Evol,delta,A,B,J2,E = ',
c     $        iPas_Evol,delta,A,B_inf,J22_SE,Ener_SE,' !!!' !##########
c      end if

      if (delta.LE.0.0d0) 
     $     then                 ! peut-etre une orbite circulaire...
         if (delta.GE.Tol_circul*A*A) then
           Apo = -2.0d0*J22_SE/A
         else                   ! Il y a une bulle ...
            write(0,*)'!!!  delta,A,B,J2,E = ',
     $           delta,A,B_inf,J22_SE,Ener_SE,' !!!'
            if (lStop_err_PeriApo) then
               call Meurt('Resoud_Rapo')
            else
               iErr_PeriApo = iNegDelta_err_PeriApo
               Ese_err_PeriApo = Ener_SE
               Jse_err_PeriApo = sqrt(2.0d0*J22_SE)
               A_err_PeriApo = A
               B_err_PeriApo = B_inf
            end if
         end if
      
      else
        Apo = (-A+sqrt(delta))/DeuxBE
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine write_err_PeriApo(iFile)
c----------------------------------------------------------------------
      include 'PeriApo_SE2_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFile
c......................................................................
c
      if (iErr_PeriApo.NE.iOK_err_PeriApo) then
         write(iFile,100) iErr_PeriApo
 100     format(
     $        ' Error in peri/apo-centre determination',/,
     $        '     iErr_PeriApo = ',I2,' : ',$)
         if (iErr_PeriApo.EQ.iRsupGTPeri_err_PeriApo) then
            write(iFile,*) 'iRsupGTPeri_err_PeriApo'
         else if (iErr_PeriApo.EQ.iRinfLTPeri_err_PeriApo) then
            write(iFile,*) 'iRinfLTPeri_err_PeriApo'
         else if (iErr_PeriApo.EQ.iRinfGTApo_err_PeriApo) then
            write(iFile,*) 'iRinfGTApo_err_PeriApo'
         else if (iErr_PeriApo.EQ.iRsupLTApo_err_PeriApo) then
            write(iFile,*) 'iRsupLTApo_err_PeriApo'
         else if (iErr_PeriApo.EQ.iRinfNotFound_err_PeriApo) then
            write(iFile,*) 'iRinfNotFound_err_PeriApo'
         else if (iErr_PeriApo.EQ.iRsupNotFound_err_PeriApo) then
            write(iFile,*) 'iRsupNotFound_err_PeriApo'
         else if (iErr_PeriApo.EQ.iRsupMaxBissec_err_PeriApo) then
            write(iFile,*) 'iRsupMaxBissec_err_PeriApo'
         else if (iErr_PeriApo.EQ.iNegDelta_err_PeriApo) then
            write(iFile,*) 'iNegDelta_err_PeriApo'
         end if
         write(iFile,200) 
     $        Mse_err_PeriApo, Ese_err_PeriApo, Jse_err_PeriApo,
     $        Rloc_err_PeriApo, Qloc_err_PeriApo,
     $        dQdRloc_err_PeriApo,
     $        A_err_PeriApo, B_err_PeriApo
 200     format(
     $        'Mse_err_PeriApo : ',E12.5,/,
     $        'Ese_err_PeriApo : ',E12.5,/,
     $        'Jse_err_PeriApo : ',E12.5,/,
     $        'Rloc_err_PeriApo : ',E12.5,/,
     $        'Qloc_err_PeriApo : ',E12.5,/,
     $        'dQdRloc_err_PeriApo : ',E12.5,/,
     $        'A_err_PeriApo : ',E12.5,/,
     $        'B_err_PeriApo : ',E12.5
     $        )
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine write_QdQ_err_PeriApo
c----------------------------------------------------------------------
      include 'PeriApo_SE2_inc.f'
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer i,iSE
      double precision Q_loc,dQdR_loc,
     $     A_parcours, B_parcours
c......................................................................

      open(11,file='_erreur_CalcPeri_SE_ext2.asc')
      write(11,900)
 900  format('# 1:R 2:Q 3:dQdR')
      do i=0,1000
         R_loc=R_inf+0.001d0*(R_sup-R_inf)*i
         call Calc_Q_dQ(Q_loc,dQdR_loc)
         write(11,910) R_loc,Q_loc,dQdR_loc
 910     format(3(E20.12,' '))
      end do
      close(11)
      open(11,file='_erreur_CalcPeri_SE_ext2_AB.asc')
      write(11,920)
 920  format('# 1:R 2:A 3:B')
      B_parcours=-B0
      A_parcours=0.0d0
      do i=1,NbSE-1
         iSE=iListeRay(i)
         A_parcours=A_parcours-M_SE(iSE)
         B_parcours=B_parcours+M_SE(iSE)/R_SE(iSE)
         write(11,'(3(E20.12," "))') R_SE(iSE),A_parcours,B_parcours
      end do
      close(11)
c
         end      
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      block data init_PeriApo_SE2
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      include 'PeriApo_SE2_inc.f'

      data iErr_PeriApo /iOK_err_PeriApo/ ! No error
      data lStop_err_PeriApo /.TRUE./ ! By default, one stops on errors; this is safer
c
      end
