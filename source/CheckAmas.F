#if (_CHECK_AMAS_STAND_ALONE_==1)
      program CheckAmas_SA
c======================================================================
c     programme calculant, par quadrature, les periodes orbitales
c     de toutes les SE d'un amas dont le nom de fichier est donne en 
c     parametre sur la ligne de commande
c     fournit egalement les distance peri- et apocentriques
c======================================================================
      implicit none
c
c var locales :
c ^^^^^^^^^^^^^
      character*256 arg
      integer iarg,narg
      integer iNbErrors, iNbErrors_serious
c
c fonction :
c ^^^^^^^^^^
      integer iargc, LongChaine
c......................................................................
c
c---- Lecture des options eventuelles
c
      iarg = 0
      narg = iargc()
      call getarg(iarg+1,arg)
      do while(arg(1:1).EQ.'-')
         iarg = iarg+1
         if (arg(1:LongChaine(arg)).EQ.'-xxx') then
            continue
         else
            write(0,*) '!!! option ',arg(1:LongChaine(arg)),
     $           ' non comprise !!!'
            call exit(1)
         end if
         if (iarg.LT.narg) call getarg(iarg+1,arg)
      end do
c
c---- charger l'amas en memoire
c
      call Charge_Amas(arg)
c
c---- carry out consistency checks
c
      iNbErrors = 0
      call CheckAmas(6,iNbErrors,iNbErrors_serious)
      if (iNbErrors.GT.0) then
         write (0,210) iNbErrors,iNbErrors_Serious
 210     format('!!! ',I6,' Error(s) found in cluster data',
     $        ' including ',I6,' serious one(s) !!!')
         if (iNbErrors_Serious.GT.0) then
            call exit(2)
         else
            call exit(1)
         end if
      else
         write (0,*) '> No error found in cluster data'
         call exit(0)
      end if
c
      end

#endif /* (_CHECK_AMAS_STAND_ALONE_==1) */

c......................................................................
      subroutine CheckAmas(iFile,iNbErrors,iNbErrors_serious)
c......................................................................
c     Sort cluster before using!!

      implicit none

      include 'CstesMath.f'
      include 'Param_Amas.f'
      include 'VarSE.f'
      include 'EvolAmas_Pertes_inc.f'
      include 'EvolAmas_TN_common.f'
      include 'PeriApo_SE2_diags_param.f'
      include 'PeriApo_SE2_diags_common.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFile, iNbErrors, iNbErrors_serious
c
c var locales :
c ^^^^^^^^^^^^^
      integer iRang,iRang_dum,iSE, itick,inext_tick,ierr
      double precision A,B,Peri,Apo,E,R_tid
c
c function :
c ^^^^^^^^^^
      double precision R_tronc_mar
c......................................................................
c
      write(0,90)
 90   format(/,
     $     '##############################',/,
     $     '## Checking cluster data... ##',/,
     $     '##############################',//,
     $     '0    10   20   30   40   50   60   70   80   90  100%',/,
     $     '|----|----|----|----|----|----|----|----|----|----|',/,
     $     $)
      inext_tick=0
      itick=0
      ierr=0
c
c---- boucle de verification
c
      iNbErrors = 0
      iNbErrors_serious = 0
      R_tid = R_tronc_mar()

      do iRang = 1,NbSE

         if (iRang.GE.inext_tick) then
            if (ierr.EQ.0) then
               write(0,95)
 95            format('^',$)
            else if (ierr.LT.10) then
               write(0,96) ierr
 96            format((I1),$)
            else
               write(0,97)
 97            format('!',$)
            end if
            itick=itick+1
            inext_tick = 0.02d0*NbSE*itick
            ierr=0
         end if

         iSE = iListeRay(iRang)
            
         if (R_SE(iSE).LT.R_perte) then
c
c---- Check that quantities are positive
c
            if (
     $           M_SE(iSE).LE.0.0d0 .OR.
     $           R_SE(iSE).LE.0.0d0 .OR.
     $           T_SE(iSE).LE.0.0d0 .OR.
     $           J_SE(iSE).LE.0.0d0
     $           ) then
               write(iFile,100) iSE,
     $              M_SE(iSE),R_SE(iSE),T_SE(iSE),J_SE(iSE)
 100           format('!!! particule iSE = ',I7,
     $              ' : M_SE,R_SE,T_SE,J_SE = ',4(E12.5,' '),'!!!')
               iNbErrors = iNbErrors+1
               iNbErrors_serious = iNbErrors_serious+1
               goto 999
            end if
                                ! extraction de la SE
            call Extraire_SE_Arbre(iSE,iRang_dum,A,B)
            A = A + M_SE(iSE)
            B0 = B0 - M_SE(iSE)/R_SE(iSE)
            M_amas = M_amas - M_SE(iSE)
                                ! Calcul du Peri et de l'apocentre
            E = T_SE(iSE)+(A-FracMasseInt*M_SE(iSE))/R_SE(iSE)+B
c
c---- Check if particle is bound to cluster
c
            if (E.GT.0.0d0) then
               write(iFile,300) iSE, E
 300           format('!!! particule iSE = ',I7,' non liee avec E = ',
     $              E12.5,' !!!')
               write(iFile,400)
     $              M_SE(iSE),R_SE(iSE),T_SE(iSE),J_SE(iSE),
     $              sqrt(2*T_SE(iSE)),J_SE(iSE)/R_SE(iSE)
 400           format('  M_SE,R_SE,T_SE,J_SE = ',4(E23.16,' '),/,
     $              '   V,Vtg = ',2(E23.16,' '))
               iNbErrors = iNbErrors+1
               ierr = ierr+1
               goto 99
            end if
c
c---- Check if one can compute peri and apocentre distances
c
            iErr_PeriApo = iOK_err_PeriApo

            lStop_err_PeriApo = .FALSE.
            call CalcPeri_SE_ext2(M_SE(iSE),E,J_SE(iSE),Peri)
            if (iErr_PeriApo.NE.iOK_err_PeriApo) then
               write(iFile,305) iSE
 305           format('!!! particule iSE = ',I7,
     $              '; erreur lors du calcul de R_peri !!!')
               call write_err_PeriApo(iFile)
               write(iFile,400)
     $              M_SE(iSE),R_SE(iSE),T_SE(iSE),J_SE(iSE),
     $              sqrt(2*T_SE(iSE)),J_SE(iSE)/R_SE(iSE)
               iNbErrors = iNbErrors+1
               iNbErrors_serious = iNbErrors_serious+1
               ierr = ierr+1
               goto 99
            endif
               
            lStop_err_PeriApo = .FALSE.
            call CalcApo_SE_ext2 (M_SE(iSE),E,J_SE(iSE),Apo )
            if (iErr_PeriApo.NE.iOK_err_PeriApo) then
               write(iFile,307) iSE
 307           format('!!! particule iSE = ',I7,
     $              '; erreur lors du calcul de R_apo !!!')
               call write_err_PeriApo(iFile)
               write(iFile,400)
     $              M_SE(iSE),R_SE(iSE),T_SE(iSE),J_SE(iSE),
     $              sqrt(2*T_SE(iSE)),J_SE(iSE)/R_SE(iSE)
               iNbErrors = iNbErrors+1
               iNbErrors_serious = iNbErrors_serious+1
               ierr = ierr+1
               goto 99
            endif

            if (R_SE(iSE).GT.Apo .OR.R_SE(iSE).LT.Peri) then
               write(iFile,310) iSE, R_SE(iSE), Peri, Apo
 310           format('!!! particule iSE = ',I7,
     $              ' en position interdite avec R,Peri,Apo = ',
     $              3(E12.5,' '),'!!!')
               write(iFile,400)
     $              M_SE(iSE),R_SE(iSE),T_SE(iSE),J_SE(iSE),
     $              sqrt(2*T_SE(iSE)),J_SE(iSE)/R_SE(iSE)
               iNbErrors = iNbErrors+1
               iNbErrors_serious = iNbErrors_serious+1
               ierr = ierr+1
               goto 99
            end if
c
c---- Check that particle doesn't go beyond tidal radius
c
            if (Apo.GT.R_tid) then 
               write(iFile,320) iSE, Apo, R_tid
 320           format('!!! particule iSE = ',I7,
     $              ' franchit le rayon de maree; Apo, R_tid = ',
     $              2(E12.5,' '),'!!!')
               write(iFile,400)
     $              M_SE(iSE),R_SE(iSE),T_SE(iSE),J_SE(iSE),
     $              sqrt(2*T_SE(iSE)),J_SE(iSE)/R_SE(iSE)
               iNbErrors = iNbErrors+1
               ierr = ierr+1
               goto 99
            end if

 99         continue
                                ! Re-introduction de la SE
            call Inserer_SE_Arbre(iSE,iRang_dum,A,B)
            B0 = B0 + M_SE(iSE)/R_SE(iSE)
            M_amas = M_amas + M_SE(iSE)
         else                   ! Cas des SE "perdues"
            Peri=1.0d30
            Apo=1.0d30
         end if
 999     continue
      end do
      write(0,*) ''

c
c---- Try to find each particle from its rank
c
      write(0,92)
 92   format(/,
     $     '####################################',/,
     $     '## Accessing particles by rank... ##',/,
     $     '####################################',//,
     $     '0    10   20   30   40   50   60   70   80   90  100%',/,
     $     '|----|----|----|----|----|----|----|----|----|----|',/,
     $     $)
      inext_tick=0
      itick=0
      ierr=0

      do iRang = 1,NbSE

         if (iRang.GE.inext_tick) then
            if (ierr.EQ.0) then
               write(0,95)
            else if (ierr.LT.10) then
               write(0,96) ierr
            else
               write(0,97)
            end if
            itick=itick+1
            inext_tick = 0.02d0*NbSE*itick
            ierr=0
         end if
         
         call Trouver_SE_Rang(iRang,iSE,A,B)
         ! Should change error treatment in Trouver_SE_Rang
         ! the current version will just ptint a message and stop in case of error
      end do
      write(0,*) ''
c
      end
c
c--------------------------------------------------------------------
      subroutine Check_Amas_And_Report
c----------------------------------------------------------------------
      implicit none
c
c var locales :
c ^^^^^^^^^^^^^
      integer iNbErrors,iNbErrors_Serious
      integer iFich_Prov
      parameter (iFich_Prov=11)
c......................................................................
c
      iNbErrors = 0
      open(iFich_Prov,file='ErrorsCheckAmas.asc')
      call CheckAmas(iFich_Prov,iNbErrors,iNbErrors_Serious)
      if (iNbErrors.GT.0) then
         close(iFich_Prov)
         write (0,210) iNbErrors,iNbErrors_Serious
 210     format('!!! ',I6,' Error(s) found in cluster data',
     $        ' including ',I6,' serious one(s)',
     $        ' (see ErrorsCheckAmas.asc) !!!')
         if (iNbErrors_Serious.GT.0) call Meurt('Check_Amas_And_Report')
      else
         close(iFich_Prov,STATUS='DELETE')
      end if
c
      end
c      
