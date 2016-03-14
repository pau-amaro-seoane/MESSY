c======================================================================
c     routines d'ecriture et lecture employant le format XDR
c======================================================================

c----------------------------------------------------------------------
      subroutine XDR_Ecrit_i(ixdr,NomVar,iVar) ! ecrit un entier
c----------------------------------------------------------------------
      include 'XDR_Ecrit_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer ixdr              ! Numero du fichier XDR
      character*(*) NomVar      ! Nom de la variable a ecrire
      integer iVar              ! Variable a ecrire
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(ixdr,NomVar)
      ixdr_err = ixdrstring(ixdr,NomType_i)
      ixdr_err = ixdrint(ixdr,1)
      ixdr_err = ixdrint(ixdr,iVar)
c
      end
c
c----------------------------------------------------------------------
      subroutine XDR_Ecrit_ti(ixdr,NomVar,iDimVar,itVar) ! ecrit un tableua d'entiers
c----------------------------------------------------------------------
      include 'XDR_Ecrit_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer ixdr              ! Numero du fichier XDR
      character*(*) NomVar      ! Nom de la variable a ecrire
      integer iDimVar           ! la dimension du tableau
      integer itVar(iDimVar)    ! Variable a ecrire
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(ixdr,NomVar)
      ixdr_err = ixdrstring(ixdr,NomType_ti)
      ixdr_err = ixdrint(ixdr,iDimVar)
      ixdr_err = ixdrimat(ixdr,iDimVar,itVar)
c
      end
c
c----------------------------------------------------------------------
      subroutine XDR_Ecrit_r(ixdr,NomVar,rVar) ! ecrit un reel simple precision
c----------------------------------------------------------------------
      include 'XDR_Ecrit_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer ixdr              ! Numero du fichier XDR
      character*(*) NomVar      ! Nom de la variable a ecrire
      real rVar                 ! Variable a ecrire
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(ixdr,NomVar)
      ixdr_err = ixdrstring(ixdr,NomType_r)
      ixdr_err = ixdrint(ixdr,1)
      ixdr_err = ixdrreal(ixdr,rVar)
c
      end
c
c----------------------------------------------------------------------
      subroutine XDR_Ecrit_tr(ixdr,NomVar,iDimVar,rtVar) ! ecrit un tableau de reels simple precision
c----------------------------------------------------------------------
      include 'XDR_Ecrit_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer ixdr              ! Numero du fichier XDR
      character*(*) NomVar      ! Nom de la variable a ecrire
      integer iDimVar           ! la dimension du tableau
      real rtVar(iDimVar)       ! Variable a ecrire
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(ixdr,NomVar)
      ixdr_err = ixdrstring(ixdr,NomType_tr)
      ixdr_err = ixdrint(ixdr,iDimVar)
      ixdr_err = ixdrrmat(ixdr,iDimVar,rtVar)
c
      end
c
c----------------------------------------------------------------------
      subroutine XDR_Ecrit_d(ixdr,NomVar,dVar) ! ecrit un reel double precision
c----------------------------------------------------------------------
      include 'XDR_Ecrit_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer ixdr              ! Numero du fichier XDR
      character*(*) NomVar      ! Nom de la variable a ecrire
      double precision dVar     ! Variable a ecrire
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(ixdr,NomVar)
      ixdr_err = ixdrstring(ixdr,NomType_d)
      ixdr_err = ixdrint(ixdr,1)
      ixdr_err = ixdrdouble(ixdr,dVar)
c
      end
c
c----------------------------------------------------------------------
      subroutine XDR_Ecrit_td(ixdr,NomVar,iDimVar,dtVar) ! ecrit un tableau de reels double precision
c----------------------------------------------------------------------
      include 'XDR_Ecrit_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer ixdr              ! Numero du fichier XDR
      character*(*) NomVar      ! Nom de la variable a ecrire
      integer iDimVar           ! la dimension du tableau
      double precision dtVar(iDimVar) ! Variable a ecrire
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(ixdr,NomVar)
      ixdr_err = ixdrstring(ixdr,NomType_td)
      ixdr_err = ixdrint(ixdr,iDimVar)
      ixdr_err = ixdrdmat(ixdr,iDimVar,dtVar)
c
      end
c
c----------------------------------------------------------------------
      subroutine XDR_Ecrit_s(ixdr,NomVar,sVar) ! ecrit une chaine de caracteres
c----------------------------------------------------------------------
      include 'XDR_Ecrit_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer ixdr              ! Numero du fichier XDR
      character*(*) NomVar      ! Nom de la variable a ecrire
      character*(*) sVar        ! Variable a ecrire
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(ixdr,NomVar)
      ixdr_err = ixdrstring(ixdr,NomType_s)
      ixdr_err = ixdrstring(ixdr,sVar)
c
      end
c
