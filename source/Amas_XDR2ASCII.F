c======================================================================
c     convert an XDR file for cluster structure (MESSY's *AMAS.xdr format)
c     to /rdb format. Get the name of xdr file as parameter on comand
c     line. Writes to standard output
c======================================================================
c
      implicit none
      include 'fxdr.inc'
c
c var locales :
c ^^^^^^^^^^^^^
      character*256 NomFich
      integer iXDR, ixdr_err
c
c fonctions :
c ^^^^^^^^^^^
      integer iargc
c......................................................................
c
      if (iargc().LT.1) then
         write(0,*) '> Name of input file (MESSY''s AMAS.xdr format) :'
         read(6,'(A)') NomFich
      else 
         call getarg(1,NomFich)
      end if

      iXDR = initxdr(NomFich,'r',.TRUE.)
      call LitAmas_XDR(iXDR)
      ixdr_err = ixdrclose(iXDR)

      call EcritAmas_ASC(6,.FALSE.)
c
      end
c
      
