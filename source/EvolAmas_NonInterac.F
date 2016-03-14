#define _VERB_NONINTERAC_ 1
c----------------------------------------------------------------------
      subroutine Init_NonInterac
c----------------------------------------------------------------------
      include 'EvolAmas_NonInterac_inc.f'
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer*1 Type1, Type2
c.......................................................................
c
      do Type1=1,Nb_TypesEtoile
         do Type2=1,Nb_TypesEtoile
            lRelaxesWith(Type1,Type2) = .TRUE.
         end do
      end do
      lNonInterac = .FALSE.
      lNonInteracAsymmetrical = .FALSE.
c
      end
c
c----------------------------------------------------------------------
      subroutine Read_NonInterac_Prescriptions(iFile)
c----------------------------------------------------------------------
      include 'EvolAmas_NonInterac_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFile
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      character*32 string
      integer*1 Type1, Type2
c.......................................................................
c
      call Init_NonInterac
c
c---- read prescriptions for non-interacting combinations
c
      do while(.TRUE.)
         read(iFile,'(A)',END=9) string
         if (string(1:1).NE.'#' .AND. string.NE.'') then
            read(string,*) Type1, Type2
            if (Type1.LT.1 .OR. Type1.GT.Nb_TypesEtoile) then
               write(0,*) '!!! Type ',Type1,
     $              ' out of range of stellar types !!!'
               call Meurt('Read_NonInterac_Prescriptions')
            end if
            if (Type2.LT.1 .OR. Type2.GT.Nb_TypesEtoile) then
               write(0,*) '!!! Type ',Type2,
     $              ' out of range of stellar types !!!'
               call Meurt('Read_NonInterac_Prescriptions')
            end if
#if (_VERB_NONINTERAC_>0)
            write(0,*) '@@ NonInterac: Type ',Type1,
     $           ' does not relax with Type ',Type2
#endif
            lRelaxesWith(Type1,Type2) = .FALSE.
            lNonInterac = .TRUE.
         end if
      end do

 9    continue
c
c---- look for asymmetrical prescription (energy non-conserving)
c
      lNonInteracAsymmetrical = .FALSE.
      do Type1=1,Nb_TypesEtoile-1
         do Type2=Type1+1,Nb_TypesEtoile
            lNonInteracAsymmetrical = lNonInteracAsymmetrical .OR.
     $           (lRelaxesWith(Type1,Type2) .NEQV.
     $            lRelaxesWith(Type2,Type1))
         end do
      end do
#if (_VERB_NONINTERAC_>0)
      if (lNonInteracAsymmetrical) then
         write(0,*) '@@ NonInterac: Non symmetrical interactions, ',
     $        'energy won''t be conserved !'
      end if
#endif
c
      end
c
c----------------------------------------------------------------------
      block data init_NonInterac_common
c----------------------------------------------------------------------
      include 'EvolAmas_NonInterac_inc.f'

      integer Ntot
      parameter (Ntot=Nb_TypesEtoile*Nb_TypesEtoile)
      data lRelaxesWith /Ntot*.TRUE./
      data lNonInterac /.FALSE./
      data lNonInteracAsymmetrical /.FALSE./

      end
