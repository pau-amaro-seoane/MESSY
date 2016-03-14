c======================================================================
c     Routines de traitement de chaines
c======================================================================
c
c$$$ccccc
c$$$c ESSAIS
c$$$      character*10 c1,c2,c3
c$$$      read(*,*) c1
c$$$c     read(*,*) c2
c$$$      c2='bbb ccc'
c$$$      write (*,*) c1 // c2
c$$$      call ConcatChaines2(c1,c2)
c$$$      write (*,*) c1
c$$$c
c$$$      end
c$$$cccccc
c
c----------------------------------------------------------------------
      subroutine SauteBlancs(Chaine,i)
c----------------------------------------------------------------------
c
      implicit none
      character*(*) Chaine
      integer i,l
c
      l=LEN(Chaine)
      do while(i.LE.l .AND. Chaine(i:i).EQ.' ')
         i=i+1
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine SauteNonBlancs(Chaine,i)
c----------------------------------------------------------------------
c
      implicit none
      character*(*) Chaine
      integer i,l
c
      l=LEN(Chaine)
      do while(i.LE.l .AND. Chaine(i:i).NE.' ')
         i=i+1
      end do
c
      end
c
c------------------------------------------------------------
      subroutine ConcatChaines(c1,c2,c3)
c------------------------------------------------------------
c     renvoie dans c3 c1 et c2 bout a bout
c
      character*(*) c1,c2,c3
c
      integer i,max,l1
c
      c3 = ' '
      l1 = LEN(c1)
      i = l1+1
 10   i = i-1
      if (c1(i:i).EQ.' ') then
         goto 10
      endif
      c3(1:i) = c1(1:i)
      max = MIN(LEN(c2),LEN(c3)-i)
      c3(i+1:i+max) = c2(1:max)
c
      end
c
c------------------------------------------------------------
      subroutine ConcatChaines2(c1,c2)
c------------------------------------------------------------
c     renvoie dans c1 c1 et c2 bout a bout
c 
      implicit none
      character*(*) c1,c2
c
      integer i,j,l1
c
      l1 = LEN(c1)
      do i=l1,1,-1
         if (c1(i:i) .NE. ' ') goto 1
      end do
 1    continue
      do j=i+1,MIN(l1,i+LEN(c2))
         c1(j:j)=c2(j-i:j-i)
      end do
c
      end
c
c--------------------------------------------------
      subroutine EnleverBlancs(Ligne)
c--------------------------------------------------
      implicit none
c
      character*(*) Ligne
c
      integer i,j
c..................................................
c
      i=0
 1    i=i+1
      if (Ligne(i:i).EQ.' ') goto 1
c      
      do j=i,LEN(Ligne)
         Ligne(j-i+1:j-i+1)=Ligne(j:j)
      end do
c
      end
c
c--------------------------------------------------
      subroutine MettresBlancs(Ligne,Long_totale)
c--------------------------------------------------
      implicit none
c
      character*(*) Ligne
      integer Long_totale
c
      integer i
c..................................................
c
      do i=LEN(Ligne),Long_totale
         Ligne(i:i)=' '
      end do
c
      end
c
c----------------------------------------------------------------------
      integer function LongChaine(Chaine)
c----------------------------------------------------------------------
      implicit none
      character*(*) Chaine
      character*1 nul_
c......................................................................
      nul_=CHAR(0)
      LongChaine=len(Chaine)
      if (LongChaine.EQ.0) return
      do while (Chaine(LongChaine:LongChaine).EQ.' ' .OR.
     $     Chaine(LongChaine:LongChaine).EQ.nul_ )
         LongChaine = LongChaine-1
         if (LongChaine.EQ.0) return
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine OutputFileHeader(iFile,Header)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^
      integer iFile
      character*(*) Header
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer i, j, l, npos, LongChaine
      character*16 form_str
c......................................................................
c
      write(iFile,100)
 100  format('#',$)

      i=1
      j=0
      l=LongChaine(Header)

      do while(i.LT.l)
         do while(i.LE.l .AND. Header(i:i) .EQ. ' ')
            i = i + 1
         end do

         if (i.LE.l) then
            j = j+1
            npos = int(log10(dble(j)))+1
            write(form_str,90) npos
 90         format('(" ",I',I1,',": ",$)')
            write(iFile,form_str) j
            do while(i.LE.l .AND. Header(i:i) .NE. ' ')
               write(iFile,110) Header(i:i)
 110           format(A1,$)
               i = i + 1
            end do
         end if
      end do
      write(iFile,*) ''
c
      end
c
c----------------------------------------------------------------------
      subroutine GetVariableValue(iFile,CommentChars,EqualChars,
     $     VarName,VarValue,iStatus)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
                                ! INPUTS:
      integer iFile
      character*(*)
     $     CommentChars,        ! the character(s) at the beginning of a line (if any)
     $     EqualChars,          ! the character(s) to separate variable name from value (if any)
     $                          ! OUTPUTS:
     $     VarName,VarValue     ! Name and value of variable
      integer iStatus           ! 0 if OK
                                ! 1 if CommentChars not found
                                ! 2 if VarName not found
                                ! 3 if EqualChars not found
                                ! 4 if VarValue not found
                                ! 5 if VarName chain too short
                                ! 6 if VarValue chain too short
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      character*256 OneLine
c
c functions :
c ^^^^^^^^^^^
      integer iCharInString
c......................................................................
c
      end
c
c----------------------------------------------------------------------
      integer function iCharInString(Char,String)
c----------------------------------------------------------------------
      integer none
c
c arguments :
c ^^^^^^^^^^^
      character Char
      character*(*) String
c
c local varaibles :
c ^^^^^^^^^^^^^^^^^
      integer l
c......................................................................
c
      l = len(String)
      iCharInString = 0
      do i=1,l
         if (Char.EQ.String(i:i)) then
            iCharInString = i
            return
         end if
      end do
c
      end
c
