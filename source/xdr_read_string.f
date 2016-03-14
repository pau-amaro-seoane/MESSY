c----------------------------------------------------------------------
      integer function ixdr_read_string(ixdr,str)
c----------------------------------------------------------------------
      include 'fxdr.inc'
      integer ixdr,l, length_str
      character*(*) str
c......................................................................
c
      str=''

      ixdr_read_string = ixdrstring(ixdr,str)
      if (str.NE.'') str=str(1:length_str(str))
      !write(0,'("|",A,"|")') str(1:length_str(str))
c
      end

c----------------------------------------------------------------------
      integer function length_str(str)
c----------------------------------------------------------------------
      character*(*) str
      character*1 nul_
c......................................................................
      nul_=CHAR(0)
      length_str=len(str)
      if (length_str.GT.0) then
         do while (str(length_str:length_str).EQ.' ' .OR.
     $        str(length_str:length_str).EQ.nul_ )
            length_str = length_str-1
         end do
      end if
c
      end
