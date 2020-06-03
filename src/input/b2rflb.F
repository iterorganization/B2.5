*-----------------------------------------------------------------------
*.specification

      subroutine b2rflb (nrid, lbl1)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nrid
*   ..output arguments (unspecified on entry)
      character lbl1*60

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2RFLB reads from a formatted file a label for one run of one
*     of the programs that make up B2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine reads from unit (nrid) a label to identify the
*     present run of one of the programs that make up B2. Unit nrid
*     must be connected to a formatted file.
*     The input file must be positioned properly: the first record must
*     begin with '*label  ' or '*LABEL  ' (without the quotes), and the
*     second record must contain the character*60 label in free format
*     (enclosed in quotes). Only those two records are read.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nrid - integer, input.
*     nrid specifies a unit number for formatted input.
*
*     lbl1 - character*60, output.
*     lbl1 specifies the label read from input.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      character id*8
*   ..procedures
      logical streql
      external subini, subend, xertst, xerrab, streql

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2rflb')
*   ..test nrid
      call xertst (1.le.nrid, 'faulty argument nrid')
*   ..read header line
      read (nrid,'(a8)',err=91) id
      call xertst (streql(id,'*label'),
     &  'faulty input id at *label, id='//id)
*   ..read label
      read (nrid,*,err=92) lbl1

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read header *label')
   92 call xerrab ('error trying to read after '//id)

*-----------------------------------------------------------------------
*.end b2rflb

      end
