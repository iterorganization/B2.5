      subroutine b2xxid (prname, lbl0)
      use b2mod_types
      implicit none
      character prname*(*), lbl0*60
*     ------------------------------------------------------------------
*     B2XXID returns an identification of the current version of b2.
*     The argument prname should be the name of the main program; e.g.,
*     prname = 'b2mn'.
*     ------------------------------------------------------------------
      character id0*12, id1*12, prnam0*8
      parameter (id0='B2.5 bjb99/7', id1=' SOLPS-ITER ')
      character jobnam*26
      external jobnam
*     id0 identifies the core version of b2; its value is set by bjb.
*     id1 identifies the local version of b2; its character*12
*     value may be set by the user.
*     The function jobnam returns some character*26 identification of
*     the present job. The form of this identification may vary between
*     computer systems; it will typically be constructed as some
*     combination of date, time, jobnumber, and machine identification.
*     ------------------------------------------------------------------
      prnam0 = prname
      lbl0 = prnam0//' '//id0//id1//' '//jobnam()
      return
*     ------------------------------------------------------------------
      end
