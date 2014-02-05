      subroutine gnopen (ioun, chname)
      use b2mod_types
      implicit none
      integer ioun
      character chname*(*)
*     ------------------------------------------------------------------
*     GNOPEN opens a graphics file for use with GKS.
*     This routine is system dependent and not completely general.
*     A copy of the argument chname is edited to conform to system
*     requirements on file names. See the code for details.
*     ------------------------------------------------------------------
      character fname*16
      intrinsic len
      external xertst, chcase, gopwk, gwsrep
*     ------------------------------------------------------------------
*   ..construct fname
*     (max length is len(fname) characters; make lowercase)
      fname = chname
      call chcase (0, fname)
      call xertst (len(chname).le.len(fname),
     &  'gnopen--faulty argument chname, too long: '//fname)
*   ..specify the filename to GKS
*     (inquire GKS implementation level)
*$$$  call gqlvks (gkserr, gkslev)
*     (gkslev: 0=L0a, 1=L0b, 2=L0c, 3=L1a, ..., 8=L2c)
*     (specify file name)
*$$$  if (1.le.gkslev) then
*$$$   call gscnid (ioun, fname)
*$$$  endif
*!!   Note: The routine gscnid (set connection identifier) is perhaps
*     not an official GKS routine. It occurs in GRAFPAK_GKS but not in
*     NCAR_GKS. NCAR_GKS allows just a single graphics file, which by
*     default receives the name "gmeta"; the environment variable
*     NCARG_GKS_OUTPUT may be set to specify a different file name.
*     Alternatively, in NCARG_GKS-3.2 the filename may be specified by
*      call gesc (-1391, 1, fname, 1, 1, chdum)
*     where chdum is a dummy character*80 variable.
*   ..open the GKS workstation
      call gopwk (ioun, ioun, 1)
*   ..set workstation properties
      call gwsrep (ioun, 0)
      return
*     ------------------------------------------------------------------
      end subroutine gnopen
