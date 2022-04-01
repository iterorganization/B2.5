      function eval (x)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   eval,x
*     ------------------------------------------------------------------
*     EVAL returns its real*8 argument through the function name.  It may
*     be used to force a compiler to store an intermediate operand.
*     ------------------------------------------------------------------
      eval = x
      return
*     ------------------------------------------------------------------
      end function eval
