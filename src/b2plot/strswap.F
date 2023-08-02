
      subroutine strswap(s1,s2)
      use b2mod_types
      character*(*) s1,s2
      character*256 s3
      s3=s1
      s1=s2
      s2=s3
      return
      end subroutine strswap
