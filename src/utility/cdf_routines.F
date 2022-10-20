      subroutine check_cdf_status(status)
      implicit none
      integer status
#ifndef NO_CDF
#     include <netcdf.inc>
      if(status.ne.NF_NOERR) write(*,*) 'CDFERR: ',
     1 trim(nf_strerror(status))
#else
      status=0
#endif
      return
      end subroutine check_cdf_status
