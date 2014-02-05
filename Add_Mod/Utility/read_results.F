      subroutine read_results
      use b2mod_optim
      use b2mod_exp_file_names
      implicit none
	integer iy

*   ..for inner lower
      fnm = 'mod_'//trim(fnm_ne_inner_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_ne_inner_lower(iy), ne_inner_lower_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_te_inner_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_te_inner_lower(iy), te_inner_lower_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_ti_inner_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_ti_inner_lower(iy), ti_inner_lower_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_js_inner_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_js_inner_lower(iy), js_inner_lower_mod(iy)
	 enddo
	close(70)

*   ..for inner upper
      fnm = 'mod_'//trim(fnm_ne_inner_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_ne_inner_upper(iy), ne_inner_upper_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_te_inner_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_te_inner_upper(iy), te_inner_upper_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_ti_inner_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_ti_inner_upper(iy), ti_inner_upper_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_js_inner_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_js_inner_upper(iy), js_inner_upper_mod(iy)
	 enddo
	close(70)

*   ..for outer midpane
      fnm = 'mod_'//trim(fnm_ne_outer_mdpl)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_ne_outer_mdpl(iy), ne_outer_mdpl_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_te_outer_mdpl)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_te_outer_mdpl(iy), te_outer_mdpl_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_ti_outer_mdpl)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_ti_outer_mdpl(iy), ti_outer_mdpl_mod(iy)
	 enddo
	close(70)

*   ..for outer lower
      fnm = 'mod_'//trim(fnm_ne_outer_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_ne_outer_lower(iy), ne_outer_lower_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_te_outer_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_te_outer_lower(iy), te_outer_lower_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_ti_outer_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_ti_outer_lower(iy), ti_outer_lower_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_js_outer_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_js_outer_lower(iy), js_outer_lower_mod(iy)
	 enddo
	close(70)

*   ..for outer upper
      fnm = 'mod_'//trim(fnm_ne_outer_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_ne_outer_upper(iy), ne_outer_upper_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_te_outer_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_te_outer_upper(iy), te_outer_upper_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_ti_outer_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_ti_outer_upper(iy), ti_outer_upper_mod(iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_js_outer_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  read(70) cy_js_outer_upper(iy), js_outer_upper_mod(iy)
	 enddo
	close(70)

      fnm = 'mod_'//trim(fnm_power)//'.dat'
      open(70,file=fnm,form='unformatted')
	  read(70) power_mod
	close(70)

	end subroutine read_results
