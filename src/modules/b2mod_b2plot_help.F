      module b2mod_b2plot_help
      use b2mod_types
      use b2mod_constants

      integer nlabels,nfns,nsparameters,nconstants
      integer nwlabels,nwfns,ntrgfns
      parameter(nlabels=247,nfns=192,nsparameters=27,
     *     nconstants=10,nwfns=61,nwlabels=85,ntrgfns=82)

      integer, save ::
     * dimensions(nlabels), wdimensions(nwlabels) 
      real (kind=R8), save ::
     * constant_value(nconstants)

      character*8, save ::
     *  labels(nlabels), fns(nfns),
     *  wlabels(nwlabels), wfns(nwfns), trgfns(ntrgfns),
     *  sparameters(nsparameters), constant_name(nconstants)
      character*64, save ::
     *  descriptions(nlabels), wdescriptions(nwlabels)

      data labels(1:50) /
     1  'pit' , 'bx'  , 'by'  , 'bz'  , 'bb'  ,
     2  'b'   , 'na'  , 'ua'  , 'vv'  , 'ww'  ,
     3  'te'  , 'ti'  , 'pr'  , 'ne'  , 'po'  ,
     4  'dfhe', 'dfhi', 'dhet', 'dhit', 'dfhj',
     5  'dfch', 'dfna', 'dfmo', 'dfhm', 'dfhp',
     6  'zeff', 'up'  , 'fnax', 'fnay', 'fmox',
     7  'fmoy', 'fhex', 'fhey', 'fhetx','fhety',
     8  'fhix', 'fhiy', 'fhitx','fhity','fhmx', 
     9  'fhmy', 'fhpx', 'fhpy', 'fhjx', 'fhjy', 
     *  'fhtx', 'fhty', 'cvsx', 'cvsy', 'ebld'/
      data labels(51:100) /
     1  'ebls', 'ebal', 'pbal', 'hx'  , 'hy'  , 
     2  'hy1' , 'sx'  , 'sy'  , 'vol' , 'mach', 
     3  'ehx' , 'ehy' , 'dab2', 'dmb2', 'ahal', 
     4  'mhal', 'rqae', 'rsai', 'rsan', 'rsam',
     5  'rran', 'rrai', 'rram', 'rcxn', 'rcxi',
     6  'rcxm', 'b2br', 'fetx', 'fety', 'tab2',
     7  'tmb2', 'rfla', 'rflm', 'refa', 'refm',
     8  'pfla', 'pflm', 'pefa', 'pefm', 'mch2',
     9  'mi'  , 'zi'  , 'zn'  , 'fmx' , 'fmy' ,
     *  'p*'  , 'conn', 'cr'  , 'cz'  , 'na0' /
      data labels(101:150) /
     1  'te0' , 'ti0' , 'ua0' , 'po0' , 'smpr',
     2  'smpt', 'smfr', 'pot' , 'poti', 'pot*',
     3  'smp' , 'nreg', 'smo' , 'smq' , 'smaf',
     4  'smag', 'smav', 'rco' , 'rmo' , 'rmt' ,
     5  'ree' , 'rei' , 'rpo' , 'cpa' , 'cua' ,
     6  'cut' , 'cte' , 'cti' , 'cpo' , 'ixlt',
     7  'ixeq', 'ixgt', 'iylt', 'iyeq', 'iygt',
     8  'islt', 'iseq', 'isgt', 'q'   , 's'   ,
     9  'alph', 'd'   , 'dp'  , 'kye' , 'kyi' , 
     *  'kyi0', 'vlax', 'vlay', 'visc', 'bvis'/
      data labels(151:200) /
     1  'sig' , 'alf' , 'fchx', 'fchy', 'ud'  ,
     2  'vd'  , 'vadia.x','vadia.y','wadia.x','wadia.y',
     3  'vexb.x','vexb.y','cd'  , 'dd'  , 'niso',
     4  'ix'  , 'iy'  , 'iix' , 'iiy' , 'lix' ,
     5  'liy' , 'rix' , 'riy' , 'bix' , 'biy' ,
     6  'tix' , 'tiy' , 'rlcl', 'gdcl', 'dib2',
     7  'tib2', 'brna', 'brmo', 'brhe', 'brhi',
     8  'brch', 'bcna', 'bcmo', 'bche', 'bchi', 
     9  'bcch', 'dvue', 'dvua', 'exbe', 'exba',
     *  'visa', 'joul', 'fraa', 'str' , 'ke'  / 
      data labels(201:247) /
     1  'facdrift', 'fac_ExB' , 'fac_vis', 'xymap'  , 'fnax32' ,
     2  'fnay32'  , 'fnax52'  , 'fnay52' , 'fnix32' , 'fniy32' ,
     3  'fnix52'  , 'fniy52'  , 'fnex32' , 'fney32' , 'fnex52' ,
     4  'fney52'  , 'fchx32'  , 'fchy32' , 'fchx52' , 'fchy52' ,
     5  'exna'    , 'exmo'    , 'exhe'   , 'exhi'   , 'exch'   ,
     6  'b2ra'    , 'mach4'   , 'fle'    , 'fli'    , 'nflfnax',
     7  'nflfnay' , 'nflfhix' , 'nflfhiy', 'flvisc' , 'sna'    ,
     8  'she'     , 'shi'     , 'sch'    , 'crx'    , 'cry'    ,
     9  'nh'      , 'qz0'     , 'qz1'    , 'sxperp' , 'ue'     ,
     *  'fchp'    , 'za2'     /
      data descriptions(1:10) /
     1    'Pitch',
     2    'bx (T) (poloidal)',
     3    'by (T)',
     4    'bz (T) (toroidal)',
     5    'bb (T)',
     6    'B (T) (components)',
     7    'Ion density (m^-3)',
     8    'Parallel velocity (ua) (m.s^-1)',
     9    'Radial velocity (vv) (m.s^-1)',
     *    'Toroidal velocity (ww) (m.s^-1)'/
      data descriptions(11:20) /
     1    'Electron temperature (eV)',
     2    'Ion temperature (eV)',
     3    'Pressure (Pa)',
     4    'Electron density (m^-3)',
     5    'phi (V)',
     6    'diff.fhe (W)',
     7    'diff.fhi (W)',
     8    'diff.fhet (W)',
     9    'diff.fhit (W)',
     *    'diff.fhj (W)'/
      data descriptions(21:30) /
     1    'diff.fch (A)',
     2    'div.fna (particles/s)',
     3    'div.fmo (N.s)',
     4    'div.fhm (W)',
     5    'div.fhp (W)',
     6    'Zeff',
     7    'Poloidal velocity (up) (m.s^-1)',
     8    'Poloidal particle flux (s^-1)',
     9    'Radial particle flux (s^-1)',
     *    'Poloidal momentum flux (N)'/
      data descriptions(31:40) /
     1    'Radial momentum flux (N)',
     2    'Poloidal electron energy flux (W)',
     3    'Radial electron energy flux (W)',
     4    'Poloidal electron energy flux [fhetx] (W)',
     5    'Radial electron energy flux [fhety] (W)',
     6    'Poloidal ion energy flux (W)',
     7    'Radial ion energy flux (W)',
     8    'Poloidal ion energy flux [fhitx] (W)',
     9    'Radial ion energy flux [fhity] (W)',
     *    'Poloidal kinetic energy flux (W)'/
      data descriptions(41:50) /
     1    'Radial kinetic energy flux (W)',
     2    'Poloidal potential energy flux (W)',
     3    'Radial potential energy flux (W)',
     4    'Poloidal electrostatic energy flux (W)',
     5    'Radial electrostatic energy flux (W)',
     6    'Poloidal total energy flux (W)',
     7    'Radial total energy flux (W)',
     8    'cvsa_cl_x',
     9    'cvsa_cl_y',
     *    'Energy balance (-divergence of fluxes) (W)'/
      data descriptions(51:60) /
     1    'Energy balance (sources) (W)',
     2    'Energy balance (W)',
     3    'Particle balance (particles)',
     4    'hx (m)',
     5    'hy (m)',
     6    'hy1 (m) (hy corrected)',
     7    'sx (m^2)',
     8    'sy (m^2)',
     9    'vol (m^3)',
     *    'Mach number'/
      data descriptions(61:70) /
     1    'ehx',
     2    'ehy',
     3    'Atomic density (m^-3)',
     4    'Molecular density (m^-3)',
     5    'Atomic H-alpha emissivity (photons.m^-3.s^-1)',
     6    'Molecular H-alpha emissivity (photons.m^-3.s^-1)',
     7    'Electron cooling rate (W) [rqahe]',
     8    'Ionisation ion energy losses (W) [rsahi]',
     9    'Ionisation particle losses (s^-1) [rsana]',
     *    'Ionisation momentum losses (m.s^-2) [rsamo]'/
      data descriptions(71:80) /
     1    'Recombination particle losses (s^-1) [rrana]',
     2    'Recombination ion energy losses (W) [rrahi]',
     3    'Recombination momentum losses (m.s^-2) [rramo]',
     4    'Particle neutral losses (s^-1) [rcxna]',
     5    'Ion energy neutral losses (W) [rcxhi]',
     6    'Ion momentum neutral losses (m.s^-2) [rcxmo]',
     7    'Bremsstrahlung losses (W)',
     8    'Total poloidal internal energy flux (W)',
     9    'Total radial internal energy flux (W)',
     *    'Neutral atom temperature (eV)'/
      data descriptions(81:90) /
     1    'Neutral molecule temperature (eV)',
     2    'Radial atomic flux (s^-1.m^-2)',
     3    'Radial molecular flux (s^-1.m^-2)',
     4    'Radial atomic energy flux (W.m^-2)',
     5    'Radial molecular energy flux (W.m^-2)',
     6    'Poloidal atomic flux (s^-1.m^-2)',
     7    'Poloidal molecular flux (s^-1.m^-2)',
     8    'Poloidal atomic energy flux (W.m^-2)',
     9    'Poloidal molecular energy flux (W.m^-2)',
     *    'Normalised parallel velocity'/
      data descriptions(91:100) /
     1    'Ion mass (kg)',
     2    'Ion charge',
     3    'Nuclear charge',
     4    'Poloidal mass flux (kg/s)',
     5    'Radial mass flux (kg/s)',
     6    'Static pressure ()',
     7    'Connection length (m)',
     8    'CR (m)',
     9    'CZ (m)',
     *    'na at previous time step (m^-3)'/
      data descriptions(101:110) /
     1    'te at previous time step (eV)',
     2    'ti at previous time step (eV)',
     3    'ua at previous time step (m.s^-1)',
     4    'po at previous time step (V)',
     5    'Electrical field force',
     6    'Thermal force',
     7    'Friction force',
     8    'Potential energy (eV) (cum. ionisation energy)',
     9    'Ionisation energy (eV)',
     *    'Ionisation energy (eV) [shifted down one]'/
      data descriptions(111:120) /
     1    'Self pressure force',
     2    'Region',
     3    'smo',
     4    'smq',
     5    'smaf',
     6    'smag',
     7    'smav',
     8    'RESCO',
     9    'RESMO',
     *    'RESMT'/
      data descriptions(121:130) /
     1    'RESEE',
     2    'RESEI',
     3    'RESPO',
     4    'CORPA',
     5    'CORUA',
     6    'CORUT',
     7    'CORTE',
     8    'CORTI',
     9    'CORPO',
     *    'ixlt'/
      data descriptions(131:140) /
     1    'ixeq',
     2    'ixgt',
     3    'iylt',
     4    'iyeq',
     5    'iygt',
     6    'islt',
     7    'iseq',
     8    'isgt',
     9    'q',
     *    'Shear'/
      data descriptions(141:150) /
     1'alpha (ratio of pressure gradient to critical ballooning value)',
     2    'D (m^2.s^-1)',
     3    'Dpa (m^2.s^-1)',
     4    'kye (m^2.s^-1)',
     5    'kyi (m^2.s^-1)',
     6    'kyi0 (m^2.s^-1)',
     7    'Poloidal pinch velocity (m.s^-1)',
     8    'Radial pinch velocity (m.s^-1)',
     9    'Viscosity (m.kg^-1.s^-1)',
     *    'Bas'' viscosity piece ()'/
      data descriptions(151:160) /
     1    'Anomalous conductivity (mho.m^-1)',
     2    'Thermo-electric coefficient (mho.m^-1)',
     3    'Poloidal current (A)',
     4    'Radial current (A)',
     5    'pol comp of the drift velocity [uadia.x] (m.s^-1)',
     6    'rad comp of the drift velocity [uadia.y] (m.s^-1)',
     7    'pol comp of div part of diamag velocity [vadia.x] (m.s^-1)',
     8    'rad comp of div part of diamag velocity [vadia.y] (m.s^-1)',
     9    'pol comp of the total diamag velocity [wadia.x] (m.s^-1)',
     *    'rad comp of the total diamag velocity [wadia.y] (m.s^-1)'/
      data descriptions(161:170) /
     1    'pol comp of the ExB velocity [vaecrb.x] (m.s^-1)',
     2    'rad comp of the ExB velocity [vaecrb.y] (m.s^-1)',
     3    'fchdia.x',
     4    'fchdia.y',
     5    'Isolating region',
     6    'ix (normalised)',
     7    'iy (normalised)',
     8    'ix',
     9    'iy',
     *    'leftix'/
      data descriptions(171:180) /
     1    'leftiy',
     2    'rightix',
     3    'rightiy',
     4    'bottomix',
     5    'bottomiy',
     6    'topix',
     7    'topiy',
     8    'real cells',
     9    'guard cells',
     *    'Test ion density (m^-3)'/
      data descriptions(181:190) /
     1    'Test ion temperature (eV)',
     2    'b2stbr_sna',
     3    'b2stbr_smo',
     4    'b2stbr_she',
     5    'b2stbr_shi',
     6    'b2stbr_sch',
     7    'b2stbc_sna',
     8    'b2stbc_smo',
     9    'b2stbc_she',
     *    'b2stbc_shi'/
      data descriptions(191:200) /
     1    'b2stbc_sch',
     2    'b2sihs_divue (phm0)',
     3    'b2sihs_divua (phm1)',
     4    'b2sihs_exbe (phm5)',
     5    'b2sihs_exba (phm6)',
     6    'b2sihs_visa (phm2)',
     7    'b2sihs_joule (phm3)',
     8    'b2sihs_fraa (phm4)',
     9    'b2sihs_str (style)',
     *    'Kinetic energy (J)'/
      data descriptions(201:210) /
     1    'facdrift',
     2    'fac_ExB',
     3    'fac_vis',
     4    'xymap',
     5    'fnax_32 (s^-1)',
     6    'fnay_32 (s^-1)',
     7    'fnax_52 (s^-1)',
     8    'fnay_52 (s^-1)',
     9    'fnix_32 (s^-1)',
     *    'fniy_32 (s^-1)'/
      data descriptions(211:220) /
     1    'fnix_52 (s^-1)',
     2    'fniy_52 (s^-1)',
     3    'fnex_32 (s^-1)',
     4    'fney_32 (s^-1)',
     5    'fnex_52 (s^-1)',
     6    'fney_52 (s^-1)',
     7    'fchx_32 (A)',
     8    'fchy_32 (A)',
     9    'fchx_52 (A)',
     *    'fchy_52 (A)'/
      data descriptions(221:230) /
     1    'ext_sna',
     2    'ext_smo',
     3    'ext_she',
     4    'ext_shi',
     5    'ext_sch',
     6    'Radiation rate (W)',
     7    'Mach number (SOLPS4.0)',
     8    'Electron heat flux limit factor',
     9    'Ion heat flux limit factor',
     *    'fllim0fna_x'/
      data descriptions(231:247) /
     1    'fllim0fna_y',
     2    'fllim0fhi_x',
     3    'fllim0fhi_y',
     4    'fllimvisc',
     5    'sna',
     6    'she',
     7    'shi',
     8    'sch',
     9    'Corner radial coordinate (m)',
     *    'Corner vertical coordinate (m)',
     1    'Hydrogenic atom density (m^-3)',
     2	  'qz0',
     3	  'qz1',
     4	  'sxperp (m^2)',
     5    'Parallel electron velocity (ue) (m.s^-1)',
     6    'Parallel current (A)',
     7    'Squared ion charge' /
      data dimensions(1:50) /
     1    2,      2,      2,      2,      2,
     2    3,      3,      3,      3,      3,
     3    2,      2,      2,      2,      2,
     4    2,      2,      2,      2,      2,
     5    2,      3,      3,      3,      3,
     6    2,      3,      3,      3,      3,
     7    3,      2,      2,      2,      2,
     8    2,      2,      2,      2,      3,
     9    3,      3,      3,      2,      2,      
     *    2,      2,      3,      3,      2/      
      data dimensions(51:100) /
     1    2,      2,      3,      2,      2,
     2    2,      2,      2,      2,      2,
     3    2,      2,      3,      3,      2,
     4    2,      3,      3,      3,      3,
     5    3,      3,      3,      3,      3,
     6    3,      2,      2,      2,      3,
     7    3,      3,      3,      3,      3,
     8    3,      3,      3,      3,      3,
     9    3,      3,      3,      3,      3,
     *    2,      2,      2,      2,      3/
      data dimensions(101:150) /
     1    2,      2,      3,      2,      3,
     2    3,      3,      3,      3,      3,
     3    3,      2,      3,      3,      3,
     4    3,      3,      3,      3,      2,
     5    2,      2,      2,      3,      3,
     6    2,      2,      2,      2,      2,
     7    2,      2,      2,      2,      2,
     8    3,      3,      3,      2,      2,
     9    2,      3,      3,      2,      3,
     *    2,      3,      3,      3,      3/
      data dimensions(151:200) /
     1    2,      2,      2,      2,      3,
     2    3,      3,      3,      3,      3,
     3    3,      3,      2,      2,      2,
     4    2,      2,      2,      2,      2,
     5    2,      2,      2,      2,      2,
     6    2,      2,      2,      2,      3,
     7    3,      3,      3,      2,      2,      
     8    2,      3,      3,      2,      2,
     9    2,      2,      2,      2,      2,
     *    2,      2,      2,      2,      3/
      data dimensions(201:247) /
     1    2,      2,      2,      2,      3,
     2    3,      3,      3,      2,      2,
     3    2,      2,      2,      2,      2,
     4    2,      2,      2,      2,      2,
     5    3,      3,      2,      2,      2,
     6    3,      2,      3,      2,      2, 
     7    3,      3,      3,      3,      3,
     8    2,      2,      2,      2,      2,
     9    2,      2,      2,      2,      2,
     *    2,      3/    
      data fns(1:50) /
     1 'sumx', 'sumy', 'f.x' , 'f.y' , 'intx',
     2 'inty', 'f.lx', 'f.ly', 'cont', 'surf',
     3 'line', 'both', 'phys', 'comp', 'mesh',
     4 'page', 'logf', 'linf', 'logx', 'linx', 
     5 'psx' , 'psy' , 'abs' , 'clvs', 'clvd',
     6 'clvn', 'eclv', 'fclv', 'font', 'lacl', 
     7 'lhor', 'lbon', 'lbof', 'lwid', 'grid',
     8 'clvl', 'sumz', 'spln', 'arow', 'arww', 
     9 'pvol', 'sep' , 'vesl', 'nocp', 'sepw', 
     * 'pne' , 'sq'  , 'a4p' , 'a4l' , 'same'/
      data fns(51:100) /
     1 'filo', 'spec', 'writ', 'm+'  , 'm-'  ,
     2 'm*'  , 'm/'  , 'rm+' , 'rm*' , 'rm**',
     3 'drop', 'dup' , 'swap', 'over', 'r+'  ,
     4 'r-'  , 'r*'  , 'r/'  , 'r**' , 'xsel',
     5 'ysel', 'zsel', 'mmin', 'mmax', 'rmmi',
     6 'rmma', 'new' , 'user', 'log' , 'ln'  ,
     7 'lbsz', 'ndec', 'zrng', 'aspc', 'chs' , 
     8 'sqrt', 'exp' , 'col' , 'bw'  , 'zero',
     9 'bala', 'outl', 'stdl', 'mvsk', 'ave' ,
     * 'traj', 't_#p', 'reso', 'coro', 'cora'/
      data fns(101:150) /
     1 'coni', 'shco', 'swid', 'xyrz', 't_r+',
     2 't_r-', 'mcol', 'scol', 'vcol', 'ocol',
     3 'rigz', 'noax', 'nspe', 'stow', 'xtow',
     4 'ytow', 'part', 'lbsa', 'nese', 'tese',
     5 'dase', 'ytra', 'ecol', 'uuse', 'vvse',
     6 'arox', 'aroy', 'pos' , 'ddx' , 'ddy' ,
     7 'f30' , 'f31' , 'eirc', 'debu', 'scat',
     8 'sdrop','sdup','idrop', 'idup','rdrop',
     9 'rdup', 'logt', 'lint','movie', 'omp' ,
     * 'imp' , 'yref', 'mmid','target','trgoff'/
      data fns(151:192) /
     1 'mlt'  , 'mle' , 'meq' , 'mne' , 'mgt' ,
     2 'mge'  , 'clvp', 'xpol', 'xpar', 'i+'  ,
     3 'i-'   , 'i*'  , 'i/'  , 'i**' ,'rswap',
     4 'rover','iswap','iover', 'zcol','sswap',
     5 'sover', 'itor', 'itom', 'itow', 'itot',
     6 'rtom' , 'rtow', 'rtot', 'itostr', 'rtostr',
     7 'owid' , 'vwid', 'ewid', 'iprt', 'rprt',
     8 'sprt' , 'list', 'acol', 'dcol', 'wcol',
     9 'tcol' , 'psxperp'/ 
      data wlabels /
     1 'wallarea', 'ei'      , 'wallflux', 'wallhflx', 'r1'      ,
     2 'r2'      , 'z1'      , 'z2'      , 'walltemp', 'plathick',
     3 'cooltemp', 'vaprtemp', 'north'   , 'south'   , 'west'    ,
     4 'east'    , 'nseg'    , 'iwlt'    , 'iweq'    , 'iwgt'    ,
     5 'iswlt'   , 'isweq'   , 'iswgt'   , 'ibnd'    , 'walldist',
     6 'wallngth', 'sput'    , 'sput2'   , 'sputchem', 'thermev' ,
     7 'sputres' , 'backscat', 'backscae', 'erosion' , 'deposit' ,
     8 'ero_ml'  , 'depo_ml' , 'time'    , 'imapx'   , 'imapy'   ,
     9 'kappa_a' , 'kappa_b' , 'kappa_a1', 'kappa_a2', 'kappa_a3',
     * 'kappa_a4', 'kappapa1', 'kappapa2', 'kappapa3', 'kappapa4',
     1 'tgdns'   , 'cp_b1'   , 'cp_b2'   , 'cp_b3'   , 'erel'    , 
     2 'bbdyemis', 'vapprs_a', 'vapprs_b', 'vapprs'  , 'qevap'   ,
     3 'qrad'    , 'stckcoef', 'sublheat', 'ckappa_a', 'ckappa_b', 
     4 'ckappaa1', 'ckappaa2', 'ckappaa3', 'ckappaa4', 'ckapppa1', 
     5 'ckapppa2', 'ckapppa3', 'ckapppa4', 'ctgdns'  , 'c_cp_b1' , 
     6 'c_cp_b2' , 'c_cp_b3' , 'coathick', 'lkappa'  , 'lrho'    ,
     7 'lcp'     , 'lemis'   , 'sputfrac', 'lthick'  , 'ldepth'  / 
      data wdimensions /
     1 1,       2,      2,      1,      1,
     2 1,       1,      1,      1,      1,
     3 1,       1,      1,      1,      1,
     4 1,       1,      1,      1,      1,
     5 2,       2,      2,      1,      1,
     6 1,       2,      2,      2,      1,
     7 2,       2,      2,      2,      2,
     8 2,       2,      1,      1,      1,
     9 1,       1,      1,      1,      1,
     * 1,       1,      1,      1,      1,
     1 1,       1,      1,      1,      1,
     2 1,       1,      1,      1,      1,
     3 1,       1,      1,      1,      1,
     4 1,       1,      1,      1,      1,
     5 1,       1,      1,      1,      1,
     6 1,       1,      1,      1,      1,
     7 1,       1,      2,      1,      1/
      data wdescriptions(1:10) /
     1 'Wall area (m^2)',
     2 'Incident energy (eV)',
     3 'Wall particle flux (s-1)',
     4 'Wall heat flux (W/m^2)',
     5 'Wall start R position (m)',
     6 'Wall end R position (m)',
     7 'Wall start Z position (m)',
     8 'Wall end Z position (m)',
     9 'Wall surface temperature (K)',
     * 'Plate thickness (m)'/
      data wdescriptions(11:20) /
     1 'Backplate temperature (K)',
     2 'Vaporisation temperature (K)',
     3 'North boundary',
     4 'South boundary',
     5 'West boundary',
     6 'East boundary',
     7 'Wall segment',
     8 'iwlt',
     9 'iweq',
     * 'iwgt'/
      data wdescriptions(21:30) /
     1 'iswlt',
     2 'isweq',
     3 'iswgt',
     4 'ibnd',
     5 'Distance along wall (m)',
     6 'Wall length (m)',
     7 'Physical sputtering yield',
     8 'Physical sputtering energy yield',
     9 'Chemical sputtering yield',
     * 'Thermal Evaporation Fluxes (/m2/s)'/
      data wdescriptions(31:40) /
     1 'RES sputtering yield',
     2 'Backscattering fraction',
     3 'Backscattering energy fraction',
     4 'Gross erosion',
     5 'Gross deposition',
     6 'Monolayer erosion',
     7 'Monolayer deposition',
     8 'Plate time (s)',
     9 'imapx',
     * 'imapy'/
      data wdescriptions(41:50) /
     1 'Kappa A factor (1/AT+B) (m/W)',
     2 'Kappa B factor (1/AT+B) (m.K/W)',
     3 'Kappa A1 polynomial factor (W/m/K)',
     4 'Kappa A2 polynomial factor (W/m/K^2)',
     5 'Kappa A3 polynomial factor (W/m/K^3)',
     6 'Kappa A4 polynomial factor (W/m/K^4)',
     7 'Kappa_perp A1 polynomial factor (W/m/K)',
     8 'Kappa_perp A2 polynomial factor (W/m/K^2)',
     9 'Kappa_perp A3 polynomial factor (W/m/K^3)',
     * 'Kappa_perp A4 polynomial factor (W/m/K^4)'/
      data wdescriptions(51:60) /
     1 'Bulk mass density (kg/m^3)',
     2 'Heat capacity B1 factor (J/g/K)',
     3 'Heat capacity B2 factor (J/g/K^2)',
     4 'Heat capacity B3 factor (J.K/g)',
     5 'Chemical sputtering factor erel',
     6 'Blackbody emissivity',
     7 'Vapor pressure A factor (K)',
     8 'Vapor pressure B factor',
     9 'Vapor pressure (Pa)',
     * 'Energy lost by vaporisation (J/m^2/s)'/
      data wdescriptions(61:80) /
     1 'Energy lost by blackbody radiation (J/m^2/s)',
     2 'Sticking coefficient',
     3 'Heat of sublimation (J/kg)',
     4 'Coating Kappa A factor (1/AT+B) (m/W)',
     5 'Coating Kappa B factor (1/AT+B) (m.K/W)',
     6 'Coating Kappa A1 polynomial factor (W/m/K)',
     7 'Coating Kappa A2 polynomial factor (W/m/K^2)',
     8 'Coating Kappa A3 polynomial factor (W/m/K^3)',
     9 'Coating Kappa A4 polynomial factor (W/m/K^4)',
     * 'Coating Kappa_perp A1 polynomial factor (W/m/K)',
     1 'Coating Kappa_perp A2 polynomial factor (W/m/K^2)',
     2 'Coating Kappa_perp A3 polynomial factor (W/m/K^3)',
     3 'Coating Kappa_perp A4 polynomial factor (W/m/K^4)',
     4 'Coating mass density (kg/m^3)',
     5 'Coating heat capacity B1 factor (J/g/K)',
     6 'Coating heat capacity B2 factor (J/g/K^2)',
     7 'Coating heat capacity B3 factor (J.K/g)',
     8 'Coating thickness (m)',
     9 'Layer heat conductivity (W/m/K)',
     * 'Layer mass density (kg/m^3)'/
      data wdescriptions(81:85) /
     1 'Layer heat capacity (J/kg/K)',
     2 'Layer blackbody emissivity',
     3 'Layer fraction',
     4 'Layer thickness (m)',
     5 'Layer depth (m)'/
      data wfns /
     1 'wall' , 'wall_ank', 'logw', 'linw', 'logxw', 'linxw', 'plwl', 
     2 'sumw' , 'intw','intw*', 'pswa','wclvs','wclvd','wclvn', 'wclvp',
     3 'wwid' , 'wabs', 'ddw' ,'wsumz', 'pwne', 'wm+' , 'wm-' , 'wm*',
     4 'wm/'  , 'rwm+', 'rwm-', 'rwm*', 'rwm/',
     5 'rwm**','wdrop','wdup','wswap','wover',
     6 'wsel' ,'wzsel','wmmin','wmmax','rwmmi',
     7 'rwmma','wmlt', 'wmle', 'wmeq', 'wmne',
     8 'wmgt' , 'wmge', 'wlog', 'wln' ,'wsqrt', 
     9 'wexp' ,'wzrng', 'wchs', 'wtos', 'wtot', 'wpos', 
     * 'wsaw' ,'wallf','outpw', 'wprt', 'wlld', 
     1 'woff' ,'wlist'/
      data trgfns /
     1 'logtrg'  , 'lintrg'  , 'tclvs'   , 'tclvd'   , 'tclvn'   ,
     2 'trg.w'   , 'trg.z'   , 'sumtrgw' , 'sumtrgz' , 'inttrgw' , 
     3 'inttrgz' , 'trgddw'  , 'trgddp'  , 'trgddz'  , 'rtrgm+'  , 
     4 'rtrgm-'  , 'rtrgm*'  , 'rtrgm/'  , 'rtrgm**' , 'trgdrop' , 
     5 'trgtemp' , 'trgwsel' , 'rtrgmmi' , 'rtrgmma' , 'trglt'   , 
     6 'trgle'   , 'trgeq'   , 'trgne'   , 'trggt'   , 'trgge'   , 
     7 'trglog'  , 'trgln'   , 'trgsqrt' , 'trgexp'  , 'trgchs'  , 
     8 'trgtow'  , 'trgpos'  , 'outptrg' , 'trgprt'  , 'trgl'    , 
     9 'trgsurf' , 'kappa_ab', 'kappa'   , 'dkappadt', 'kappap'  , 
     * 'kappapdt', 'cp'      , 'trgwidth', 'trgdepth', 'trgdup'  , 
     1 'trgswap' , 'trgover' , 'trgmmin' , 'trgmmax' , 'trgr1'   , 
     2 'trgr2'   , 'trgz1'   , 'trgz2'   , 'trgnorth', 'trgsouth', 
     3 'trgwest' , 'trgeast' , 'trgseg'  , 'itrglt'  , 'itrgeq'  , 
     4 'itrggt'  , 'trgm+'   , 'trgm-'   , 'trgm*'   , 'trgm/'   ,
     5 'trgvol'  , 'tclvp'   , 'inttrgw*', 'ztrglt'  , 'ztrgeq'  ,
     6 'ztrggt'  , 'trgabs'  , 'trgiz'   , 'logxtrg' , 'linxtrg' ,
     7 'trglist' , 'rho'/
      data sparameters /
     1 'extl', 'cltb', 'chor', 'chpl', 'cmd' ,
     2 'cmdcdf','chbx','chvl', 'flux', 'chmp',
     3 'cam' , 'movi', 'sonn', 'saw' , 'oson', 'uedgegrd',
     4 'expf', 'outp', 'pdev', 'mprt', 'amju',
     5 'wallblbl', 'trgblbl', 'glbl' , 'wlbl', 'tlbl',
     6 'tranprof'/
      data constant_name /
     1 'me'    ,'mp' , 'ev'  , 'qe'   , 'mu0',
     2 'eps0'  ,'pi' , 'c'   , 'avogr', 'kbolt'/
      data constant_value /
     1  me,      mp,    ev,     qe,      mu0,  
     2  eps0,    pi,    c,      avogr,   kbolt/

      end module b2mod_b2plot_help
