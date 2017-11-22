B2.5 IDS documentation
============================================

Required python2.7 packages
----------------------------

  - Sphinx-Fortran:

    $ pip install --user sphinx-fortran

  - VACUMM:
  Download .tar file *from http://www.ifremer.fr/vacumm/user.install.installs.html* and

    $ tar xzf vacumm-3.4.0.tar.gz
    $ cd vacumm-3.4.0
    $ python setup.py build
    $ python setup.py install --user

Generating HTML Fortran Documentation (using Sphinx 1.5.6)
----------------------------------------------------------

    $ module load imas/3.13.0/ual/3.6.3
    $ make clean
    $ make html
    $ firefox $SOLPSTOP/modules/B2.5/src/documentation/build/index.html




