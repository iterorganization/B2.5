B2.5 IDS documentation
======================

Generating and viewing documentation (HTML)
-------------------------------------------

    $ module load doxygen/1.8.8
    $ make doc
    $ firefox source/Doxygen/html/index.html

Generating and viewing documentation (PDF)
-------------------------------------------

    $ module load doxygen/1.8.8 texlive/2015
    $ make doc
    $ cd source/Doxygen/latex
    $ make
    $ firefox refman.pdf

or

    $ evince refman.pdf
