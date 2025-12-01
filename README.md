# B2.5

The **B2** code solves a system of fluid equations that models
a multi-species plasma in a two-dimensional geometry. Its main
application is simulation of the edge plasma in tokamaks.

**B2** is composed of several independent programs that communicate
through external files. The main sequence involves the programs
`b2ag` (prepares a geometry and magnetic field), `b2ah` (prepares
a table of default physics parameters), `b2ai` (prepares an
initial plasma state), `b2mn` (principal calculation), and `b2yp`
(postprocessing).

Additional programs may prepare tables of rate coefficients for
atomic processes, tables of plasma-wall interaction processes,
external source terms, and other physics, or may perform further
postprocessing.
