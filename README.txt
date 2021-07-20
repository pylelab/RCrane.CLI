RCrane.CLI, version 1.0.2

OVERVIEW

RCrane.CLI takes input of a PDB file containing phosphates, bases, and C1'
atoms, and uses these coordinates to predict and build backbone conformers.
This is done using the methodology described in:
    Keating KS and Pyle AM.  Semiautomated model building for RNA
      crystallography using a directed rotameric approach.  Proc Natl Acad
      Sci USA, 107: 8177-8182 (2010).
Note that all publications resulting from use of this program must
acknowledge this manuscript.

Please also note that RCrane.CLI is a command line program and is not
intended for interactive use.  Crystallographers wishing to use this
methodology for structure building should use RCrane, which is a Coot plugin
that allows for faster and easier building of RNA into electron density.
RCrane is available from http://pylelab.org/software/.

RCrane.CLI is copyright 2010, Kevin Keating, and is licensed under the
Educational Community License, Version 2.0.


INSTALLATION INSTRUCTIONS

RCrane.CLI requires Perl 5, version 5.8 or newer.  If you are using Linux or
Macintosh OS X, Perl should already be installed on your computer.  If you are
using Windows, you may install ActivePerl
(http://www.activestate.com/activeperl/) or Strawberry Perl
(http://strawberryperl.com/).

Once Perl and the necessary modules are installed, RCrane.CLI may be
installed by simply unzipping the program into any desired directory
(ex. /usr/local/xtal/rcrane.cli/ or C:\Program Files\RCrane.CLI\).  Users
may wish to place this directory in their path (or create a symbolic link in
the appropriate directory).


RUNNING RCRANE.CLI

RCrane.CLI must be run from the command line.  Input for RCrane.CLI is a PDB
file containing phosphates, bases, and C1' atoms.  (Note that RCrane.CLI does
not help with the placement of these atoms.  The Coot plugin will help a
crystallographer in placing these atoms and will be released in summer 2010.)
If RCrane.CLI is installed in one of the locations described above and the
necessary atoms are placed in struc.pdb, the program may be invoked as follows:
    /usr/local/xtal/rcrane.cli/rcrane.cli.pl struc.pdb
        or
    "C:\Program Files\RCrane.CLI\rcrane.cli.pl" struc.pdb

After the program has finished running, four output files will be generated:
struc.built.pdb
    A PDB file containing a molecule constructed using the most likely
    conformers
struc.data.csv
    A table containing all data used for the conformers predictions:
    theta', eta', base-phosphate perpendicular distances, C1'-C1'
    distances, and P-P distances
struc.probs.csv
    A table containing the probabilities for each conformer at each suite
    of the structure
struc.rots.csv
    A table containing the three most likely conformers for each suite,
    along with their sugar puckers. The bottom of this table also contains
    a conformer string for the most likely conformers.

Additional options are described in MANUAL.txt or by running rcrane.cli.pl -m.


CHANGELOG

version 1.0.2
Remove dependency on perl Math modules.

version 1.0.1
Changed name from CONDOR.CLI to RCrane.CLI
Updated PNAS citation

version 1.0
Initial release
