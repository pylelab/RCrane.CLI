#!/usr/bin/perl

#module to calculate torsions given a set of atomic coordinates

# Copyright 2010 Kevin Keating Licensed under the
# Educational Community License, Version 2.0 (the "License"); you may
# not use this file except in compliance with the License. You may
# obtain a copy of the License at
#
# http://www.osedu.org/licenses/ECL-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an "AS IS"
# BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
# or implied. See the License for the specific language governing
# permissions and limitations under the License.

package MeasureTorsions;

use strict;
use English '-no_match_vars';
use Carp;

use StrucCalc qw(calcTorsion);

use Exporter;
our @ISA       = ('Exporter');
our @EXPORT    = qw(calcAlpha calcBeta calcGamma calcDelta calcEpsilon calcZeta calcChi calcXi);
our @EXPORT_OK = qw(calcAlpha calcBeta calcGamma calcDelta calcEpsilon calcZeta calcChi calcXi);
our %EXPORT_TAGS = (backbone => [qw(calcAlpha calcBeta calcGamma calcDelta calcEpsilon calcZeta)], sugar => [qw(calcChi calcXi)]);


sub calcAlpha {
    #calculate alpha
    #ARGUMENTS:
    #   $atoms - a hash of atomic coordinates
    #            the O3' of the previous nucleotide should be named O3'-1
    #RETURNS:
    #   the alpha torsion in degrees
    #ERRORS RAISED:
    #   will croak if the required atoms are not present
    
    my $atoms = shift;
    
    unless (defined $atoms->{"O3'-1"} and defined $atoms->{"P"} and defined $atoms->{"O5'"} and defined $atoms->{"C5'"}) {
        croak "Missing atom for calculating alpha";
    }
    
    return calcTorsion($atoms->{"C5'"}, $atoms->{"O5'"}, $atoms->{"P"}, $atoms->{"O3'-1"});
}

sub calcBeta {
    #calculate beta
    #ARGUMENTS:
    #   $atoms - a hash of atomic coordinates
    #RETURNS:
    #   the beta torsion in degrees
    #ERRORS RAISED:
    #   will croak if the required atoms are not present
    
    my $atoms = shift;
    
    unless (defined $atoms->{"P"} and defined $atoms->{"O5'"} and defined $atoms->{"C5'"} and defined $atoms->{"C4'"}) {
        croak "Missing atom for calculating beta";
    }
    
    return calcTorsion($atoms->{"P"}, $atoms->{"O5'"}, $atoms->{"C5'"}, $atoms->{"C4'"});
}

sub calcGamma {
    #calculate beta
    #ARGUMENTS:
    #   $atoms - a hash of atomic coordinates
    #RETURNS:
    #   the gamma torsion in degrees
    #ERRORS RAISED:
    #   will croak if the required atoms are not present
    
    my $atoms = shift;
    
    unless (defined $atoms->{"O5'"} and defined $atoms->{"C5'"} and defined $atoms->{"C4'"} and defined $atoms->{"C3'"}) {
        croak "Missing atom for calculating gamma";
    }
    
    return calcTorsion($atoms->{"O5'"}, $atoms->{"C5'"}, $atoms->{"C4'"}, $atoms->{"C3'"});
}

sub calcDelta {
    #calculate beta
    #ARGUMENTS:
    #   $atoms - a hash of atomic coordinates
    #RETURNS:
    #   the delta torsion in degrees
    #ERRORS RAISED:
    #   will croak if the required atoms are not present
    my $atoms = shift;
    
    unless (defined $atoms->{"C5'"} and defined $atoms->{"C4'"} and defined $atoms->{"C3'"} and defined $atoms->{"O3'"}) {
        croak "Missing atom for calculating delta";
    }
    
    return calcTorsion($atoms->{"C5'"}, $atoms->{"C4'"}, $atoms->{"C3'"}, $atoms->{"O3'"});
}

sub calcEpsilon {
    #calculate beta
    #ARGUMENTS:
    #   $atoms - a hash of atomic coordinates
    #            the P of the next nucleotide should be named P+1
    #RETURNS:
    #   the epsilon torsion in degrees
    #ERRORS RAISED:
    #   will croak if the required atoms are not present
    
    my $atoms = shift;
    
    unless (defined $atoms->{"C4'"} and defined $atoms->{"C3'"} and defined $atoms->{"O3'"} and defined $atoms->{"P+1"}) {
        croak "Missing atom for calculating epsilon";
    }
    
    return calcTorsion($atoms->{"C4'"}, $atoms->{"C3'"}, $atoms->{"O3'"}, $atoms->{"P+1"});
}

sub calcZeta {
    #calculate beta
    #ARGUMENTS:
    #   $atoms - a hash of atomic coordinates
    #            the P and O5' of the next nucleotide should be named P+1 and O5'+1
    #RETURNS:
    #   the zeta torsion in degrees
    #ERRORS RAISED:
    #   will croak if the required atoms are not present
    
    my $atoms = shift;
    
    unless (defined $atoms->{"C3'"} and defined $atoms->{"O3'"} and defined $atoms->{"P+1"} and defined $atoms->{"O5'+1"}) {
        croak "Missing atom for calculating zeta";
    }
    
    return calcTorsion($atoms->{"C3'"}, $atoms->{"O3'"}, $atoms->{"P+1"}, $atoms->{"O5'+1"});
}


sub calcChi {
    #calculate chi
    #ARGUMENTS:
    #   $atoms - a hash of atomic coordinates
    #RETURNS:
    #   the chi torsion in degrees
    #ERRORS RAISED:
    #   will croak if the required atoms are not present
    
    my $atoms   = shift;
    
    my ($baseAtom1, $baseAtom2);
    if (defined $atoms->{"N9"}) {
        $baseAtom1 = "N9";
        $baseAtom2 = "C4";
    } else {
        $baseAtom1 = "N1";
        $baseAtom2 = "C2";
    }
    
    unless (defined $atoms->{"O4'"} and defined $atoms->{"C1'"} and defined $atoms->{$baseAtom1} and defined $atoms->{$baseAtom2}) {
        croak "Missing atom for calculating chi";
    }
    
    return calcTorsion($atoms->{"O4'"}, $atoms->{"C1'"}, $atoms->{$baseAtom1}, $atoms->{$baseAtom2});
}

sub calcXi {
    #calculate xi (the angle (not torsion) between the sugar and the base)
    #ARGUMENTS:
    #   $atoms - a hash of atomic coordinates
    #RETURNS:
    #   the xi torsion in degrees
    #ERRORS RAISED:
    #   will croak if the required atoms are not present
    
    my $atoms   = shift;
    
    my $baseAtom;
    if (defined $atoms->{"N9"}) {
        $baseAtom = "N9";
    } else {
        $baseAtom = "N1";
    }
    
    unless (defined $atoms->{"C4'"} and defined $atoms->{"O4'"} and defined $atoms->{"C1'"} and defined $atoms->{$baseAtom}) {
        croak "Missing atom for calculating xi (the sugar/base angle)";
    }
    
    return calcTorsion($atoms->{"C4'"}, $atoms->{"O4'"}, $atoms->{"C1'"}, $atoms->{$baseAtom});
}


1;
