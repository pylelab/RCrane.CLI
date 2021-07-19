#!/usr/bin/perl

#module to build a sugar onto a base
#This module will only build sugars in a default orientation, which means that the sugars will not
#necessarily be positioned appropriately.  Use FullSuiteMinimizer.pm for that.

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

package BuildInitSugar;

use strict;
use English '-no_match_vars';
use Math::Trig qw(deg2rad acos);
use Carp;

use PDB;
use StrucCalc qw(minus crossProd magnitude scalarProd plus dotProd calcTorsion);

use constant STARTINGCHI => -150; #the chi to build the sugars with


sub new {
    #create and initialize a BuildInitSugar object
    #ARGUMENTS:
    #   $c3pStruc - the filename of a PDB file containing a sample C3' endo pucker sugar
    #   $c2pStruc - the filename of a PDB file containing a sample C2' endo pucker sugar
    #RETURNS:
    #   $self     - a blessed BuildInitSugar object

    my $class       = shift;
    my $c3pStruc    = shift;
    my $c2pStruc    = shift;
    
    #read in the sugar structures
    my $c3pPdb   = new PDB($c3pStruc);
    my $c3pAtoms = $c3pPdb->firstRes->atoms;
    
    my $c2pPdb   = new PDB($c2pStruc);
    my $c2pAtoms = $c2pPdb->firstRes->atoms;
    
    #translate the sugars so that the C1' atom is at zero (so that we don't have to do it every time we align a sugar)
    for my $curSugar ($c3pAtoms, $c2pAtoms) {
        my $c1p = $curSugar->{"C1'"};
        for my $curAtom (values(%{$curSugar})) {
            $curAtom = minus($curAtom, $c1p);
        }
    }
    
    #create the object
    my $self = { c3pStrucFilename   => $c3pStruc,
                 c2pStrucFilename   => $c2pStruc,
                 c3pAtoms           => $c3pAtoms,
                 c2pAtoms           => $c2pAtoms,
               };
    bless $self, $class;
    
    return $self;
}

sub buildInitSugar {
    #build a sugar (and other backbone atoms) onto a base
    #ARGUMENTS:
    #   $baseAtoms - a hash of atom coordinates for base atoms
    #   $pucker    - the desired pucker of the sugar (either C3' or C2')
    #RETURNS:
    #   a hash of atom coordinates for base and sugar atoms
    
    my $self      = shift;
    my $baseAtoms = shift;
    my $pucker    = shift;
    
    #figure out which sugar pucker structure to use
    $pucker = lc($pucker);
    my $sugarAtoms;
    if      ($pucker eq "c3p" or $pucker eq "c3" or $pucker eq "c3'" or $pucker == 3) {
        $sugarAtoms = $self->{c3pAtoms};
    } elsif ($pucker eq "c2p" or $pucker eq "c2" or $pucker eq "c2'" or $pucker == 2) {
        $sugarAtoms = $self->{c2pAtoms};
    } else {
        croak "buildInitSugar called with unrecognized sugar pucker ($pucker)";
    }
    
    #move the sugar so its attached to the base
    my $alignedSugar = $self->alignSugar($baseAtoms, $sugarAtoms);
    
    #return the new sugar atoms along with all atoms that we were given
    return {(%{$alignedSugar}, %{$baseAtoms})}
}


sub alignSugar {
    #align a sugar to a base
    #ARGUMENTS:
    #   $baseAtoms      - a hash of atom coordinates for the base atoms
    #   $sugarAtoms     - a hash of atom coordinates for the sugar atoms
    #RETURNS:
    #   $alignedStruc - a hash of the aligned rotamer, in the form of {atom name => [x,y,z]}
    
    my $self            = shift;
    my $baseAtoms       = shift;
    my $sugarAtoms      = shift;
    
    #figure out which N and C atoms to use (if the nucleotide is a purine or a pyrimidine)
    my ($Natom, $Catom);
    if (defined $baseAtoms->{"N9"}) {
        $Natom = "N9";
        $Catom = "C4";
    } else {
        $Natom = "N1";
        $Catom = "C2";
    }
    
    #rotate the sugar so the glycosidic bond is at the appropriate angle
    #calculate an axis for the rotation
    my $translatedBaseN = minus($baseAtoms->{$Natom}, $baseAtoms->{"C1'"});
    my $sugarN = $sugarAtoms->{$Natom};
    
    my $axis  = crossProd(   $sugarN, $translatedBaseN);
    my $angle = calcTorsion($translatedBaseN, $axis, [0,0,0], $sugarN);
    
    unless ($angle == 0 or magnitude($axis) == 0) {
        $sugarAtoms = rotateAtoms($sugarAtoms, $angle, $axis);
    }
    
    
    #rotate the sugar so that chi is appropriate
    my $translatedBaseC = minus($baseAtoms->{$Catom}, $baseAtoms->{"C1'"});
    my $curChi = calcTorsion($translatedBaseC, $translatedBaseN, [0,0,0], $sugarAtoms->{"O4'"});
    $sugarAtoms = rotateAtoms($sugarAtoms, $curChi - STARTINGCHI, $translatedBaseN);
    
    
    #translate the sugar to the C1' atom of the base
    for my $curAtom (values(%{$sugarAtoms})) {
        $curAtom = plus($curAtom, $baseAtoms->{"C1'"});
    }
    
    
    #remove unneccessary atoms from $sugarAtoms
    delete(@{$sugarAtoms}{qw(N1 N9 C2 C4)});
    
    return $sugarAtoms;
}


sub rotateAtoms {
    #rotate atoms about a specified axis
    #Note: this function will *not* translate atoms to the origin before rotation
    #ARGUMENTS:
    #   $atoms - a hash of atomic coordinates for all atoms to rotate
    #   $angle - the angle to rotate by
    #   $axis  - the axis to rotate through
    #RETURNS:
    #   $rotatedAtoms - a hash of atomic coordinates for rotated atoms
    
    my $atoms = shift;
    my $angle = shift;
    my $axis = shift;
    
    #normalize the axis
    $axis = scalarProd( 1/magnitude($axis), $axis);
    
    $angle = deg2rad($angle);
    my $cosTheta = cos($angle);
    my $sinTheta = sin($angle);
    
    my ($u,$v,$w) = @{$axis};
    my $rotatedAtoms = {};
    
    #go through each atom and rotate it
    for my $curAtom (keys %{$atoms}) {
        my ($x,$y,$z) = @{$atoms->{$curAtom}};
        my $a = $u*$x + $v*$y + $w*$z;
        my $newX = $a*$u + ($x-$a*$u)*$cosTheta + ($v*$z-$w*$y)*$sinTheta;
        my $newY = $a*$v + ($y-$a*$v)*$cosTheta + ($w*$x-$u*$z)*$sinTheta;
        my $newZ = $a*$w + ($z-$a*$w)*$cosTheta + ($u*$y-$v*$x)*$sinTheta;
        $rotatedAtoms->{$curAtom} = [$newX, $newY, $newZ];
    }
    
    return $rotatedAtoms;
}

1;