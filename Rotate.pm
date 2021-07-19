#!/usr/bin/perl

#a module for rotating atoms and the sugar

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

#   these functions would probably be better (and likely faster) to do with a good matrix library
#   but the math is easy enough that implementing them in Perl is far simpler

package Rotate;

use strict;
use English '-no_match_vars';
use Math::Trig;
#use Data::Dumper::Simple;

use StrucCalc qw(minus magnitude scalarProd plus calcDist calcTorsion);

use Exporter;
our @ISA       = ('Exporter');
our @EXPORT    = qw(rotateSugar rotateAtoms);
our @EXPORT_OK = qw(rotateSugar rotateAtoms);




sub rotateSugar {
    #rotate the sugar and return all coordinates
    #ARGUMENTS:
    #   $atoms      - a hash containing atomic coordinates
    #   $angleChi   - the number of degrees to rotate chi by
    #                 note that this number is relative to the current chi torsion
    #                 (i.e. this is not the angle to rotate chi *to*)
    #OPTIONAL ARGUMENTS:
    #   $onlyC3C4   - how many atoms to rotate (this is useful for keeping runtime down if you don't need all atoms rotated)
    #                 if zero (or not given), then all sugar atoms will be rotated
    #                 if 1, then only the C3', C4', and O5' atoms will be rotated
    #                 if 2, then only the C1', C2', O2', C3', C4', and O4' atoms will be rotated
    #   $angleBase  - the number of degrees to rotate xi (the sugar/base angle) by
    #                 as with chi, this is relative to the current xi angle
    
    my $atoms = shift;
    my $angleChi = deg2rad(shift);
    my $onlyC3C4 = shift;  #whether to rotate all the atoms, or just C3' and C4'
    my $angleBase = deg2rad(shift);
    
    my $atomList;
    if ($onlyC3C4 == 2) {
        $atomList = [qw(C1' C2' O2' C3' C4' O4')];
    } elsif ($onlyC3C4) {
        #$atomList = [qw(C3' C4' O5' O3')];
        $atomList = [qw(C3' C4' O5')];
    } else {
        $atomList = [qw(C1' C2' O2' C3' O3' C4' O4' C5' O5')];
    }
    
    my $C1ploc = $atoms->{"C1'"};
    
    my $translatedAtoms = {};
    #translate all the points so that C1' is at the origin
    for my $curAtom (@{$atomList}) {
        $translatedAtoms->{$curAtom} = minus($atoms->{$curAtom}, $C1ploc);
    }
    
    
    #calculate axis for rotation about chi (a unit vector along the C1'-N9/1 glycosidic bond)
    my $baseAtom = $atoms->{"N9"};
    if (not defined $baseAtom) {
        $baseAtom = $atoms->{"N1"};
    }
    my $axis = minus($C1ploc, $baseAtom);
    $axis = scalarProd($axis, 1/magnitude($axis));
    
    #perform the rotation about chi
    my $cosTheta = cos($angleChi);
    my $sinTheta = sin($angleChi);
    my ($u,$v,$w) = @{$axis};
    my $rotatedAtoms = {};
    
    for my $curAtom (@{$atomList}) {
        my ($x,$y,$z) = @{$translatedAtoms->{$curAtom}};
        my $a = $u*$x + $v*$y + $w*$z;
        my $newX = $a*$u + ($x-$a*$u)*$cosTheta + ($v*$z-$w*$y)*$sinTheta;
        my $newY = $a*$v + ($y-$a*$v)*$cosTheta + ($w*$x-$u*$z)*$sinTheta;
        my $newZ = $a*$w + ($z-$a*$w)*$cosTheta + ($u*$y-$v*$x)*$sinTheta;
        $rotatedAtoms->{$curAtom} = [$newX, $newY, $newZ];
    }
    
    
    #if a rotation of the angle (not torsion) of the sugar relative to the base is desired
    if ($angleBase) {
        #calculate the axis for the rotation (a unit vector along the C1'-O4' bond)
        $axis = minus($rotatedAtoms->{"C1'"}, $rotatedAtoms->{"O4'"});
        $axis = scalarProd($axis, 1/magnitude($axis));
        
        #perform the rotation about chi
        my $cosTheta = cos($angleBase);
        my $sinTheta = sin($angleBase);
        my ($u,$v,$w) = @{$axis};
        
        for my $curAtom (@{$atomList}) {
            my ($x,$y,$z) = @{$rotatedAtoms->{$curAtom}};
            my $a = $u*$x + $v*$y + $w*$z;
            my $newX = $a*$u + ($x-$a*$u)*$cosTheta + ($v*$z-$w*$y)*$sinTheta;
            my $newY = $a*$v + ($y-$a*$v)*$cosTheta + ($w*$x-$u*$z)*$sinTheta;
            my $newZ = $a*$w + ($z-$a*$w)*$cosTheta + ($u*$y-$v*$x)*$sinTheta;
            $rotatedAtoms->{$curAtom} = [$newX, $newY, $newZ];
        }
    }
    
    #translate all the points back
    my $newAtoms = {};
    for my $curAtom (@{$atomList}) {
        $newAtoms->{$curAtom} = plus($rotatedAtoms->{$curAtom}, $C1ploc);
    }
    
    return $newAtoms;
}


sub rotateAtoms {
    #rotate atoms about a given axis
    #ARGUMENTS:
    #   $atoms          - a hash of atom name->[x,y,z]
    #   $atomsToRotate  - an array containing the names of all atoms to rotate
    #   $axisAtoms      - an array containing the names of the two atoms that make up the rotation axis
    #   $angle          - the angle to rotate by in degrees
    #                     note that this angle is relative to the current torsion.  It is NOT an absolute torsion
    #RETURNS:
    #   a hash containing coordinates for all rotated atoms in the form of atom name->[x,y,z]
    
    my $atoms         = shift;
    my $atomsToRotate = shift;
    my $axisAtoms     = shift;
    my $angle         = shift;
    
    #convert the angle to radians
    $angle = deg2rad($angle);

    #one end of the rotation axis needs to be translated to (0,0,0) for the rotation and then translated back
    #so remember it's coordinates
    my $zeroLoc = $atoms->{$axisAtoms->[0]};
        
    #calculate a unit vector along rotation axis
    my $axis = minus($zeroLoc, $atoms->{$axisAtoms->[1]});
    $axis = scalarProd($axis, 1/magnitude($axis));
    
    
    my $translatedAtoms = {};
    #translate all the points so that $zeroLoc is at the origin
    for my $curAtom (@{$atomsToRotate}) {
        $translatedAtoms->{$curAtom} = minus($atoms->{$curAtom}, $zeroLoc);
    }
    
    #perform the rotation
    my $cosTheta = cos($angle);
    my $sinTheta = sin($angle);
    my ($u,$v,$w) = @{$axis};
    my $rotatedAtoms = {};
    
    for my $curAtom (@{$atomsToRotate}) {
        my ($x,$y,$z) = @{$translatedAtoms->{$curAtom}};
        my $a = $u*$x + $v*$y + $w*$z;
        my $newX = $a*$u + ($x-$a*$u)*$cosTheta + ($v*$z-$w*$y)*$sinTheta;
        my $newY = $a*$v + ($y-$a*$v)*$cosTheta + ($w*$x-$u*$z)*$sinTheta;
        my $newZ = $a*$w + ($z-$a*$w)*$cosTheta + ($u*$y-$v*$x)*$sinTheta;
        $rotatedAtoms->{$curAtom} = [$newX, $newY, $newZ];
    }
    
    
    #translate all the points back
    my $newAtoms = {(%{$atoms})};
    for my $curAtom (@{$atomsToRotate}) {
        $newAtoms->{$curAtom} = plus($rotatedAtoms->{$curAtom}, $zeroLoc);
    }
    
    return $newAtoms;
}