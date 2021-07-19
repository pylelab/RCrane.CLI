#!/usr/bin/perl

#Functions to build phosphoryl oxygens
#if you have an O3' atom from the previous nucleotide, use buildPhosOxy
#if you don't, use buildInitPhosOxy

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

package BuildPhosOxy;

use strict;
use English '-no_match_vars';
use Exporter;
use Math::Trig;
use Data::Dumper;

use StrucCalc qw(minus dotProd plus scalarProd magnitude crossProd);

use constant { PHOSBONDLENGTH => 1.485,
               PHOSBONDANGLE  => 119.6,
               INITPHOSANGLE  => 108     };   #the angle between the O5' and the phosphoryl oxygens in the first phosphate of a segment

our @ISA       = ('Exporter');
our @EXPORT    = qw(buildPhosOxy buildInitPhosOxy);


sub buildPhosOxy {
    #build phosphoryl oxygens for the current nucleotide
    #this function requires the coordinates for the O3' atom from the previous nucleotide
    #ARGUMENTS:
    #   $curResAtoms  - a hash of coordinates for current residue atoms
    #OPTIONAL ARGUMENTS:
    #   $prevResAtoms - a hash of coordinates for previous residue atoms
    #                   if this is not provided, $curResAtoms must contain an O3'-1 entry with coordinates for the previous O3' atom
    #RETURNS:
    #   #curResAtoms  - a hash of coordinates for current residue atoms, including OP1 and OP2
    
    my $curResAtoms  = shift;
    my $prevResAtoms = shift;
    
    $curResAtoms = {(%{$curResAtoms})}; #do a shallow copy
    #we don't modify $prevResAtoms, so we don't need to bother with a shallow copy of that
    
    my $O3;
    if (defined $prevResAtoms->{"O3'"}) {
        $O3 = $prevResAtoms->{"O3'"};
    } else {
        $O3 = $curResAtoms->{"O3'-1"};
    }
    
    my $P  = $curResAtoms->{"P"};
    my $O5 = $curResAtoms->{"O5'"};
    
    #calculate a line from O5' to O3'
    my $norm = minus($O5, $O3);
    
    #calculate the intersection of a plane (with normal $norm and point $P) and a line (from O5' to O3')
    #using formula from http://local.wasp.uwa.edu.au/~pbourke/geometry/planeline/
    my $i = dotProd($norm, minus($P, $O3)) / dotProd($norm, minus($O5,$O3));
        #       $norm dot ($P - $O3)
        # i = ------------------------
        #       $norm dot ($O5 - $O3)
    
    
    my $interPoint = plus($O3, scalarProd($i, minus($O5, $O3)));
        #$interPoint = $O3 + $u($O5 - $O3)
    
    #print Dumper($interPoint);
    
    
    #move $interPoint so that the distance from $P to $interPoint is 1.485 (the length of the P-O1P bond)
    #we also reflect the point about P
    my $PIline = minus($P, $interPoint); #here's is where the reflection occurs, because we do $P-$interPoint instead of $interPoint-$P
    my $scaledPoint = scalarProd (1/magnitude($PIline) * PHOSBONDLENGTH, $PIline);
    #to get the new point location, we would do $P + $scaledPoint
    #but we need to rotate the point first before translating it back
    
    
    #rotate this new point by 59.8 and -59.8 degrees to determine the phosphoryl oxygen locations
    #we rotate about the axis defined by $norm
    
    my $angle = deg2rad(PHOSBONDANGLE / 2);
    my ($x, $y, $z) = @{$scaledPoint};
    
    my $unitNorm = scalarProd( 1/magnitude($norm), $norm);
    my ($u,$v,$w) = @{$unitNorm};
    
    my $newCoords = [];
    for my $theta ($angle, -$angle) {
        my $cosTheta = cos($theta);
        my $sinTheta = sin($theta);
        
        #perform the rotation, and then add $P to the coordinates
        my $a = $u*$x + $v*$y + $w*$z;
        my $newX = $a*$u + ($x-$a*$u)*$cosTheta + ($v*$z-$w*$y)*$sinTheta + $P->[0];
        my $newY = $a*$v + ($y-$a*$v)*$cosTheta + ($w*$x-$u*$z)*$sinTheta + $P->[1];
        my $newZ = $a*$w + ($z-$a*$w)*$cosTheta + ($u*$y-$v*$x)*$sinTheta + $P->[2];
        
        push(@{$newCoords}, [$newX, $newY, $newZ]);
    }
    
    $curResAtoms->{"OP1"} = $newCoords->[0];
    $curResAtoms->{"OP2"} = $newCoords->[1];
    
    return $curResAtoms;
}

sub buildInitPhosOxy {
    #build phosphoryl oxygens using only the O3' or O5'atom
    #this is useful for the first or last nucleotide after/before a chain break
    #ARGUMENTS:
    #   $curAtoms  - a hash of coordinates for current residue atoms
    #OPTIONAL ARGUMENTS:
    #   $prevAtoms - a hash of coordinates for the previous residue atoms
    #                if $prevAtoms is present, the oxygens will be placed using the O3' and C3' coordinates from $prevAtoms
    #                otherwise, the O5' and C5' coordinates from $curAtoms will be used
    #                so only provide $prevAtoms if you are building phosphoryl oxygens for the last phosphate before a chain break
    #RETURNS:
    #   #atoms  - a hash of coordinates for current residue atoms, including OP1 and OP2
    
    my $curAtoms  = shift;
    my $prevAtoms = shift;
    
    $curAtoms = {(%{$curAtoms})}; #do a shallow copy so we don't modify $atoms
    
    #fetch the atoms that we need
    my $P  = $curAtoms->{"P"};
    
    my ($O, $C);
    if (defined $prevAtoms) {
        $O = $prevAtoms->{"O3'"};
        $C = $prevAtoms->{"C3'"};
    } else {
        $O = $curAtoms->{"O5'"};
        $C = $curAtoms->{"C5'"};
    }
    
    #place atom along the P-O5' bond that is the appropriate distance from P
    my $phosOxy = minus($O, $P);
    $phosOxy = scalarProd( 1/magnitude($phosOxy) * PHOSBONDLENGTH, $phosOxy);
    
    
    #define plane with C5'-O5'-P
    my $norm = crossProd( minus($C,$O), minus($P,$O));
    $norm = scalarProd( 1/magnitude($norm), $norm);
    
    
    #rotate dummy atom in plane about P by INITPHOSANGLE (which is the appropriate O5'-O1P angle)
    #do it in a block since we declare lots of temporary values that we don't want to accidentally access later
    {
        my ($x,$y,$z) = @{$phosOxy};
        my ($u,$v,$w) = @{$norm};
        my $theta = 0-deg2rad(INITPHOSANGLE); #we use the negative angle so that we rotate away from the C5'
        
        my $cosTheta = cos($theta);
        my $sinTheta = sin($theta);
        
        #perform the rotation, and then add $P to the coordinates
        my $a = $u*$x + $v*$y + $w*$z;
        $phosOxy->[0] = $a*$u + ($x-$a*$u)*$cosTheta + ($v*$z-$w*$y)*$sinTheta;
        $phosOxy->[1] = $a*$v + ($y-$a*$v)*$cosTheta + ($w*$x-$u*$z)*$sinTheta;
        $phosOxy->[2] = $a*$w + ($z-$a*$w)*$cosTheta + ($u*$y-$v*$x)*$sinTheta;
    }
    
    
    #rotate dummy atom about O5'-P axis by 59.8 and -59.8 degrees
    $norm = minus($O, $P);
    $norm = scalarProd( 1/magnitude($norm), $norm);
    
    my ($x,$y,$z) = @{$phosOxy};
    my ($u,$v,$w) = @{$norm};
    my $angle = deg2rad(PHOSBONDANGLE / 2);
    
    my $newCoords = [];
    for my $theta ($angle, -$angle) {
        my $cosTheta = cos($theta);
        my $sinTheta = sin($theta);
        
        #perform the rotation, and then add $P to the coordinates
        my $a = $u*$x + $v*$y + $w*$z;
        my $newX = $a*$u + ($x-$a*$u)*$cosTheta + ($v*$z-$w*$y)*$sinTheta + $P->[0];
        my $newY = $a*$v + ($y-$a*$v)*$cosTheta + ($w*$x-$u*$z)*$sinTheta + $P->[1];
        my $newZ = $a*$w + ($z-$a*$w)*$cosTheta + ($u*$y-$v*$x)*$sinTheta + $P->[2];
        
        push(@{$newCoords}, [$newX, $newY, $newZ]);
    }
    
    $curAtoms->{"OP1"} = $newCoords->[0];
    $curAtoms->{"OP2"} = $newCoords->[1];
    
    return $curAtoms;
}

1;