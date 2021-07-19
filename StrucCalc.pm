#!/usr/bin/perl

# This package is for calculations performed on the structure (such as distances and torsions)
# (not for calculating structures, which would be considerably more difficult)
# all coordinates handled by this package are assumed to be ref. to arrays containing three numbers
# (i.e. x, y, and z coordinates)

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

package StrucCalc;

use strict;
use English '-no_match_vars';
use Exporter;
use Math::Trig qw(rad2deg acos);
use Data::Dumper;

our @ISA       = ('Exporter');
our @EXPORT    = qw(calcDist calcTorsion calcDistToLine);
our @EXPORT_OK = qw(calcDist calcTorsion calcDistToLine minus crossProd plus scalarProd magnitude dotProd angle);


sub calcDist {
    #calculate the distance between two points
    #ARGUMENTS:
    #   $a - the first point
    #   $b - the second point
    #RETURNS:
    #   the distance between the two points
    
    my $a = shift;
    my $b = shift;
    
    return magnitude(minus($a, $b));
}

sub calcTorsion {
    #calculate the torsional angle between four points
    #ARGUMENTS:
    #   $atom1   - point 1
    #   $atom2   - point 2
    #   $atom3   - point 3
    #   $atom4   - point 4
    #RETURNS
    #   $torsion - the angle between points 1 and 4 about the axis defined by points 2 and 3
    
    my $atom1 = shift;
    my $atom2 = shift;
    my $atom3 = shift;
    my $atom4 = shift;
    
    my $vector1 = minus($atom2, $atom1);
    my $vector2 = minus($atom3, $atom2);
    my $vector3 = minus($atom4, $atom3);
    
    #torsion angle = atan2 (|b2|b1 . (b2 x b3), (b1 x b2) . (b2 x b3)
    #source: http://en.wikipedia.org/wiki/Dihedral_angle
    
    my $y1 = scalarProd( magnitude($vector2), $vector1);
    my $y2 = crossProd($vector2, $vector3);
    my $y  = dotProd($y1, $y2);
    
    my $x1 = crossProd($vector1, $vector2);
    my $x2 = $y2;
    my $x  = dotProd($x1, $x2);
    
    my $torsion = atan2($y, $x);
    $torsion = rad2deg($torsion);
    $torsion += 360 if $torsion < 0; #put the torsion in the [0,360) range
    
    return $torsion;
}
    

sub calcDistToPlane {
    #not yet implemented
}

sub calcDistToLine {
    #calculate the distance from a point to a line, where the line is defined by two points on the line
    #ARGUMENTS:
    #   $pointP - the point not on the line
    #   $pointC - point 1 that defined the line
    #   $pointN - point 2 that defined the line
    #RETURNS:
    #   $dist   - the distance from $pointP to a line through $pointC and $pointN
    
    my $pointP  = shift; #the point not on the line
    my $pointC  = shift;
    my $pointN  = shift;
    
    #the distance from point P to line CN is
    #               |(C-N) x (N-P)|
    #   distance = -----------------
    #                    |C-N|
    
    my $cMinusN = minus($pointC, $pointN);
    my $nMinusP = minus($pointN, $pointP);
    
    my $crossProd = crossProd($cMinusN, $nMinusP);
    
    my $dist = magnitude($crossProd) / magnitude($cMinusN);
    
    return $dist;
}

sub magnitude {
    #calculate the magnitude of a vector
    #ARGUMENTS:
    #   $a  - the vector
    #RETURNS:
    #   the magnitude of $a (the distance from (0,0,0) to $a)
    
    my $a = shift;
    
    return sqrt( $a->[0]**2 + $a->[1]**2 + $a->[2]**2);
}

sub crossProd {
    #calculate the cross product of two vectors
    #ARGUMENTS:
    #   $a - vector 1
    #   $b - vector 2
    #RETURNS:
    #   $a cross $b
    
    my $a = shift;
    my $b = shift;
    
    my ($ax, $ay, $az) = @{$a};
    my ($bx, $by, $bz) = @{$b};
    
    return [ $ay*$bz - $az*$by,
             $az*$bx - $ax*$bz,
             $ax*$by - $ay*$bx ];
    
}

sub dotProd {
    #calculate the dot product of two vectors
    #ARGUMENTS:
    #   $a - vector 1
    #   $b - vector 2
    #RETURNS:
    #   $a . $b
    
    my $a = shift;
    my $b = shift;
    
    return ($a->[0]*$b->[0] + $a->[1]*$b->[1] + $a->[2]*$b->[2]);
}

sub scalarProd {
    #calculate the scalar product of a vector and a scalar
    #ARGUMENTS:
    #   Note that the order of arguments for this function is unimportant
    #   $scalar - the scalar
    #   $vector - the vector
    #RETURNS:
    #   the scalar product of $scalar x $vector
    
    my $scalar = shift;
    my $vector = shift;
    
    #flip the variables if necessary
    if (ref $scalar eq 'ARRAY') {
        ($scalar, $vector) = ($vector, $scalar);
    }
    
    return [ $scalar * $vector->[0], $scalar * $vector->[1], $scalar * $vector->[2] ]
}

sub minus {
    #subtract two vectors
    #ARGUMENTS:
    #   $a - vector 1
    #   $b - vector 2
    #RETURNS:
    #   $a - $b
    
    my $a = shift;
    my $b = shift;
    
    return [ $a->[0] - $b->[0],  $a->[1] - $b->[1], $a->[2] - $b->[2] ];
}

sub plus {
    #add two vectors
    #ARGUMENTS:
    #   $a - vector 1
    #   $b - vector 2
    #RETURNS:
    #   $a + $b
    
    my $a = shift;
    my $b = shift;
    
    return [ $a->[0] + $b->[0],  $a->[1] + $b->[1], $a->[2] + $b->[2] ];
}

sub angle {
    #calculate the angle between three points
    #ARGUMENTS:
    #   $a - point 1
    #   $b - point 2
    #   $c - point 3
    #RETURNS:
    #   angle <abc (i.e. the angle about point b
    
    my $a = shift;
    my $b = shift;
    my $c = shift;
    
    #calculate the two vectors
    my $vec1 = minus($a, $b);
    my $vec2 = minus($c, $b);
    
    #normalize the two vectors
    $vec1 = scalarProd(1/magnitude($vec1), $vec1);
    $vec2 = scalarProd(1/magnitude($vec2), $vec2);
    
    return rad2deg(acos(dotProd($vec1, $vec2)));
    
}
    
1;