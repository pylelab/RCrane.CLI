#!/usr/bin/perl

#a module for calculating the minimum distance between two angles, keeping the appropriate sign

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


package AngleDist;

use Exporter;
our @ISA       = ('Exporter');
our @EXPORT    = qw(angleDist);
our @EXPORT_OK = qw(angleDist fmodpos);


use strict;
use English '-no_match_vars';
use POSIX 'fmod';

sub angleDist {
    #calculate the minimum distance between two angles, keeping the appropriate sign
    #ARGUMENTS:
    #   $num1 - the first measurement
    #   $num2 - the second measurement
    #RETURNS:
    #   the distance between these two angles in the range [-180, 180)
    
    my $num1 = shift;
    my $num2 = shift;
    
    return fmodpos($num1 - $num2 + 180, 360) - 180;
}

sub fmodpos {
    #calculate the modulus in the range of 0 to N
    #ARGUMENTS:
    #   $num - the number to take the modulus of
    #   $mod - the modulus base
    #RETURNS:
    #   $num mod $mod
    
    my $num = shift;
    my $mod = shift;
    
    my $res = fmod($num, $mod);
    if ($res < 0) {
        $res += $mod;
    }
    
    return $res;
}