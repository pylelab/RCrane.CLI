#!/usr/bin/perl

#a package containing two useful variables:
#   $puckerList - a hash of all rotamers and their starting and ending puckers
#   $rotList    - a list of all rotamers

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

package PuckerList;

use strict;
use English '-no_match_vars';

use Exporter;
our @ISA       = ('Exporter');
our @EXPORT    = qw($puckerList);
our @EXPORT_OK = qw($puckerList $rotList);

our $puckerList = { '1a' => [3,3],
                    '1m' => [3,3],
                    '1L' => [3,3],
                    '&a' => [3,3],
                    '7a' => [3,3],
                    '3a' => [3,3],
                    '9a' => [3,3],
                    '1g' => [3,3],
                    '7d' => [3,3],
                    '3d' => [3,3],
                    '5d' => [3,3],
                    '1e' => [3,3],
                    '1c' => [3,3],
                    '1f' => [3,3],
                    '5j' => [3,3],
                    
                    '1b' => [3,2],
                    '1[' => [3,2],
                    '3b' => [3,2],
                    '1z' => [3,2],
                    '5z' => [3,2],
                    '7p' => [3,2],
                    '1t' => [3,2],
                    '5q' => [3,2],
                    '1o' => [3,2],
                    '7r' => [3,2],
                    
                    '2a' => [2,3],
                    '4a' => [2,3],
                    '0a' => [2,3],
                    '#a' => [2,3],
                    '4g' => [2,3],
                    '6g' => [2,3],
                    '8d' => [2,3],
                    '4d' => [2,3],
                    '6d' => [2,3],
                    '2h' => [2,3],
                    '4n' => [2,3],
                    '0i' => [2,3],
                    '6n' => [2,3],
                    '6j' => [2,3],
                    
                    '2[' => [2,2],
                    '4b' => [2,2],
                    '0b' => [2,2],
                    '4p' => [2,2],
                    '6p' => [2,2],
                    '4s' => [2,2],
                    '2o' => [2,2],
                    
                    '5n' => [3,3],
                    '5p' => [3,2],
                    '5r' => [3,2],
                    '2g' => [2,3],
                    '0k' => [2,3],
                    '2z' => [2,2],
                    '2u' => [2,2],

                 };

our $rotList = [keys %{$puckerList}];

1;