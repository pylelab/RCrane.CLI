#!/usr/bin/perl

#a function to remove trailing and leading spaces from a string

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

package Strip;

use Exporter;

our @ISA       = ('Exporter');
our @EXPORT    = ('strip');


sub strip {
    #remove trailing and leading spaces from a string
    #ARGUMENTS:
    #   a string
    #RETURNS:
    #   the string with trailing and leading spaces removed
    
    $_[0] =~ /^\s*(.*?)\s*$/ or return "";
    return $1;
}

1;