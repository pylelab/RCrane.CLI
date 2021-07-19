#!/usr/bin/perl

# A class for accessing data about a suite.
# Note that the data is actually stored in the associated Residue objects

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

package Suite;

use strict;
use English '-no_match_vars';
use StrucCalc;


sub new {
    #create a Suite object
    #ARGUMENTS:
    #   $startingRes - the Residue object containing the starting nucleotide of the suite
    #   $endingRes   - the Residue object containing the ending nucleotide of the suite
    
    my $class       = shift;
    my $startingRes = shift;
    my $endingRes   = shift;
    
    my $self = { startingRes => $startingRes,
                 endingRes   => $endingRes     };
    bless $self, $class;    
    return $self;
}


sub sugarDist {
    #get the sugar-sugar distance (C1'-C1' distance) of the suite
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the sugar-sugar distance
    
    my $self = shift;
    
    return $self->endingRes->startingSugarDist;
}


sub fullNumber {
    #get the two number version of the suite number (such as 111-112)
    #see number for the one number version of the suite number
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the two number version of the suite number

    my $self = shift;
    return ($self->{startingRes}->number . "-" . $self->{endingRes}->number);
}


sub theta {
    #get the theta value of the suite
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the theta value

    my $self = shift;
    return $self->startingRes->theta;
}

sub eta {
    #get the eta value of the suite
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the eta value

    my $self = shift;    
    return $self->endingRes->eta;
}


sub nextConnectedSuite {
    #get the next suite, but only if it is connected to this one (i.e. if they share a sugar)
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a Suite object containing the next connected suite
    #   undef if there is no connected next suite
    
    my $self = shift;
    
    if ($self->connectedToNext) {
        return new Suite($self->{endingRes}, $self->{endingRes}->nextRes);
    } else {
        return undef;
    }
}

sub prevConnectedSuite {
    #get the previous suite, but only if it is connected to this one (i.e. if they share a sugar)
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a Suite object containing the previous connected suite
    #   undef if there is no connected previous suite
    
    my $self = shift;
    
    if ($self->connectedToPrev) {
        return new Suite($self->{startingRes}->prevRes, $self->{startingRes});
    } else {
        return undef;
    }
}

sub connectedToNext {
    #is this suite connected to the next suite in the chain? (i.e. do they share a sugar?)
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   true if this suite is connected to the next one, false if it isn't
    
    my $self = shift;
    
    if ($self->{endingRes}->connectedToNext and $self->{endingRes}->nextRes->atoms->{"P"} and $self->{endingRes}->nextRes->atoms->{"C1'"}) {
        return 1;
    } else {
        return 0;
    }
}

sub connectedToPrev {
    #is this suite connected to the previous suite in the chain? (i.e. do they share a sugar?)
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   true if this suite is connected to the previous one, false if it isn't
    
    my $self = shift;
    
    if ($self->{startingRes}->connectedToPrev and $self->{startingRes}->atoms->{"P"} and $self->{startingRes}->prevRes->atoms->{"C1'"}) {
        return 1;
    } else {
        return 0;
    }
}


sub nextSuite {
    #find the next suite, regardless of whether or not it is connected to this one
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a Suite object containing the next suite
    #   undef if there are no more suites in the chain
    
    my $self = shift;
    
    my $curRes = $self->endingRes;
    
    while (defined $curRes) {
        if ($curRes->connectedToNext and $curRes->atoms->{"C1'"} and $curRes->nextRes->atoms->{"P"} and $curRes->nextRes->atoms->{"C1'"}) {
            return new Suite($curRes, $curRes->nextRes);
        }
        $curRes = $curRes->nextRes;
    }
    
    return undef;
}

sub prevSuite {
    #find the prev suite, regardless of whether or not it is connected to this one
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a Suite object containing the previous suite
    #   undef if there are no prior suites in the chain
    
    my $self = shift;
    
    my $curRes = $self->startingRes;
    
    while (defined $curRes) {
        if ($curRes->connectedToPrev and $curRes->atoms->{"C1'"} and $curRes->atoms->{"P"} and $curRes->prevRes->atoms->{"C1'"}) {
            return new Suite($curRes->prevRes, $curRes);
        }
        $curRes = $curRes->prevRes;
    }
    
    return undef;
}


sub structure {
    #get the Structure object that this suite is part of
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the Structure object that this suite is part of
    
    my $self = shift;
    return $self->{startingRes}->structure;
}

sub number {
    #get the one number version of the suite number (i.e. the number of the ending residue of the suite)
    #see fullNumber for the two number version of the suite number
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the one number version of the suite number

    my $self = shift;
    return $self->{endingRes}->number;
}


sub atom {
    #get atomic coordinates for a given atom
    #note that this function does not return an lvalue (i.e. you cannot modify coordinates through this function)
    #ARGUMENTS:
    #   $atomName - the name of the atom to get coordinates for
    #RETURNS:
    #   coordinates of the atom
    
    my $self     = shift;
    my $atomName = shift;
    
    if (substr($atomName, -2) eq "-1") {
        return $self->{startingRes}->atoms->{substr($atomName, 0, -2)};
    } elsif ($atomName eq "03'") {
        return $self->{startingRes}->atoms->{$atomName};
    } else {
        return $self->{endingRes}->atoms->{$atomName};
    }
}


sub startingRes {   $ARG[0]->{startingRes};   }
sub endingRes   {   $ARG[0]->{endingRes};     }


1;