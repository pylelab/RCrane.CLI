#!/usr/bin/perl

# A class for storing data about a chain of nucleotides.

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

package Chain;

use strict;
use English '-no_match_vars';
use Residue;
use Suite;


sub new {
    #create a new chain object
    #ARGUMENTS:
    #   $structure  - what structure this residue is part of
    #   $id         - what the chain ID letter is
    #RETURNS:
    #   $self       - a blessed Chain object
    
    my $class     = shift;
    my $structure = shift; #what structure this residue is part of
    my $id        = shift; #what the chain ID letter is
    
    
    my $self = { resHash    => {},           #a hash of residue number -> location in the resList/suiteList array
                 resList    => [],           #an array containing the Residue objects
                 structure  => $structure,
                 id         => $id
                };
    
    bless $self, $class;
    
    return $self;
}

sub addRes {
    #add a residue to the end of the chain
    #ARGUMENTS:
    #   this function can be passed either a Residue object, or the data required to create a Residue object
    #   Residue-object mode:
    #       $res    - the Residue object to be added to the chain
    #   non-Residue-object mode
    #       $resNum - the residue number (which may include an insertion code)
    #       $type   - what type of residue this is (i.e. A, C, G, or U)
    #       $atoms  - a ref to hash of atom names and coordintes: { atomName => [x, y, z]}
    #RETURNS:
    #   none
    #EFFECTS:
    #   adds a new residue onto the end of the chain
    
    my $self = shift;
    
    my $res; #a Residue object representing the new residue
    
    if (@ARG == 1 and ref($ARG[0]) eq "Residue") {
        #if we are passed a Residue object, then use that
        $res  = shift;
        $res->{structure} = $self->{structure};
        $res->{chain}     = $self->{id};
        $res->{seqNum}    = $self->numRes;
    } else {
        #otherwise, assume the rest of the arguments are details about the residue
        $res = new Residue($self->{structure}, $self, $self->numRes, @ARG);
    }
    
    #create an entry in the residue hash
    $self->{resHash}->{$res->number} = $self->numRes;
    
    #add the new residue object to the end of the seqList array
    push(@{$self->{resList}}, $res);
}

sub numRes {
    #return the number of residues in this chain
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   the number of residues in this chain
    
    my $self = shift;
    
    return scalar(@{$self->{resList}});
}

sub firstRes {
    #get the first residue of the chain
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   the first residue of the chain
    
    my $self = shift;
    return $self->{resList}->[0];
}

sub lastRes {
    #get the last residue of the chain
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   the last residue of the chain
    
    my $self = shift;
    return $self->{resList}->[$self->numRes - 1];
}

sub res {
    #get the residue with the given number (which may include an insertion code)
    #ARGUMENTS:
    #   $resNum - the residue number
    #RETURNS:
    #   the requested residue, or undef if the residue does not exist
    
    my $self   = shift;
    my $resNum = shift;
    
    return $self->{resList}->[$self->{resHash}->{$resNum}];    
}

sub hasRes {
    #check existence of a specified residue
    #ARGUMENTS:
    #   $resNum - the residue number
    #RETURNS:
    #   true if the requested residue exists, false otherwise
    
    my $self   = shift;
    my $resNum = shift;
    
    if (defined($self->{resHash}->{$resNum})) {
        return 1;
    } else {
        return 0;
    }
}
    

sub suite {
    #return the suite with the given number (which may include an insertion code)
    #ARGUMENTS:
    #   $suiteNum - the suite number (should be a single number - i.e. 124, not 123-124)
    #RETURNS:
    #   the requested suite, or undef if the residue does not exist
    
    my $self  = shift;
    my $suiteNum = shift;
    
    
    my $seqNum = $self->{resHash}->{$suiteNum};
    
    if ($seqNum == 0) {
        return undef;
    } elsif (not defined $self->{resList}->[$seqNum]) {
        return undef;
    } elsif (not $self->{resList}->[$seqNum]->connectedToPrev) {
        return undef;
    } else {
        return new Suite ($self->{resList}->[$seqNum - 1], $self->{resList}->[$seqNum]);
    }
}


sub firstSuite {
    #return the first suite of this chain
    #this may not necessarily correspond to the first two nucleotides if they aren't connected
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   the first suite of the chain
    
    my $self  = shift;
    
    for (my $i = 1; $i < $self->numRes; $i++) {
        if ($self->{resList}->[$i]->connectedToPrev) {
            return new Suite ($self->{resList}->[$i - 1], $self->{resList}->[$i]);
        }
    }
    
    #if this chain doesn't have any suites
    return undef;
}

sub lastSuite {
    #return the last suite of this chain
    #this may not necessarily correspond to the last two nucleotides if they aren't connected
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   the last residue of the chain
    
    my $self  = shift;
    
    for (my $i = ($self->numRes - 1); $i >= 1; $i--) {
        if ($self->{resList}->[$i]->connectedToPrev) {
            return new Suite ($self->{resList}->[$i - 1], $self->{resList}->[$i]);
        }
    }
    
    #if this chain doesn't have any suites
    return undef;
}

sub structure   {   $ARG[0]->{structure};   }
sub id          {   $ARG[0]->{id};          }


1;