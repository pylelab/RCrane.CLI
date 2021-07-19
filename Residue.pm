#!/usr/bin/perl

# A class for storing data about a residue.  Note that this class assumes that
# the residue will not change once it has been initialized, so measurements (eta,
# theta, pperp) will not be updated if the atomic coordinates are changed

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

package Residue;

use strict;
use English '-no_match_vars';
use StrucCalc;


my $sugarAtom = "C1'"; #which sugar atom to use for pseudotorsion and distance calculations (C1' or C4')
use constant NOTCALC => "notcalc"; # a flag to indicate that the value has not yet been calculated
use constant LOOSECONNECTIVITY => 1; #if true, be very generous when deciding if residues are connected


sub new {
    #create a Residue object
    #OPTIONAL ARGUMENTS:
    #   $structure  - a reference to the PDB object that this residue is part of
    #   $number     - what residue number this is
    #   $type       - what type of residue this is (i.e. A, C, G, or U)
    #   $atoms      - a ref to hash of atom names and coordintes: { atomName => [x, y, z]}
    #RETURNS:
    #   a Residue object initialized with any provided data

    my $class     = shift;
    my $structure = shift; #what structure this residue is part of
    my $chain     = shift; #what chain this residue is part of
    my $seqNum    = shift; #what number residue this is counting sequentially
    my $resNum    = shift; #what number residue this is following the counting of the PDB file
                           #(note that this may include an insertion code)
    my $type      = shift;
    my $atoms     = (shift or {});
    
    
    my $self = { atoms                  => $atoms,
                 number                 => $resNum,
                 type                   => $type,
                 structure              => $structure,
                 chain                  => $chain,
                 seqNum                 => $seqNum,
                 etaVal                 => NOTCALC,
                 thetaVal               => NOTCALC,
                 pperpVal               => NOTCALC,
                 phosDistVal            => NOTCALC,
                 startingSugarDistVal   => NOTCALC,
                 connectedToPrev        => NOTCALC };
    bless $self, $class;
    
    return $self;
}


sub eta {
    #get the eta value of the residue
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the eta value
    #NOTE:
    #   The first time this function is called, the value will be calculated
    #   and stored so that it does not need to be calculated again.

    my $self = shift;
    
    if ($self->{etaVal} eq NOTCALC) {
        #make sure that all the necessary atoms exist to calculate eta
        if (not $self->connectedToPrev or
            not $self->connectedToNext or
            not defined $self->prevRes->atoms->{$sugarAtom} or
            not defined $self->atoms->{"P"} or
            not defined $self->atoms->{$sugarAtom} or 
            not defined $self->nextRes->atoms->{"P"}) {
            
            #if they don't, put an undef for eta
            $self->{etaVal} = undef;
        } else {
            #calculate and store eta
            $self->{etaVal} = calcTorsion( $self->prevRes->atoms->{$sugarAtom},
                                           $self->atoms->{"P"},
                                           $self->atoms->{$sugarAtom},
                                           $self->nextRes->atoms->{"P"});
        }
    }
    
    return $self->{etaVal};
}


sub theta {
    #get the theta value of the residue
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the theta value
    #NOTE:
    #   The first time this function is called, the value will be calculated
    #   and stored so that it does not need to be calculated again.
    
    my $self = shift;
    
    if ($self->{thetaVal} eq NOTCALC) {
        
        #make sure that all the necessary atoms exist to calculate theta
        if (not $self->connectedToNext or
            not defined $self->atoms->{"P"} or
            not defined $self->atoms->{$sugarAtom} or 
            not defined $self->nextRes->atoms->{"P"} or
            not defined $self->nextRes->atoms->{$sugarAtom}) {
            
            #if they don't, put an undef for eta
            $self->{thetaVal} = undef;
        } else {
            #calculate and store eta
            $self->{thetaVal} = calcTorsion( $self->atoms->{"P"},
                                             $self->atoms->{$sugarAtom},
                                             $self->nextRes->atoms->{"P"},
                                             $self->nextRes->atoms->{$sugarAtom});
        }
    }
    
    return $self->{thetaVal};
}


sub pperp {
    #get the P-perp value of the residue (the distance from the glycosidic bond to the 3' phosphate)
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the P-perp value
    #NOTE:
    #   The first time this function is called, the value will be calculated
    #   and stored so that it does not need to be calculated again.
    my $self = shift;
    
    if ($self->{pperpVal} eq NOTCALC) {
        #figure out the name of the nitrogen involved in the glycosidic bond
        my $nitrogen;
        if (exists $self->atoms->{"N9"}) {
            $nitrogen = $self->atoms->{"N9"};
        } elsif (exists $self->atoms->{"N1"}) {
            $nitrogen = $self->atoms->{"N1"};
        }
        
        #make sure that all the necessary atoms exist to calculate p-perp
        if (not defined $nitrogen or
            not defined $self->atoms->{"C1'"} or
            not $self->connectedToNext or
            not defined $self->nextRes->atoms->{"P"}) {
            
            $self->{pperpVal} = undef;
        } else {
            $self->{pperpVal} = calcDistToLine( $self->nextRes->atoms->{"P"},
                                                $self->atoms->{"C1'"},
                                                $nitrogen                      );
        }
    }
    
    return $self->{pperpVal};
}

sub phosDist {
    #get the phosphate-phosphate distance of this residue
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the phosphate-phosphate distance
    #NOTE:
    #   The first time this function is called, the value will be calculated
    #   and stored so that it does not need to be calculated again.
    
    my $self = shift;
    
    if ($self->{phosDistVal} eq NOTCALC) {
        if (not $self->connectedToNext or
            not defined $self->atoms->{"P"} or
            not defined $self->nextRes->atoms->{"P"}) {
            
            $self->{phosDistVal} = undef;
        } else {
            $self->{phosDistVal} = calcDist($self->atoms->{"P"},
                                            $self->nextRes->atoms->{"P"});
        }                                
    }
    
    return $self->{phosDistVal};
}

sub nextRes {
    #get the next residue
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a Residue object containing the next residue
    
    my $self = shift;
    return $self->{chain}->{resList}->[$self->{seqNum}+1];
}

sub prevRes {
    #get the previous residue
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a Residue object containing the previous residue
    
    my $self = shift;
    return $self->{chain}->{resList}->[$self->{seqNum}-1];
}

sub connectedToNext {
    #is this residue covalently connected to the next residue in the chain?
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   true if this residue is connected to the next one, false if it isn't
    
    my $self = shift;
    my $nextRes = $self->nextRes;
    
    unless (defined $nextRes) {
        return undef;
    } else {
        return $nextRes->connectedToPrev;
    }
}


sub connectedToPrev {
    #is this residue covalently connected to the previous residue in the chain?
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   true if this residue is connected to the previous one, false if it isn't
    
    my $self = shift;
    my $prevRes = $self->prevRes;
    
    if ($self->{connectedToPrev} eq NOTCALC) {
        my $prevRes = $self->prevRes;
        if ((defined $prevRes) and (defined $self->atoms->{"P"}) and (defined $prevRes->atoms->{"O3'"}) and not (defined $self->structure and $self->structure->{pseudoatom})) {
            #check the P-O3' bond distance
            my $bondDist = calcDist($self->atoms->{"P"}, $prevRes->atoms->{"O3'"});
            if ($bondDist < 3) {
                $self->{connectedToPrev} = 1;
            } else {
                $self->{connectedToPrev} = 0;
            }
        } elsif ((defined $prevRes) and (defined $self->atoms->{"P"}) and (defined $prevRes->atoms->{"C1'"})) {
            #check the P-C1' pseudobond distance
            my $bondDist = calcDist($self->atoms->{"P"}, $prevRes->atoms->{"C1'"});
            if ($bondDist < 6) {
                $self->{connectedToPrev} = 1;
            } else {
                $self->{connectedToPrev} = 0;
            }
        } else {
            $self->{connectedToPrev} = 0;
        }
    }
    
    return $self->{connectedToPrev};
}

sub nextConnectedRes {
    #get the next residue, but only if it is connected to this residue
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a Residue object containing the next connected residue
    #   undef if there is no connected next residue
    
    my $self = shift;
    if ($self->connectedToNext) {
        return $self->{chain}->{resList}->[$self->{seqNum}+1];
    } else {
        return undef;
    }
}

sub prevConnectedRes {
    #get the previous residue, but only if it is connected to this residue
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a Residue object containing the previous connected residue
    #   undef if there is no connected previous residue
    
    my $self = shift;
    if ($self->connectedToPrev) {
        return $self->{chain}->{resList}->[$self->{seqNum}-1];
    } else {
        return undef;
    }
}


sub startingSuite {
    #get the suite that contains the start of this residue
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a Suite object containing the starting suite
    
    my $self = shift;
    if ($self->connectedToPrev) {
        return $self->{chain}->{suiteList}->[$self->{seqNum}];
    } else {
        return undef;
    }
}

sub endingSuite {
    #get the suite that contains the end of this residue
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a Suite object containing the ending suite
    
    my $self = shift;
    if (self->connectedToNext) {
        return $self->{chain}->{suiteList}->[$self->{seqNum} + 1];
    } else {
        return undef;
    }
}

sub startingSugarDist {
    #get the distance from the $sugarAtom (C1' or C4') of this nucleotide to the $sugarAtom of the previous nucleotide
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the C1'-C1' (or C4'-C4') distance, or undef if there is no previous connected residue
    
    my $self = shift;
    
    if ($self->{startingSugarDistVal} eq NOTCALC) {
        if ($self->connectedToPrev and
            defined $self->atoms->{$sugarAtom} and
            defined $self->prevRes->atoms->{$sugarAtom}) {
            
            $self->{startingSugarDistVal} = calcDist($self->atoms->{$sugarAtom},
                                             $self->prevRes->atoms->{$sugarAtom});
            
        } else {
            $self->{startingSugarDistVal} = undef;
        }                                
    }
    
    return $self->{startingSugarDistVal};
}

sub endingSugarDist {
    #get the distance from the $sugarAtom (C1' or C4') of this nucleotide to the $sugarAtom of the next nucleotide
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   the C1'-C1' (or C4'-C4') distance, or undef if there is no next connected residue
    
    my $self = shift;
    
    if ($self->connectedToNext) {
        return $self->nextRes->startingSugarDist;
    } else {
        return undef;
    }
}

#functions to access data about the residue
sub structure   {   $ARG[0]->{structure};   }
sub number      {   $ARG[0]->{number};      }
sub type        {   $ARG[0]->{type};        }
sub name        {   $ARG[0]->{type};        }
sub atoms       {   $ARG[0]->{atoms};       }
sub resNum      {   $ARG[0]->{number};      }
sub chain       {   $ARG[0]->{chain};       }

1;