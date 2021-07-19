#!/usr/bin/perl
#
# A class to read coordinate data from a PDB file and access it as either residues (nucleotides) or suites
# EXAMPLE:
#   my $pdb = new PDB("test.pdb");
#   my $curSuite = $pdb->firstSuite;
#   while ($curSuite) {
#       print "The eta value of suite " . $curSuite->fullNumber . " is " . $curSuite->theta . "\n";
#       $curSuite = $curSuite->nextSuite;
#   }
#
# RESTRICTIONS:
#   All alternate locations are ignored

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

package PDB;

use strict;
use English '-no_match_vars';
use Chain;
use Residue;
use Suite;
use Strip;


my $pseudoatomHash = {map {$ARG => 1} qw(N1 N2 O2 C2 N3 N4 O4 C4 C5 N6 O6 C6 N7 C8 N9 C1' P)};
    #a hash of atoms to read in pseudoatom-only mode


sub new {
    #create and initialize a PDB object
    #ARGUMENTS:
    #   $filename       - what filename to read
    #OPTIONAL ARGUMENTS:
    #   $pseudoatom     - if true, will read in only base, C1', and P atoms
    #                     defaults to false
    #   $manualConnect  - if true, consecutive residues are assumed to be connected unless there is a BREAK record between them
    #                     defaults to $pseudoatom
    #RETURNS:
    #   a PDB object initialized with data from $filename, if provided

    my $class           = shift;
    my $filename        = shift;
    my $pseudoatom      = shift;
    my $manualConnect   = shift;
    
    my $self = { filename       => undef,
                 pseudoatom     => undef,
                 manualConnect  => undef,
                 chainHash      => {},
                 chainList      => []      };
    bless $self, $class;
    
    $self->read($filename, $pseudoatom, $manualConnect);
    
    return $self;
}

sub read {
    #read the contents of a PDB file
    #ARGUMENTS:
    #   $filename   - what filename to read
    #OPTIONAL ARGUMENTS:
    #   $pseudoatom     - if true, will read in only base, C1', and P atoms
    #                     defaults to false
    #   $manualConnect  - if true, consecutive residues are assumed to be connected unless there is a BREAK record between them
    #                     defaults to $pseudoatom
    #RETURNS:
    #   None
    #EFFECTS:
    #   $self will now contain all coordinate data from the PDB file

    my $self            = shift;
    my $filename        = shift;
    my $pseudoatom      = shift;
    my $manualConnect   = shift;
    
    if (not defined $manualConnect) {
        $manualConnect = $pseudoatom;
    }
    
    $self->{filename}   = $filename;
    $self->{pseudoatom} = $pseudoatom;
    $self->{manualConnect} = $manualConnect;
    
    open(IN, $filename) or die "Could not open $filename for reading\n";
    
    #read in each line from the PDB file and create the necessary residue objects
    my $prevResObj; #the previous residue object that we read in
        #we need to remember this in case we find a BREAK record
    while (my $curline = <IN>) {
        my $recordType = strip(substr($curline, 0, 6));
        
        if ($recordType eq "BREAK") {
            $prevResObj->{connectedToPrev} = 0;
            next;
        }
        
        next unless ($recordType eq "ATOM" or $recordType eq "HETATM");
        
        my $atomName    = strip(substr($curline, 12, 4));
        my $resName     = strip(substr($curline, 17, 3));
        my $chain       = substr($curline, 21, 1);
        my $resNum      = substr($curline, 22, 4) + 0;
        my $insCode     = strip(substr($curline, 26, 1));
        my $x           = substr($curline, 30, 8);
        my $y           = substr($curline, 38, 8);
        my $z           = substr($curline, 46, 8);
        
        $resNum .= $insCode;
        
        #convert the atom names to those in PDB file format 3.0
        $atomName =~ s/\*/'/g; #apostrophes instead of asterisk
                               #since we should only be getting RNA in this PDB file, don't worry about accidentally converting extra asterisks
        $atomName = "OP1" if $atomName eq "O1P";
        $atomName = "OP2" if $atomName eq "O2P";
        
        #if we're in pseudoatom only mode and this isn't a pseudoatom, skip it
        next if ($pseudoatom and not $pseudoatomHash->{$atomName});
        
        #if this is the first residue of a chain, create a new Chain object
        if (not defined $self->{chainHash}->{$chain}) {
            $self->{chainHash}->{$chain} = scalar @{$self->{chainList}};
            push(@{$self->{chainList}}, new Chain($self, $chain));
        }
        
        my $curChainObj = $self->{chainList}->[$self->{chainHash}->{$chain}];
        
        #if this is the first atom of this residue, create a new Residue object
        unless ($curChainObj->hasRes($resNum)) {
            $curChainObj->addRes($resNum, $resName);
            
            #if we are assigning all connectivity manually, assume that this residue is connected to the previous one
            #unless it's the first residue
            #(if we find a BREAK card, then we'll change it so that it's not connected
            if ($manualConnect) {
                if ($curChainObj->numRes == 1) {
                    $curChainObj->res($resNum)->{connectedToPrev} = 0;
                } else {
                    $curChainObj->res($resNum)->{connectedToPrev} = 1;
                }
            }
        }
        
        #warn if this atom already exists
        if (defined $curChainObj->res($resNum)->atoms->{$atomName}) {
            warn "Duplicate atom found: file $filename, chain $chain, residue $resNum ($resName), atom $atomName\n";
        } else {
            #otherwise, store the atomic coordinates
            $curChainObj->res($resNum)->atoms->{$atomName} = [$x, $y, $z];
        }
        
        $prevResObj = $curChainObj->res($resNum)
    }
    close(IN);
    
}


sub firstChain {
    #get the first chain of the file
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   a Chain object containing the first chain
    
    my $self = shift;
    
    return $self->{chainList}->[0];
}

sub lastChain {
    #get the last chain of the file
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   a Chain object containing the last chain
    
    my $self = shift;
    
    return $self->{chainList}->[-1];
}

sub chain {
    #retreive the specified chain
    #ARGUMENTS:
    #   $chain  - the ID of the requested chain
    #RETURNS:
    #   a Chain object containing the specified chain, or undef if the chain does not exist
    
    my $self = shift;
    my $chain = shift;
    
    return $self->{chainList}->[$self->{chainHash}->{$chain}];
}

sub firstRes {
    #get the first residue of the first chain
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   a Residue object containing the first residue
    
    my $self = shift;
    
    return $self->firstChain->firstRes;
}

sub firstSuite {
    #get the first suite of the first chain
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   a Suite object containing the first suite
    
    my $self = shift;
    
    return $self->firstChain->firstSuite;
}

sub lastRes {
    #get the last residue of the last chain
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   a Residue object containing the last residue
    
    my $self = shift;
    
    return $self->lastChain->lastRes;
}

sub lastSuite {
    #get the last suite of the last chain
    #ARGUMENTS:
    #   none
    #RETURNS:
    #   a Suite object containing the last suite
    
    my $self = shift;
    
    return $self->lastChain->lastSuite;
}

1;