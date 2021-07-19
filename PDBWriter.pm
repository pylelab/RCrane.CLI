#!/usr/bin/perl

# A class for outputting a PDB file
# EXAMPLE:
#   my $pdb = new PDB("test.pdb");
#   my $res = $pdb->firstResidue;
#   my $pdbWriter = new PDBWriter("new.pdb");
#   $pdbWriter->printRes($res);
#   $pdbWriter->close;

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

package PDBWriter;

use strict;
use English '-no_match_vars';
use Data::Dumper;

use constant SUGARATOMS     => ("C4'", "O4'", "C1'", "C2'", "O2'", "C3'");
use constant NONSUGARATOMS  => ("P", "OP1", "OP2", "O5'", "C5'");
use constant BASEATOMS      => qw(N1 N2 O2 C2 N3 N4 O4 C4 C5 N6 O6 C6 N7 C8 N9);

sub new {
    #create and initialize a PDBWriter object
    #ARGUMENTS:
    #   $fileName - the file name to write the PDB to
    #OPTIONAL ARGUMENTS:
    #   $startingSuiteNum   - the number of the first suite (all suites will be numbered consecutively if no Suite object is given to printSuite)
    #                         if not given, 2 is used (making the first residue 1)
    #                         Note that this variable is ignored if Suite objects are given to the printSuite method
    #RETURNS:
    #   a PDBWriter object

    my $class            = shift;
    my $fileName         = shift;
    my $startingSuiteNum = shift;
    
    #decrement $startingSuiteNum, since it will be incremented when printSuite is called
    if (defined $startingSuiteNum) {
        $startingSuiteNum--;
    } else {
        $startingSuiteNum = 1;
    }
    
    my $OUT;
    open($OUT, ">", $fileName) or die "Could not open $fileName\n";
    
    my $self = { fileName   => $fileName,
                 OUT        => $OUT,
                 suiteNum   => $startingSuiteNum,
                 atomNum    => 1      };
    bless $self, $class;
    
    return $self;
}

sub close {
    #close the PDB file
    #ARUGMENTS:
    #   None
    #RETURNS:
    #   None
    #SIDE EFFECTS:
    #   prints TER to the pdb file and then closes the output stream
    
    my $self = shift;
    
    print {$self->{OUT}} "TER\n";
    close($self->{OUT});
}

sub printBlank {
    #prints a blank line to the PDB file
    #ARUGMENTS:
    #   None
    #RETURNS:
    #   None
    #SIDE EFFECTS:
    #   prints a blank line to the PDB file
    
    my $self = shift;
    
    print {$self->{OUT}} "\n";
}

sub printSuite {
    #prints the backbone of a given suite to the PDB file
    #ARGUMENTS:
    #   $backbone           - the backbone atoms to print in the format of {atom name => [x,y,x]}
    #OPTIONAL ARGUMENTS:
    #   $base               - a Suite object containing the base atoms to print (for both bases)
    #                         this object will also be used to name and number the residues
    #                         if not given, no base atoms will be printed, residues will be named UNK,
    #                         and the residues will be numbered based on the previous suite printed
    #   $printStartingSugar - whether or not to print the starting sugar of the suite
    #                         (since the ending sugar of the previous suite may have been printed)
    #                         defaults to true
    #RETURNS:
    #   None
    #SIDE EFFECTS:
    #   prints the given suite to the PDB file
    
    my $self                = shift;
    my $backbone            = shift;
    my $bases               = shift;
    my $printStartingSugar  = shift;
    $printStartingSugar = 1 unless defined($printStartingSugar);
    
    my ($startingResName, $endingResName, $suiteNum);
    if (defined $bases) {
        #if a Suite object is given, use the residue names and numbers from that object
        $startingResName    = $bases->startingRes->name;
        $endingResName      = $bases->endingRes->name;
        
        $suiteNum = $bases->number;
        $self->{suiteNum} = $suiteNum
    } else {
        #if no Suite object is given, name the residues UNK and base the residue numbers on what was stored in $self->{suiteNum}
        $startingResName    = "UNK";
        $endingResName      = "UNK";
        
        $self->{suiteNum}++;
        $suiteNum           = $self->{suiteNum}
    }
    
    #if desired, print out the first sugar
    if ($printStartingSugar) {
        for my $curAtom (SUGARATOMS) {
            $self->printAtom($curAtom, $startingResName, $suiteNum-1, $backbone->{"$curAtom-1"});
        }
    }
    
    #print out the first base
    if (defined $bases) {
        for my $curAtom (BASEATOMS) {
            if (exists $bases->atoms->{"$curAtom-1"}) {
                $self->printAtom($curAtom, $startingResName, $suiteNum-1, $bases->atoms->{"$curAtom-1"});
            }
        }
    }
    
    #print out the backbone atoms in between the sugars
    $self->printAtom("O3'", $startingResName, $suiteNum-1, $backbone->{"O3'-1"});
    for my $curAtom (NONSUGARATOMS) {
        $self->printAtom($curAtom, $endingResName, $suiteNum, $backbone->{"$curAtom"});
    }
    
    #print out the second sugar
    for my $curAtom (SUGARATOMS) {
        $self->printAtom($curAtom, $endingResName, $suiteNum, $backbone->{"$curAtom"});
    }
    
    #print out the second base
    if (defined $bases) {
        for my $curAtom (BASEATOMS) {
            if (exists $bases->atoms->{"$curAtom"}) {
                $self->printAtom($curAtom, $endingResName, $suiteNum, $bases->atoms->{"$curAtom"});
            }
        }
    }
}


sub printRes {
    #prints a given residue to the PDB file
    #ARGUMENTS:
    #   this function can be given either a Residue object or a hash of atomic coordinates
    #   Residue-object mode:
    #       $res    - a Residue object to print
    #   Hash mode:
    #       $atoms      - a hash of atomic coordinates
    #       $resName    - (optional) the name of the residue (i.e. A, C, G, or U)
    #                     if not given, defaults to UNK
    #       $resNum     - the residue number
    #                     if not given, uses and increments $self->{suiteNum}
    #       $chain      - the chain ID
    #                     if not given, uses a blank space
    
    my $self    = shift;
    
    my ($atoms, $resName, $resNum, $chain);
    if (ref $ARG[0] eq "Residue") {
        my $res  = shift;
        $atoms   = $res->atoms;
        $resName = $res->name;
        $resNum  = $res->number;
        $chain   = $res->chain->id;
    } else {
        $atoms   = shift;
        $resName = (shift or "UNK");
        $resNum  = shift;
        $chain   = (shift or " ");
    }
    
    unless (defined $resNum) {
         $resNum = $self->{suiteNum};
         $self->{suiteNum}++;
    }
    
    for my $curAtom (NONSUGARATOMS, SUGARATOMS, BASEATOMS, "O3'") {
        if (defined $atoms->{$curAtom}) {
            $self->printAtom($curAtom, $resName, $resNum, $atoms->{$curAtom}, $chain);
        }
    }
}


sub printAtom {
    #prints an atom to the PDB file
    #ARGUMENTS:
    #   $atomName   - the name of the atom
    #   $resName    - the name of the residue
    #   $resNum     - the residue number of the atom
    #   $coords     - a ref to array of the [x,y,z] coordinates of the atom
    #RETURNS:
    #   None
    #SIDE EFFECTS:
    #   prints atom to the PDB file
    
    my $self     = shift;
    my $atomName = shift;
    my $resName  = shift;
    my $resNum   = shift;
    my $coords   = shift;
    my $chain    = shift;
    
    my ($resNumNum, $insCode);
    if ($resNum =~ /^(\d+)([A-Z]?)$/i) {
        $resNumNum = $1;
        $insCode   = $2;
    }
    
    print {$self->{OUT}} sprintf("ATOM  %5i %-4s %-3s %1s%4i%1s   %8.3f%8.3f%8.3f  1.00 10.00           %s  \n", $self->{atomNum}, " " . $atomName, $resName, $chain, $resNumNum, $insCode, @{$coords}, substr($atomName,0,1));
    $self->{atomNum}++;
}

1;