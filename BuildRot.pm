#!/usr/bin/perl

#module to built initial coordinates for nucleotides
#these initial coordinates will not necessarily be in chemically reasonable positions
#use FullSuiteMin.pm for that

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

package BuildRot;

use strict;
use English '-no_match_vars';
use Data::Dumper;

use Rotate qw(rotateAtoms);
use MeasureTorsions qw(calcZeta);
use PuckerList;
use BuildPhosOxy;
use BuildInitSugar;

use StrucCalc qw(calcDistToLine minus magnitude plus scalarProd crossProd);

use constant { PO5pIDEAL            => 1.593,   #ideal bond length for the P-O5' bond
               PANGLEIDEAL          => 104.0,   #ideal bond angle for the angle about the P atom
};

sub new {
    #create and initialize a BuildRot object
    #ARGUMENTS:
    #   $dihedData - filename containing data on the dihedral angles for each rotamer
    #   $c3pStruc  - filename for a PDB file containing a C3' endo pucker sugar
    #   $c2pStruc  - filename for a PDB file containing a C2' endo pucker sugar
    #RETURNS:
    #   $self      - a blessed BuildRot object

    my $class       = shift;
    my $dihedData   = shift;
    my $c3pStruc    = shift;
    my $c2pStruc    = shift;
    
    my $self = { dihedDataFilename      => $dihedData,
                 c3pStrucFilename       => $c3pStruc,
                 c2pStrucFilename       => $c2pStruc,
                 sugarBuilder           => undef,
                 torsionMinimizer       => undef,
                 ABGdata                => {},
                 DEZData                => {},
               };
    bless $self, $class;
    
    $self->{sugarBuilder}     = new BuildInitSugar($c3pStruc, $c2pStruc);
    $self->readDihedData;
    
    return $self;
}

sub buildNt {
    #build initial atom locations for a nucleotide
    #actually, we're not quite building a full nucleotide - we normally build only the O5' from the next nucleotide, not the current one
    #(unless the current one doesn't have an O5', in which case we also build it too)
    #ARGUMENTS:
    #   $rot     - the rotamer of the *next* suite
    #   $curNuc  - coordinates for atoms in the current nucleotide (at least base atoms must be provided)
    #   $nextNuc - coordinates for atoms in the next nucleotide (at least the phosphate atom must be provided)
    #              (if you do not have coordinates for the next phosphate, use buildLastNt instead)
    #RETURNS:
    #   $curNuc  - coordinates for atoms in the current nucleotide with backbone atoms built
    #   $nextNuc - coordinates for atoms in the next nucleotide with the O5' atom built
    
    
    my $self = shift;
    my $rot  = shift;
    my $curNuc   = shift;
    my $nextNuc  = shift;
    
    #it would probably be faster (run-time-wise) to add on the O5'+1 later, since we can directly calculate an
    #optimal location for it given the C3', O3' and P locations (as we're not placing the C5'+1 yet and
    #therefore can trivially satisfy ideal bond length, angle, and torsion)
    #but for now, we add the O5'+1 here since it will (hopefully) make it easier to modify the minimization
    #later if we decide to add on more constraints
    
    
    #figure out the pucker of this nucleotide
    my $pucker = $puckerList->{$rot}->[0];
    
    #build the sugar and backbone of the current nucleotide
    my $builtNuc = $self->sugarBuilder->buildInitSugar($curNuc, $pucker);
    if (defined ($curNuc->{"O5'"})) {
        #if we already had an O5', then keep it, otherwise use the default position
        $builtNuc->{"O5'"} = $curNuc->{"O5'"};
    }
    $curNuc = $builtNuc;
    
    
    #place the O5' atom of the next nucleotide
    #we use the ideal zeta torsion for the current rotamer, along with the ideal bond length and andgle
    #this won't be particularly accurate given that C5' and O5' haven't been appropriately positioned yet
    #but the O5' location will be appropriately adjusted during the minimization
    
    #calculate a point along the P-O3' bond that is PO5pIDEAL Angstroms away from the phosphate
    $curNuc->{"P+1"} = $nextNuc->{"P"};
    my $pO3pVector = minus($curNuc->{"O3'"}, $curNuc->{"P+1"});
    my $pO3pDist   = magnitude($pO3pVector);
    my $pO5pVector = scalarProd($pO3pVector, PO5pIDEAL / $pO3pDist);
    $curNuc->{"O5'+1"} = plus($curNuc->{"P+1"}, $pO5pVector);
    
    #rotate O5' about the P
    #first, calculate an arbitrary axis to do the rotation about
    #the axis just needs to be through the phosphate and perpendicular to the P-O3' bond
    my $axisVector = crossProd( minus($curNuc->{"P+1"}, $curNuc->{"O3'"}), minus($curNuc->{"P+1"}, $curNuc->{"C3'"}));
    $curNuc->{"axis"} = plus($curNuc->{"P+1"}, $axisVector);
    #$curNuc->{"O5'+1"} = $curNuc->{"axis"};
    my $rotatedAtoms = rotateAtoms($curNuc, ["O5'+1"], ["P+1", "axis"], PANGLEIDEAL);
    $curNuc->{"O5'+1"} = $rotatedAtoms->{"O5'+1"};
    
    #rotate O5' about zeta
    #first, fetch the ideal zeta value for the current rotamer
    my $idealZeta = $self->{DEZdata}->{$rot}->[2];
    $rotatedAtoms = rotateAtoms($curNuc, ["O5'+1"], ["O3'", "P+1"],  calcZeta($curNuc) - $idealZeta);
    $curNuc->{"O5'+1"} = $rotatedAtoms->{"O5'+1"};
    
    #copy the appropriate atoms from $curNuc to $nextNuc
    $nextNuc->{"P"}   = $curNuc->{"P+1"};
    $nextNuc->{"O5'"} = $curNuc->{"O5'+1"};
    
    #delete the temporary atoms from $curNuc
    delete $curNuc->{"axis"};
    delete $curNuc->{"P+1"};
    delete $curNuc->{"O5'+1"};
    
    return ($curNuc, $nextNuc);
}


sub buildLastNt {
    #build initial atom locations for the last nucleotide of a segment
    #this differs from buildNt because we don't build the next O5' (and we take the current rotamer as an argument instead of the last one)
    #ARGUMENTS:
    #   $rot     - the rotamer of the *current* suite
    #   $curNuc  - coordinates for atoms in the current nucleotide (at least base atoms must be provided)
    #RETURNS:
    #   $curNuc  - coordinates for atoms in the current nucleotide with backbone atoms built
    
    my $self     = shift;
    my $rot      = shift;
    my $curNuc   = shift;
    
    #figure out the pucker of this nucleotide
    my $pucker = $puckerList->{$rot}->[1];
    
    #build the sugar and backbone of the current nucleotide
    my $builtNuc = $self->sugarBuilder->buildInitSugar($curNuc, $pucker);
    if (defined ($curNuc->{"O5'"})) {
        #if we already had an O5', then keep it, otherwise use the default position
        $builtNuc->{"O5'"} = $curNuc->{"O5'"};
    }
    
    return $builtNuc;

}


sub readDihedData {
    #read in the data on torsion means and standard deviations
    #ARGUMENTS:
    #   none (filename to read is taken from $self->{dihedDataFilename})
    #RETURNS:
    #   none
    #EFFECTS:
    #   reads data into $self->{DEZdata} and $self->{ABGdata}
    #NOTE:
    #   this module only needs information on the zeta torsion
    
    
    my $self = shift;
    
    open(IN, $self->{dihedDataFilename}) or die "Could not open " . $self->{dihedDataFilename} . " for reading\n";
    <IN>; #skip the header line
    
    while (my $curline = <IN>) {
        my @curdata = split(",", $curline);
        
        my ($rot, undef, undef, $deltaMean, $epMean, $zetaMean, $alphaMean, $betaMean, $gammaMean, undef, undef, undef, undef, $deltaSD, $epSD, $zetaSD, $alphaSD, $betaSD, $gammaSD) = @curdata;
        
        $self->{DEZdata}->{$rot} = [$deltaMean, $epMean,   $zetaMean,  $deltaSD, $epSD,   $zetaSD];
        $self->{ABGdata}->{$rot} = [$alphaMean, $betaMean, $gammaMean, $alphaSD, $betaSD, $gammaSD];
    }
    
    close(IN);
}

sub sugarBuilder       {    return $ARG[0]->{sugarBuilder};      }

1;