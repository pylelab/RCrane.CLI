#!/usr/bin/perl

#given a structure with phosphates and base atoms, calculate and build a rotameric backbone
#for more information, type rcrane.cli.pl -h
#
#this program requires the following CPAN modules:
#   Math::Interpolate
#   Math::Amoeba

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

use strict;
use English '-no_match_vars';
use Getopt::Long qw(:config bundling ignore_case_always);
use File::Basename;
use Pod::Usage;
use FindBin;                     #Find the directory where this script is located so that modules can be imported
use lib "$FindBin::RealBin";
use Pod::Usage;                  #used to print help
use Data::Dumper;

use PDB;
use PseudoPredictor;
use RotamerHMM;
use BuildRot;
use PDBWriter;
use FullSuiteMinimizer;
use BuildPhosOxy;
use PuckerList;

my $delim = ","; #the field deliter to use for the output CSV files

#parse the input arguments
my ($input, $outputPrefix, $excelOutput, $numStrucs, $singleStruc, $rotString, $noStruc, $noTables, $verbose) = parseArgs();


#####################################################################
#  Open input and output files and initialize all predictors
#####################################################################

#read in the input file
my $pdb = new PDB($input, 1, 0); #read in only the base, C1', and P atoms of the input

#initialize the predictor
my $pseudoPredic = new PseudoPredictor( THETAETA      => "$FindBin::RealBin/data/thetaEtaClusts.csv",
                                        PUCKER        => "$FindBin::RealBin/data/smoothedPuckerDist.csv",
                                        SUGARDIST     => "$FindBin::RealBin/data/sugarDists.csv",
                                        STARTPHOSDIST => "$FindBin::RealBin/data/startingPDists.csv",
                                        ENDPHOSDIST   => "$FindBin::RealBin/data/endingPDists.csv",
                                        THETAONLY     => "$FindBin::RealBin/data/thetaOnly.csv",
                                        ETAONLY       => "$FindBin::RealBin/data/etaOnly.csv", );

my $rotList = [sort $pseudoPredic->getRotList()];

my $bestPath;
my $predictedProbs = [];

if ($rotString) {
    #if the user gave us a rotamer string to use, make sure that it's valid
    
    #strip out all spaces from $rotString
    $rotString =~ tr/ //d;
    
    #convert $rotString into a list of rotamers
    $bestPath = [unpack("(a2)*", $rotString)];
    
    #make sure the rotamer string only contains valid rotamers
    my $rotHash = {map {$ARG => 1} @{$rotList}};                    #create a hash of {rotamer} => True
    my $invalidRots = [grep((not $rotHash->{$ARG}), @{$bestPath})]; #pull out all invalid rotamers from $bestPath
    if (@{$invalidRots}) {
        die "Error: Specified rotamer string contains invalid rotamers: " . join(", ", @{$invalidRots}) . "\n";
    }
    
    #create a list of suites in the PDB file
    my $curSuite = $pdb->firstSuite; #read in the first suite
    while ($curSuite) {
        push (@{$predictedProbs}, [$curSuite, undef]);
        $curSuite = $curSuite->nextSuite;
    }
    
    #check how long the PDB file is (in suites) so we can make sure our rotamer string is the appropriate length
    unless (@{$bestPath} == @{$predictedProbs}) {
        die "Error: Specified rotamer string is not appropriate length.  The input PDB contains\n    " . scalar(@{$predictedProbs}) . " suites and the rotamer string contains " . scalar(@{$bestPath}) . " rotamers.\n";
    }
    
    #make sure that the rotamers don't have any conflicting sugar puckers
    for (my $i = 1; $i < @{$bestPath}; $i++) {
        if ($predictedProbs->[$i]->[0]->connectedToPrev and ($puckerList->{$bestPath->[$i-1]}->[1] != $puckerList->{$bestPath->[$i]}->[0])) {
            die "Error: Specified rotamer string contains conflicting puckers.  Suite " . $predictedProbs->[$i-1]->[0]->fullNumber . " (" . $bestPath->[$i-1] . ") ends with a " . $puckerList->{$bestPath->[$i-1]}->[1] . "' endo\n    pucker and suite " . $predictedProbs->[$i]->[0]->fullNumber . " (" . $bestPath->[$i] . ") starts with a " . $puckerList->{$bestPath->[$i]}->[0] . "' endo pucker.\n";
        }
    }
    
    
    print "Using input rotamer string: " . join("", @{$bestPath}) . "\n";
    
} else {
    #otherwise, predict the most likely rotamer sequence
    
    #####################################################################
    #  Go through each suite of the input file and generate predictions
    #####################################################################
    
    my $suiteStrings = [];           #a list of the suite strings for each predicted structure
    my $curSuite = $pdb->firstSuite; #read in the first suite
    #print Dumper($curSuite);
    
    while ($curSuite) {
        #print "In loop - suite " . $curSuite->fullNumber . "\n";
        
        #################################################################
        #  Make sure we can calculate all of the information needed for the prediction
        #################################################################
        
        #if we can't figure out the pseudotorsions or the pucker, then warn that prediction accuracy may be decreased
        if (not defined $curSuite->theta) {
            warn "Warning: Theta could not be calculated for suite "            . $curSuite->fullNumber . ".  Prediction accuracy may be reduced for this suite.\n";
        }
        if (not defined $curSuite->eta) {
            warn "Warning: Eta could not be calculated for suite "              . $curSuite->fullNumber . ".  Prediction accuracy may be reduced for this suite.\n";
        }
        if (not defined $curSuite->startingRes->pperp) {
            warn "Warning: Starting P-perp could not be calculated for suite "  . $curSuite->fullNumber . ".  Prediction accuracy may be reduced for this suite.\n";
        }
        if (not defined $curSuite->endingRes->pperp) {
            warn "Warning: Ending P-perp could not be calculated for suite "    . $curSuite->fullNumber . ".  Prediction accuracy may be reduced for this suite.\n";
        }
        if (not defined $curSuite->sugarDist) {
            warn "Warning: Sugar (C1'-C1') distance could not be calculated for suite "       . $curSuite->fullNumber . "\n"; }
        if (not defined $curSuite->startingRes->phosDist) {
            warn "Warning: Starting phosphate distance could not be calculated for suite "    . $curSuite->fullNumber . "\n"; }
        if (not defined $curSuite->endingRes->phosDist) {
            warn "Warning: Ending phosphate distance could not be calculated for suite "      . $curSuite->fullNumber . "\n"; }
        
        
        #################################################################
        #  Make the prediction and store the results
        #################################################################
        
        
        #calculate the probabilities
        my $curProbs = $pseudoPredic->calcProb( THETA          => $curSuite->theta,
                                                ETA            => $curSuite->eta,
                                                STARTPPERP     => $curSuite->startingRes->pperp,
                                                ENDPPERP       => $curSuite->endingRes->pperp,
                                                SUGARDIST      => $curSuite->sugarDist,
                                                STARTPHOSDIST  => $curSuite->startingRes->phosDist,
                                                ENDPHOSDIST    => $curSuite->endingRes->phosDist);
        
        
        push (@{$predictedProbs}, [$curSuite, $curProbs]);
    } continue {
        $curSuite = $curSuite->nextSuite;   #go to the next suite
    }
    
    #print Dumper($predictedProbs);
    
    #################################################################
    #  Calculate the best sequence of suites using an HMM
    #  (this ensures that consecutive suites have compatable sugar puckers)
    #################################################################
    $bestPath = rotamerHMM($predictedProbs, $pseudoPredic);
    print "Predicted rotamer string:\n\t" . join("", @{$bestPath}) . "\n";
    
    unless ($noTables) {
        printProbsTable("$outputPrefix.probs.csv", $predictedProbs, $rotList);
        printTestData  ("$outputPrefix.data.csv",  $predictedProbs);
        printBestRots  ("$outputPrefix.rots.csv",  $predictedProbs, $pseudoPredic, $numStrucs, $bestPath);
    }
}

unless ($noStruc) {
    
    #################################################################
    #  Build the structure using the HMM-calculated path
    #################################################################
    
    #initialize BuildRot objects
    my $buildRot = new BuildRot("$FindBin::RealBin/data/dihedData.csv",
                                "$FindBin::RealBin/data/c3p.pdb",
                                "$FindBin::RealBin/data/c2p.pdb"         );
    
    #initialize the full suite minimizer
    initializeFullSuiteMinimizer("$FindBin::RealBin/data/dihedData.csv", $verbose);
    
    #this loop assumes that the only atoms that currently exist in the PDB object are base atoms, C1', and P
    #because the pseudo-atom only mode in the PDB module will only read in those atoms
    
    my ($builtPhosLoc, $nextBuiltPhosLoc);
    for my $loopVar (@{$predictedProbs}) {
        my $curSuite = $loopVar->[0];
        my $curRot   = shift(@{$bestPath});
        my $nextRot  = $bestPath->[0];
        my $builtInitSugar; #whether we have built the initial sugar
        
        #make sure that the initial sugar is present
        unless ($curSuite->startingRes->atoms->{"O4'"}) {
            #if it isn't, build it
            print "Building initial nt " . $curSuite->startingRes->number . "\n";
            ($curSuite->startingRes->{atoms}, $curSuite->endingRes->{atoms}) = $buildRot->buildNt($curRot, $curSuite->startingRes->atoms, $curSuite->endingRes->atoms);
            
            $builtInitSugar = 1; #remeber that we've built the init sugar
        }
        
        #remember where the next phosphate is before we do the minimization,
        #so that we can use it's original position as a constraint for the next minimization run
        if ($builtInitSugar) {
            #if we're right after a chain break, then the phosphate hasn't been moved yet
            undef $builtPhosLoc;
        } else {
            $builtPhosLoc = $nextBuiltPhosLoc;
        }
        if (defined $curSuite->endingRes->nextRes) {
            $nextBuiltPhosLoc = $curSuite->endingRes->nextRes->atoms->{"P"};
        }
        
        #print Dumper($builtPhosLoc);
        
        if ($curSuite->connectedToNext and defined($nextRot)) {
            #if we have a next suite, then we build the nucleotide normally
            
            print "Building nt " . $curSuite->endingRes->number . "\n";
            #build initial coordinates for the nucleotide
            ($curSuite->endingRes->{atoms}, $curSuite->endingRes->nextRes->{atoms}) = $buildRot->buildNt($nextRot, $curSuite->endingRes->atoms, $curSuite->endingRes->nextRes->atoms);
            
            #perform the minimization
            ($curSuite->startingRes->{atoms}, $curSuite->endingRes->{atoms}, $curSuite->endingRes->nextRes->{atoms}) = fullNtMinimize($curSuite->startingRes->atoms, $curSuite->endingRes->atoms, $curSuite->endingRes->nextRes->atoms, $curRot, $bestPath->[0], $builtPhosLoc, $verbose);
            
        } else {
            #if this is the last nucleotide of a chain, then we have to modify the minimization slightly since
            #we don't have a next rotamer nor a next O5' (and possibly no next phosphate, although then we wouldn't
            #have made a rotamer prediction on the current suite)
            
            print "Building final nt " . $curSuite->endingRes->number . "\n";
            
            #build initial coordinates for the nucleotide
            $curSuite->endingRes->{atoms} = $buildRot->buildLastNt($curRot, $curSuite->endingRes->atoms);
            
            #if there is a final phosphate present, then minimize it's location
            my $nextResAtoms;
            if ($curSuite->endingRes->connectedToNext) {
                $nextResAtoms = $curSuite->endingRes->nextRes->{atoms};
            }
            
            ($curSuite->startingRes->{atoms}, $curSuite->endingRes->{atoms}, $nextResAtoms) = fullNtMinimize($curSuite->startingRes->atoms, $curSuite->endingRes->atoms, $nextResAtoms, $curRot, undef, $builtPhosLoc, $verbose);
            
            if ($curSuite->endingRes->connectedToNext) {
                #if there is a final phosphate, build phosphate oxygens on it
                $nextResAtoms = buildInitPhosOxy($nextResAtoms, $curSuite->endingRes->atoms);
                
                $curSuite->endingRes->nextRes->{atoms} = $nextResAtoms;
            }
        }
        
        if ($builtInitSugar and $curSuite->startingRes->atoms->{"P"}) {
            #if we're at the start of the chain and there is an initial phosphate,
            #minimize the initial C5' and O5' coordinates and then build the initial phosphate oxygens
            
            print "Minimizing initial nt " . $curSuite->startingRes->number . "\n" if $verbose;
            
            #minimize the initial C5' and O5' coordinates
            $curSuite->startingRes->{atoms} = minInitAtoms($curSuite->startingRes->atoms, $verbose);
            
            #build the initial phoshate oxygens
            $curSuite->startingRes->{atoms} = buildInitPhosOxy($curSuite->startingRes->atoms);
            
        }
        
        
        #build the phosphoryl oxygens
        $curSuite->endingRes->{atoms} = buildPhosOxy($curSuite->endingRes->atoms, $curSuite->startingRes->atoms)
    }
    
    #print out the built structure
    my $pdbWriter = new PDBWriter("$outputPrefix.built.pdb");
    
    my $curRes = $pdb->firstRes;
    while ($curRes) {
        $pdbWriter->printRes($curRes);
        $curRes = $curRes->nextRes;
    }
    
    $pdbWriter->close;
    
    #report on how long the program took to run
    my $userTime = times;
    $userTime = int($userTime + 0.5); #round the time to the nearest second
    my $hours = int( $userTime / 3600);
    my $min   = int(($userTime % 3600) / 60);
    my $sec   = $userTime % 60;
    
    my $runtime;
    if ($hours) {
        $runtime = "$hours hour";
        $runtime .= "s" if $hours > 1;
        $runtime .= " ";
    }
    if ($hours or $min) {
        $runtime .= "$min min ";
    }
    $runtime .= "$sec sec";
    
    print "rcrane.cli.pl completed.  (Run time: $runtime)\n"; 
}


sub parseArgs {
    #parse the input arguments
    #ARGUMENTS:
    #   None (data is taken from @ARGV)
    #RETURNS:
    #   $input      - the name of the input file
    #   $output     - what to prefix all output filenames with
    #   $excel      - whether the output CSV files should be formatted for Excel
    #   $numStrucs  - the number of structures to generate
    #SIDE EFFECTS:
    #   empties @ARGV
    
    my ($input, $output, $excel, $numStrucs, $singleStruc, $rotString, $noStruc, $noTables, $verbose, $help, $man);
    
    my $progName = fileparse($0);
    
    #set the defaults
    $excel     = 0;
    $numStrucs = 3;
    $singleStruc = 0;
    
    my $success = GetOptions ('i|input=s'       => \$input,
                              'o|output=s'      => \$output,
                              'x|excel'         => \$excel,
                              'n|numStrucs=i'   => \$numStrucs,
                              '1|singleStruc'   => \$singleStruc,
                              'r|rot=s'         => \$rotString,
                              's|noStruc'       => \$noStruc,
                              't|noTables'      => \$noTables,
                              'v|verbose+'      => \$verbose,
                              'h|?|help'        => \$help,
                              'm|man'           => \$man);
    
    #print help if the user requested it
    if ($man) {
        pod2usage(-verbose => 2);
    } elsif ($help) {
        pod2usage();
    }
    
    #make sure GetOpt was able to parse the command line arguments
    unless ($success) {
        die "Type $progName -h for help.\n";
    }
    
    #if the user didn't specify the input file with the -i option, see if it they gave it on the command line without an option
    unless (defined $input) {
        $input = shift(@ARGV);
    }
    
    #if there is still something in @ARGV, then die since we don't know what to do with it
    if (@ARGV) {
        die "Unknown argument: " . shift(@ARGV) . "\nType $progName -h for help.\n";
    }
    
    #if the input file still isn't defined, then die
    unless (defined $input) {
        die "No input file given.  Type $progName -h for help.\n";
    }
    
    #make sure that the specified input file exists
    unless (-e $input) {
        if (-e "$input.pdb") {
            $input .= ".pdb";
        } elsif (-e "$input.PDB") {
            $input .= ".PDB";
        } elsif (-e "$input.ent") {
            $input .= ".ent";
        } elsif (-e "$input.ENT") {
            $input .= ".ENT";
        } else {
            die "Could not find input file $input.  Type $progName -h for help.\n";
        }
    }
    
    #see if the user gave an output prefix.  If not, use the prefix of the input filename
    unless ($output) {
        $output = fileparse($input, qr/\.(?:pdb|ent)$/i);
    }
    
    return ($input, $output, $excel, $numStrucs, $singleStruc, $rotString, $noStruc, $noTables, $verbose);
}
    

sub suiteNumOut {
    #if $excelOutput is set, format $suiteNum to be printed so that Excel will properly read the CSV file
    #ARGUMENTS:
    #   $suiteNum       - the suite number to print (formatted as "2-3")
    #VARIABLES TAKEN FROM THE CALLING SCOPE:
    #   $excelOutput    - if this is true, then the suite number will be formatted for Excel
    #                     otherwise, the suite number will be returned as is
    #RETURNS:
    #   $suiteNum formatted as desired
    #NOTE:
    #       When opening a CSV file, Microsoft Excel decides that everything that looks like a date must be a date.
    #   As a result, it will not only format the cell of data as date format, but will also replace the contents
    #   of the cell with Excel's internal representation of dates (an integer representing the number of days since
    #   the 0th of January, 1900).  As such, there is no way to recover the data once Excel has "helpfully" applied
    #   it's formatting to the cell.  As far as I know, there is no way to turn off this automatic formatting
    #   without having to go through the full Import External Data process on the CSV files or renaming them to
    #   a .txt extension.  The same thing occurs in OpenOffice Calc, although that at least allows you to specify
    #   that the column is text before it opens the CSV file.
    #       Because the full suite numbers look like dates to Excel (for example, "2-3", which represents suite 3,
    #   gets formatted as the 3rd of February of the current year), we must protect these suite numbers so Excel
    #   doesn't reformat them.  To do this, we turn the suite number into a formula: ="2-3".  However, other, less
    #   "helpful" programs probably aren't stupid enough to destroy data just because it looks like a date, so the
    #   default of this program is to print out just the suite number, and only protect the data from Excel if
    #   the -x command line option is given.
    
    my $suiteNum = shift;
    
    if ($excelOutput) {
        return '="' . $suiteNum . '"';
    } else {
        return $suiteNum;
    }
}

sub printProbsTable {
    #print a table of likelihoods for each rotamer at each suite
    #ARGUMENTS:
    #   $output         - the filename to output to
    #   $predictedProbs - a list of [$suite, $probs] for each suite, where $suite is a Suite object
    #                     and $probs is a hash of probabilities for each rotamer
    #   $rotList        - a list of all rotamers
    #RETURNS:
    #   none
    #EFFECTS:
    #   prints out a table of likelihoods to $output
    #   will die if $output cannot be written to
    
    my $output          = shift;
    my $predictedProbs  = shift;
    my $rotList         = shift;
    
    
    open(OUT, ">", $output) or die "Could not open $output for writing\n";
    
    #print the header
    print OUT join($delim, "Suite", @{$rotList}) . "\n";

    foreach my $curData (@{$predictedProbs}) {
        my ($curSuite, $curProbs) = @{$curData};
        print OUT join( $delim, suiteNumOut($curSuite->fullNumber), map($curProbs->{$ARG}, @{$rotList})) . "\n";
    }
    
    close(OUT);
}


sub printBestRots {
    #print a table of the best rotamer for each suite
    #ARGUMENTS:
    #   $output         - the filename to output to
    #   $predictedProbs - a list of [$suite, $probs] for each suite, where $suite is a Suite object
    #                     and $probs is a hash of probabilities for each rotamer
    #   $pseudoPredic   - an initialized PseudoPredictor object (used to get the puckers for each rotamer)
    #   $numStrucs      - the number of rotamers to print out for each suite
    #   $bestPath       - a list of the best (pucker-compatible) rotamers for each suite
    #RETURNS:
    #   none
    #EFFECTS:
    #   prints out a table to $output
    #   will die if $output cannot be written to
    
    my $output          = shift;
    my $predictedProbs  = shift;
    my $pseudoPredic    = shift;
    my $numStrucs       = shift;
    my $bestPath        = shift;
    
    open(OUT, ">", $output) or die "Could not open $output for writing\n";
    
    #print the header
    print OUT "Suite${delim}Predicted Structure";
    print OUT "$delim ${delim}Best Rotamer" . $delim x 2;
    for my $curStruc (2..$numStrucs) {
        print OUT $delim x 2 . "Rotamer $curStruc" . $delim x 2;
    }
    print OUT "\n";
    
    print OUT $delim x 2;
    for my $curStruc (1..$numStrucs) {
        print OUT join($delim, "Rotamer", "Score", "Puckers") . $delim x 2;
    }
    print OUT "\n";

    for (my $i = 0; $i < @{$predictedProbs}; $i++) {
        my ($curSuite, $curProbs) = @{$predictedProbs->[$i]};
        my $bestPathRot = $bestPath->[$i];
        
        my $curBestRots = [sort { $curProbs->{$b} <=> $curProbs->{$a} } keys(%{$curProbs})];
        
        #print out the formatted rotamer report
        print OUT $delim;
        for my $curStruc (0..($numStrucs-1)) { 
            print OUT $delim x 3 . $pseudoPredic->getStartPucker($curBestRots->[$curStruc]) . $delim;
        }
        print OUT "\n";
        
        print OUT suiteNumOut($curSuite->fullNumber) . $delim . $bestPathRot . $delim;
        for my $curStruc (0..($numStrucs-1)) {
            my $curRot = $curBestRots->[$curStruc];
            print OUT join($delim, $curRot, $curProbs->{$curRot}) . $delim x 3;
        }
        print OUT "\n";
        
        print OUT $delim;
        for my $curStruc (0..($numStrucs-1)) {
            print OUT $delim x 3 . $pseudoPredic->getEndPucker($curBestRots->[$curStruc]) . $delim;
        }
        print OUT "\n";
    }
    
    print OUT "\nBest rotamer string:$delim" . join("", @{$bestPath}) . "\n";


    close(OUT);
}

sub printTestData {
    #print out all of the data being used for predictions
    #ARGUMENTS:
    #   $output         - the filename to output to
    #   $predictedProbs - a list of [$suite, $probs] for each suite, where $suite is a Suite object
    #                     and $probs is a hash of probabilities for each rotamer
    #RETURNS:
    #   none
    #EFFECTS:
    #   prints out a table of data to $output
    #   will die if $output cannot be written to
    
    my $output          = shift;
    my $predictedProbs  = shift;
    
    open(OUT, ">", $output) or die "Could not open $output for writing\n";
    
    #print the header
    print OUT join($delim, qw(Suite Theta' Eta'), "Starting P-perp", "Ending P-perp", "Sugar (C1'-C1') Distance", "Starting Phosphate Distance", "Ending Phosphate Distance") . "\n";

    foreach my $curData (@{$predictedProbs}) {
        my ($curSuite, $curProbs) = @{$curData};
        
        print OUT join( $delim, suiteNumOut($curSuite->fullNumber),
                                     $curSuite->theta,
                                     $curSuite->eta,
                                     $curSuite->startingRes->pperp,
                                     $curSuite->endingRes->pperp,
                                     $curSuite->sugarDist,
                                     $curSuite->startingRes->phosDist,
                                     $curSuite->endingRes->phosDist) . "\n";
    }
    
    close(OUT);
}


#####################################################################
#####################################################################
########                 POD documentation                   ########
#####################################################################
#####################################################################

=head1 NAME

RCrane.CLI - Calculate and build the RNA backbone given phosphate and base atoms

=head1 SYNOPSIS

rcrane.cli.pl [OPTIONS] PDBFILE

 
Options:
   
   -r=STRING   Build the given conformer string
   
   -i=NAME     Alternate method for providing the input PDB file
   -o=NAME     Prefix for the output filenames
   
   -s          Do not build the structure
   -t          Do not output data tables
   -x          Format the output tables for Excel

   -h          Print a brief help message and exit
   -m          Print the full manual

=head1 DESCRIPTION

RCrane.CLI is designed to aid in crystallographic model building.  Specifically, it helps in building the RNA backbone into electron density.  To use RCrane.CLI, first build phosphate, base, and C1' atoms and save these coordinates to a PDB file.  Using this PDB file as input, RCrane.CLI will automatically construct the RNA backbone and output a PDB containing the constructed molecule.

RCrane.CLI builds the backbone using a combination of pseudotorsions and the all-atom consensus backbone conformer library.  The phosphates and bases are used to calculate pseudotorsions and determine sugar puckers, and these measurements are then used to predict and build the appropriate conformers.  For more information, see Keating KS and Pyle AM.  Semi-automated model building for RNA crystallography using a directed rotameric approach.  Proc Natl Acad Sci U S A, In press (2010).

By default, RCrane.CLI will output a PDB structure and several csv data tables.  See the C<-t> and C<-s> options for ways to modify this.  For an input of struc.pdb, the default output files will be

=over 2

=item F<struc.built.pdb>

A PDB file containing a molecule constructed using the most likely conformers

=item F<struc.data.csv>

A table containing all data used for the conformers predictions: theta', eta', base-phosphate perpendicular distances, C1'-C1' distances, and P-P distances

=item F<struc.probs.csv>

A table containing the probabilities for each conformer at each suite of the structure

=item F<struc.rots.csv>

A table containing the three most likely conformers for each suite, along with their sugar puckers.  The bottom of this table also contains a conformer string for the most likely conformers.

=back

=head1 EXAMPLES

=over 2

=item rcrane.cli.pl F<struc.pdb>

Build and predict conformers for file F<struc.pdb>.  Output files F<struc.built.pdb>, F<struc.data.csv>, F<struc.probs.csv>, and F<struc.rots.csv> will be created
    
=item rcrane.cli.pl -t F<struc.pdb>

Build and predict conformers for file F<struc.pdb>.  Only F<struc.built.pdb> will be output.
    
=item rcrane.cli.pl -r 1a1a1a7r6g1L1f1m1a F<struc.pdb>

Build F<struc.pdb> using the specified suites.  Note that F<struc.pdb> must contain exactly 9 suites.

=back

=head1 OPTIONS

=over 8

=item B<-r, --rot=STRING>

Build the given conformer string instead of predicting conformers.  Note that the given conformer string must contain the same number of suites as the input PDB file.  Implies C<-t>.

=item B<-i, --input=FILE>

Alternate method for providing the input PDB file.  If this option is used, then the PDB filename does not need to be given at the end of the arguments.

=item B<-o, --output=PREFIX>

Prefix for the output filenames.  If this option is used, output files will be named F<prefix.built.pdb>, F<prefix.data.csv>, F<prefix.probs.csv>, and F<prefix.rots.csv>.

=item B<-s, --nostruc>

Do not build the structure.  Appropriate conformers will be predicted and output, but not built.

=item B<-t, --notables>

Do not output the tables containing data about the conformer predictions.

=item B<-x, --excel>

Format the output F<.csv> data tables for Excel.  If this option is not used and the tables are opened in Excel (or OpenOffice.org if the importing defaults are used), then suite names may be incorrectly interpreted as dates - for example, suite 2-3 will be displayed as February 3rd of the current year.  This option will protect all suite names using formulas.

=item B<-h, -?, --help>

Print a brief help message and exit.

=item B<-m, --man>

Print the full manual and exit.

=back

=head1 AUTHOR

RCrane.CLI programmed by Kevin Keating, Yale University, 2010.  Kevin.Keating@yale.edu

All publications resulting from use of this program must acknowledge:
Keating KS and Pyle AM.  Semiautomated model building for RNA crystallography using a directed rotameric approach.  Proc Natl Acad Sci USA, 107: 8177-8182 (2010).

RCrane.CLI is copyright 2010, Kevin Keating, and is licensed under the Educational Community License, Version 2.0

=cut