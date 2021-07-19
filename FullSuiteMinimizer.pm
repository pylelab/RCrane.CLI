#!/usr/bin/perl

#module to minimize nucleotide structure in terms of bond lengths, angles, and torsions
#actual minimization is handled by the Math::Amoeba module

#this module requires the following CPAN modules:
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

package FullSuiteMinimizer;

use strict;
use English '-no_match_vars';
#use Math::Amoeba qw(MinimiseND);
use Data::Dumper;
use Exporter;
use Storable qw(dclone);
use PuckerList;

use MeasureTorsions;
use Rotate qw(rotateSugar);
use StrucCalc qw(plus magnitude minus angle calcDist);

our @ISA       = ('Exporter');
our @EXPORT    = qw(initializeFullSuiteMinimizer fullNtMinimize minInitAtoms);


my $rotData; #data on ideal torsion values for all rotamers - initialized in initializeFullSuiteMinimizer

#Math::Amoeba minimizes over all arguments to a function, so we keep the things that should stay constant in module variables
#(there are a few other ways to do this that might be slightly cleaner, such as generating functions to pass to Math::Amoeba,
#but this is simpler and works well enough)
my ($modulePrevResAtoms, $moduleCurResAtoms, $moduleNextResAtoms); #atomic coordinates for the previous, current, and next residue
my ($moduleCurRot, $moduleNextRot);     #the rotamer for the current and next suite
my ($moduleCurPhos, $moduleNextPhos);   #the coordinates for the current and next phosphate atom

my $numEvals = 0;                       #the number of times we have evaluated the objective function (evalNt)


use constant { IGNOREXI          => 0,  #whether to ignore the angle between the base and the sugar during the minimization
               IGNOREPHOSDEV     => 0,  #whether to ignore the penalty for moving the phosphate during mimization
               IGNORENEXTPHOSDEV => 0,  #whether to ignore the penalty for moving the phosphate of the next residue during mimization
};

#ideal bond lengths and angles taken from CNS 1.2 dna-rep_rep.param
use constant { C3pO3pIDEAL          => 1.423,       #ideal bond length for the C3'-O3' bond
               O3pPIDEAL            => 1.607,       #ideal bond length for the O3'-P bond
               PO5pIDEAL            => 1.593,       #ideal bond length for the P-O5' bond
               O5pC5pIDEAL          => 1.425,       #ideal bond length for the O5'-C5' bond
               C5pC4pIDEAL          => 1.510,       #ideal bond length for the C5'-C4' bond
               
               BONDSTDDEV           => 0.06,        #standard deviation for all bond lengths
                                                    #lower this value to make the bond length constraints stronger
               
               C3pANGLEIDEAL        => 110.5,       #ideal bond angle for the angle about the C3' atom
               O3pANGLEIDEAL        => 119.7,       #ideal bond angle for the angle about the O3' atom
               PANGLEIDEAL          => 104.0,       #ideal bond angle for the angle about the P atom
               O5pANGLEIDEAL        => 120.9,       #ideal bond angle for the angle about the O5' atom
               C5pANGLEIDEAL        => 110.2,       #ideal bond angle for the angle about the C5' atom
               C4pANGLEIDEAL        => 115.5,       #ideal bond angle for the angle about the C4' atom
               
               ANGLESTDDEV          => 4,           #standard deviation for all bond angles
                                                    #lower this value to make the bond angle constraints stronger
               
               TORSIONSTDDEVMOD     => (1/3),       #standard deviation modifier for all bond torsions
                                                    #standard deviations for any specific torsion in a given rotamer are multiplied by this value
                                                    #lower this value to make the bond torsion constraints stronger
               
               #Sugar/base angle (Xi) information calculated (by me) from RNA05 data set
               C3SUGARBASEIDEAL     => 241.0,       #the ideal sugar/base angle (xi) for C3' endo sugar puckers
               #C3SUGARBASEDEV       => 5.59,       #standard deviation for the sugar/base angle (xi) for C3' endo sugar puckers
               
               C2SUGARBASEIDEAL     => 214.2,       #the ideal sugar/base angle (xi) for C2' endo sugar puckers
               #C2SUGARBASEDEV       => 7.77,       #standard deviation for the sugar/base angle (xi) for C2' endo sugar puckers
               
               C3SUGARBASEDEV       => 2,           #these values are used instead of the ones above to make the sugar/base angle constraints stronger
               C2SUGARBASEDEV       => 2,           #lower these values to make the sugar/base angle constraints even stronger
               
               PHOSSTDDEV           => 0.1,         #standard deviation (in Angstrom) for moving the phosphate atoms
               NEXTPHOSSTDDEV       => 0.1,         #lower these values to make the phosphate constraints stronger
                                                    #(i.e. make it harder to move the phosphates during minimization)
};

#collect the sugar/base angle data into hashes
use constant SUGARBASEIDEAL => { 3 => C3SUGARBASEIDEAL, 2 => C2SUGARBASEIDEAL};
use constant SUGARBASEDEV   => { 3 => C3SUGARBASEDEV,   2 => C2SUGARBASEDEV};


#fullNtMinimize will perform a downhill simplex minimization run and restart the run until the minimum converges
#these constants determine the definition of convergence
use constant { MINSTOPPERCENTAGE => 0.999,  #the minimization will stop when the previous objective function
                                            #minimum is at least this percentage of the new objective function minimum
               MINSTOPCOUNT      => 25,     #the maximum number of restarts to allow, regardless of whether convergence has occurred
               SYNSUGARCUTOFF    => 200,    #if the minimization run starting with an anti sugar configuration doesn't
                                            #give an objective function minimum lower than this, then a new minimization
                                            #run will be started with a syn sugar
             };


#### Math:Amoeba ####

use Carp;
use constant TINY => 1e-16;

my ($ALPHA,$BETA,$GAMMA)=(1.0,0.5,2.0);

sub MinimiseND {
	my ($guesses,$scales,$func,$tol,$itmax, $verbose)=@_;
	my @p=ConstructVertices($guesses,$scales);
	my @y=EvaluateVertices(\@p,$func);
	return Amoeba(\@p,\@y,$func,$tol,$itmax, $verbose);
}

sub ConstructVertices {
	# given 2 vector references constructs an amoeba
	# returning the vertices
	my ($vector,$ofs)=@_;
	my $n=$#{$vector};
	my @vector=@{$vector};
	my (@p,@y,$i);

	$p[0]=[]; @{$p[0]}=@{$vector};
	for($i=0; $i<=$n; $i++) {
		my $v=[]; @{$v}=@{$vector};
		$v->[$i]+=$ofs->[$i];
		$p[$i+1]=$v;
	}
	return @p;
}

sub EvaluateVertices {
	# evaluates functions for all vertices of the amoeba
	my ($p,$func)=@_;
	my ($i,@y);
	for($i=0; $i<=$#{$p}; $i++) {
		$y[$i]=&$func(@{$p->[$i]});
	}
	return @y;
}

sub Amoeba {

    my ($p,$y,$func,$ftol,$itmax, $verbose)=@_;
	
    my $n=$#{$p}; # no points
    
	# Default parameters
	$verbose = (defined($verbose)) ? $verbose : 1;
	if (!$itmax) { $itmax=200; }
    if (!$ftol) { $ftol=1e-6; }

	# Member variables
    my ($i,$j);
    my $iter=0;
    my ($ilo, $inhi, $ihi);

	my ($pbar, $pr, $pe, $pc, $ypr, $ype, $ypc);

	# To control the recalculation of centroid
	my $recalc = 1;
	my $ihi_o;

	# Loop until any of stopping conditions hit 
	while (1)
	{
    	($ilo, $inhi, $ihi) = _FindMarkers($y);

		# Stopping conditions	
		my $rtol = 2*abs($y->[$ihi]-$y->[$ilo])/(abs($y->[$ihi])+abs($y->[$ilo])+TINY);
		if ($rtol<$ftol) { last; } 
		if ($iter++>$itmax) {
		  	carp "Amoeba exceeded maximum iterations\n" if ($verbose); 
		  	last;
		}

		# Determine the Centroid
		if ($recalc) {
			$pbar = _CalcCentroid($p, $ihi);
		} else {
			_AdjustCentroid($pbar, $p, $ihi_o, $ihi);
		}

		# Reset the re-calculation flag, and remember the current highest
		$recalc = 0;

		# Determine the reflection point, evaluate its value
		$pr = _CalcReflection($pbar, $p->[$ihi], $ALPHA);
		$ypr = &$func(@$pr);

		 # if it gives a better value than best point, try an
		 # additional extrapolation by a factor gamma, accept best
		if ($ypr < $y->[$ilo]) {

			$pe = _CalcReflection($pbar, $pr, -$GAMMA);
		    $ype=&$func(@$pe);
		    if ($ype < $y->[$ilo]) { 
				$p->[$ihi] = $pe; $y->[$ihi] = $ype; 
			}
		    else { 
				$p->[$ihi] = $pr; $y->[$ihi] = $ypr; 
			}
		}
		# if reflected point worse than 2nd highest
		elsif ($ypr >= $y->[$inhi]) {

			# if it is better than highest, replace it
		    if ($ypr < $y->[$ihi] ) {
 		        $p->[$ihi] = $pr; $y->[$ihi] = $ypr; 
		    }

		    # look for an intermediate lower point by performing a
		    # contraction of the simplex along one dimension
		    $pc = _CalcReflection($pbar, $p->[$ihi], -$BETA);
		    $ypc = &$func(@$pc);
		    
			# if contraction gives an improvement, accept it
			if ($ypc < $y->[$ihi]) { 		        
				$p->[$ihi] = $pc; $y->[$ihi] = $ypc;
		    }
			# otherwise cant seem to remove high point
			# so contract around lo (best) point
		    else {
				for($i=0; $i<=$n; $i++) {
					if ($i!=$ilo) {
						$p->[$i] = _CalcReflection($p->[$ilo], $p->[$i], -$BETA);
						$y->[$i] = &$func(@{$p->[$i]});
					}
		        }
				$recalc = 1;
		    }
		}
		# if reflected point is in-between lowest and 2nd highest 
		else {
		    $p->[$ihi] = $pr; $y->[$ihi] = $ypr;
		}
		
		# Remember the replacing position and its value
		$ihi_o = $ihi;
	}

	return ($p->[$ilo],$y->[$ilo]);
}

# Helper function - find the lowest, 2nd highest and highest position
sub _FindMarkers
{
	my $y = shift;
    
	my ($ilo, $inhi, $ihi);
	my ($i, $n);

	$n = @$y - 1;
	
	$ilo=0;
	if ($y->[0]>$y->[1]) {
	    $ihi=0; $inhi=1;
	}
	else {
	    $ihi=1; $inhi=0;
	}

	for($i = 0; $i <= $n; $i++) {
	    if ($y->[$i] < $y->[$ilo]) { $ilo = $i; }
	    if ($y->[$i] > $y->[$ihi]) { $inhi = $ihi; $ihi = $i; }
	    elsif ($y->[$i] > $y->[$inhi] && $ihi != $i) { $inhi = $i; }
	}

	return ($ilo, $inhi, $ihi);
}

# Determine the centroid (except the highest point)		
sub _CalcCentroid
{
	my ($p, $ihi) = @_;
	my ($i, $j, $n);
	
	$n = @$p - 1; 

	my $pbar = [];
	for($j=0; $j<$n; $j++) {
		for($i=0; $i<=$n; $i++) {
		    if ($i!=$ihi) {
				$pbar->[$j] += $p->[$i][$j];
	        }
	    }
		$pbar->[$j] /= $n;
	}

	return $pbar;
}

# Adjust the centroid only
sub _AdjustCentroid
{
	my ($pbar, $p, $ihi_o, $ihi) = @_;
	my ($j, $n);

	$n = @$pbar;
	
	if ($ihi_o != $ihi) {
		for($j=0; $j<$n; $j++) {
			$pbar->[$j] += ($p->[$ihi_o][$j] - $p->[$ihi][$j]) / $n;
		}
	}
}	

# Determine the reflection point
sub _CalcReflection
{
	my ($p1, $p2, $scale) = @_;
	my $j;

	my $n = @$p1;

	my $pr = [];
	for($j=0; $j<$n; $j++) {
	    $pr->[$j] = $p1->[$j] + $scale*($p1->[$j]-$p2->[$j]);
	}

	return $pr;
}

#### Math:Amoeba end ####


sub initializeFullSuiteMinimizer {
    #read in the data on torsion means and standard deviations
    #ARGUMENTS:
    #   $filename - a csv file containing data about torsion means and standard deviations for each rotamer
    #OPTIONAL ARGUMENTS:
    #   $verbose  - whether to print the random number seed to standard out
    #RETURNS:
    #   none
    #SIDE-EFFECTS:
    #   reads the data from $filename into the $rotData hash
    #   if $verbose is true, the random number seed will be printed to standard out
    
    my $filename = shift;
    my $verbose  = shift;
    
    open(IN, $filename) or die "Could not open $filename for reading\n";
    <IN>; #skip the header line
    
    while (my $curline = <IN>) {
        my @curdata = split(",", $curline);
        
        my ($rot, undef, undef, $prevDeltaMean, $epMean, $zetaMean, $alphaMean, $betaMean, $gammaMean, $curDeltaMean, undef, undef, undef, $prevDeltaSD, $epSD, $zetaSD, $alphaSD, $betaSD, $gammaSD, $curDeltaSD) = @curdata;
        
        $rotData->{$rot} = [$prevDeltaMean, $epMean, $zetaMean, $alphaMean, $betaMean, $gammaMean, $curDeltaMean, $prevDeltaSD, $epSD, $zetaSD, $alphaSD, $betaSD, $gammaSD, $curDeltaSD];
    }
    
    close(IN);
    
    if ($verbose) {
        #seed the random number generator and report on the seed used
        #this is only for debugging purposes (so that we can recreate a run if we want to)
        #the default seed is probably a better source of entropy
        #but we're not really concerned with *how* random the numbers are, so this is more than good enough
        my $seed = int(time ^ $$);
        $seed = 1234800781;        #if you want to recreate a run, uncomment this line and put the appropriate seed here
        print "Seeding random number generator with $seed\n";
        srand($seed);
    }

}


sub modifyNt {
    #calculate the coordinates given the minimization values
    #ARGUMENTS:
    #   $prevResAtoms                  - a hash containing atomic coordinates for the previous residue
    #   $curResAtoms                   - a hash containing atomic coordinates for the current residue
    #   $nextResAtoms                  - a hash containing atomic coordinates for the next residue
    #                                    leave as undef if there is no next residue (i.e. if this is the last residue of a segment)
    #   $chiRotation                   - how much to rotate chi by
    #   $sugarBaseAngle                - how much to rotate the sugar/base angle (xi) by
    #   $prevO3x, $prevO3y, $prevO3z   - x,y,z coordinates for how much to move the previous O3' atom
    #   $Px, $Py, $Pz                  - x,y,z coordinates for how much to move the phosphate atom
    #   $O5x, $O5y, $O5z               - x,y,z coordinates for how much to move the O5' atom
    #   $C5x, $C5y, $C5z               - x,y,z coordinates for how much to move the C5' atom
    #   $O3x, $O3y, $O3z               - x,y,z coordinates for how much to move the O3' atom
    #   $nextPx, $nextPy, $nextPz      - x,y,z coordinates for how much to move the next phosphate atom
    #   $nextO5x, $nextO5y, $nextO5z   - x,y,z coordinates for how much to move the next O5' atom
    #RETURNS:
    #   $prevResAtoms                  - a hash containing the modified atomic coordinates for the previous residue
    #   $curResAtoms                   - a hash containing the modified atomic coordinates for the current residue
    #   $nextResAtoms                  - a hash containing the modified atomic coordinates for the next residue
    
    
    my $prevResAtoms = shift;
    my $curResAtoms  = shift;
    my $nextResAtoms = shift;
    my $chiRotation    = shift;
    my $sugarBaseAngle = shift;
    my ($prevO3x, $prevO3y, $prevO3z) = (shift, shift, shift);
    my ($Px, $Py, $Pz)    = (shift, shift, shift);
    my ($O5x, $O5y, $O5z) = (shift, shift, shift);
    my ($C5x, $C5y, $C5z) = (shift, shift, shift);
    my ($O3x, $O3y, $O3z) = (shift, shift, shift);
    my ($nextPx, $nextPy, $nextPz) = (shift, shift, shift);
    my ($nextO5x, $nextO5y, $nextO5z) = (shift, shift, shift);
    
    #rotate the sugar atoms
    my $rotatedSugar = rotateSugar($curResAtoms, $chiRotation, 2, $sugarBaseAngle);
    $curResAtoms = {(%{$curResAtoms}, %{$rotatedSugar})};
    
    #do a shallow copy of $prevResAtoms and $nextResAtoms so that we don't modify the input variable (i.e. so that there aren't side-effects of this function)
    #the previous line (where we combine $curResAtoms and $rotatedSugar) already does a shallow copy of $curResAtoms
    $prevResAtoms = { (%{$prevResAtoms}) };
    $nextResAtoms = { (%{$nextResAtoms}) } if defined $nextResAtoms;
    
    #modify the coordinates
    $prevResAtoms->{"O3'"} = plus($prevResAtoms->{"O3'"}, [$prevO3x, $prevO3y, $prevO3z]);
    $curResAtoms ->{"P"}   = plus($curResAtoms ->{"P"},   [$Px, $Py, $Pz]);
    $curResAtoms ->{"O5'"} = plus($curResAtoms ->{"O5'"}, [$O5x, $O5y, $O5z]);
    $curResAtoms ->{"C5'"} = plus($curResAtoms ->{"C5'"}, [$C5x, $C5y, $C5z]);
    $curResAtoms ->{"O3'"} = plus($curResAtoms ->{"O3'"}, [$O3x, $O3y, $O3z]);
    
    if (defined $nextResAtoms) {
        $nextResAtoms->{"P"}   = plus($nextResAtoms->{"P"},   [$nextPx, $nextPy, $nextPz]) ;
        if (exists($nextResAtoms->{"O5'"})) {
            $nextResAtoms->{"O5'"} = plus($nextResAtoms->{"O5'"}, [$nextO5x, $nextO5y, $nextO5z]);
        }
    }
    
    return($prevResAtoms, $curResAtoms, $nextResAtoms);
}

sub evalNt {
    #Calculate the objective function value (i.e. how good of a nucleotide structure is this)
    #ARGUMENTS:
    #   $prevResAtoms                  - a hash containing atomic coordinates for the previous residue
    #   $curResAtoms                   - a hash containing atomic coordinates for the current residue
    #   $nextResAtoms                  - a hash containing atomic coordinates for the next residue, if any
    #                                    should be undef otherwise
    #   $curRot                        - the rotamer for the current suite
    #   $nextRot                       - the rotamer for the next suite (if there is one)
    #                                    should be undef otherwise
    #   $chiRotation                   - how much to rotate chi by
    #   $sugarBaseAngle                - how much to rotate the sugar/base angle (xi) by
    #   $prevO3x, $prevO3y, $prevO3z   - x,y,z coordinates for how much to move the previous O3' atom
    #   $Px, $Py, $Pz                  - x,y,z coordinates for how much to move the phosphate atom
    #   $O5x, $O5y, $O5z               - x,y,z coordinates for how much to move the O5' atom
    #   $C5x, $C5y, $C5z               - x,y,z coordinates for how much to move the C5' atom
    #   $O3x, $O3y, $O3z               - x,y,z coordinates for how much to move the O3' atom
    #   $nextPx, $nextPy, $nextPz      - x,y,z coordinates for how much to move the next phosphate atom, if present
    #                                    should be undef otherwise
    #   $nextO5x, $nextO5y, $nextO5z   - x,y,z coordinates for how much to move the next O5' atom
    #                                    should be undef otherwise
    #OPTIONAL ARGUMENTS
    #   $verbose                       - how much information (if any) about this structure to print to standard out
    #                                     if not given, no information is printed
    #   $resetEvalCount                - whether or not to reset the count of how many times this function has been run
    #                                     if not given, defaults to $verbose
    #RETURNS:
    #   $totalDev                       - the value of the objective function (a lower number means a better nucleotide structure)
    #SIDE-EFFECTS:
    #   if $verbose if true, prints information to standard out
    #   if $resetEvalCount is true, sets the module variable $numEvals to 0
    
    my $prevResAtoms = shift;
    my $curResAtoms  = shift;
    my $nextResAtoms = shift;
    my $curRot   = shift;
    my $nextRot  = shift;
    my $chiRotation    = shift;
    my $sugarBaseAngle = shift;
    my ($prevO3x, $prevO3y, $prevO3z) = (shift, shift, shift);
    my ($Px, $Py, $Pz)    = (shift, shift, shift);
    my ($O5x, $O5y, $O5z) = (shift, shift, shift);
    my ($C5x, $C5y, $C5z) = (shift, shift, shift);
    my ($O3x, $O3y, $O3z) = (shift, shift, shift);
    my ($nextPx, $nextPy, $nextPz) = (shift, shift, shift);
    my ($nextO5x, $nextO5y, $nextO5z) = (shift, shift, shift);
    my $verbose = shift;
    my $resetEvalCount = shift;
    
    $resetEvalCount = $verbose unless defined $resetEvalCount;
    
    my $totalDev = 0; #the value of the objective function
    
    #modify the coordinates
    ($prevResAtoms, $curResAtoms, $nextResAtoms) = modifyNt($prevResAtoms, $curResAtoms, $nextResAtoms, $chiRotation, $sugarBaseAngle, $prevO3x, $prevO3y, $prevO3z, $Px, $Py, $Pz, $O5x, $O5y, $O5z, $C5x, $C5y, $C5z, $O3x, $O3y, $O3z, $nextPx, $nextPy, $nextPz, $nextO5x, $nextO5y, $nextO5z);
    
    #calculate torsions
    $prevResAtoms->{"P+1"}   = $curResAtoms ->{"P"};
    $prevResAtoms->{"O5'+1"} = $curResAtoms ->{"O5'"};
    $curResAtoms ->{"O3'-1"} = $prevResAtoms->{"O3'"};
    if (defined $nextResAtoms) {
        $curResAtoms ->{"P+1"}   = $nextResAtoms->{"P"};
        $curResAtoms ->{"O5'+1"} = $nextResAtoms->{"O5'"};
    }
    
    my ($prevDelta, $prevEpsilon,$prevZeta, $alpha, $beta, $gamma, $delta, $epsilon, $zeta);
    $prevDelta   = calcDelta($prevResAtoms);
    $prevEpsilon = calcEpsilon($prevResAtoms);
    $prevZeta    = calcZeta($prevResAtoms);
    $alpha       = calcAlpha($curResAtoms);
    $beta        = calcBeta($curResAtoms);
    $gamma       = calcGamma($curResAtoms);
    $delta       = calcDelta($curResAtoms);
    if (defined $nextRot) {
        $epsilon     = calcEpsilon($curResAtoms);
        $zeta        = calcZeta($curResAtoms);
    }
    
    my ($prevDeltaMean, $prevEpMean, $prevZetaMean, $alphaMean, $betaMean, $gammaMean, $deltaMean, $prevDeltaSD, $prevEpSD, $prevZetaSD, $alphaSD, $betaSD, $gammaSD, $deltaSD) = @{$rotData->{$curRot}};
    my ($epMean, $zetaMean, $epSD, $zetaSD);
    if (defined $nextRot) {
        $epMean   = $rotData->{$nextRot}->[1];
        $zetaMean = $rotData->{$nextRot}->[2];
        $epSD     = $rotData->{$nextRot}->[8];
        $zetaSD   = $rotData->{$nextRot}->[9];
    }
    
    
    $totalDev += (($prevDelta     - $prevDeltaMean)/($prevDeltaSD*TORSIONSTDDEVMOD))**2;
    $totalDev += (($prevEpsilon   - $prevEpMean)   /($prevEpSD   *TORSIONSTDDEVMOD))**2;
    $totalDev += (($prevZeta      - $prevZetaMean) /($prevZetaSD *TORSIONSTDDEVMOD))**2;
    $totalDev += (($alpha         - $alphaMean)    /($alphaSD    *TORSIONSTDDEVMOD))**2;
    $totalDev += (($beta          - $betaMean)     /($betaSD     *TORSIONSTDDEVMOD))**2;
    $totalDev += (($gamma         - $gammaMean)    /($gammaSD    *TORSIONSTDDEVMOD))**2;
    $totalDev += (($delta         - $deltaMean)    /($deltaSD    *TORSIONSTDDEVMOD))**2;
    if (defined $nextRot) {
        $totalDev += (($epsilon       - $epMean)       /($epSD       *TORSIONSTDDEVMOD))**2;
        $totalDev += (($zeta          - $zetaMean)     /($zetaSD     *TORSIONSTDDEVMOD))**2;
    }
    
    
    
    #calculate bond lengths
    my ($prevc3po3pBond, $prevo3pPBond, $Po5pBond, $o5pc5pBond, $c5pc4pBond, $c3po3pBond, $o3pPBond, $nextPo5pBond);
    $prevc3po3pBond  = magnitude(minus($prevResAtoms->{"C3'"}, $prevResAtoms->{"O3'"} ));
    $prevo3pPBond    = magnitude(minus($prevResAtoms->{"O3'"}, $curResAtoms ->{"P"}   ));
    $Po5pBond        = magnitude(minus($curResAtoms ->{"P"},   $curResAtoms ->{"O5'"} ));
    $o5pc5pBond      = magnitude(minus($curResAtoms ->{"O5'"}, $curResAtoms ->{"C5'"} ));
    $c5pc4pBond      = magnitude(minus($curResAtoms ->{"C5'"}, $curResAtoms ->{"C4'"} ));
    $c3po3pBond      = magnitude(minus($curResAtoms ->{"C3'"}, $curResAtoms ->{"O3'"} ));
    if (defined $nextResAtoms and exists $nextResAtoms->{"P"}) {
        $o3pPBond        = magnitude(minus($curResAtoms ->{"O3'"}, $nextResAtoms->{"P"}   ));
        if (exists $nextResAtoms->{"O5'"}) {
            $nextPo5pBond    = magnitude(minus($nextResAtoms->{"P"},   $nextResAtoms->{"O5'"} ));
        }
    }
    
    #sum up the bond length deviation
    $totalDev += (($prevc3po3pBond - C3pO3pIDEAL)/BONDSTDDEV)**2;
    $totalDev += (($prevo3pPBond   - O3pPIDEAL)  /BONDSTDDEV)**2;
    $totalDev += (($Po5pBond       - PO5pIDEAL)  /BONDSTDDEV)**2;
    $totalDev += (($o5pc5pBond     - O5pC5pIDEAL)/BONDSTDDEV)**2;
    $totalDev += (($c5pc4pBond     - C5pC4pIDEAL)/BONDSTDDEV)**2;
    $totalDev += (($c3po3pBond     - C3pO3pIDEAL)/BONDSTDDEV)**2;
    $totalDev += (($o3pPBond       - O3pPIDEAL)  /BONDSTDDEV)**2 if defined $o3pPBond;
    $totalDev += (($nextPo5pBond   - PO5pIDEAL)  /BONDSTDDEV)**2 if defined $nextPo5pBond;
    
    
    
    #calculate bond angles
    my ($prevc3pAngle, $prevo3pAngle, $PAngle, $o5pAngle, $c5pAngle, $c4pAngle, $c3pAngle, $o3pAngle, $nextPAngle);
    $prevc3pAngle = angle($prevResAtoms->{"C4'"}, $prevResAtoms->{"C3'"}, $prevResAtoms->{"O3'"} );
    $prevo3pAngle = angle($prevResAtoms->{"C3'"}, $prevResAtoms->{"O3'"}, $curResAtoms ->{"P"}   );
    $PAngle       = angle($prevResAtoms->{"O3'"}, $curResAtoms ->{"P"},   $curResAtoms ->{"O5'"} );
    $o5pAngle     = angle($curResAtoms ->{"P"},   $curResAtoms ->{"O5'"}, $curResAtoms ->{"C5'"} );
    $c5pAngle     = angle($curResAtoms ->{"O5'"}, $curResAtoms ->{"C5'"}, $curResAtoms ->{"C4'"} );
    $c4pAngle     = angle($curResAtoms ->{"C5'"}, $curResAtoms ->{"C4'"}, $curResAtoms ->{"C3'"} );
    $c3pAngle     = angle($curResAtoms ->{"C4'"}, $curResAtoms ->{"C3'"}, $curResAtoms ->{"O3'"} );
    if (defined $nextResAtoms and exists $nextResAtoms->{"P"}) {
        $o3pAngle     = angle($curResAtoms ->{"C3'"}, $curResAtoms ->{"O3'"}, $nextResAtoms->{"P"}   );
        if (exists $nextResAtoms->{"O5'"}) {
            $nextPAngle   = angle($curResAtoms ->{"O3'"}, $nextResAtoms->{"P"},   $nextResAtoms->{"O5'"} );
        }
    }
    
    #sum up the angle deviations
    $totalDev += (($prevc3pAngle - C3pANGLEIDEAL)/ANGLESTDDEV)**2;
    $totalDev += (($prevo3pAngle - O3pANGLEIDEAL)/ANGLESTDDEV)**2;
    $totalDev += (($PAngle       - PANGLEIDEAL)  /ANGLESTDDEV)**2;
    $totalDev += (($o5pAngle     - O5pANGLEIDEAL)/ANGLESTDDEV)**2;
    $totalDev += (($c5pAngle     - C5pANGLEIDEAL)/ANGLESTDDEV)**2;
    $totalDev += (($c4pAngle     - C4pANGLEIDEAL)/ANGLESTDDEV)**2;
    $totalDev += (($c3pAngle     - C3pANGLEIDEAL)/ANGLESTDDEV)**2;
    $totalDev += (($o3pAngle     - O3pANGLEIDEAL)/ANGLESTDDEV)**2 if defined $o3pAngle;
    $totalDev += (($nextPAngle   - PANGLEIDEAL)  /ANGLESTDDEV)**2 if defined $nextPAngle;
    
    #calculate the base/sugar angle (assuming we have the necessary atoms, otherwise, ignore this term)
    my $sugarBaseAngleReal = "notcalc";
    if (not IGNOREXI and ($curResAtoms->{"N9"} or $curResAtoms->{"N1"})) {
        $sugarBaseAngleReal = calcXi($curResAtoms);
        $totalDev += (($sugarBaseAngleReal - (SUGARBASEIDEAL->{$puckerList->{$curRot}->[1]})) / (SUGARBASEDEV->{$puckerList->{$curRot}->[1]})) **2;
    }
    
    #calculate the movement of the phosphate
    my ($curPhosDist, $nextPhosDist);
    $curPhosDist  = calcDist($curResAtoms->{"P"}, $moduleCurPhos);
    $nextPhosDist = calcDist($nextResAtoms->{"P"}, $moduleNextPhos) if defined $nextResAtoms;
    unless (IGNOREPHOSDEV) {
        $totalDev += ($curPhosDist / PHOSSTDDEV)**2;
    }
    if ((not IGNORENEXTPHOSDEV) and (defined $nextResAtoms)) {
        $totalDev += ($nextPhosDist / NEXTPHOSSTDDEV)**2;
    }
    
    
    #if desired, print out information about the end point of the minimization
    if ($verbose) {
        print "\tNum of func evals: $numEvals\n";
        
        if ($verbose >= 2) {
            print "\tTorsions:\n";
            print "\t\tprevDelta:   $prevDelta\t"   . (($prevDelta     - $prevDeltaMean)/($prevDeltaSD*TORSIONSTDDEVMOD)) . "\n";
            print "\t\tprevEpsilon: $prevEpsilon\t" . (($prevEpsilon   - $prevEpMean)   /($prevEpSD   *TORSIONSTDDEVMOD)) . "\n";
            print "\t\tprevZeta:    $prevZeta\t"    . (($prevZeta      - $prevZetaMean) /($prevZetaSD *TORSIONSTDDEVMOD)) . "\n";
            print "\t\talpha:       $alpha\t"       . (($alpha         - $alphaMean)    /($alphaSD    *TORSIONSTDDEVMOD)) . "\n";
            print "\t\tbeta:        $beta\t"        . (($beta          - $betaMean)     /($betaSD     *TORSIONSTDDEVMOD)) . "\n";
            print "\t\tgamma:       $gamma\t"       . (($gamma         - $gammaMean)    /($gammaSD    *TORSIONSTDDEVMOD)) . "\n";
            print "\t\tdelta:       $delta\t"       . (($delta         - $deltaMean)    /($deltaSD    *TORSIONSTDDEVMOD)) . "\n";
            if (defined $nextRot) {
                print "\t\tepsilon:     $epsilon\t"     . (($epsilon       - $epMean)       /($epSD       *TORSIONSTDDEVMOD)) . "\n";
                print "\t\tzeta:        $zeta\t"        . (($zeta          - $zetaMean)     /($zetaSD     *TORSIONSTDDEVMOD)) . "\n";
            }
            
            print "\tBond lengths:\n";
            print "\t\tpC3'-O3' bond: $prevc3po3pBond\t" . (($prevc3po3pBond - C3pO3pIDEAL)/BONDSTDDEV) . "\n";
            print "\t\tpO3'-P bond:   $prevo3pPBond\t"   . (($prevo3pPBond   - O3pPIDEAL)  /BONDSTDDEV) . "\n";
            print "\t\tP-O5' bond:    $Po5pBond\t"       . (($Po5pBond       - PO5pIDEAL)  /BONDSTDDEV) . "\n";
            print "\t\tO5'-C5' bond:  $o5pc5pBond\t"     . (($o5pc5pBond     - O5pC5pIDEAL)/BONDSTDDEV) . "\n";
            print "\t\tC5'-C4' bond:  $c5pc4pBond\t"     . (($c5pc4pBond     - C5pC4pIDEAL)/BONDSTDDEV) . "\n";
            print "\t\tC3'-O3' bond:  $c3po3pBond\t"     . (($c3po3pBond     - C3pO3pIDEAL)/BONDSTDDEV) . "\n";
            print "\t\tO3'-P' bond:   $o3pPBond\t"       . (($o3pPBond       - O3pPIDEAL)  /BONDSTDDEV) . "\n" if defined $o3pPBond;
            print "\t\tnP-O5' bond:   $nextPo5pBond\t"   . (($nextPo5pBond   - PO5pIDEAL)  /BONDSTDDEV) . "\n" if defined $nextPo5pBond;
                                    
            print "\tBond angles\n";
            print "\t\tpC3' angle:  $prevc3pAngle\t" . (($prevc3pAngle - C3pANGLEIDEAL)/ANGLESTDDEV) . "\n";
            print "\t\tpO3' angle:  $prevo3pAngle\t" . (($prevo3pAngle - O3pANGLEIDEAL)/ANGLESTDDEV) . "\n";
            print "\t\tP angle:     $PAngle\t"       . (($PAngle       - PANGLEIDEAL)  /ANGLESTDDEV) . "\n";
            print "\t\tO5' angle:   $o5pAngle\t"     . (($o5pAngle     - O5pANGLEIDEAL)/ANGLESTDDEV) . "\n";
            print "\t\tC5' angle:   $c5pAngle\t"     . (($c5pAngle     - C5pANGLEIDEAL)/ANGLESTDDEV) . "\n";
            print "\t\tC4' angle:   $c4pAngle\t"     . (($c4pAngle     - C4pANGLEIDEAL)/ANGLESTDDEV) . "\n";
            print "\t\tC3' angle:   $c3pAngle\t"     . (($c3pAngle     - C3pANGLEIDEAL)/ANGLESTDDEV) . "\n";
            print "\t\tO3' angle:   $o3pAngle\t"     . (($o3pAngle     - O3pANGLEIDEAL)/ANGLESTDDEV) . "\n" if defined $o3pAngle;
            print "\t\tnP angle:    $nextPAngle\t"   . (($nextPAngle   - PANGLEIDEAL)  /ANGLESTDDEV) . "\n" if defined $nextPAngle;
                                  
            print "\tSugar and phosphates:\n";
            print "\t\tsugarBaseAngle: $sugarBaseAngleReal\t" . (($sugarBaseAngleReal - (SUGARBASEIDEAL->{$puckerList->{$curRot}->[1]})) / (SUGARBASEDEV->{$puckerList->{$curRot}->[1]})) . "\n";
            print "\t\tcurPhosDev:     $curPhosDist\t"        . ($curPhosDist  / PHOSSTDDEV)     . "\n";
            print "\t\tnextPhosDev:    $nextPhosDist\t"       . ($nextPhosDist / NEXTPHOSSTDDEV) . "\n" if defined $nextPhosDist;
        }
        
        print "\tmin value: $totalDev\n\n";
        
        print "\n" if $verbose >= 2;
        
    }
    
    $numEvals++;
    $numEvals = 0 if $resetEvalCount;
    
    return $totalDev;
}

sub evalNtCaller {
    #call evalSuite with the module variable values for $atoms and $rot
    #(this allows Math::Amoeba to evaluate evalSuite)
    
    return evalNt($modulePrevResAtoms, $moduleCurResAtoms, $moduleNextResAtoms, $moduleCurRot, $moduleNextRot, @ARG);
}

sub fullNtMinimizeSingle {
    #do a single minimization run (i.e. a single run of downhill simplex without restarts)
    #ARGUMENTS:
    #   $prevResAtoms       - a hash containing atomic coordinates for the previous residue
    #   $curResAtoms        - a hash containing atomic coordinates for the current residue
    #OPTIONAL ARGUMENTS:
    #   $nextResAtoms       - a hash containing atomic coordinates for the next residue
    #                         if undef, no next residue coordinates will be minimized
    #   $verbose            - how much information (if any) about this structure to print to standard out
    #                         if not given, no information is printed 
    #RETURNS:
    #   $prevResAtoms       - a hash containing the modified atomic coordinates for the previous residue
    #   $curResAtoms        - a hash containing the modified atomic coordinates for the current residue
    #   $nextResAtoms       - a hash containing the modified atomic coordinates for the next residue
    #   $minObjectiveVal    - the minimized value of the objective function
    #                         (a lower number means a better nucleotide structure)
    
    my $prevResAtoms = shift;
    my $curResAtoms  = shift;
    my $nextResAtoms = shift;
    my $verbose      = shift; #optional
    
    #put the atom coordinates into module variables
    #so that evalSuiteCaller can access them (which allows them to be used with MinimiseND)
    $modulePrevResAtoms = {(%{$prevResAtoms})}; #shallow copy
    $moduleCurResAtoms  = {(%{$curResAtoms})};  #shallow copy
    if (defined $nextResAtoms) {
        $moduleNextResAtoms = {(%{$nextResAtoms})}; #shallow copy
    } else {
        undef $moduleNextResAtoms;
    }
    
    #randomize the signs for the standard deviations (to help us jump out of local minima)
    my $minimizeSDs = randomizeSign([10,  3,   1,    1,    1,    1,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,    1,    1]);
    
    #calculate the minimumized suite structure
                                    #             chi, xi,  pO3x, pO3y, pO3z, Px, Py, Pz, O5x, O5y, O5z, C5x, C5y, C5z, O3x, O3y, O3z, nPx, nPy, nPz, nO5x, nO5y, nO5z
    my ($minVals, $minObjectiveVal) = MinimiseND([0,   0,   0,    0,    0,    0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,    0],
                                                 $minimizeSDs,
                                                 \&evalNtCaller, 0.001, 5000, 0);
    
    
    if ($verbose) {
        #if desired, do another run of the evaluation function so that we can print out the verbose information
        evalNtCaller(@{$minVals}, $verbose);
    }
    
    #calculate the new coordinates
    ($prevResAtoms, $curResAtoms, $nextResAtoms) = modifyNt($prevResAtoms, $curResAtoms, $nextResAtoms, @{$minVals});
    
    return($prevResAtoms, $curResAtoms, $nextResAtoms, $minObjectiveVal);
}

sub fullNtMinimize {
    #minimize the structure of a nucleotide (relative to bond lengths, angles, and torsions)
    #this function will perform a downhill simplex minimization run and restart the run until the minimum converges
    #it will also try a minimization run starting with a syn sugar if the run with an anti sugar didn't produce a good result
    #ARGUMENTS:
    #   $prevResAtoms       - a hash containing atomic coordinates for the previous residue
    #   $curResAtoms        - a hash containing atomic coordinates for the current residue
    #   $nextResAtoms       - a hash containing atomic coordinates for the next residue, if any
    #                         should be undef otherwise
    #   $curRot             - the rotamer for the current suite
    #   $nextRot            - the rotamer for the next suite, if present
    #                         should be undef otherwise
    #OPTIONAL ARGUMENTS:
    #   $phosLoc            - the coordinates of the phosphate to use for calculating how much it's moved
    #                         this should be the location of the phosphate before our previous minimization run moved it
    #                         if not given, the phosphate location in $curResAtoms is used
    #   $verbose            - how much information (if any) about this structure to print to standard out
    #                         if not given, no information is printed 
    #RETURNS:
    #   $prevResAtoms       - a hash containing the modified atomic coordinates for the previous residue
    #   $curResAtoms        - a hash containing the modified atomic coordinates for the current residue
    #   $nextResAtoms       - a hash containing the modified atomic coordinates for the next residue
    #   $curMinVal          - the minimized value of the objective function
    #                         (a lower number means a better nucleotide structure)
    
    my $prevResAtoms = shift;
    my $curResAtoms  = shift;
    my $nextResAtoms = shift;
    my $curRot       = shift;
    my $nextRot      = shift;
    my $phosLoc      = shift; #optional
    my $verbose      = shift; #optional
    
    #put the rotamers and phosphate locations into module variables
    #so that evalSuiteCaller can access them (which allows them to be used with MinimiseND)
    $moduleCurRot  = $curRot;
    $moduleNextRot = $nextRot;
    
    if (defined $nextResAtoms) {
        $moduleNextPhos = $nextResAtoms->{"P"};
    } else {
        $moduleNextPhos = undef;
    }
    
    #if we were given a phosphate location (presumably its location before we did the previous round of minimization), use that
    if (defined $phosLoc) {
        $moduleCurPhos = $phosLoc;
    } else {
        #otherwise use the current phosphate location
        $moduleCurPhos = $curResAtoms ->{"P"};
    }
    
    #remember the atomic coordinates so that we can restart the minimization with a syn sugar if we need to
    my $synPrevResAtoms = dclone($prevResAtoms);
    my $synCurResAtoms  = dclone($curResAtoms);
    my $synNextResAtoms;
    $synNextResAtoms = dclone($nextResAtoms) if defined $nextResAtoms;
    
    #do the first minimization run
    my ($curMinVal, $prevMinVal, $numMinRuns);
    ($prevResAtoms, $curResAtoms, $nextResAtoms, $prevMinVal) = fullNtMinimizeSingle($prevResAtoms, $curResAtoms, $nextResAtoms, $verbose);
    
    #keep restarting the minimization until the objective function value stops going down or we've done MINSTOPCOUNT minimization runs
    for ($numMinRuns = 0; $numMinRuns < MINSTOPCOUNT; $numMinRuns++) {
        ($prevResAtoms, $curResAtoms, $nextResAtoms, $curMinVal) = fullNtMinimizeSingle($prevResAtoms, $curResAtoms, $nextResAtoms, $verbose);
        last if (($curMinVal / $prevMinVal) > MINSTOPPERCENTAGE);
        $prevMinVal = $curMinVal;
    }
    
    #if this minimization run doesn't give us a great answer, then try again starting from a syn sugar
    #(first, make sure that we have the atoms required to calculate chi, otherwise we can't rotate the sugar to syn)
    if ($curMinVal > SYNSUGARCUTOFF) {
        unless (($synCurResAtoms->{"N9"} and $synCurResAtoms->{"C4"}) or ($synCurResAtoms->{"N1"} and $synCurResAtoms->{"C2"})) {
            print "Cannot rebuild nt with a syn sugar - no base present\n" if $verbose;
        } else {
            print "Rebuilding nt with a syn sugar\n" if $verbose;
            my $antiMinVal = $curMinVal;
            
            #rotate the sugar to a syn value
            my $curChi = calcChi($synCurResAtoms);
            my $synSugar = rotateSugar($synCurResAtoms,  70 - $curChi);
            $synCurResAtoms = { %{$synCurResAtoms}, %{$synSugar} };
            
            #do another minimization run
            ($synPrevResAtoms, $synCurResAtoms, $synNextResAtoms, $prevMinVal) = fullNtMinimizeSingle($synPrevResAtoms, $synCurResAtoms, $synNextResAtoms, $verbose);
            for ($numMinRuns = 0; $numMinRuns < MINSTOPCOUNT; $numMinRuns++) {
                ($synPrevResAtoms, $synCurResAtoms, $synNextResAtoms, $curMinVal) = fullNtMinimizeSingle($synPrevResAtoms, $synCurResAtoms, $synNextResAtoms, $verbose);
                last if ($curMinVal < 0.1 or (($curMinVal / $prevMinVal) > MINSTOPPERCENTAGE));
                $prevMinVal = $curMinVal;
            }
            
            #if this worked better than the anti minimization, then use these coordinates
            if ($curMinVal < $antiMinVal) {
                print "Using syn structure\n" if $verbose;
                $prevResAtoms = $synPrevResAtoms;
                $curResAtoms  = $synCurResAtoms;
                $nextResAtoms = $synNextResAtoms;
            } else {
                print "Using anti structure\n" if $verbose;
                #otherwise, set $curMinVal back to what it was (since we're about to return its value)
                $curMinVal = $antiMinVal;
            }
        }
    }
    
    return($prevResAtoms, $curResAtoms, $nextResAtoms, $curMinVal);
}

sub randomizeSign {
    #randomize the sign for all entries in an array
    #(used to randomize the sign for the initial test values used during minimization)
    #ARGUMENTS:
    #   $array - a list of numbers
    #RETURNS:
    #   $array - the array with the signs randomized
    
    my $array = shift;
    
    for my $i (@{$array}) {
        $i *= 1 - (2 * int(rand(2)));
    }
    
    return $array;
}


sub minInitAtoms {
    #minimize the structure of the initial P, O5', and C5' in a chain (relative to bond lengths and angles)
    #this function will perform a downhill simplex minimization run and restart the run until the minimum converges
    #ARGUMENTS:
    #   $atoms       - a hash containing atomic coordinates
    #   $verbose     - how much information (if any) about this structure to print to standard out
    #                  if not given, no information is printed 
    #RETURNS:
    #   $atoms       - a hash containing the modified atomic coordinates
    
    my $atoms = shift;
    my $verbose = shift;
    
    #put the phosphate location into a module variable
    #so that evalInitAtomsCaller can access them (which allows them to be used with MinimiseND)
    $moduleCurPhos = $atoms ->{"P"};
    
    my ($curMinVal, $prevMinVal, $numMinRuns);
    ($atoms, $prevMinVal) = minInitAtomsSingle($atoms, $verbose);
    
    #keep restarting the minimization until the objective function value stops going down or we've done MINSTOPCOUNT minimization runs
    for ($numMinRuns = 0; $numMinRuns < MINSTOPCOUNT; $numMinRuns++) {
        ($atoms, $curMinVal) = minInitAtomsSingle($atoms, $verbose);
        last if ($curMinVal < 0.1 or (($curMinVal / $prevMinVal) > MINSTOPPERCENTAGE));
        $prevMinVal = $curMinVal;
    }
    
    if (wantarray) {
        return ($atoms, $curMinVal);
    } else {
        return $atoms;
    }
}


sub minInitAtomsSingle {
    #do a single minimization run (i.e. a single run of downhill simplex without restarts)
    #ARGUMENTS:
    #   $atoms       - a hash containing atomic coordinates
    #   $verbose     - how much information (if any) about this structure to print to standard out
    #                  if not given, no information is printed 
    #RETURNS:
    #   $atoms       - a hash containing the modified atomic coordinates
    
    my $atoms = shift;
    my $verbose = shift;
    
    #put the atom coordinates into module variables
    #so that evalSuiteCaller can access them (which allows them to be used with MinimiseND)
    $moduleCurResAtoms  = {(%{$atoms})};  #shallow copy
    
    #randomize the signs for the standard deviations (to help us jump out of local minima)
                                    #Px, Py, Pz, O5x, O5y, O5z, C5x, C5y, C5z
    my $minimizeSDs = randomizeSign([1,  1,  1,   1,   1,   1,   1,   1,   1 ]);
    
                                                 #Px, Py, Pz, O5x, O5y, O5z, C5x, C5y, C5z
    my ($minVals, $minObjectiveVal) = MinimiseND([0,  0,  0,   0,   0,   0,   0,   0,   0],
                                                 $minimizeSDs, \&evalInitAtomsCaller, 0.001, 5000, 0);
    
    if ($verbose) {
        #if desired, do another run of the evaluation function so that we can print out the verbose information
        evalInitAtomsCaller(@{$minVals}, $verbose);
    }
    
    #calculate the new coordinates
    $atoms = modifyInitAtoms($atoms, @{$minVals});
    
    return ($atoms, $minObjectiveVal);
}


sub evalInitAtoms {
    #Calculate the objective function value (i.e. how good of a nucleotide structure is this)
    #ARGUMENTS:
    #   $atoms                  - a hash containing atomic coordinates
    #   $Px, $Py, $Pz           - x,y,z coordinates for how much to move the phosphate atom
    #   $O5x, $O5y, $O5z        - x,y,z coordinates for how much to move the O5' atom
    #   $C5x, $C5y, $C5z        - x,y,z coordinates for how much to move the C5' atom
    #OPTIONAL ARGUMENTS
    #   $verbose                - how much information (if any) about this structure to print to standard out
    #                             if not given, no information is printed
    #   $resetEvalCount         - whether or not to reset the count of how many times this function has been run
    #                             if not given, defaults to $verbose
    #RETURNS:
    #   $totalDev               - the value of the objective function (a lower number means a better nucleotide structure)
    
    my $atoms  = shift;
    my ($Px, $Py, $Pz)    = (shift, shift, shift);
    my ($O5x, $O5y, $O5z) = (shift, shift, shift);
    my ($C5x, $C5y, $C5z) = (shift, shift, shift);
    my $verbose = shift;
    my $resetEvalCount = shift;
    
    $resetEvalCount = $verbose unless defined $resetEvalCount;
    
    my $totalDev = 0; #the value of the objective function
    
    #modify the coordinates
    $atoms = modifyInitAtoms($atoms, $Px, $Py, $Pz, $O5x, $O5y, $O5z, $C5x, $C5y, $C5z);
    
    #calculate bond lengths
    my ($Po5pBond, $o5pc5pBond, $c5pc4pBond);
    $Po5pBond        = magnitude(minus($atoms ->{"P"},   $atoms ->{"O5'"} ));
    $o5pc5pBond      = magnitude(minus($atoms ->{"O5'"}, $atoms ->{"C5'"} ));
    $c5pc4pBond      = magnitude(minus($atoms ->{"C5'"}, $atoms ->{"C4'"} ));
    
    #sum up the bond length deviation
    $totalDev += (($Po5pBond       - PO5pIDEAL)  /BONDSTDDEV)**2;
    $totalDev += (($o5pc5pBond     - O5pC5pIDEAL)/BONDSTDDEV)**2;
    $totalDev += (($c5pc4pBond     - C5pC4pIDEAL)/BONDSTDDEV)**2;
    
    #calculate bond angles
    my ($o5pAngle, $c5pAngle, $c4pAngle);
    $o5pAngle     = angle($atoms ->{"P"},   $atoms ->{"O5'"}, $atoms ->{"C5'"} );
    $c5pAngle     = angle($atoms ->{"O5'"}, $atoms ->{"C5'"}, $atoms ->{"C4'"} );
    $c4pAngle     = angle($atoms ->{"C5'"}, $atoms ->{"C4'"}, $atoms ->{"C3'"} );
    
    #sum up the angle deviations
    $totalDev += (($o5pAngle     - O5pANGLEIDEAL)/ANGLESTDDEV)**2;
    $totalDev += (($c5pAngle     - C5pANGLEIDEAL)/ANGLESTDDEV)**2;
    $totalDev += (($c4pAngle     - C4pANGLEIDEAL)/ANGLESTDDEV)**2;
    
    #calculate the movement of the phosphate
    my $curPhosDist  = calcDist($atoms->{"P"}, $moduleCurPhos);
    unless (IGNOREPHOSDEV) {
        $totalDev += ($curPhosDist / PHOSSTDDEV)**2;
    }
    
    #if desired, print out information about the end point of the minimization
    if ($verbose) {
        print "\tNum of func evals: $numEvals\n";
        
        if ($verbose >= 2) {
            print "\tBond lengths:\n";
            print "\t\tP-O5' bond:    $Po5pBond\t"       . (($Po5pBond       - PO5pIDEAL)  /BONDSTDDEV) . "\n";
            print "\t\tO5'-C5' bond:  $o5pc5pBond\t"     . (($o5pc5pBond     - O5pC5pIDEAL)/BONDSTDDEV) . "\n";
            print "\t\tC5'-C4' bond:  $c5pc4pBond\t"     . (($c5pc4pBond     - C5pC4pIDEAL)/BONDSTDDEV) . "\n";
            
            print "\tBond angles\n";
            print "\t\tO5' angle:   $o5pAngle\t"     . (($o5pAngle     - O5pANGLEIDEAL)/ANGLESTDDEV) . "\n";
            print "\t\tC5' angle:   $c5pAngle\t"     . (($c5pAngle     - C5pANGLEIDEAL)/ANGLESTDDEV) . "\n";
            print "\t\tC4' angle:   $c4pAngle\t"     . (($c4pAngle     - C4pANGLEIDEAL)/ANGLESTDDEV) . "\n";
            
            print "\tPhosphates:\n";
            print "\t\tcurPhosDev:     $curPhosDist\t"        . ($curPhosDist  / PHOSSTDDEV)     . "\n";
        }
            
        print "\tmin value: $totalDev\n\n";
        print "\n" if $verbose >= 2;
    }
    
    $numEvals++;
    $numEvals = 0 if $resetEvalCount;
    
    return $totalDev;
}

sub evalInitAtomsCaller {
    #call evalInitAtoms with the module variable values for $atoms
    #(this allows Math::Amoeba to evaluate evalSuite)
    
    return evalInitAtoms($moduleCurResAtoms, @ARG);
}

sub modifyInitAtoms {
    #calculate coordinates given the minimization values
    #ARGUMENTS:
    #   $atoms                  - a hash containing atomic coordinates
    #   $Px, $Py, $Pz           - x,y,z coordinates for how much to move the phosphate atom
    #   $O5x, $O5y, $O5z        - x,y,z coordinates for how much to move the O5' atom
    #   $C5x, $C5y, $C5z        - x,y,z coordinates for how much to move the C5' atom
    #RETURNS:
    #   $toms                  - a hash containing the modified atomic coordinates
    
    my $atoms  = shift;
    my ($Px, $Py, $Pz)    = (shift, shift, shift);
    my ($O5x, $O5y, $O5z) = (shift, shift, shift);
    my ($C5x, $C5y, $C5z) = (shift, shift, shift);
    
    #do a shallow copy of $atoms so that we don't modify the input variable (i.e. so that there aren't side-effects of this function)
    $atoms = { (%{$atoms}) };
    
    #modify the coordinates
    $atoms ->{"P"}   = plus($atoms ->{"P"},   [$Px, $Py, $Pz]);
    $atoms ->{"O5'"} = plus($atoms ->{"O5'"}, [$O5x, $O5y, $O5z]);
    $atoms ->{"C5'"} = plus($atoms ->{"C5'"}, [$C5x, $C5y, $C5z]);
    
    return($atoms);
}

1;
