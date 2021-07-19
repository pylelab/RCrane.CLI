#!/usr/bin/perl
# Predict rotamers from either a theta or eta value
# Example of common use:
#   my $thetaPredic = new ThetaEtaOnlyPredictor("thetaOnly.csv");
#   
#   my $prob = $thetaPredic->calcProb(164, "1o") . "\n";
#   print "Probability of rotamer 1o: $prob\n";
#
#   my $probs = $thetaEtaPredic->calcProb(164);
#   print "Probabilities for all rotamers are:\n";
#   print Dumper($probs);

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

package ThetaEtaOnlyPredictor;

use strict;
use English '-no_match_vars';
use Math::Trig; #defines pi constant and deg2rad
use POSIX 'fmod';


sub new {
    #create a ThetaEtaOnlyPredictor object
    #ARGUMENTS:
    #   $inputFile - a file containing training data for theta or eta Gaussian clusters
    #RETURNS:
    #   a ThetaEtaOnlyPredictor object initialized with the cluster data from $inputFile
    
    my $class = shift;
    my $inputFile = shift;
    
    my $self = { center   => {},
                 variance => {}, };
    bless $self, $class;
    
    $self->init($inputFile);
    
    return $self;
    
}

sub init {
    #initialize the predictor
    #ARGUMENTS:
    #   $inputFile  - a file containing training data for theta or eta Gaussian clusters
    #RETURNS:
    #   none
    #SIDE EFFECTS:
    #   reads in all Gaussian cluster data from $inputFile and stores it in $self
    
    my $self = shift;
    my $inputFile = shift;
    
    open(CLUSTS, $inputFile) or die "Could not open $inputFile for reading\n";
    <CLUSTS>; #ignore the header
    
    my $clusters = {}; #the data about all of the clusters
    while (my $curline = <CLUSTS>) {
        my @data = split(",", $curline);
        
        my $rotName     = $data[0];
        my $center      = $data[1];
        my $variance    = $data[2];
        
        $self->{center}->{$rotName} = $center;
        $self->{variance}->{$rotName} = $center;
    }
    close(CLUSTS);
}

sub calcProb {
    #calculate the probability of each rotamer for a given theta or eta value
    #ARGUMENTS:
    #   $val    - the theta or eta value of the point
    #OPTIONAL ARGUMENTS
    #   $rot    - a specific rotamer
    #RETURNS:
    #   if $rot is not given, a ref to hash of {rotamer => probability}
    #   if $rot is given, the probability that the rotamer is $rot
    
    my $self  = shift;
    my $val   = shift;
    my $rot   = shift; #optional
    
    my $probs = {};
    my $totalprob = 0;
    
    #calculate the probality for each rotamer
    foreach my $currot (keys(%{$self->{center}})) {
        my $curprob = $self->calcRawProb($val, $currot);
        $probs->{$currot} = $curprob;
        $totalprob += $curprob;
    }
    
    #scale the raw probabilities so they add up to 1
    foreach my $curprob (values(%{$probs})) {
        $curprob /= $totalprob;
    }
    
    #return either the desired probability or the hash containing all probabilities
    if (defined $rot) {
        return $probs->{$rot};
    } else {
        return $probs;
    }
}


sub calcRawProb {
    #calculate the responsibility of a given rotamer cluster for a given point based only on the theta or eta value
    #ARGUMENTS:
    #   $val    - the theta value of the point
    #   $rot    - the rotamer
    #RETURNS:
    #   the unscaled responsibility of $rot for theta or eta value $val
    
    my $self  = shift;
    my $val = shift;
    my $rot   = shift;
    
    my $curcenter   = $self->{center}->{$rot};
    my $curvariance = $self->{variance}->{$rot};
    
    my $dist = subtractCoords($val, $curcenter);
    
    #calculate the sum and product required by the formula
    my $product = sqrt(2*pi*$curvariance);;
    my $sum     = ($dist**2) / (2*$curvariance);
    
    return 1/$product * exp(-$sum);
}

sub subtractCoords {
    #calculate the minimum distance between two measures, keeping the appropriate sign
    #ARGUMENTS:
    #   $num1 - the first measurement
    #   $num2 - the second measurement
    #RETURNS:
    #   the distance between these two angles of no more than 180 degrees
    
    my $num1 = shift;
    my $num2 = shift;
    
    return fmodpos($num1 - $num2 + 180, 360) - 180;
}

sub fmodpos {
    #calculate the modulus in the range of 0 to N
    #ARGUMENTS:
    #   $num - the number to take the modulus of
    #   $mod - the modulus base
    #RETURNS:
    #   $num mod $mod
    
    my $num = shift;
    my $mod = shift;
    
    my $res = fmod($num, $mod);
    if ($res < 0) {
        $res += $mod;
    }
    
    return $res;
}


1;