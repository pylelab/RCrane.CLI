#!/usr/bin/perl
# Predict rotamers from (theta, eta) values
# Example of common use:
#   my $thetaEtaPredic = new ThetaEtaPredictor("thetaEtaClusts.csv");
#   
#   my $prob = $thetaEtaPredic->calcProb(164, 109, "1o") . "\n";
#   print "Probability of rotamer 1o: $prob\n";
#
#   my $probs = $thetaEtaPredic->calcProb(164, 109);
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

package ThetaEtaPredictor;

use strict;
use English '-no_match_vars';
use Math::Trig; #defines pi constant and deg2rad
use POSIX 'fmod';

use constant { STARTPUCKER => 0,
               ENDPUCKER   => 1,
               CENTER      => 2,
               VARIANCE    => 3,
               ROTATION    => 4  };

sub new {
    #create a ThetaEtaPredictor object
    #ARGUMENTS:
    #   $inputFile - a file containing training data for theta-eta Gaussian clusters
    #RETURNS:
    #   a ThetaEtaPredictor object initialized with the cluster data from $inputFile
    
    my $class = shift;
    my $inputFile = shift;
    
    my $self = { clusters => undef };
    bless $self, $class;
    
    $self->init($inputFile);
    
    return $self;
    
}

sub init {
    #initialize the predictor
    #ARGUMENTS:
    #   $inputFile  - a file containing training data for theta-eta Gaussian clusters
    #RETURNS:
    #   none
    #SIDE EFFECTS:
    #   reads in all Gaussian cluster data from $inputFile and stores it in $self->{clusters}
    
    my $self = shift;
    my $inputFile = shift;
    
    open(CLUSTS, $inputFile) or die "Could not open $inputFile for reading\n";
    <CLUSTS>; #ignore the header
    
    my $clusters = {}; #the data about all of the clusters
    while (my $curline = <CLUSTS>) {
        my @data = split(",", $curline);
        
        my $rotName     = $data[0];
        my $startPucker = $data[1];
        my $endPucker   = $data[2];
        my $center      = [@data[3..4]];
        my $variance    = [@data[5..6]];
        my $rotation    = $data[7] or 0;
        
        $clusters->{$rotName} = [$startPucker, $endPucker, $center, $variance, $rotation];
    }
    close(CLUSTS);
    
    $self->{clusters} = $clusters;
}

sub calcProb {
    #calculate the probability of each rotamer for a given (theta, eta) value
    #ARGUMENTS:
    #   $theta  - the theta value of the point
    #   $eta    - the eta value of the point
    #OPTIONAL ARGUMENTS
    #   $rot    - a specific rotamer
    #RETURNS:
    #   if $rot is not given, a ref to hash of {rotamer => probability}
    #   if $rot is given, the probability that the rotamer is $rot
    
    my $self  = shift;
    my $theta = shift;
    my $eta   = shift;
    my $rot   = shift; #optional
    
    my $probs = {};
    my $totalprob = 0;
    
    #calculate the probality for each rotamer
    foreach my $currot (keys(%{$self->{clusters}})) {
        my $curprob = $self->calcRawProb($theta, $eta, $currot);
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


sub getStartPucker {
    #get the starting pucker associated with a given rotamer
    #ARGUMENTS:
    #   $rot    - the rotamer
    #RETURNS:
    #   the starting pucker for $rot
    
    my $self = shift;
    my $rot  = shift;
    
    return $self->{clusters}->{$rot}->[STARTPUCKER];
}

sub getEndPucker {
    #get the ending pucker associated with a given rotamer
    #ARGUMENTS:
    #   $rot    - the rotamer
    #RETURNS:
    #   the ending pucker for $rot
    
    my $self = shift;
    my $rot  = shift;
    
    return $self->{clusters}->{$rot}->[ENDPUCKER];
}

sub getPuckers {
    #get the starting and ending pucker associated with a given rotamer
    #ARGUMENTS:
    #   $rot    - the rotamer
    #RETURNS:
    #   a list of starting and ending pucker for $rot
    
    my $self = shift;
    my $rot  = shift;
    
    my $rotData = $self->{clusters}->{$rot};
    
    return ($rotData->[STARTPUCKER], $rotData->[ENDPUCKER]);
}

sub getRotList {
    #get the list of all rotamers
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a list (*not* a ref to list) of all rotamers in the input file
    
    my $self = shift;
    return keys(%{$self->{clusters}});
}

sub calcRawProb {
    #calculate the responsibility of a given rotamer cluster for a given point based only on (theta, eta) values
    #ARGUMENTS:
    #   $theta  - the theta value of the point
    #   $eta    - the eta value of the point
    #   $rot    - the rotamer
    #RETURNS:
    #   the unscaled responsibility of $rot for ($theta, $eta)
    
    my $self  = shift;
    my $theta = shift;
    my $eta   = shift;
    my $rot   = shift;
    
    my $curcenter   = $self->{clusters}->{$rot}->[CENTER];
    my $curvariance = $self->{clusters}->{$rot}->[VARIANCE];
    my $curangle    = $self->{clusters}->{$rot}->[ROTATION];
    my $curpoint    = [$theta, $eta];
    
    #first, translate the point
    my $x = subtractCoords($curpoint->[0], $curcenter->[0]);
    my $y = subtractCoords($curpoint->[1], $curcenter->[1]);
    
    #then rotate the point by the opposite of the angle (I think)
    $curangle = deg2rad($curangle);
    $curpoint->[0] = $x * cos($curangle) - $y * sin($curangle);
    $curpoint->[1] = $x * sin($curangle) + $y * cos($curangle);
    
    #calculate the sum and product required by the formula
    my $product = 1;
    my $sum = 0;
    for (my $curaxis = 0; $curaxis < 2; $curaxis++) {
        $product *= (sqrt(2*pi) * sqrt($curvariance->[$curaxis]));
        $sum += ( ($curpoint->[$curaxis]**2) / (2*$curvariance->[$curaxis]));
    }
    
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