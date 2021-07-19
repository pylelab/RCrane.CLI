#!/usr/bin/perl

# Predict rotamers from low resolution information
# Note that any of all of the prediction methods may be used
# Predictions will be made with all available information
# Example of common use:
#   my $pseudoPredic = new PseudoPredictor( THETAETA      => "thetaEtaClusts.csv",
#                                           PUCKER        => "smoothedPuckerDist.csv",
#                                           SUGARDIST     => "sugarDists.csv",
#                                           STARTPHOSDIST => "endingPDists.csv",
#                                           ENDPHOSDIST   => "endingPDists.csv" );
#
#   my $prob  = $pseudoPredic->calcProb( THETA          => 180,
#                                        ETA            => 180,
#                                        STARTPPERP     => 4.7,
#                                        ENDPPERP       => 4.7,
#                                        SUGARDIST      => 6,
#                                        STARTPHOSDIST  => 6,
#                                        ENDPHOSDIST    => 6
#                                        ROTAMER        => '&a');
#   print "Probability of rotamer &a: $prob\n;
#
#   my $probs = $pseudoPredic->calcProb( THETA          => 180,
#                                        ETA            => 180,
#                                        STARTPPERP     => 4.7,
#                                        ENDPPERP       => 4.7,
#                                        SUGARDIST      => 6,
#                                        STARTPHOSDIST  => 6,
#                                        ENDPHOSDIST    => 6);
#   print "Probability for all rotamers: " . Dumper($probs);

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

package PseudoPredictor;

use strict;
use English '-no_match_vars';
use SmoothProb;
use ThetaEtaPredictor;
use ThetaEtaOnlyPredictor;


sub new {
    #create a PseudoPredictor object
    #ARGUMENTS (AS HASH): one or more of the following must be given
    #   THETAETA      => A file contianing training data for theta eta Gaussian clusters
    #   PUCKER        => A file contianing training data for P-perp distances
    #   SUGARDIST     => A file contianing training data for sugar-sugar (C1'-C1') distances
    #   STARTPHOSDIST => A file contianing training data for starting phosphate-phosphate distances
    #   ENDPHOSDIST   => A file contianing training data for ending phosphate-phosphate distances
    #   THETAONLY     => A file contianing training data for theta-only Gaussian clusters
    #   ETAONLY       => A file contianing training data for theta-only Gaussian clusters
    #RETURNS:
    #   a PseudoPredictor object initialized with all provided input data
    
    my $class = shift;
    my $inputFiles = {@ARG};
    
    my $self = { inputFiles          => $inputFiles,
                 thetaEtaPredic      => undef,
                 puckerPredic        => undef,
                 sugarDistPredic     => undef,
                 startPhosDistPredic => undef,
                 endPhosDistPredic   => undef,
                 thetaOnlyPredic     => undef,
                 etaOnlyPredic       => undef   };
    bless $self, $class;
    
    #read in training data in the input files
    $self->init();
    
    return $self;
}

sub init {
    #initialize the predictor
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   None
    #SIDE EFFECTS:
    #   reads in all data from $self->{inputFiles} and generates the appropriate predictors
    
    my $self = shift;
    
    if (defined $self->{inputFiles}->{"THETAETA"}) {
        $self->{thetaEtaPredic} = new ThetaEtaPredictor( $self->{inputFiles}->{"THETAETA"} );
    }
    
    if (defined $self->{inputFiles}->{"PUCKER"}) {
        $self->{puckerPredic} = new SmoothProb( $self->{inputFiles}->{"PUCKER"} );
    }
    
    if (defined $self->{inputFiles}->{"SUGARDIST"}) {
        $self->{sugarDistPredic} = new SmoothProb( $self->{inputFiles}->{"SUGARDIST"} );
    }
    
    if (defined $self->{inputFiles}->{"STARTPHOSDIST"}) {
        $self->{startPhosDistPredic} = new SmoothProb( $self->{inputFiles}->{"STARTPHOSDIST"} );
    }
    
    if (defined $self->{inputFiles}->{"ENDPHOSDIST"}) {
        $self->{endPhosDistPredic} = new SmoothProb( $self->{inputFiles}->{"ENDPHOSDIST"} );
    }
    
    if (defined $self->{inputFiles}->{"THETAONLY"}) {
        $self->{thetaPredic} = new ThetaEtaOnlyPredictor( $self->{inputFiles}->{"THETAONLY"} );
    }
    
    if (defined $self->{inputFiles}->{"ETAONLY"}) {
        $self->{etaPredic} = new ThetaEtaOnlyPredictor( $self->{inputFiles}->{"ETAONLY"} );
    }
}

sub calcProb {
    #calculate the probability of each rotamer from all given information
    #ARGUMENTS (AS HASH): one or more of the following (other than ROTAMER) must be given
    #   THETAETA       => a ref to list of [theta, eta] values
    #   THETA          => the theta value
    #   ETA            => the eta value
    #                  Note: give *either* THETAETA or THETA and ETA
    #   STARTPPERP     => the P-perp value for the starting sugar
    #   ENDPPERP       => the P-perp value for the ending sugar
    #   SUGARDIST      => the sugar-sugar (C1'-C1' distance)
    #   STARTPHOSDIST  => the starting phosphate-phosphate distance
    #   ENDPHOSDIST    => the ending phosphate-phosphate distance
    #   ROTAMER        => the rotamer to calculate the probability of
    #RETURNS:
    #   if ROTAMER is not given, a ref to hash of {rotamer => probability}
    #   if ROTAMER is given, the probability that the rotamer is $rot
    
    my $self  = shift;
    my $inputs = {@ARG};
    
    #parse the input arguments
    my ($theta, $eta, $startPPerp, $endPPerp, $sugarDist, $startPhosDist, $endPhosDist, $rot);
    ($theta, $eta) = @{$inputs->{THETAETA}} if (defined $inputs->{THETAETA});
    $theta         = $inputs->{THETA}       if (defined $inputs->{THETA});
    $eta           = $inputs->{ETA}         if (defined $inputs->{ETA});
    $startPPerp    = $inputs->{STARTPPERP};
    $endPPerp      = $inputs->{ENDPPERP};
    $sugarDist     = $inputs->{SUGARDIST};
    $startPhosDist = $inputs->{STARTPHOSDIST};
    $endPhosDist   = $inputs->{ENDPHOSDIST};
    $rot           = $inputs->{ROTAMER};
    
    #calculate the all probabilites for the given information
    my ($thetaEtaProbs, $startPuckerProbs, $endPuckerProbs, $sugarDistProbs, $startPhosDistProbs, $endPhosDistProbs);
    if (defined $theta and defined $eta) {
        $thetaEtaProbs  = $self->{thetaEtaPredic}     ->calcProb($theta, $eta)   if (defined $self->{thetaEtaPredic});
    } elsif (defined $theta and defined $self->{thetaOnlyPredic}) {
        $thetaEtaProbs  = $self->{thetaOnlyPredic}    ->calcProb($theta);
    } elsif (defined $eta and defined $self->{etaOnlyPredic}) {
        $thetaEtaProbs  = $self->{etaOnlyPredic}      ->calcProb($eta);
    }
    
    $startPuckerProbs   = $self->{puckerPredic}       ->calcProb($startPPerp)    if (defined $self->{puckerPredic}        and defined $startPPerp);
    $endPuckerProbs     = $self->{puckerPredic}       ->calcProb($endPPerp)      if (defined $self->{puckerPredic}        and defined $endPPerp);
    $sugarDistProbs     = $self->{sugarDistPredic}    ->calcProb($sugarDist)     if (defined $self->{sugarDistPredic}     and defined $sugarDist);
    $startPhosDistProbs = $self->{startPhosDistPredic}->calcProb($startPhosDist) if (defined $self->{startPhosDistPredic} and defined $startPhosDist);
    $endPhosDistProbs   = $self->{endPhosDistPredic}  ->calcProb($endPhosDist)   if (defined $self->{endPhosDistPredic}   and defined $endPhosDist);
    
    my $probs = {};
    my $totalprob = 0;
    
    #calculate the probability for each rotamer
    foreach my $currot ($self->getRotList()) {
        my $curprob = 1;
        $curprob *= $thetaEtaProbs     ->{$currot} if defined $thetaEtaProbs;
        $curprob *= $sugarDistProbs    ->{$currot} if defined $sugarDistProbs;
        $curprob *= $startPhosDistProbs->{$currot} if defined $startPhosDistProbs;
        $curprob *= $endPhosDistProbs  ->{$currot} if defined $endPhosDistProbs;
        
        $curprob *= $startPuckerProbs  ->{$self->getStartPucker($currot)} if defined $startPuckerProbs;
        $curprob *= $endPuckerProbs    ->{$self->getEndPucker($currot)}   if defined $endPuckerProbs;
        
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
    
    return $self->{thetaEtaPredic}->getStartPucker($rot);
}

sub getEndPucker {
    #get the ending pucker associated with a given rotamer
    #ARGUMENTS:
    #   $rot    - the rotamer
    #RETURNS:
    #   the ending pucker for $rot
    
    my $self = shift;
    my $rot  = shift;
    
    return $self->{thetaEtaPredic}->getEndPucker($rot);
}

sub getPuckers {
    #get the starting and ending pucker associated with a given rotamer
    #ARGUMENTS:
    #   $rot    - the rotamer
    #RETURNS:
    #   a list of starting and ending pucker for $rot
    
    my $self = shift;
    my $rot  = shift;
    
    my $thetaEtaPredic = $self->{thetaEtaPredic};
    
    return ($thetaEtaPredic->getStartPucker($rot), $thetaEtaPredic->getEndPucker($rot));
}

sub getRotList {
    #get the list of all rotamers
    #ARGUMENTS:
    #   None
    #RETURNS:
    #   a list (*not* a ref to list) of all rotamers in the input file
    
    my $self = shift;
    return $self->{thetaEtaPredic}->getRotList();
}

1;