#!/usr/bin/perl

# Calculate the probabilities of various models for a given distance measure
# Must be initialized with smoothed data for each model
# Example of common use:
#   my $puckerDistPredic = new smoothProb("smoothedPuckerDist.csv");
#   print "Probality of C3' pucker: " . $smoother->calcProb(3.7, "C3") . "\n";
#   print "Probability of all puckers: " . Dumper($smoother->calcProb(3.7)) . "\n";

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

package SmoothProb;

use strict;
use English '-no_match_vars';
use Strip;
use Math::Interpolate 'linear_interpolate';#, 'constant_interpolate';


sub new {
    #create a SmoothProb object
    #OPTIONAL ARGUMENTS:
    #   $inputFile - a file containing smoothed training data
    #   $delim - the field delimeter used in $inputfile
    #            comma delimeted is assumed if $delim is not given
    #RETURNS:
    #   a SmoothProb object
    #   if $inputFile is given, the SmoothProb object will have been initialized with the data from $inputFile
    
    my $class = shift;
    my $inputFile = shift;
    my $delim = shift;
    
    my $self = { dist => [],
                 probs => {}   };
    bless $self, $class;
    
    if (defined($inputFile)) {
        $self->init($inputFile, $delim);
    }
    
    return $self;
}

sub init {
    #initialize a SmoothProb object
    #ARGUMENTS:
    #   $inputFile - a file containing smoothed training data
    #OPTIONAL ARGUMENTS:
    #   $delim - the field delimeter used in $inputfile
    #            comma delimeted is assumed if $delim is not given
    #RETURNS:
    #   None
    #SIDE EFFECTS:
    #   stores the distances and probabilites from $inputFile in $self->{dist} and $self-{probs}
    
    my $self = shift;
    my $inputFile = shift;
    my $delim = (shift or ",");
    
    open (IN, $inputFile) or die "Could not open $inputFile for reading\n";
    
    #read in the header line, which contains the labels for our probabilities
    my @labels = split($delim, <IN>);
    shift(@labels); #we don't care about the header for the first column, since that's just distance
    @labels = map(strip($ARG), @labels); #get rid of any leading and trailing spaces
    $self->{probs} = {map {$ARG=>[]} @labels};
    
    #read in each line and store all of the probabilities
    while (my $curline = <IN>) {
        chomp($curline);
        my @data = split($delim, $curline);
        push(@{$self->{dist}}, shift(@data));
        
        foreach my $curlabel (@labels) {
            push(@{$self->{probs}->{$curlabel}}, shift(@data));
        }
    }
    
    close(IN);
}

sub calcProb {
    #calculate the probabilities for a given distance
    #ARGUMENTS:
    #   $dist   - the distance
    #OPTIONAL ARGUMENTS:
    #   $label  - what model to calculate the probability for
    #             must match a label in the first line of the input file
    #RETURNS:
    #   if $label is given, the probability that the model $label is correct
    #   if $label is not given, a ref to hash containing probabilities for all models
    
    my $self = shift;
    my $dist = shift;
    my $label = shift; #optional
    
    #go through each label and determine a raw probability (from interpolation)
    my $probs = {};
    my $totalprob = 0;
    foreach my $curlabel (keys(%{$self->{probs}})) {
        my $curprob = linear_interpolate($dist, $self->{dist}, $self->{probs}->{$curlabel});
        $probs->{$curlabel} = $curprob;
        $totalprob += $curprob;
    }
    
    #scale the raw probabilities so they add up to 1
    foreach my $curprob (values(%{$probs})) {
        $curprob /= $totalprob;
    }
    
    #return either the desired probability or the hash containing all probabilities
    if (defined $label) {
        return $probs->{$label};
    } else {
        return $probs;
    }
}


1;