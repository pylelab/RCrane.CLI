#!/usr/bin/perl

# use a Hidden Markov Model to determine the most likely series of rotamers

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

package RotamerHMM;

use strict;
use English '-no_match_vars';
use Exporter;

use PuckerList qw($puckerList $rotList);

our @ISA       = ('Exporter');
our @EXPORT    = qw(rotamerHMM);

use constant NEGINF => -999999999;  #a stand-in for negative infinity
    #this value used for the log of 0

#this module could be sped up by preprocessing the predicted probabilities
#to remove everything but the best scoring rotamer for each starting and ending pucker combination
#(i.e. there would then only be four possible rotamers for each position)
#That probably won't save all that much time, and this way I can add in real transition probabilities if I wanted

sub rotamerHMM {
    #given rotamer likelihoods for all suites, predict the most likely rotamer string using a Hidden Markov Model
    #ARGUMENTS:
    #   $predictedProbs - a list in the form containing [suite, rotamer probabilities] pairs, where suite is an object
    #                     defined by Suite.pm and rotamer probabilities is a hash ref of {rotamer}->likelihood
    #                     (as is created by PseudoPredictor.pm)
    #RETURNS:
    #   $bestPath       - a array ref of rotamers indicating the best rotamer for each suite
    
    my $predictedProbs  = shift;
    
    my $pathProbs = []; #the log probability of having followed a given path (the delta or phi array)
    my $path      = []; #the path followed for $pathProbs (the psi array)
    
    #initialize the $pathProbs array
    foreach my $curRot (@{$rotList}) {
        $pathProbs->[0]->{$curRot} = ln($predictedProbs->[0]->[1]->{$curRot});
    }
    
    my $curPos = 1;
    #go through each suite and determine the most likely path
    while (defined $predictedProbs->[$curPos]) {
        
        #check to see if this suite is connected to the previous one
        my $connected = $predictedProbs->[$curPos]->[0]->connectedToPrev;
        
        foreach my $curRot (@{$rotList}) {
            
            #figure out what the best previous rotamer is for ending up at the current rotamer
            my $bestPrevRot;
            my $bestScore = NEGINF;
            foreach my $prevRot (@{$rotList}) {
                
                my $curScore = $pathProbs->[$curPos - 1]->{$prevRot} + transitionProb($prevRot, $curRot, $connected);
                if ($curScore > $bestScore) {
                    $bestScore   = $curScore;
                    $bestPrevRot = $prevRot;
                }
            }
            
            $pathProbs->[$curPos]->{$curRot} = $bestScore + ln($predictedProbs->[$curPos]->[1]->{$curRot});
            $path     ->[$curPos]->{$curRot} = $bestPrevRot;
        }
        
        $curPos++;
    }
    
    #figure out what the best path was
    $curPos--; #the last position that was defined
    my $bestPath = [];
    
    #figure out the best last position
    my $bestLastRot;
    my $bestScore = NEGINF;
    foreach my $curRot (@{$rotList}) {
        if ($pathProbs->[$curPos]->{$curRot} > $bestScore) {
            $bestScore   = $pathProbs->[$curPos]->{$curRot};
            $bestLastRot = $curRot;
        }
    }
    $bestPath->[$curPos] = $bestLastRot;
    
    #follow the path back
    while ($curPos > 0) {
        $bestPath->[$curPos-1] = $path->[$curPos]->{$bestPath->[$curPos]};
        $curPos--;
    }
    
    return $bestPath;
}

sub transitionProb {
    #retrieve the natural log of transition probability between two rotamers (which is either ln(1) or ln(0))
    #ARGUMENTS:
    #   $startingRot - the rotamer we are transitioning from
    #   $endingRot   - the rotamer we are transitioning to
    #OPTIONAL ARGUMENTS:
    #   $connected   - whether this suite is connected to the previous one
    #                  if not defined, assumed to be true
    #RETURNS:
    #   0 (which is ln(1)) if the two rotamers have compatible puckers
    #   NEGING (which is ln(0)) if they don't
    my $startingRot  = shift;
    my $endingRot    = shift;
    my $connected    = shift;
    
    $connected = 1 if not defined $connected;
    
    if (not $connected) {
        return 0; #ln(1)
    } elsif ($puckerList->{$startingRot}->[1] eq $puckerList->{$endingRot}->[0]) {
        return 0; #ln(1)
    } else {
        return NEGINF; #ln(0)
    }
}

sub ln {
    #calculate the natural log, and return NEGINF for ln(0)
    #ARGUMENTS
    #   $a  - the value to calculate the natural log of
    #RETURNS:
    #   log($a), or NEGINF if $a is 0
    
    my $a = shift;
    if ($a == 0) {
        return NEGINF;
    } else {
        return log($a);
    }
}