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
#use Math::Interpolate 'linear_interpolate';#, 'constant_interpolate';
#use Math::IntervalSearch qw(interval_search);
use Carp;

sub cluck { warn Carp::longmess @_ }

sub LessThan {
  $_[0] < $_[1];
}

sub LessThanEqualTo {
  $_[0] <= $_[1];
}

# This holds the result from the last interval search.
my $last_interval_result = undef;

sub interval_search {
  if ( @_ > 4 ) {
    cluck "interval called with too many parameters";
    return;
  }

  # Get the input arguments.
  my $x           = shift;
  my $sequenceRef = shift;

  return unless defined($x);
  return unless defined($sequenceRef);
  return unless ref($sequenceRef);

  # Check the input arguments for any code references and use them.
  my $LessThan        = \&LessThan;
  my $LessThanEqualTo = \&LessThanEqualTo;
  @_ and defined(ref($_[0])) and ref($_[0]) eq 'CODE' and
    $LessThan = shift;
  @_ and defined(ref($_[0])) and ref($_[0]) eq 'CODE' and
    $LessThanEqualTo = shift;

  # Get the number of points in the data.
  my $num = @$sequenceRef;

  # Return -1 if there's no data.
  if ( $num <= 0 ) {
    $last_interval_result = 0;
    return -1;
  }

  # Use the result from the last time through the subroutine, if it
  # exists.  Force the result into the range required by the array
  # size.
  $last_interval_result = 0 unless defined($last_interval_result);

  # Which side of the data point is x on if there's only one point?
  if ( $num == 1 ) {
    $last_interval_result = 0;
    if ( &$LessThan($x, $sequenceRef->[0]) ) {
      return -1;
    }
    else {
      return 0;
    }
  }

  # Is the point less than the smallest point in the sequence?
  if ( &$LessThan($x, $sequenceRef->[0]) ) {
    $last_interval_result = 0;
    return -1;
  }

  # Is the point greater than the largest point in the sequence?
  if ( &$LessThanEqualTo($sequenceRef->[$num-1], $x) ) {
    return $last_interval_result = $num - 1;
  }

  # Use the result from the last run as a start for this run.
  if ( $last_interval_result > $num-1 ) {
    $last_interval_result = $num - 2;
  }
  my $ilo = $last_interval_result;
  my $ihi = $ilo + 1;

  # Is the new upper ihi beyond the extent of the sequence?
  if ( $ihi >= $num ) {
    $ihi = $num - 1;
    $ilo = $ihi - 1;
  }

  # If x < sequence(ilo), then decrease ilo to capture x.
  if ( &$LessThan($x, $sequenceRef->[$ilo]) ) {
    my $istep = 1;
    for (;;) {
      $ihi = $ilo;
      $ilo = $ihi - $istep;
      if ( $ilo <= 0 ) {
	$ilo = 0;
	last;
      }
      if ( &$LessThanEqualTo($sequenceRef->[$ilo], $x) ) {
	last;
      }
      $istep *= 2;
    }
  }

  # If x >= sequence(ihi), then increase ihi to capture x.
  if ( &$LessThanEqualTo($sequenceRef->[$ihi], $x) ) {
    my $istep = 1;
    for (;;) {
      $ilo = $ihi;
      $ihi = $ilo + $istep;
      if ( $ihi >= $num-1 ) {
	$ihi = $num - 1;
	last;
      }
      if ( &$LessThan($x, $sequenceRef->[$ihi]) ) {
	last;
      }
      $istep *= 2;
    }
  }

  # Now sequence(ilo) <= x < sequence(ihi).  Narrow the interval.
  for (;;) {
    # Find the middle point of the sequence.
    my $middle = int(($ilo + $ihi)/2);

    # The division above was integer, so if ihi = ilo+1, then
    # middle=ilo, which tests if x has been trapped.
    if ( $middle == $ilo ) {
      $last_interval_result = $ilo;
      return $ilo;
    }
    if ( &$LessThan($x, $sequenceRef->[$middle]) ) {
      $ihi = $middle;
    }
    else {
      $ilo = $middle;
    }
  }
}

sub linear_interpolate {
  my $x = shift;
  return unless defined($x);

  my $X = shift;
  return unless defined($X);
  return unless ref($X);

  my $Y = shift;
  return unless defined($Y);
  return unless ref($Y);

  my $num_x = @$X;
  my $num_y = @$Y;
  return unless $num_x == $num_y;

  # Find where the point to be interpolated lies in the input sequence.
  # If the point lies outside, then coerce the index value to be legal for
  # the routine to work.  Remember, this is only an interpreter, not an
  # extrapolator.
  my $j = interval_search($x, $X);
  if ( $j < 0 ) {
    $j = 0;
  }
  elsif ( $j >= $num_x - 1 ) {
    $j = $num_x - 2;
  }
  my $k = $j + 1;

  # Calculate the linear slope between the two points.
  my $dy = ($Y->[$k] - $Y->[$j]) / ($X->[$k] - $X->[$j]);

  # Use the straight line between the two points to interpolate.
  my $y  = $dy*($x - $X->[$j]) + $Y->[$j];

  return wantarray ? ($y, $dy) : $y;
}

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
