package Microarray::Reporter;

use 5.006;
use strict;
use warnings;
our $VERSION = '1.2';

require Microarray;

# an array_reporter contains a number of spot objects
# the reporter objects are identified by a biologically relevant id
# the reporter summarises the averaged spot data
{ package array_reporter;

	sub new {
		my $class = shift;
		my $self = { _id => shift, _spots => [] };
		bless $self, $class;
		return $self;
	}
	sub get_reporter_ratios {
		my $self = shift;
		my $hAnalysed_Data = { };
		$hAnalysed_Data->{ M_mean_of_ratios } = $self->mean_log_ratios;
		$hAnalysed_Data->{ M_ratio_of_means } = $self->log_ratio_means;
		$hAnalysed_Data->{ ch1_mean } = $self->mean_ch1;
		$hAnalysed_Data->{ ch2_mean } = $self->mean_ch2;
		return $hAnalysed_Data;
	}	
	# genetic_data() will set some relevant value(s)
	# usually, a sub-class such as bac_reporter or gene_reporter 
	# will have its own methods for getting and setting the genetic data
	sub genetic_data {
		my $self = shift;
		@_	?	$self->{ _genetic_data } = shift
			:	$self->{ _genetic_data };
	}
	# reporter_id() will be an alias for the genetic reporter id
	sub reporter_id {
		my $self = shift;
		$self->{ _id };
	}
	sub add_reporter_spot {
		my $self = shift;
		my $aSpots = $self->get_reporter_spots;
		push (@$aSpots, shift);	
	}
	sub get_reporter_spots {
		my $self = shift;
		$self->{ _spots };
	}
	sub get_reporter_replicates {
		my $self = shift;
		my $aSpots = $self->get_reporter_spots;
		if (@$aSpots){
			my $replicates = @$aSpots;
			return $replicates;
		} else {
			return;
		}
	}
	sub good_spot {
		my $self = shift;
		$self->{ _spots_passed_qc }++;
	}
	sub spots_passed_qc {
		my $self = shift;
		$self->{ _spots_passed_qc };
	}
	sub all_ch1 {
		my $self = shift;
		unless (defined $self->{ _all_ch1 }){
			$self->{ _all_ch1 } = [];
		}
		if (@_){
			my $aCh1_Signals = $self->{ _all_ch1 };
			push (@$aCh1_Signals, shift);
		} else {
			$self->{ _all_ch1 };
		}
	}
	sub all_ch2 {
		my $self = shift;
		unless (defined $self->{ _all_ch2 }){
			$self->{ _all_ch2 } = [];
		}
		if (@_){
			my $aCh2_Signals = $self->{ _all_ch2 };
			push (@$aCh2_Signals, shift);
		} else {
			$self->{ _all_ch2 };
		}
	}
	sub all_ratios {
		my $self = shift;
		unless (defined $self->{ _all_ratios }){
			$self->{ _all_ratios } = [];
		}
		if (@_){
			my $aRatios = $self->{ _all_ratios };
			push (@$aRatios, shift);
		} else {
			$self->{ _all_ratios };
		}
	}
	sub mean_ch1 {
		my $self = shift;
		stats_or_return($self->all_ch1);
	}
	sub mean_ch2 {
		my $self = shift;
		stats_or_return($self->all_ch2);
	}
	sub mean_ratios {
		my $self = shift;
		stats_or_return($self->all_ratios);
	}
	sub stats_or_return {
		my $aValues = shift;
		if (@$aValues == 1){
			return $$aValues[0];
		} else {
			return mean_values($aValues);
		}
	}
	sub ratio_means {
		my $self = shift;
		if ($self->mean_ch2 > 0){
			return ($self->mean_ch1) / ($self->mean_ch2);
		} else {
			return;
		}
	}
	sub log_ratio_means {
		my $self = shift;
		return calculate_log2($self->ratio_means);
	}
	sub mean_log_ratios {
		my $self = shift;
		my $aLog_Ratios = calculate_log2($self->all_ratios);
		return mean_values($aLog_Ratios);
	}
	sub mean_values {
		use Statistics::Descriptive;
		my $aValues = shift;
		my $data = Statistics::Descriptive::Full->new();
		$data->add_data(@$aValues);
		return $data->mean();
	}
	sub calculate_log2 {
		my $value = shift;
		if (ref $value){
			for my $val (@$value){
				next if ($val<=0);
				$val = log ($val) / log(2);
			}
			return $value;
		} else {
			return if ($value<=0);
			return log ($value) / log(2);
		}
	}	
}

1;


__END__

=head1 NAME

Microarray::Reporter - A Perl module for creating and manipulating microarray reporter objects

=head1 SYNOPSIS

	use Microarray;

	my $reporter = array_reporter->new('reporter 1');
	$reporter->add_reporter_spot($spot);

=head1 DESCRIPTION

Microarray::Reporter is an object-oriented Perl module for creating and manipulating microarray reporter objects. It serves as a container into which you place spot objects that are replicates of the same genetic reporter, and returns average information about those spots. 

=head1 METHODS

=over

=item B<reporter_id>

Name of the reporter

=item B<genetic_data>

An object containing relevant genetic data. 

=item B<get_reporter_spots>

Returns a list of spot objects attributed to a reporter

=item B<get_reporter_replicates>

Returns the number of spots attributed to a reporter

=item B<spots_passed_qc>

Returns the number of spots that passed QC criteria and are included in the reporter data

=item B<mean_ch1> and B<mean_ch2>

Mean signal of all spots representing a reporter

=item B<mean_ratios> and B<mean_log_ratios>

Calculates the ratio (or log2 ratio) between the two signal channels for each replicate, and returns the mean of those values

=item B<ratio_means> and B<log_ratio_means>

Calculates the mean of the replicate signals for each channel, and returns their ratio (or log2 ratio)

=back

=head1 SEE ALSO

L<Microarray|Microarray>, L<Microarray::Spot|Microarray::Spot>

=head1 AUTHOR

Christopher Jones, Translational Research Laboratories, Institute for Women's Health, University College London.

L<http://www.instituteforwomenshealth.ucl.ac.uk/trl>

c.jones@ucl.ac.uk

=head1 COPYRIGHT AND LICENSE

Copyright 2006 by Christopher Jones, University College London

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 

=cut
