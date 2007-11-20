package Microarray::File::Data::Agilent;

use 5.006;
use strict;
use warnings;
our $VERSION = '0.18';


{ package agilent_file;

	require Microarray::File::Data;

	our @ISA = qw( data_file );

	# setter for { _spot_data }, { _data_fields } and { _header_info }
	sub sort_data {
		my $self = shift;
		my $aaData = shift;
		$self->set_header_info($aaData);
		$self->set_spot_objects($aaData);	# all the numbers
	}
	# information about the scan
	sub set_header_info {
		my $self = shift;
		my $aaData = shift;
		my $hHeader_Info = { };
		my $hFeparams_Info = { };
		my $hStats_Info = { };
		my $which_info;
		ROW:while (my $aRow = shift @$aaData){
			next ROW if (($$aRow[0] eq 'TYPE')||($$aRow[0] eq '*'));
			if ($$aRow[0] eq 'FEATURES'){
				$self->set_data_fields($aRow);
				last ROW;
			} elsif ($$aRow[0] eq 'FEPARAMS'){
				for (my $i=1; $i<@$aRow; $i++){
					$hFeparams_Info->{ $i } = $$aRow[$i];
				}
				$which_info = 'feparams';
				next ROW;
			} elsif ($$aRow[0] eq 'STATS'){
				for (my $i=1; $i<@$aRow; $i++){
					$hStats_Info->{ $i } = $$aRow[$i];
				}
				$which_info = 'stats';
				next ROW;
			} elsif ($$aRow[0] eq 'DATA'){
				if ($which_info eq 'stats'){
					for (my $i=1; $i<@$aRow; $i++){
						$hHeader_Info->{ $hStats_Info->{ $i } } = $$aRow[$i];
					}
				} else {
					for (my $i=1; $i<@$aRow; $i++){
						$hHeader_Info->{ $hFeparams_Info->{ $i } } = $$aRow[$i];
					}
				}
			}
		}
		$self->{ _header_info } = $hHeader_Info;
	}
	sub data_file_fields {	# minimum required fields
		[	'spot_index','feature_id','synonym_id',
			'channel1_signal','channel2_signal',
			'channel1_snr','channel2_snr',
			'channel1_quality','channel2_quality',
			'channel1_sat','channel2_sat',
			'spot_pixels','flag_id',
			'spot_row','spot_col','x_pos','y_pos','bg_pixels', 
			'ch1_mean_f','ch1_median_f','ch1_sd_f','ch1_mean_b',
			'ch1_median_b','ch1_sd_b','ch2_mean_f',
			'ch2_median_f','ch2_sd_f','ch2_mean_b','ch2_median_b',
			'ch2_sd_b',	];
	}
	sub pixel_size {
		my $self = shift;
		unless ($self->pixelx_size == $self->pixely_size){
			# worry about this if it ever becomes an issue!
		}
		$self->pixelx_size;
	}
	sub pixelx_size {
		my $self = shift;
		$self->get_header_info('Scan_MicronsPerPixelX');
	}
	sub pixely_size {
		my $self = shift;
		$self->get_header_info('Scan_MicronsPerPixelY');
	}
	sub channel1_name {
		'r';
	}
	sub channel2_name {
		'g';
	}
	sub channel_name {
		my $self = shift;
		my $ch = shift;
		my $method = 'channel'.$ch.'_name';
		$self->$method;
	}
	sub num_channels {
		my $self = shift;
		$self->get_header_info('Scan_NumChannels');
	}
	sub slide_barcode {
		my $self = shift;
		$self->get_header_info('FeatureExtractor_Barcode');
	}
	sub gal_file {
		my $self = shift;
		$self->get_header_info('FeatureExtractor_DesignFileName');
	}
	sub analysis_software {
		'Agilent';
	}
	sub scanner {
		my $self = shift;
		$self->get_header_info('Scan_ScannerName');
	}
	### data file fields ###	
	sub spot_index {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('FeatureNum') ];
	}
	sub spot_row {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('Row') ];
	}
	sub spot_col {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('Col') ];
	}
	sub x_pos {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('PositionX') ];
	}
	sub y_pos {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('PositionY') ];
	}
	sub feature_id {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('GeneName') ];
	}
	sub synonym_id {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('SystematicName') ];
	}
	sub ch1_spot_pixels {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('rNumPix') ];
	}
	sub ch2_spot_pixels {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('gNumPix') ];
	}
	sub spot_pixels {
		my $self = shift;
		my $aData = shift;
		my $ch1_pix = $self->ch1_spot_pixels($aData);
		my $ch2_pix = $self->ch2_spot_pixels($aData);
		if ($ch2_pix > $ch1_pix){
			return $ch2_pix;
		} else {
			return $ch1_pix;
		}
	}	
	sub ch1_bg_pixels {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('rBGNumPix') ];
	}
	sub ch2_bg_pixels {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('gBGNumPix') ];
	}
	sub bg_pixels {
		my $self = shift;
		my $aData = shift;
		my $ch1_pix = $self->ch1_bg_pixels($aData);
		my $ch2_pix = $self->ch2_bg_pixels($aData);
		if ($ch2_pix > $ch1_pix){
			return $ch2_pix;
		} else {
			return $ch1_pix;
		}
	}
	sub flag_id {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('IsManualFlag') ];
	}
	sub ch1_mean_f {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('rMeanSignal') ];
	}
	sub ch1_median_f {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('rMedianSignal') ];
	}
	sub ch1_sd_f {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('rPixSDev') ];
	}
	sub ch1_mean_b {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('rBGMeanSignal') ];
	}
	sub ch1_median_b {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('rBGMedianSignal') ];
	}
	sub ch1_sd_b {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('rBGPixSDev') ];
	}
	sub channel1_quality {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('rPValFeatEqBG') ];
	}
	sub channel1_sat {
		my $self = shift;
		my $aData = shift;
		my $num_pix = $self->ch1_spot_pixels($aData);
		my $num_sat_pix = $self->channel1_sat_pix($aData);
		return ($num_sat_pix/$num_pix) * 100;
	}
	sub channel1_sat_pix {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('rNumSatPix') ];
	}
	sub ch2_mean_f {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('gMeanSignal') ];
	}
	sub ch2_median_f {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('gMedianSignal') ];
	}
	sub ch2_sd_f {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('gPixSDev') ];
	}
	sub ch2_mean_b {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('gBGMeanSignal') ];
	}
	sub ch2_median_b {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('gBGMedianSignal') ];
	}
	sub ch2_sd_b {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('gBGPixSDev') ];
	}
	sub channel2_quality {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('gPValFeatEqBG') ];
	}
	sub channel2_sat {
		my $self = shift;
		my $aData = shift;
		my $num_pix = $self->ch2_spot_pixels($aData);
		my $num_sat_pix = $self->channel2_sat_pix($aData);
		return ($num_sat_pix/$num_pix) * 100;
	}
	sub channel2_sat_pix {
		my $self = shift;
		my $aData = shift;
		$$aData[ $self->get_column_id('gNumSatPix') ];
	}

	### calculated fields ###
	sub channel1_signal {
		my $self = shift;
		my $aData = shift;
#		$self->ch1_mean_f($aData) - $self->ch1_median_b($aData);
		$$aData[ $self->get_column_id('rProcessedSignal') ];
	}
	sub channel2_signal {
		my $self = shift;
		my $aData = shift;
#		$self->ch2_mean_f($aData) - $self->ch2_median_b($aData);
		$$aData[ $self->get_column_id('gProcessedSignal') ];
	}
	sub channel_signal {
		my $self = shift;
		my $aData = shift;
		my $ch = shift;
		my $method = 'channel'.$ch.'_signal';
		$self->$method($aData);
	}
	###Êquantarray calculates SNR as (median signal/background SD) ###
	sub channel1_snr {
		my $self = shift;
		my $aData = shift;
		my $median = $self->ch1_median_f($aData);
		my $sd = $self->ch1_sd_b($aData);
		if ($median && $sd) {
			return $median / $sd;
		} elsif ($median) {
			return $median;
		} else {
			return 0;
		}
#		$self->ch1_median_f($aData) / $self->ch1_sd_b($aData);
#		$$aData[ $self->get_column_id('Ch1 SignalNoiseRatio') ];
	}
	sub channel2_snr {
		my $self = shift;
		my $aData = shift;
		my $median = $self->ch2_median_f($aData);
		my $sd = $self->ch2_sd_b($aData);
		if ($median && $sd) {
			return $median / $sd;
		} elsif ($median) {
			return $median;
		} else {
			return 0;
		}
#		$self->ch2_median_f($aData) / $self->ch2_sd_b($aData);
#		$$aData[ $self->get_column_id('Ch2 SignalNoiseRatio') ];
	}
	sub channel_snr {
		my $self = shift;
		my $aData = shift;
		my $ch = shift;
		my $method = 'channel'.$ch.'_snr';
		$self->$method($aData);
	}
	sub channel_quality {
		my $self = shift;
		my $aData = shift;
		my $ch = shift;
		my $method = 'channel'.$ch.'_quality';
		$self->$method($aData);
	}
	sub channel_sat {
		my $self = shift;
		my $aData = shift;
		my $ch = shift;
		my $method = 'channel'.$ch.'_sat';
		$self->$method($aData);
	}
	sub bad_flags {
		{ '1'=>'1' };
	}
}

1;

__END__

=head1 NAME

Microarray::File::Data::Agilent - A Perl module for managing Agilent microarray files

=head1 SYNOPSIS

	use Microarray::File::Data::Agilent;

	my $data_file = agilent_file->new("/file.csv");

=head1 DESCRIPTION

Microarray::File::Data::Agilent is an object-oriented Perl module for managing microarray files created by Agilent Microarray Scanner software. It inherits from L<Microarray::File|Microarray::File>, and maps data fields in an Agilent data file to those used by L<Microarray::File::Data|Microarray::File::Data>.

=head1 METHODS

To be added

=head1 SEE ALSO

L<Microarray|Microarray>, L<Microarray::File|Microarray::File>, L<Microarray::File::Data|Microarray::File::Data>

=head1 AUTHOR

Christopher Jones, Translational Research Laboratories, Institute for Women's Health, University College London.

L<http://www.instituteforwomenshealth.ucl.ac.uk/trl>

c.jones@ucl.ac.uk

=head1 COPYRIGHT AND LICENSE

Copyright 2006 by Christopher Jones, University College London

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 

=cut
