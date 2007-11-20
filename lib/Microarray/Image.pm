package Microarray::Image;

use 5.006;
use strict;
use warnings;
our $VERSION = '0.112';

use GD::Image;
use Microarray::File;

{ package plot;
  
	sub new {
		my $class = shift;
		my $self  = { };
		if (@_){
			$self->{ _data_object } = shift;
			bless $self, $class;
			$self->set_data;
		} else {
			bless $self, $class;
		}
		return $self;
	}
	sub gd_object {
		my $self = shift;
		$self->{ _gd_object };
	}
	sub parse_args {
		my $self = shift;
		my %hArgs = @_;
		while(my($arg,$val) = each %hArgs){
			if ($self->can($arg)){
				$self->$arg($val);
			} else {
				die "Microrray::Image ERROR; No parameter '$arg' is defined\n";
			}
		}
	}
	sub set_data {
		my $self = shift;
		if (@_) {
			$self->{ _data_object } = shift;
		} 
		die "Microarray::Image ERROR: No data object provided\n" unless $self->data_object;
		$self->sort_data;
		## polymorphic method to process data such as MA RI log2 etc
		$self->process_data;
	}
	## from data object the background adjusted intensity of ch1 and ch2 and x,y coordinates 
	## are set in the plot image object
	sub sort_data {
		my $self   = shift;
		my $oData  = $self->data_object;
		my $spot_count = $oData->spot_count;
		my $aCh1   = [];
		my $aCh2   = [];
		my $aXcoords = [];
		my $aYcoords = [];		
		for (my $i=0; $i<$spot_count; $i++){
			my $ch1 = $oData->channel1_signal($i);
			my $ch2 = $oData->channel2_signal($i);
			my $x_pos = $oData->x_pos($i);
			my $y_pos = $oData->y_pos($i);
			next if (($ch1 <= 0)||($ch2<=0));
			push(@$aCh1, $ch1);
			push(@$aCh2, $ch2);
			push(@$aXcoords, $x_pos);
			push(@$aYcoords, $y_pos);
		}
		$self->{ _ch1_values } = $aCh1;
		$self->{ _ch2_values } = $aCh2;
		$self->{ _x_coords }   = $aXcoords;
		$self->{ _y_coords }   = $aYcoords;
	}
	sub process_data {
		my $self = shift;
		my $aCh1 = $self->ch1_values;
		my $aCh2 = $self->ch2_values;
		## if a process_data method is not set in plot class, simply use the raw intensity data
		$self->{ _x_values } = $aCh1;
		$self->{ _y_values } = $aCh2;
	}
	## create a GD image and draw on it the desired plot
	sub make_plot {
		my $self  = shift;
		if (@_) {
			$self->parse_args(@_);
		}
		## Get the x and y coordiantes of the plot, dynamicaly set by the size of the dataset
		my ($x, $y) = $self->plot_dimensions;
		## normalise the plot data according to the plot dimensions
		$self->plot_values;
		$self->{ _gd_object } = GD::Image->new($x,$y);	

		$self->plot_outline;
		$self->plot_spots;
		$self->return_image;
	}
	sub return_image {
		my $self = shift;
		my $image = $self->gd_object;
		$image->png;
	}
	sub plot_dimensions {
		my $self = shift;
		my $scale = $self->scale;
		my $aY   = $self->y_values;
		my $aX   = $self->x_values;
		my ($y_min,$y_max,$y_range) = $self->data_range($aY);
		my ($x_min,$x_max,$x_range) = $self->data_range($aX);
		my $x = ($x_range * $scale);
		my $y = ($y_range * $scale);
		my $x_margin = 100;
		my $y_margin = 100;
		$self->{ _x_length } = $x;
		$self->{ _y_length } = $y;
		$self->{ _x_margin } = $x_margin;
		$self->{ _y_margin } = $y_margin;
		$self->{ _middle } = ($y + $y_margin)/2 ;
		return(($x + $x_margin), ($y + $y_margin));
	}
	sub plot_values {
		my $self = shift;
		my $aX = $self->x_values;
		my $aY = $self->y_values;
		## if a plot_values method is not set in plot class, simply plot the raw data
		$self->{ _plotx_values } = $aX;
		$self->{ _ploty_values } = $aY;
	}
	sub set_plot_background {
		my $self = shift;
		my $image = $self->gd_object;
		## first allocated colour is set to background
		my $grey = $image->colorAllocate(180,180,180);
	}
	#sub set_plot_axis {
	#
	#}
	#sub set_plot_labels {
	#
	#}
	#sub set_plot_title {
	#
	#}
	# plot the outline of the diagram, ready for spots to be added
	sub plot_outline {
		my $self = shift;
		my $image = $self->gd_object;
		my $scale = $self->scale;
		$self->set_plot_background;
		my $aY   = $self->ploty_values;
		my $aX   = $self->plotx_values;
		my ($x_min,$x_max,$x_range) = $self->data_range($aX);
		my ($y_min,$y_max,$y_range) = $self->data_range($aY);
		my $middle = $self->middle;
		my $x_length = $self->x_length;
		my $y_length = $self->y_length;
		my $x_margin = $self->x_margin;
		my $y_margin = $self->y_margin;
		# get colours from the GD colour table 
		my $black = $image->colorExact(0,0,0);
		## x axis
		$image->filledRectangle(50,$y_length+50,$x_length+50,$y_length+50,$black);
		## y axis
		$image->filledRectangle(50,50,50,$y_length+50,$black);
		## x axis label
		my $max_label = int($x_max);
		my $min_label = int($x_min);
		$image->string($image->gdGiantFont,($x_length+50)/2,$y_length+80,"(0.5)*Log2(R*G)",$black);
		$image->string($image->gdGiantFont,(($x_min-$x_min)*$scale)+50,$y_length+70,"$min_label",$black);
		$image->string($image->gdGiantFont,(($x_max-$x_min)*$scale)+50,$y_length+70,"$max_label",$black);
		## y axis label
		$image->stringUp($image->gdGiantFont,10,$middle,"Log2(R/G)",$black);
		$image->stringUp($image->gdGiantFont,30,$middle - 25,"0",$black);
	}
	## given a GD image object and X and Y data arrays this method will plot on the 
	## image each X,Y data point in black
	sub plot_spots {
		my $self = shift;
		my $image = $self->gd_object;
		my $aX   = $self->plotx_values;
		my $aY   = $self->ploty_values;
		my $black = $image->colorExact(0,0,0);
		for (my $i=0; $i<@$aX; $i++){
			my $x = $aX->[$i];
			my $y = $aY->[$i];
			next unless ($x && $y);
			$image->filledEllipse($x,$y,3,3,$black);
		}
	}
	## finds the minimum, maximum and range of an array of data 
	sub data_range {
		use Statistics::Descriptive;
		my $self = shift;
		my $aData = shift;
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@$aData);
		my $min = $stat->min();
		my $max = $stat->max();
		my $range = $stat->sample_range();
		return($min,$max,$range);
	}
	## set graduated colours in the GD colour table, for use in the plot
	## from red-yellow-green 
	sub make_colour_grads {
		my $self  = shift;
		my $image = $self->gd_object;
		my $count = 0;
		for (my $i = 0; $i<=255; $i+=5){ 	
			$image->colorAllocate($i,255,0);
			$image->colorAllocate(255,$i,0);
		}	
		$image->colorAllocate(0,0,0);		
	}
	sub print_image {
		my $self = shift;
		my $image = $self->gd_object;
		if (@_) {
			$self->{ _print_location} = shift;
			$self->{ _image_name} = shift;
		}
		my $name  = $self->get_image_name;  
		my $print_location = $self->print_location;
		open FH, $print_location.$name.".png";
		binmode(FH);
		print FH $image->png;
		close(FH);
	}
	sub get_image_name {
		my $self = shift;
		unless (defined $self->{ _image_name }){
			$self->set_image_name;
		}
		$self->{ _image_name };
	}
	sub set_image_name {
		my $self = shift;
		## need to concat plot type and data source(file name/db id)
		my $image_name = 'image';
		$self->{ _image_name } = $image_name;
	}
	sub print_location {
		my $self = shift;
		if (@_)	{
			$self->{ _print_location } = shift;
		} else {
			if (defined $self->{ _print_location }){
				$self->{ _print_location };
			} else {
				die "Microarray::Image ERROR; No print destination defined\n";
			}
		}
	}
	sub file_path {
		my $self = shift;
		$self->{ _file_path };
	}
	sub microarray_file {
		my $self = shift;
		$self->{ _microarray_file };
	}
	sub ch1_values {
		my $self = shift;
		$self->{ _ch1_values };
	}
	sub ch2_values {
		my $self = shift;
		$self->{ _ch2_values };
	}
	sub middle {
		my $self = shift;
		$self->{ _middle };
	}
	sub scale {
		my $self = shift;
		if (@_){
			$self->{ _scale } = shift;
		} else {
			if (defined $self->{ _scale }){
				$self->{ _scale };
			} else {
				$self->default_scale;
			}
		}
	}
	sub default_scale {
		100; 
	}
	sub x_length {
		my $self = shift;
		$self->{ _x_length };
	}
	sub y_length {
		my $self = shift;
		$self->{ _y_length };
	}
	sub x_margin {
		my $self = shift;
		$self->{ _x_margin };
	}
	sub y_margin {
		my $self = shift;
		$self->{ _y_margin };
	}
	sub x_values {
		my $self = shift;
		$self->{ _x_values };
	}
	sub y_values {
		my $self = shift;
		$self->{ _y_values };
	}
	sub plotx_values {
		my $self = shift;
		$self->{ _plotx_values };
	}
	sub ploty_values {
		my $self = shift;
		$self->{ _ploty_values };
	}
	sub ratio_values {
		my $self = shift;
		$self->{ _ratio_values };
	}
	sub x_coords {
		my $self = shift;
		$self->{ _x_coords };
	}
	sub y_coords {
		my $self = shift;
		$self->{ _y_coords };
	}
	sub data_object {
		my $self = shift;
		@_	?	$self->{ _data_object } = shift
			:	$self->{ _data_object };
	}
	sub data {
		my $self = shift;
		$self->{ _data };
	}
}

{ package ma_plot;
  
	our @ISA = qw( plot );

	sub process_data {
		my $self = shift;
		my $aCh1 = $self->ch1_values;
		my $aCh2 = $self->ch2_values;
		my $aMvalues = [];
		my $aAvalues = [];
		for (my $i=0; $i<@$aCh1; $i++){
			my $ch1 = $aCh1->[$i];
			my $ch2 = $aCh2->[$i];
			next unless ($ch1 && $ch2);
			my $m = $self->calc_m($ch1, $ch2);
			my $a = $self->calc_a($ch1, $ch2);
			next unless ($m && $a);
			push(@$aMvalues, $m);
			push(@$aAvalues, $a);  
		}
		$self->{ _y_values } = $aMvalues;
		$self->{ _x_values } = $aAvalues;
	}
	sub plot_values {
		my $self = shift;
		my $aX = $self->x_values;
		my $aY = $self->y_values;
		my $scale  = $self->scale;
		my $middle = $self->middle;
		my ($y_min,$y_max,$y_range) = $self->data_range($aY);
		my ($x_min,$x_max,$x_range) = $self->data_range($aX);
		my $aXadjusted = [];
		my $aYadjusted = [];
		for (my $i=0; $i<@$aX; $i++){
			my $x = $aX->[$i];
			my $y = $aY->[$i];
			next unless ($x && $y);
			$x = (($x - $x_min) * $scale) + 50;
			$y = ($middle - ($y * $scale));
			push(@$aXadjusted, $x);
			push(@$aYadjusted, $y);  
		}
		$self->{ _plotx_values } = $aXadjusted;
		$self->{ _ploty_values } = $aYadjusted;
	}
	sub calc_m {
		my $self = shift;
		my $ch1  = shift;
		my $ch2  = shift;
		if (($ch1 == 0 )||($ch2 == 0)){
			return();  
		}
		return(log($ch1/$ch2)/log(2));
	}
	sub calc_a {
		my $self = shift;
		my $ch1  = shift;
		my $ch2  = shift;
		if (($ch1 == 0 )||($ch2 == 0)){
			return();  
		}
		return(log(0.5*($ch1*$ch2))/log(2));
	}
}

{ package ri_plot;

	our @ISA = qw( ma_plot );

	sub process_data {
		my $self = shift;
		my $aCh1 = $self->ch1_values;
		my $aCh2 = $self->ch2_values;
		my $aRvalues = [];
		my $aIvalues = [];
		for (my $i=0; $i<@$aCh1; $i++){
			my $ch1 = $aCh1->[$i];
			my $ch2 = $aCh2->[$i];
			next unless ($ch1 && $ch2);
			my $r_val = $self->calc_r($ch1, $ch2);
			my $i_val = $self->calc_i($ch1, $ch2);  
			next unless ($r_val && $i_val);  
			push(@$aRvalues, $r_val);
			push(@$aIvalues, $i_val);
		}
		$self->{ _y_values } = $aRvalues;
		$self->{ _x_values } = $aIvalues;
	}
	sub calc_r {
		my $self = shift;
		my $ch1  = shift;
		my $ch2  = shift;
		if (($ch1 == 0 )||($ch2 == 0)){
			return();  
		}
		return(log($ch1/$ch2)/log(2));
	}
	sub calc_i {
		my $self = shift;
		my $ch1  = shift;
		my $ch2  = shift;
		if (($ch1 == 0 )||($ch2 == 0)){
			return();  
		}
		return(log($ch1*$ch2)/log(2));
	}
	# plot the outline of the diagram, ready for spots to be added
	sub plot_outline {
		my $self = shift;
		my $image = $self->gd_object;
		my $scale = $self->scale;
		## first allocated colour is set to background
		my $grey = $image->colorAllocate(180,180,180);
		my $aM   = $self->ploty_values;
		my $aA   = $self->plotx_values;
		my ($m_min,$m_max,$m_range) = $self->data_range($aM);
		my ($a_min,$a_max,$a_range) = $self->data_range($aA);
		my $middle = $self->middle;
		my $x_length = $self->x_length;
		my $y_length = $self->y_length;
		my $x_margin = $self->x_margin;
		my $y_margin = $self->y_margin;
		# get colours from the GD colour table 
		my $black = $image->colorExact(0,0,0);
		## x axis
		$image->filledRectangle(50,$y_length+50,$x_length+50,$y_length+50,$black);
		## y axis
		$image->filledRectangle(50,50,50,$y_length+50,$black);
		## x axis label
		my $max_label = int($a_max);
		my $min_label = int($a_min);
		$image->string($image->gdGiantFont,($x_length+50)/2,$y_length+80,"Log2(R*G)",$black);
		$image->string($image->gdGiantFont,(($a_min-$a_min)*$scale)+50,$y_length+70,"$min_label",$black);
		$image->string($image->gdGiantFont,(($a_max-$a_min)*$scale)+50,$y_length+70,"$max_label",$black);
		## y axis label
		$image->stringUp($image->gdGiantFont,10,$middle,"Log2(R/G)",$black);
		$image->stringUp($image->gdGiantFont,30,$middle - 25,"0",$black);
	}
}
{ package intensity_scatter;

	our @ISA = qw( plot );

	sub plot_dimensions {
		my $self = shift;
		my $scale = $self->scale;
		my $aX    = $self->x_values;
		my $aY    = $self->y_values;
		my ($x_min,$x_max,$x_range) = $self->data_range($aX);
		my ($y_min,$y_max,$y_range) = $self->data_range($aY);
		my $x = ($x_max / $scale);
		my $y = ($y_max / $scale);
		my $x_margin = 100;
		my $y_margin = 100;
		$self->{ _x_length } = $x;
		$self->{ _y_length } = $y;
		$self->{ _x_margin } = $x_margin;
		$self->{ _y_margin } = $y_margin;
		$self->{ _middle } = ($y + $y_margin)/2 ;
		return(($x + $x_margin), ($y + $y_margin));
	}
	sub plot_values {
		my $self = shift;
		my $aX = $self->x_values;
		my $aY = $self->y_values;
		my $scale  = $self->scale;
		my $x_length = $self->x_length;
		my $y_length = $self->y_length;
		my $x_margin = $self->x_margin;
		my $y_margin = $self->y_margin;
		my ($y_min,$y_max,$y_range) = $self->data_range($aY);
		my ($x_min,$x_max,$x_range) = $self->data_range($aX);
		my $aXadjusted = [];
		my $aYadjusted = [];
		for (my $i=0; $i<@$aX; $i++){
			my $x = $aX->[$i];
			my $y = $aY->[$i];
			next unless ($x && $y);
			$x = (($x/$scale))+50;
			$y = ($y_length-($y/$scale))+50;
			push(@$aXadjusted, $x);
			push(@$aYadjusted, $y);  
		}
		$self->{ _plotx_values } = $aXadjusted;
		$self->{ _ploty_values } = $aYadjusted;
	}
	# plot the outline of the diagram, ready for spots to be added
	sub plot_outline {
		my $self = shift;
		my $image = $self->gd_object;
		my $scale = $self->scale;
		## first allocated colour is set to background
		my $grey = $image->colorAllocate(180,180,180);
		my $aM   = $self->ch1_values;
		my $aA   = $self->ch2_values;
		my ($m_min,$m_max,$m_range) = $self->data_range($aM);
		my ($a_min,$a_max,$a_range) = $self->data_range($aA);
		my $middle = $self->middle;
		my $x_length = $self->x_length;
		my $y_length = $self->y_length;
		my $x_margin = $self->x_margin;
		my $y_margin = $self->y_margin;
		# get colours from the GD colour table 
		my $black = $image->colorExact(0,0,0);
		## x axis
		$image->line(50,$y_length+50,$x_length+50,$y_length+50,$black);
		## y axis
		$image->line(50,50,50,$y_length+50,$black);
		## x axis label
		my $x_max_label = int($a_max);
		my $y_max_label = int($m_max);
		$image->string($image->gdGiantFont,($x_length+50)/2,$y_length+80,"Channel 1 Intensity",$black);
		$image->string($image->gdGiantFont,50,$y_length+50,"0",$black);
		$image->string($image->gdGiantFont,$x_length+50,$y_length+50,"$x_max_label",$black);
		## y axis label
		$image->stringUp($image->gdGiantFont,10,$middle,"Channel 2 Intensity",$black);
		$image->stringUp($image->gdGiantFont,30,$y_length+50,"0",$black);
		$image->stringUp($image->gdGiantFont,30,50,"$y_max_label",$black);
	}
}

{ package heatmap;
  
	our @ISA = qw( plot );

	sub plot_dimensions {
		my $self  = shift;
		my $aY    = $self->y_coords;
		my $aX    = $self->x_coords;
		my ($y_min,$y_max,$y_range) = $self->data_range($aY);
		my ($x_min,$x_max,$x_range) = $self->data_range($aX);		
		my $scale = $self->calculate_scale($x_range,$y_range);  
		my $x = int((($x_max + $x_min) / $scale) +1);	# include margin, and round up 
		my $y = int((($y_max + $y_min) / $scale) +1);
		$self->{ _x_length } = $x;
		$self->{ _y_length } = $y;		
		return($x, $y);
	}
	sub plot_values {
		return;
	}
	# we dynamically set the scale to best match the array layout.
	#Êcan get the spot layout from the data file, if available 
	# or alternatively, the user can set the number of x/y spots
	# or alternatively, the default scale will be used 
	sub calculate_scale {
		my $self = shift;
		my ($x_range,$y_range) = @_;
		my $x_spots = $self->x_spots;
		my $y_spots = $self->y_spots;
		if ($x_spots && $y_spots){
			my $x_scale = $x_range/$x_spots;
			my $y_scale = $y_range/$y_spots;
			if ($x_scale < $y_scale){
				$self->scale(int($x_scale));
			} else {
				$self->scale(int($y_scale));
			}
		} 
		return $self->scale;
	}
	sub x_spots {
		my $self = shift;		
		if (@_){
			$self->{ _x_spots } = shift;
		} else {
			unless (defined $self->{ _x_spots }){
				my $oData  = $self->data_object;
				if ($oData->can('array_columns') && $oData->can('spot_columns')){
					$self->{ _x_spots } = $oData->array_columns * ($oData->spot_columns + 1);
				}
			}
			$self->{ _x_spots };
		}
	}
	sub y_spots {
		my $self = shift;		
		if (@_){
			$self->{ _y_spots } = shift;
		} else {
			unless (defined $self->{ _y_spots }){
				my $oData  = $self->data_object;
				if ($oData->can('array_rows') && $oData->can('spot_rows')){
					$self->{ _y_spots } = $oData->array_rows * ($oData->spot_rows + 1);
				}
			}
			$self->{ _y_spots };
		}
	}
	sub default_scale {
		50;
	}
}

{ package log2_heatmap;

	our @ISA = qw( heatmap );
  
	sub process_data {
		my $self = shift;
		my $aCh1 = $self->ch1_values;
		my $aCh2 = $self->ch2_values;
		my $aRatio = [];
		for (my $i=0; $i<@$aCh1; $i++){
			my $ch1 = $aCh1->[$i];
			my $ch2 = $aCh2->[$i];
			next unless ($ch1 && $ch2);
			push(@$aRatio, log($ch1/$ch2)/log(2));
		}
		$self->{ _x_values } = $aCh1;
		$self->{ _y_values } = $aCh2;
		$self->{ _ratio_values } = $aRatio;
	}
	sub plot_outline {
		my $self = shift;
		my $image = $self->gd_object;
		## first allocated colour is set to background
		my $grey = $image->colorAllocate(180,180,180);
	}
	sub plot_spots {
		my $self = shift;
		my $image = $self->gd_object;
		my $aXcoords = $self->x_coords;
		my $aYcoords = $self->y_coords;
		my $aRatio   = $self->ratio_values;
		my $scale    = $self->scale;
		$self->make_colour_grads;
		for (my $i=0; $i<@$aXcoords; $i++){
			my $x = $aXcoords->[$i];
			my $y = $aYcoords->[$i];
			my $ratio = $aRatio->[$i];
			next unless ($x && $y && $ratio);
			$x = int(($x / $scale) + 1);
			$y = int(($y / $scale) + 1);
			my $colour = $self->get_colour($ratio);
			$image->setPixel($x,$y,$colour);
		}
	}
	sub get_colour {
		my $self  = shift;
		my $ratio = shift;
		my $image = $self->gd_object;
		my $colour;
		if ($ratio <= -1.1){
			$colour = $image->colorExact(255,0,0);
		} elsif ($ratio >= 1.1){
			$colour = $image->colorExact(0,255,0); 
		} elsif ((0.1 > $ratio)&&($ratio > -0.1)) {
			$colour = $image->colorExact(255,255,0);
		} elsif ($ratio >= 0.1) {
			my $red_hue = 255 - (255 * ($ratio - 0.1));
			$colour = $image->colorClosest($red_hue,255,0);		# reducing red, closer to green
		} else {
			my $green_hue = 255 + (255 * ($ratio + 0.1));
			$colour = $image->colorClosest(255,$green_hue,0);	# reducing green, closer to red
		}
		return($colour);
	}  
	sub make_colour_grads {
		my $self  = shift;
		my $image = $self->gd_object;
		$image->colorAllocate(255,255,0);
		for (my $i = 0; $i<255; $i+=2){ 
			$image->colorAllocate($i,255,0);	## Add red -> green = yellow
			$image->colorAllocate(255,$i,0); 	## Add green -> red = yellow
		}	
	}
	
}

{ package intensity_heatmap;

	our @ISA = qw( heatmap );

	sub plot_outline {
		my $self = shift;
		my $image = $self->gd_object;
		## first allocated colour is set to background
		my $black = $image->colorAllocate(0,0,0);
	}
	sub plot_channel {
		my $self = shift;
		@_	?	$self->{ _plot_channel } = shift
			:	$self->{ _plot_channel };
	}
	sub plot_spots {
		my $self = shift;
		my $image = $self->gd_object;
		my $aXcoords = $self->x_coords;
		my $aYcoords = $self->y_coords;
		my $scale    = $self->scale;
		my $plot_channel;
		if ($self->plot_channel && ($self->plot_channel == 2)){
			$plot_channel = 'ch2_values';
		} else {
			$plot_channel = 'ch1_values';
		}
		my $aValues  = $self->$plot_channel;
		$self->make_colour_grads;
		for (my $i=0; $i<@$aXcoords; $i++){
			my $x = $aXcoords->[$i];
			my $y = $aYcoords->[$i];
			my $value = $aValues->[$i];
			next unless ($x && $y && $value);
			$x = int(($x / $scale) + 1);
			$y = int(($y / $scale) + 1);
			my $colour = $self->get_colour($value);
			$image->setPixel($x,$y,$colour);
		}
	}
	sub get_colour {
		my $self  = shift;
		my $value = shift;
		my $image = $self->gd_object;
		my $colour;
		if ($value == 0) {  ## if the value is 0 colour black
			$colour = $image->colorExact(0,0,0);
		} elsif ($value <= 13000) {  ## colour towards blue
			my $blue_hue = $value / 50.9;
			$colour = $image->colorClosest(0,0,$blue_hue);  
		} elsif ($value <= 26000) {  ## colour towards turquoise
			my $turquoise_hue = ($value - 13000) / 50.9;
			$colour = $image->colorClosest(0,$turquoise_hue,255);  
		} elsif ($value <= 39000) {  ## colour towards green
			my $green_hue = ($value - 26000) / 50.9;
			$colour = $image->colorClosest(0,255,255-$green_hue);  
		} elsif ($value <= 52000) {  ## colour towards yellow
			my $yellow_hue = ($value - 39000) / 50.9;
			$colour = $image->colorClosest($yellow_hue,255,0);  
		} elsif ($value < 65000) {  ## colour towards red
			my $red_hue = ($value - 52000) / 50.9;
			$colour = $image->colorClosest(255,255-$red_hue,0);  
		} elsif ($value >= 65000) {  ## if value is saturated colour white
			$colour = $image->colorExact(255,255,255);
		}
		return($colour);
	}
	# set a rainbow of graduated colours in the GD colour table, for use in the plot 
	sub make_colour_grads {
		my $self  = shift;
		my $image = $self->gd_object;
		my $count = 0;
		$image->colorAllocate(0,0,0); 
		for (my $i = 5; $i<=255; $i+=5){ 	
			$image->colorAllocate(0,0,$i);        ## Add blue up to 255 -> blue
			$image->colorAllocate(0,$i,255);      ## Add green up to 255 -> turquise
			$image->colorAllocate(0,255,255-$i);  ## Reduce blue -> green
			$image->colorAllocate($i,255,0);      ## Add red up to 255 -> yellow
			$image->colorAllocate(255,255-$i,0);  ## Reduce green -> red
		}	
	}
}

{ package cgh_plot;
  
	our @ISA = qw( plot );

	sub new {
		my $class = shift;
		my $self  = { };
		bless $self, $class;
		if (@_){
			$self->data_object(shift);
			if (@_){
				$self->clone_locns_file(shift);
			}
		} 
		return $self;
	}
	sub set_data {
		my $self = shift;
		$self->data_object(shift);
		if (@_){
			$self->clone_locns_file(shift);
		}
	}
	sub sort_data {
		my $self = shift;
		my $plot_chr = $self->plot_chromosome;
		my $oData_File = $self->data_object;
		die "Microarray::Image ERROR: No data object provided\n" unless $oData_File;
				
		if ($oData_File->isa('pipeline_processed_data') ||
			$oData_File->isa('processed_scanarray_data') ||
			$oData_File->isa('processed_bluefuse_data') ){	# already sorted by seq_start, already flip_flopped
			$self->{ _x_values } = $oData_File->all_locns;
			$self->{ _y_values } = $oData_File->all_log2_ratio;
			$self->{ _feature_names } = $oData_File->all_feature_names;
		} else {
			$oData_File->flip if ($self->flip_flop == -1);
			my $spot_count = $oData_File->spot_count;
			if ($self->embedded_locns){
				for (my $i=0; $i<$spot_count; $i++){
					if (my $embedded_locn = $oData_File->feature_id($i)){
						my ($chr,$locn) = $self->parse_embedded_locn($embedded_locn);
						next unless ($plot_chr eq $chr);
						if (my $log = $oData_File->log2_ratio($i)){
							$self->set_feature_data($oData_File->synonym_id($i),$log,$locn);
						}
					}
				}	
			} elsif (my $oClone_Positions = $self->clone_locns_file) {
				my $hClones = $oClone_Positions->clone_hash;
				for (my $i=0; $i<$spot_count; $i++){
					my $feature = $oData_File->feature_id($i);
					next unless (	(defined $$hClones{$feature}) && 
									($plot_chr eq $$hClones{$feature}{_chr}) );
					if (my $log = $oData_File->log2_ratio($i)){
						$self->set_feature_data($feature,$log,$oClone_Positions->location($feature));
					}
				}
			} else {
				die "Microarray::Image cgh_plot ERROR; No clone positions to work with\n";
			}
			$self->order_data;
		}
		$self->{ _data_sorted } = 1;
	}
	sub plot_chromosome {
		my $self = shift;
		if (@_){
			$self->{ _plot_chromosome } = shift;
		} else {
			if (defined $self->{ _plot_chromosome }){
				$self->{ _plot_chromosome };
			} else {
				die "Microarray::Image ERROR; No plot chromosome was specified\n";
			}
		}
	}
	sub data_sorted {
		my $self = shift;
		$self->{ _data_sorted };
	}
	sub order_data {
		my $self = shift;
		my $hFeatures = $self->features;
		my @aFeatures = keys %$hFeatures;
		my $aSorted_Features = [];
		if ($self->smoothing){
			@$aSorted_Features = sort { $$hFeatures{ $a }{locn} <=> $$hFeatures{ $b }{locn} } @aFeatures;
		} else {
			$aSorted_Features = \@aFeatures;
		}
		my @aLog_Ratios = ();
		my @aLocns = ();
		
		for my $feature (@$aSorted_Features){
			push(@aLocns,$self->feature_locn($feature));
			push(@aLog_Ratios,$self->feature_log($feature));
		}

		$self->{ _x_values } = [\@aLocns];
		$self->{ _y_values } = [\@aLog_Ratios];
		$self->{ _feature_names } = [$aSorted_Features];
	}
	sub flip_flop {
		my $self = shift;
		if (defined $self->{ _flip_flop }){
			$self->{ _flip_flop };
		} else {
			return 1;
		}
	}
	sub flip {
		my $self = shift;
		$self->{ _flip_flop } = -1;
	}
	sub flop {
		my $self = shift;
		$self->{ _flip_flop } = 1;
	}
	sub features {
		my $self = shift;
		unless (defined $self->{ _features }){
			$self->{ _features } = {};
		}
		$self->{ _features };
	}
	sub feature_names {
		my $self = shift;
		$self->{ _feature_names };
	}
	sub feature_chrs {
		my $self = shift;
		$self->{ _feature_chrs };
	}
	sub get_feature_data {
		my $self = shift;
		my $feature = shift;
		my $hFeatures = $self->features;
		unless (defined $hFeatures->{ $feature }){
			$hFeatures->{ $feature } = { ratios => [] };
		}
		$hFeatures->{ $feature };
	}
	sub set_feature_data {
		my $self = shift;
		my $hFeature_data = $self->get_feature_data(shift);	
		my $aRatios = $hFeature_data->{ ratios };
		push(@$aRatios,shift);
		$hFeature_data->{ locn } = shift;
	}
	sub feature_locn {
		my $self = shift;
		my $hFeature = $self->get_feature_data(shift,shift);	# second shift for genome plot - chromosome name
		return $hFeature->{locn};
	}
	sub feature_log {
		my $self = shift;
		my $hFeature = $self->get_feature_data(shift,shift);	# second shift for genome plot - chromosome name
		my $aFeat_Ratios = $hFeature->{ratios};
		if (@$aFeat_Ratios == 1){
			return $aFeat_Ratios->[0];
		} else {
			my $log_ratio;
			for my $ratio (@$aFeat_Ratios){
				$log_ratio += $ratio;
			}
			return ($log_ratio/scalar @$aFeat_Ratios);
		}
	}
	sub parse_embedded_locn {
		my $self = shift;
		my $location = shift;
		my ($chr,$start,$end) = split(/:|\.\./,$location);
		$chr =~ s/chr//;
		return ($chr,int(($start+$end)/2));
	}
	sub clone_locns_file {
		my $self = shift;
		@_	?	$self->{ _clone_locns_file } = shift
			:	$self->{ _clone_locns_file };
	}
	sub embedded_locns {
		my $self = shift;
		@_	?	$self->{ _embedded_locns } = shift
			:	$self->{ _embedded_locns };
	}
	sub plot_centromere {
		my $self = shift;
		if (@_){
			$self->{ _plot_centromere } = shift;
		} else {
			if (defined $self->{ _plot_centromere }){
				return $self->{ _plot_centromere };
			} else {
				return 1;
			}
		}
	}
	sub shift_zero {
		my $self = shift;
		@_	?	$self->{ _shift_zero } = shift
			:	$self->{ _shift_zero };
	}
	sub make_plot {
		my $self  = shift;
		if (@_){
			$self->parse_args(@_);
		}
		$self->sort_data;
		my ($x, $y) = $self->plot_dimensions;
		## set the y value the same as the genome plot
		$self->{ _gd_object } = GD::Image->new($x,$y);
		$self->set_plot_background;
		$self->make_colour_grads;
		$self->make_outline($y);
		if ($self->smoothing) {
			$self->smooth_data_by_location;
		}			
		$self->plot_spots;
		$self->return_image;
	}
	sub default_scale {
		500000;
	}
	# normalise the processed data relative to the image ready to be plotted
	sub plot_data {
		my $self = shift;
		my $aaX = $self->x_values;
		my $aaY = $self->y_values;
		my $aX  = $$aaX[0];
		my $aY  = $$aaY[0];
		my $scale    = $self->scale;
		my $middle   = $self->middle;
		my $x_margin = $self->x_margin;
		my $zero_shift = $self->shift_zero;
		my $aLocns = [];
		my $aLog2  = [];
		for (my $i=0; $i<@$aX; $i++ ){
			my $locn = $aX->[$i];
			my $log2 = $aY->[$i];
			next unless($locn && $log2);
			push(@$aLocns, int($locn/$scale));

			## divide the y axis by 3 for a +1.5/-1.5 plot
			$log2 += $zero_shift if ($zero_shift);
			push(@$aLog2, $middle - ($log2 * (450/3)));
		}
		$self->{ _plotx_values } = $aLocns;
		$self->{ _ploty_values } = $aLog2;
		return($aLocns,$aLog2);
	}
	sub plot_dimensions {
		my $self = shift;
		my $scale = $self->scale;
		my $aX   = $self->x_values;
		my ($x_min,$x_max,$x_range) = $self->data_range($aX);
		my $x = ($x_range /$scale);
		my $y = 450;
		my $x_margin = 0;
		my $y_margin = 0;
		$self->{ _x_length } = $x;
		$self->{ _y_length } = $y;
		$self->{ _x_margin } = $x_margin;
		$self->{ _y_margin } = $y_margin;
		$self->{ _middle } = ($y/2)+($y_margin/2) ;
		return(($x + $x_margin), ($y + $y_margin));
	}
	sub plot_spots {
		my $self = shift;
		my $image = $self->gd_object;
		my ($aX,$aY) = $self->plot_data;
		# get colour from the GD colour table
		for (my $i=0; $i<@$aX; $i++){
			my $x = $aX->[$i];
			my $y = $aY->[$i];
			next unless ($x && $y);
			my $colour = $self->get_colour($y);
			$image->filledEllipse($x,$y,3,3,$colour);
		}
	}
	# plot the outline of the diagram, ready for spots to be added
	sub make_outline {	
		use GD;
		my $self = shift;
		my $y = shift;

		my $image = $self->gd_object;
		my $scale = $self->scale;
		my $chr = $self->plot_chromosome;
		# get colours from the GD colour table 
		my $black 	= $image->colorExact(0,0,0);
		my $red   	= $image->colorExact(255,0,0);      
		my $green 	= $image->colorExact(0,255,0);      
		my $blue	= $image->colorExact(125,125,255);   
		# 3px wide log2 ratio lines
		$image->filledRectangle(0,150,3080,150,$red);		# +0.5
		$image->filledRectangle(0,225,3080,225,$green);	#  0.0
		$image->filledRectangle(0,300,3080,300,$red);		# -0.5
		# axis labels
		$image->string($image->gdGiantFont,10,150,'0.5',$black);
		$image->string($image->gdGiantFont,10,225,'0',$black);
		$image->string($image->gdGiantFont,10,300,'-0.5',$black);
		
		if ($self->plot_centromere){
			# dashed style for centromere lines
			$image->setStyle($blue,$blue,$blue,$blue,gdTransparent,gdTransparent);
			my $cen = int($self->chr_centromere($chr)/$scale);
			$image->line($cen,0,$cen,$y,gdStyled);
		}
	}
	sub do_smoothing {
		my $self = shift;
		@_	?	$self->{ _do_smoothing } = shift
			:	$self->{ _do_smoothing };
	}
	sub smoothing {
		my $self = shift;
		if (@_){
			$self->smooth_window(shift);
			$self->smooth_step(shift);
			$self->{ _do_smoothing }++;
		} else {
			$self->{ _do_smoothing };		
		}
	}
	sub smooth_window {
		my $self = shift;
		if (@_){
			$self->{ _smooth_window } = shift;
			$self->{ _do_smoothing }++;
		} else {
			if (defined $self->{ _smooth_window }){
				$self->{ _smooth_window };		
			} else {
				$self->default_smooth_window;
			}
		}
	}
	sub default_smooth_window {
		500000
	}
	sub smooth_step {
		my $self = shift;
		if(@_){
			$self->{ _smooth_step } = shift;
			$self->{ _do_smoothing }++;
		} else {
			if (defined $self->{ _smooth_step }){
				$self->{ _smooth_step };		
			} else {
				$self->default_smooth_step;
			}
		}
	}
	sub default_smooth_step {
		150000
	}
	sub is_smoothed {
		my $self = shift;
		$self->{ _smoothed };
	}
	# apply a moving average to the log2 data
	sub smooth_data_by_clone {
		my $self   = shift;
		my $window = $self->smooth_window;
		my $step   = $self->smooth_step;
		my $aX = $self->x_values;
		my $aY = $self->y_values;
		my $aSmooth_x = [];
		my $aSmooth_y = [];
		
		for (my $i=0; $i<(@$aX-($window-1)); $i+=$step){
			my $end = $i+($window-1);
			my @xslice = @$aX[$i..$i+($window-1)];
			my @yslice = @$aY[$i..$i+($window-1)];
			my $locn_sum = 0;
			my $log_sum = 0;
			foreach my $locn (@xslice) {
				$locn_sum += $locn;
			}
			foreach my $log (@yslice) {
				$log_sum += $log;
			}
			push(@$aSmooth_x, ($locn_sum/$window));
			push(@$aSmooth_y, ($log_sum/$window));
		}
		$self->{ _x_values } = $aSmooth_x;
		$self->{ _y_values } = $aSmooth_y;
		$self->{ _smoothed } = 1;
	}
	
	# this method uses a moving window of a specific genomic length, 
	# that moves along the genome in steps of a defined length,
	# to smooth our CGH profiles. 
	sub smooth_data_by_location {
		my $self 			= shift;
		
		# our smoothing parameters
		my $window 			= $self->smooth_window;
		my $step 			= $self->smooth_step;
		
		# all of the sorted data
		my $aAll_Locns 		= $self->x_values;
		my $aAll_Logs 		= $self->y_values;

		# arrays to hold the final smoothed data
		my @aAll_Smooth_Locns 	= ();
		my @aAll_Smooth_Logs 	= ();

		# set the chromosome we are working with - single or whole genome?
		my @aPlot_Chromosomes;
		if (my $chr = $self->plot_chromosome){
			@aPlot_Chromosomes = ($chr);
		} else {
			@aPlot_Chromosomes = (1..22,'X','Y');
		}
		
		# scroll through the individual chromosomes...
		for (my $j=0; $j<@aPlot_Chromosomes; $j++){
		
			# ...and get the sorted data for this chromosome
			my $alocns 	= $aAll_Locns->[$j];
			my $alogs 	= $aAll_Logs->[$j];
			
			# reset the window start and end location
			my $start = $alocns->[0];
			my $end = $start + $window;
			
			# arrays to hold the smoothed data for this chromosome
			my @aSmooth_Chr_Locns 	= ();
			my @aSmooth_Chr_Logs 	= ();
			
			#Êarrays to hold data in the moving window
			my @aWindow_Logs		= ();
			my @aWindow_Locns		= ();
				
			# scroll through the sorted data
			for (my $i=0; $i<@$alocns; $i++){
				my $genomic_locn = $alocns->[$i];
				my $log_value    = $alogs->[$i]; 
				
				# are we past the end of the window?
				if ($genomic_locn > $end){
					# if so, average up what's in the window...
					my ($av_locn, $av_log) = $self->moving_average(\@aWindow_Locns, \@aWindow_Logs);
					# ...add these values to the smoothed data arrays...
					push(@aSmooth_Chr_Locns, $av_locn);
					push(@aSmooth_Chr_Logs, $av_log);
					# ...move the end of the window to include the next location...
					while ($genomic_locn > $end){ 
						$end = (int($genomic_locn/100000) * 100000) + $step;
						# ...move the start up as well...
						$start = $end - $window;
						# ...and remove any data that is no longer in the window
						while ((@aWindow_Locns) && ($aWindow_Locns[0] < $start)){		# get rid of any values now out of the region 
							my $shifted1 = shift @aWindow_Locns;
							my $shifted2 = shift @aWindow_Logs;
						}		
					}
				}
				
				# either this location fell in the window to start with, 
				# or the window has now been moved to include it
				# so we add it to the window array and continue to the next location
				push (@aWindow_Locns, $genomic_locn);
				push (@aWindow_Logs, $log_value);
			
			}
			# we've finished with this chromosome, but have some values left in the last window
			my ($av_locn, $av_log) = $self->moving_average(\@aWindow_Locns, \@aWindow_Logs);
			push(@aSmooth_Chr_Locns, $av_locn);
			push(@aSmooth_Chr_Logs, $av_log);
			
			# add the smoothed data for this chromosome to our array
			# and then continue to the next chromosome
			push(@aAll_Smooth_Locns,\@aSmooth_Chr_Locns);
			push(@aAll_Smooth_Logs,\@aSmooth_Chr_Logs);
		}

		# finally, set all the smoothed data to our plotting values
		$self->{ _x_values } = \@aAll_Smooth_Locns;
		$self->{ _y_values } = \@aAll_Smooth_Logs;
		$self->{ _smoothed } = 1;
	
	}
	sub moving_average {
		my $self = shift;
		my $Locns = shift;
		my $Logs  = shift;
		my ($av_locn, $av_log2, $med_log2);
		if (@$Locns > 1){
			my $locn_stat = Statistics::Descriptive::Full->new();
			my $log_stat = Statistics::Descriptive::Full->new();
			$log_stat->add_data($Logs); 
			$locn_stat->add_data($Locns); 
			$av_locn = $locn_stat->mean;
			$av_log2 = $log_stat->mean;
			$med_log2 = $log_stat->median;
		} elsif (@$Locns == 1){
			$av_locn = @$Locns[0];
			$av_log2 = @$Logs[0];
		}
		return($av_locn, $med_log2);
	}
	# set a rainbow of graduated colours in the GD colour table, for use in the plot 
	sub set_plot_background {
		my $self = shift;
		my $image = $self->gd_object;
		## first allocated colour is set to background
		my $grey = $image->colorAllocate(180,180,180);
	}
	sub make_colour_grads {
		my $self  = shift;
		my $image = $self->gd_object;
		$image->colorAllocate(0,0,0);
		$image->colorAllocate(125,125,255);
		
		for (my $i = 0; $i<=255; $i+=3){ 
			$image->colorAllocate($i,255,0);	## Add red -> green = yellow
			$image->colorAllocate(255,$i,0); 	## Add green -> red = yellow
		}	
	}
	sub get_colour {
		my $self  = shift;
		my $ratio = shift;
		
		my $image  = $self->gd_object;

		# get colours from the GD colour table
		my $red    = $image->colorExact(255,0,0);      
		my $green  = $image->colorExact(0,255,0);      
		my $yellow = $image->colorExact(255,255,0);      
		my $colour;
		if ($ratio <= 150){
			$colour = $red;
		} elsif ($ratio >= 300){
			$colour = $green;
		} elsif ((262 > $ratio)&&($ratio > 187)) {
			$colour = $yellow;
		} elsif ($ratio >= 262) {
			my $red_hue = int(255-(6.8*($ratio-262)));				# factorial = 255/(low_yellow - green)
			$colour = $image->colorClosest($red_hue,255,0);		# reducing red, closer to green
		} else {
			my $green_hue = int(255-(6.9*(187-$ratio)));				# factorial = 255/(high_yellow - red)
			$colour = $image->colorClosest(255,$green_hue,0);	# reducing green, closer to red
		}
		return($colour);
	}
	sub image_map {
		my $self      = shift;
		my $scale     = shift;
		my $link      = shift;   
		my $aX        = $self->x_values;
		my $aY        = $self->y_values;
		my $aFeatures = $self->feature_names;
		my $middle    = shift;
		print "<MAP NAME=\"Map_one\">";
		for (my $i=0; $i<@$aX; $i++){
			my $x = $aX->[$i];
			my $y = $aY->[$i];
			my $feature = $aFeatures->[$i];
			next unless ($x && $y);
			$y = $middle - ($y * $scale);
			$x = $x * $scale;
			print "<AREA SHAPE=\"Circle\" COORDS=\"$x,$y, 5\" HREF=$link>";
		}
		print "</MAP>";
	}
	sub chr_centromere {
		my $self = shift;
		my $chr  = shift;
		my %hCentromere = (		
			1 => 124200000,  
			2 => 93400000,
			3 => 91700000,
			4 => 50900000,
			5 => 47700000,
			6 => 60500000,
			7 => 58900000,
			8 => 45200000,
			9 => 50600000,
			10 => 40300000,
			11 => 52900000,
			12 => 35400000,
			13 => 16000000,
			14 => 15600000,
			15 => 17000000,
			16 => 38200000,
			17 => 22200000,
			18 => 16100000,
			19 => 28500000,
			20 => 27100000,
			21 => 12300000,
			22 => 11800000,
			23 => 59400000,
			'X' => 59400000,
			24 => 11500000,
			'Y' => 11500000
		);
		return ($hCentromere{$chr});
	}
	
}
{ package genome_cgh_plot;

	our @ISA = qw( cgh_plot );

	sub sort_data {
		my $self = shift;
		my $oData_File = $self->data_object;
		die "Microarray::Image ERROR: No data object provided\n" unless $oData_File;
				
		if ($oData_File->isa('processed_data') ||
			$oData_File->isa('processed_bluefuse_data') || 
			$oData_File->isa('processed_scanarray_data')){	# already sorted by [chr],seq_start
			$self->{ _x_values } = $oData_File->all_locns;
			$self->{ _y_values } = $oData_File->all_log2_ratio;
			$self->{ _feature_names } = $oData_File->all_feature_names;
		} else {
			$oData_File->flip if ($self->flip_flop == -1);
			$self->{ _features } = {};	# reset the features for this chromosome
			my $spot_count = $oData_File->spot_count;
			if ($self->embedded_locns){
				for (my $i=0; $i<$spot_count; $i++){
					if (my $embedded_locn = $oData_File->feature_id($i)){
						my ($chr,$locn) = $self->parse_embedded_locn($embedded_locn);
						if (my $log = $oData_File->log2_ratio($i)){
							$self->set_feature_data($oData_File->synonym_id($i),$log,$locn,$chr);
						}
					}
				}	
			} elsif (my $oClone_Positions = $self->clone_locns_file) {
				my $hClones = $oClone_Positions->clone_hash;
				for (my $i=0; $i<$spot_count; $i++){
					my $feature = $oData_File->feature_id($i);
					next unless (defined $$hClones{$feature});
					if (my $log = $oData_File->log2_ratio($i)){
						$self->set_feature_data($feature,$log,$oClone_Positions->location($feature),$$hClones{$feature}{_chr});
					}
				}
			} else {
				die "Microarray::Image cgh_plot ERROR; No clone positions to work with\n";
			}
			$self->order_genome_data;
		}
		$self->{ _data_sorted } = 1;
	}
	sub order_genome_data {
		my $self = shift;
		my $hFeatures = $self->features;
		
		my (@aFeatures,@aLocns,@aLog_Ratios);
		
		for my $chr ((1..22,'X','Y')){
		
			my $hChr_Features = $hFeatures->{ $chr };
			my @aChr_Features = keys %$hChr_Features;
			
			my $aSorted_Chr_Features = [];
			my $aSorted_Chr_Logs = [];
			my $aSorted_Chr_Locns = [];
			
			if ($self->smoothing){
				@$aSorted_Chr_Features = sort { $$hChr_Features{ $a }{locn} <=> $$hChr_Features{ $b }{locn} } @aChr_Features;
			} else {
				$aSorted_Chr_Features = \@aChr_Features;
			}
			for my $feature (@$aSorted_Chr_Features){
				push(@$aSorted_Chr_Locns,$$hChr_Features{$feature}{ locn });
				push(@$aSorted_Chr_Logs,$self->feature_log($feature,$chr));
			}
			push (@aFeatures,$aSorted_Chr_Features);
			push (@aLocns,$aSorted_Chr_Locns);
			push (@aLog_Ratios,$aSorted_Chr_Logs);
		}
		$self->{ _x_values } = \@aLocns;
		$self->{ _y_values } = \@aLog_Ratios;
		$self->{ _feature_names } = \@aFeatures;
	}
	sub set_feature_data {
		my $self = shift;
		my $hFeature_data = $self->get_feature_data(shift,pop);	# pop chromosome name
		my $aRatios = $hFeature_data->{ ratios };
		push(@$aRatios,shift);
		$hFeature_data->{ locn } = shift;
	}
	sub get_feature_data {
		my $self = shift;
		my $feature = shift;
		my $chr = shift;
		my $hFeatures = $self->features;
		unless (defined $hFeatures->{ $chr }{ $feature }){
			$hFeatures->{ $chr }{ $feature } = { ratios => [] };
		}
		$hFeatures->{ $chr }{ $feature };
	}
	sub make_plot {
		my $self  = shift;
		if (@_){
			$self->parse_args(@_);
		}
		## chr 25 is the total length of the genome plus 25 margin
		my $x_axis = ($self->chr_offset(25)/$self->scale)+25;
		$self->{ _gd_object } = GD::Image->new($x_axis,450);	
		$self->set_plot_background;
		$self->make_colour_grads;
		$self->make_outline;

		$self->sort_data;
		if ($self->smoothing) {
			$self->smooth_data_by_location;
		}			
		$self->plot_spots;
		$self->return_image;
	}
	sub default_scale {
		## set scale for the genome plot to 2.5Mb per pixel
		2500000;
	}
	sub plot_chromosome {
		return;
	}
	sub plot_data {
		my $self = shift;
		my $scale = $self->scale;
		my $aaX = $self->x_values;
		my $aaY = $self->y_values;
		
		my $zero_shift = $self->shift_zero;
		my $aLocns = [];
		my $aLog2  = [];
		
		for my $chr ((0..23)){
			my $aX = $$aaX[$chr];
			my $aY = $$aaY[$chr];
			for (my $i=0; $i<@$aX; $i++ ){
				my $locn = $aX->[$i];
				my $log2 = $aY->[$i];			
				next unless($locn && $log2);
				my $chr_offset = $self->chr_offset($chr+1);
				push(@$aLocns, int($locn/$scale) + ($chr_offset/$scale) + 25);
				$log2 += $zero_shift if ($zero_shift);
				## multiply the log value by a quarter of the y axis to get a +2/-2 plot 
				push(@$aLog2, 225 - ($log2 * (450/3)));
			}
		}
		$self->{ _plotx_values } = $aLocns;
		$self->{ _ploty_values } = $aLog2;
		return($aLocns,$aLog2);
	}
	# Harcode the plot outline for the genome plot as dimensions do not change
	sub make_outline {
		use GD;
		my $self = shift;
		my $image = $self->gd_object;
		my $scale = $self->scale;
		# get colours from the GD colour table 
		my $black 	= $image->colorExact(0,0,0);
		my $red   	= $image->colorExact(255,0,0);      
		my $green 	= $image->colorExact(0,255,0); 
		my $blue	= $image->colorExact(125,125,255);   
		# 3px wide log2 ratio lines
		$image->filledRectangle(0,150,3080,150,$red);		# +0.5
		$image->filledRectangle(0,225,3080,225,$green);		#  0.0
		$image->filledRectangle(0,300,3080,300,$red);		# -0.5
		# axis labels
		$image->string($image->gdSmallFont,0,150,'0.5',$black);
		$image->string($image->gdSmallFont,0,225,'0',$black);
		$image->string($image->gdSmallFont,0,300,'-0.5',$black);
		# dashed style for centromere lines
		$image->setStyle($blue,$blue,$blue,$blue,gdTransparent,gdTransparent);
		
		# plot chr separator lines and chr names for each chromosome
		for my $chr ((1..24)){
			my $start 	= int($self->chr_offset($chr)/$scale);
			my $end 	= int($self->chr_offset($chr+1)/$scale);
			my $middle 	= int(($start+$end)/2);
			if ($self->plot_centromere){
				# centromere
				my $cen = int(($self->chr_offset($chr)+$self->chr_centromere($chr))/$scale);
				$image->line($cen+25,0,$cen+25,450,gdStyled);
			}
			# chr buffer
			$image->filledRectangle($start+25,0,$start+25,450,$black);
			# set chr names
			my $chr_name;
			if ($chr == 23){
				$chr_name = 'X';
			} elsif ($chr == 24){
				$chr_name = 'Y';
				# end line
				$image->line($end+25,0,$end+25,450,$black);
			} else {
				$chr_name = $chr;
			}
			# print chr name at bottom of plot
			$image->string($image->gdSmallFont,$middle+20,425,$chr_name,$black);
		}
	}
	sub chr_offset {
		my $self = shift;
		my $chr  = shift;
		my %hChromosome = (		
			# start bp  		# chr length
			1 => 0,   			# 247249719
			2 => 247249720,		# 242951149
			3 => 490200869,		# 199501827
			4 => 689702696,		# 191273063
			5 => 880975759,		# 180857866
			6 => 1061833625,	# 170899992
			7 => 1232733617,	# 158821424
			8 => 1391555041,	# 146274826
			9 => 1537829867,	# 140273252
			10 => 1678103119,	# 135374737
			11 => 1813477856,	# 134452384
			12 => 1947930240,	# 132349534
			13 => 2080279774,	# 114142980
			14 => 2194422754,	# 106368585
			15 => 2300791339,	# 100338915
			16 => 2401130254,	# 88827254
			17 => 2489957508,	# 78774742
			18 => 2568732250,	# 76117153
			19 => 2644849403,	# 63811651
			20 => 2708661054,	# 62435964
			21 => 2771097018,	# 46944323
			22 => 2818041341,	# 49691432
			23 => 2867732773,	# 154913754
			'X' => 2867732773,	# 154913754
			24 => 3022646527,	# 57772954
			'Y' => 3022646527,	# 57772954
			25 => 3080419480	# END
		);
		return ($hChromosome{$chr});
	}
}

1;

__END__

=head1 NAME

Microarray::Image - A Perl module for creating microarray data plots

=head1 SYNOPSIS

	use Microarray::Image;
	use Microarray::File::Data_File;

	my $oData_File = data_file->new($data_file);
	my $oMA_Plot = ma_plot->new($oData_File);
	my $ma_plot_png = $oMA_Plot->make_plot;	

	open (PLOT,'>ma_plot.png');
	print PLOT $ma_plot_png;
	close PLOT;

=head1 DESCRIPTION

Microarray::Image is an object-oriented Perl module for creating microarray data plots from a scan data file, using the GD module and image library. A number of different plot types are supported, including MA, RI, intensity scatter, intensity heatmap, log2 heatmap and CGH plots. Currently, only the export of PNG (Portable Network Graphics - or 'PNGs Not GIFs') images is supported.   

=head1 QC/QA PLOTS

There are several plots for viewing basic microarray data for QC/QA purposes. Most of the parameters for these plots are the same, and only the class name used to create the plot object differs from one plot to another.

=head2 Standard Data Plots

=over

=item B<ma_plot>

See the SYNOPSIS for all there is to know about how to create an MA plot. To create any of the other plot types, just append C<'ma_plot'> in the above example with one of the class names listed below. 

=item B<ri_plot>

An RI plot is basically identical to an MA plot - at least in appearance.

=item B<intensity_scatter>

This is a plot of channel 1 signal vs channel 2 signal.

=back

=head2 Heatmaps

=over 

=item B<intensity_heatmap>

An image of the slide, using a black->white rainbow colour gradient to indicate the signal intensity across the array. Uses channel 1 as the signal by default, but the channel can be changed by setting the C<plot_channel> parameter in the call to C<make_plot()>.

	my $oInt_Heatmap = intensity_heatmap->new($oData_File);
	my $int_heatmap_png = $oInt_Heatmap->make_plot(plot_channel=>2);

=item B<log2_heatmap>

An image of the slide using a red->yellow->green colour gradient to indicate the Log2 of the signal ratio across the array. 

=back

One difference between heatmaps and other plots is in their implementation of the plot scale. This is calculated dynamically in order to generate the best looking image of the array, and requires the dimensions of the array in terms of the number of spots in the x and y axes. If you are using a data file format that returns those values in its header information (such as a Scanarray file, using the Quantarray module) then the scale will be calculated automatically. If BlueFuse files are sorted such that the last data row has the highest block/spot row/column number, then again the scale can be calculated automatically. However, for GenePix files, you will have to pass these values to the make_plot() method (adding extra spots for block padding where appropriate);

	my $oLog2_Heatmap = log2_heatmap->new($oData_File);
	my $log_heatmap_png = $oLog2_Heatmap->make_plot(x_spots=>108, y_spots=>336);  

=head1 CGH PLOT

There are two types of CGH plot - a single chromosome plot (C<cgh_plot>) or a whole genome plot (C<genome_cgh_plot>). The big difference between CGH plots and the other types described above is of course that they require genomic mapping data for each feature. This can be loaded into the object using a L<C<clone_locn_file>|Microarray::File::Clone_Locn_File> object (see below) or using information embedded in the data file by setting the C<embedded_locns> flag. 

	use Microarray::Image;
	use Microarray::File::Data_File;
	use Microarray::File::Clone_Locn_File;
	
	# first make your data objects
	my $oData_File = data_file->new($data_file);
	my $oClone_File = clone_locn_file->new($clone_file);
	
	# create the plot object
	my $oGenome_Image = genome_cgh_plot->new($oData_File,$oClone_File);
	my $oChrom_Image = cgh_plot->new($oData_File,$oClone_File);
	
	# make the plot image
	# several parameters can be set when calling make_plot() 
	my $genome_png = $oGenome_Image->make_plot;
	my $chrom_png = $oChrom_Image->make_plot(plot_chromosome=>1, scale=>100000);

=head2 CGH Plot Methods

=over

=item B<new()>

Pass the data file and clone file objects at initialisation.

=item B<make_plot()>

Pass hash arguments to C<make_plot()> to set various parameters (see below). The only argument required is C<plot_chromosome>, when creating a single chromosome plot using the C<cgh_plot> class

=item B<set_data()>

The C<data_file> and C<clone_locn_file> objects do not have to be passed at initialisation, but can instead be set using the C<set_data()> method. 

=back

=head2 Plot parameters

The following parameters can be set in the call to C<make_plot()>, or separately before calling C<make_plot()>.

=over

=item B<plot_chromosome>

Set this parameter to indicate which chromosome to plot. Required for single chromosome plots using the C<cgh_plot> class. Must match the chromosome name provided by the clone positions file (or embedded data). 

=item B<has_embedded_locns>

By setting the C<has_embedded_locns> parameter to 1, the module expects the feature name to be in the 'ID' field of the data file (i.e. the C<synonym_id()> of the C<data_file> object) and the clone location to be present in the 'Name' field of the data file (i.e. the C<feature_id()> of the L<C<data_file>|Microarray::File::Data_File> object). The clone location should be of the notation 'chr1:12345..67890'. When using C<has_embedded_locns> a clone position file is not required. Disabled by default.

=item B<plot_centromere>

Set this parameter to zero to disable plotting of the centromere lines. Default is to plot the centromere locations as dashed blue lines. 

=item B<scale>

Pass an integer value to set the desired X-scale of the plot, in bp/pixel. Default for C<cgh_plot> (individual chromosome plot) is 500,000 bp per pixel; default for C<genome_cgh_plot> (whole genome plot) is 2,500,000 bp/pixel. 

=item B<smooth_data_by_location>, B<do_smoothing>

Set either of these parameters to 1 to perform data smoothing. The Log2 ratios in a window of C<$window> bp are averaged, and plotted against the central location of the window. The window moves in steps of C<$step> bp. Disabled by default. 

=item B<smooth_window>, B<smooth_step>

Set the desired window and step sizes for smoothing using these two parameters. A default window size of 500,000bp and step size of 150,000bp provide a moderate level of smoothing, removing outliers while preserving short regions of copy number change. Setting either of these parameters will invoke the smoothing process without setting C<do_smoothing>. 

=item B<shift_zero>

Set this parameter to a value by which all Log2 ratios will be adjusted. Useful to better align the plot with the zero line. 

=item B<flip>

Set this parameter to 1 in order to invert the log ratios returned by the L<C<data_file>|Microarray::File::Data_File> object.  

=back

=head1 SEE ALSO

L<Microarray|Microarray>, L<Microarray::File|Microarray::File>, L<Microarray::File::Data_File|Microarray::File::Data_File>, L<Microarray::File::Clone_Locn_File|Microarray::File::Clone_Locn_File>

=head1 PREREQUISITES

This module utilises the L<GD|GD> module, which requires installation of Thomas Boutell's GD image library (L<http://www.libgd.org>). 

=head1 AUTHOR

James Morris, Translational Research Laboratories, Institute for Women's Health, University College London.

L<http://www.instituteforwomenshealth.ucl.ac.uk/trl>

james.morris@ucl.ac.uk

=head1 COPYRIGHT AND LICENSE

Copyright 2007 by James Morris, University College London

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 

=cut
