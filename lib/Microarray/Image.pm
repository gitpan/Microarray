package Microarray::Image;

use 5.006;
use strict;
use warnings;
our $VERSION = '0.117';

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
	sub default_scale {
		25; 
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
	sub default_scale {
		25; 
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
	sub default_scale {
		200; 
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
		140;
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
  
  	use Microarray::Analysis::CGH;
  	
	our @ISA = qw( plot chromosome_cgh );

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
			$self->sort_chromosome_data;
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
	sub y_scale {
		1;
	}
	# normalise the processed data relative to the image ready to be plotted
	sub plot_data {
		my $self = shift;
		my $aaX = $self->x_values;
		my $aaY = $self->y_values;
		my $aX  = $$aaX[0];
		my $aY  = $$aaY[0];
		my $scale    = $self->scale;
		my $y_scale    = $self->y_scale;
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
			push(@$aLog2, ($middle - ($log2 * (450/3)))*$y_scale );
		}
		$self->{ _plotx_values } = $aLocns;
		$self->{ _ploty_values } = $aLog2;
		return($aLocns,$aLog2);
	}
	sub plot_dimensions {
		my $self = shift;
		my $scale = $self->scale;
		my $y_scale = $self->y_scale;

		my $chr_length = $self->chr_length($self->plot_chromosome);
		my $x = int(($chr_length/$scale) + 1);
		my $y = 450 * $y_scale;
		my $x_margin = 0;
		my $y_margin = 0;
		$self->{ _x_length } = $x;
		$self->{ _y_length } = $y;
		$self->{ _x_margin } = 0;
		$self->{ _y_margin } = 0;
		$self->{ _middle } = $y/2;
		
		return ($x,$y);
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
		my $y_scale = $self->y_scale;
		my $chr = $self->plot_chromosome;
		# get colours from the GD colour table 
		my $black 	= $image->colorExact(0,0,0);
		my $red   	= $image->colorExact(255,0,0);      
		my $green 	= $image->colorExact(0,255,0);      
		my $blue	= $image->colorExact(125,125,255);   
		# 3px wide log2 ratio lines
		$image->filledRectangle(0,(150*$y_scale),3080,(150*$y_scale),$red);		# +0.5
		$image->filledRectangle(0,(225*$y_scale),3080,(225*$y_scale),$green);	#  0.0
		$image->filledRectangle(0,(300*$y_scale),3080,(300*$y_scale),$red);		# -0.5
		# axis labels
		$image->string($image->gdGiantFont,10,(150*$y_scale),'0.5',$black);
		$image->string($image->gdGiantFont,10,(225*$y_scale),'0',$black);
		$image->string($image->gdGiantFont,10,(300*$y_scale),'-0.5',$black);
		
		if ($self->plot_centromere){
			# dashed style for centromere lines
			$image->setStyle($blue,$blue,$blue,$blue,gdTransparent,gdTransparent);
			my $cen = int($self->chr_centromere($chr)/$scale);
			$image->line($cen,0,$cen,$y,gdStyled);
		}
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
		my $y_scale = $self->y_scale;
		# get colours from the GD colour table
		my $red    = $image->colorExact(255,0,0);      
		my $green  = $image->colorExact(0,255,0);      
		my $yellow = $image->colorExact(255,255,0);      
		my $colour;
		if ($ratio <= (150*$y_scale)){
			$colour = $red;
		} elsif ($ratio >= (300*$y_scale)){
			$colour = $green;
		} elsif (((262*$y_scale) > $ratio)&&($ratio > (187*$y_scale))) {
			$colour = $yellow;
		} elsif ($ratio >= (262*$y_scale)) {
			my $red_hue = int(255-(6.8*($ratio-(262*$y_scale))));				# factorial = 255/(low_yellow - green)
			$colour = $image->colorClosest($red_hue,255,0);		# reducing red, closer to green
		} else {
			my $green_hue = int(255-(6.9*((187*$y_scale)-$ratio)));				# factorial = 255/(high_yellow - red)
			$colour = $image->colorClosest(255,$green_hue,0);	# reducing green, closer to red
		}
		return($colour);
	}
	sub image_map {
		my $self 		= shift;
		my ($aX,$aY) 	= $self->plot_data;
		
		my ($aFeatures,$feature);
		if ($self->smoothing) {
			$aFeatures = $self->x_values;		# genomic location of window centre
			$feature = 'window';
		} else {
			$aFeatures = $self->feature_names;	# bac names
			$feature = 'clone';
		}
		
		my $map_name = 'cgh_plot';
		if (@_){
			$map_name = shift;
		}
		my $map_string = "<MAP NAME=\"$map_name\">\n";
		
		for (my $i=0; $i<@$aX; $i++){
			my $x = $aX->[$i];
			my $y = int($aY->[$i]) + 225;
			next unless (($x && $y)&&($y>0)&&($x>0));
			my $element = $aFeatures->[0][$i];
			$map_string .= "<area \"[$feature]=[$element]\" shape=\"circle\" coords=\"$x,$y,3\" />\n";
		}
		$map_string .= "</MAP>\n";
		return $map_string;
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
	sub chr_length {
		my $self = shift;
		my $chr  = shift;
		my %hChromosome = (		
			# chr length
			1 => 247249719,
			2 => 242951149,
			3 => 199501827,
			4 => 191273063,
			5 => 180857866,
			6 => 170899992,
			7 => 158821424,
			8 => 146274826,
			9 => 140273252,
			10 => 135374737,
			11 => 134452384,
			12 => 132349534,
			13 => 114142980,
			14 => 106368585,
			15 => 100338915,
			16 => 88827254,
			17 => 78774742,
			18 => 76117153,
			19 => 63811651,
			20 => 62435964,
			21 => 46944323,
			22 => 49691432,
			23 => 154913754,
			'X' => 154913754,
			24 => 57772954,
			'Y' => 57772954,
			25 => 3080419480	# END
		);
		return ($hChromosome{$chr});
	}
}
{ package genome_cgh_plot;

	# must try genome_cgh first, otherwise will move through cgh_plot to chromosome_cgh
	our @ISA = qw( genome_cgh cgh_plot );	

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
			$self->sort_genome_data;
		}
		$self->{ _data_sorted } = 1;
	}
	sub make_plot {
		my $self  = shift;
		if (@_){
			$self->parse_args(@_);
		}
		## chr 25 is the total length of the genome plus 25 margin
		my $x_axis = ($self->chr_offset(25)/$self->scale)+25;
		$self->{ _gd_object } = GD::Image->new($x_axis,(450*$self->y_scale));	
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
	sub plot_dimensions {
		my $self = shift;
		my $scale = $self->scale;
		my $y_scale = $self->y_scale;

		my $genome_length = $self->chr_length(25);	# end
		my $x = int(($genome_length/$scale) + 1);
		my $y = 450 * $y_scale;
		my $x_margin = 0;
		my $y_margin = 0;
		$self->{ _x_length } = $x;
		$self->{ _y_length } = $y;
		$self->{ _x_margin } = 0;
		$self->{ _y_margin } = 0;
		$self->{ _middle } = $y/2;
		
		return ($x,$y);
	}
	sub default_scale {
		## set scale for the genome plot to 2.5Mb per pixel
		2500000;
	}
	sub y_scale {
		my $self = shift;
		my $scale = $self->scale;
		my $default_scale = $self->default_scale;
		return $default_scale/$scale;
	}
	sub plot_chromosome {
		return;
	}
	sub plot_data {
		my $self = shift;
		my $scale = $self->scale;
		my $y_scale = $self->y_scale;
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
				push(@$aLog2, (225 - ($log2 * (450/3)))*$y_scale );
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
		my $y_scale = $self->y_scale;
		# get colours from the GD colour table 
		my $black 	= $image->colorExact(0,0,0);
		my $red   	= $image->colorExact(255,0,0);      
		my $green 	= $image->colorExact(0,255,0); 
		my $blue	= $image->colorExact(125,125,255);   
		# 3px wide log2 ratio lines
		$image->filledRectangle(0,(150*$y_scale),3080,(150*$y_scale),$red);		# +0.5
		$image->filledRectangle(0,(225*$y_scale),3080,(225*$y_scale),$green);		#  0.0
		$image->filledRectangle(0,(300*$y_scale),3080,(300*$y_scale),$red);		# -0.5
		# axis labels
		$image->string($image->gdSmallFont,0,(150*$y_scale),'0.5',$black);
		$image->string($image->gdSmallFont,0,(225*$y_scale),'0',$black);
		$image->string($image->gdSmallFont,0,(300*$y_scale),'-0.5',$black);
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
				$image->line($cen+25,0,$cen+25,(450*$y_scale),gdStyled);
			}
			# chr buffer
			$image->filledRectangle($start+25,0,$start+25,(450*$y_scale),$black);
			# set chr names
			my $chr_name;
			if ($chr == 23){
				$chr_name = 'X';
			} elsif ($chr == 24){
				$chr_name = 'Y';
				# end line
				$image->line($end+25,0,$end+25,(450*$y_scale),$black);
			} else {
				$chr_name = $chr;
			}
			# print chr name at bottom of plot
			$image->string($image->gdSmallFont,$middle+20,(425*$y_scale),$chr_name,$black);
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

Pass the L<Microarray::File::Data|Microarray::File::Data> and (optional) L<Microarray::File::Clone_Locns|Microarray::File::Clone_Locns> file objects at initialisation.

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

=item B<plot_centromere>

Set this parameter to zero to disable plotting of the centromere lines. Default is to plot the centromere locations as dashed blue lines. 

=item B<scale>

Pass an integer value to set the desired X-scale of the plot, in bp/pixel. Default for C<cgh_plot> (individual chromosome plot) is 500,000 bp per pixel; default for C<genome_cgh_plot> (whole genome plot) is 2,500,000 bp/pixel. 

=item B<shift_zero>

Set this parameter to a value by which all Log2 ratios will be adjusted. Useful to better align the plot with the zero line. 

=back

=head2 Other analysis methods

The cgh_plot and genome_cgh_plot classes can use methods from the L<Microarray::Analysis::CGH|Microarray::Analysis::CGH> module. Pass analysis parameters to the make_plot() method to implement L<flip()|Microarray::Analysis::CGH/"flip">, L<has_embedded_locns()|Microarray::Analysis::CGH/"has_embedded_locns">, L<do_smoothing()|Microarray::Analysis::CGH/"do_smoothing"> etc. 

=head1 SEE ALSO

L<Microarray|Microarray>, L<Microarray::Analysis|Microarray::Analysis>, L<Microarray::Analysis::CGH|Microarray::Analysis::CGH>, L<Microarray::File|Microarray::File>, L<Microarray::File::Data_File|Microarray::File::Data_File>, L<Microarray::File::Clone_Locn_File|Microarray::File::Clone_Locn_File>

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
