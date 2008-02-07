#!/usr/bin/perl -w

use strict;

use FindBin;
use Test::More tests=>3;
use Test::Group;
use Test::Differences;
use Test::Deep;

BEGIN {
	use_ok('Microarray::Spot');
}

my ($oSpot);

test "Object creation" => sub {
	ok($oSpot = array_spot->new(1),'object creation');
	isa_ok($oSpot,'array_spot','array_spot object');
};

my @aMethods = qw( 	spot_index block_row block_col spot_row spot_col 
					x_pos y_pos spot_diameter feature_id synonym_id 
					spot_pixels bg_pixels footprint flag_id ch1_mean_f 
					ch1_median_f ch1_sd_f ch1_mean_b ch1_median_b 
					ch1_sd_b ch1_b1sd channel1_quality channel1_sat 
					ch2_mean_f ch2_median_f ch2_sd_f ch2_mean_b ch2_median_b 
					ch2_sd_b ch2_b1sd channel2_quality channel2_sat spot_status );
my @aValues = (4,2,3,1,4,2588,6736,100,'RP13-486L5','none',280,436,7,3,3306,2071,3213.69,1143,883,859.79,57.5,43.6,0,3181,1891,3237.4,790,658,603.81,61.4,50.4,0,1);
for (my $i=0; $i<@aMethods; $i++){
	my $method = $aMethods[$i];
	my $value = $aValues[$i];
	$oSpot->$method($value);
}

test "checking parameters" => sub {
	cmp_ok($oSpot->spot_index ,'==',4,'spot_index');
	cmp_ok($oSpot->block_row ,'==',2,'block_row');
	cmp_ok($oSpot->block_col ,'==',3,'block_col');
	cmp_ok($oSpot->spot_row ,'==',1,'spot_row');
	cmp_ok($oSpot->spot_col ,'==',4,'spot_col');
	cmp_ok($oSpot->x_pos ,'==',2588,'x_pos');
	cmp_ok($oSpot->y_pos ,'==',6736,'y_pos');
	cmp_ok($oSpot->spot_diameter ,'==',100,'spot_diameter');
	cmp_ok($oSpot->feature_id ,'eq','RP13-486L5','feature_id');
	cmp_ok($oSpot->synonym_id ,'eq','none','synonym_id');
	cmp_ok($oSpot->spot_pixels ,'==',280,'spot_pixels');
	cmp_ok($oSpot->bg_pixels ,'==',436,'bg_pixels');
	cmp_ok($oSpot->footprint ,'==',7,'footprint');
	cmp_ok($oSpot->flag_id ,'==',3,'flag_id');
	cmp_ok($oSpot->ch1_mean_f ,'==',3306,'ch1_mean_f');
	cmp_ok($oSpot->ch1_median_f ,'==',2071,'ch1_median_f');
	cmp_ok($oSpot->ch1_sd_f ,'==',3213.69,'ch1_sd_f');
	cmp_ok($oSpot->ch1_mean_b ,'==',1143,'ch1_mean_b');
	cmp_ok($oSpot->ch1_median_b ,'==',883,'ch1_median_b');
	cmp_ok($oSpot->ch1_sd_b ,'==',859.79,'ch1_sd_b');
	cmp_ok($oSpot->ch1_b1sd ,'==',57.5,'ch1_b1sd');
	cmp_ok($oSpot->channel1_quality ,'==',43.6,'channel1_quality');
	cmp_ok($oSpot->channel1_sat ,'==',0,'channel1_sat');
	cmp_ok($oSpot->ch2_mean_f ,'==',3181,'ch2_mean_f');
	cmp_ok($oSpot->ch2_median_f ,'==',1891,'ch2_median_f');
	cmp_ok($oSpot->ch2_sd_f ,'==',3237.4,'ch2_sd_f');
	cmp_ok($oSpot->ch2_mean_b ,'==',790,'ch2_mean_b');
	cmp_ok($oSpot->ch2_median_b ,'==',658,'ch2_median_b');
	cmp_ok($oSpot->ch2_sd_b ,'==',603.81,'ch2_sd_b');
	cmp_ok($oSpot->ch2_b1sd ,'==',61.4,'ch2_b1sd');
	cmp_ok($oSpot->channel2_quality ,'==',50.4,'channel2_quality');
	cmp_ok($oSpot->channel2_sat ,'==',0,'channel2_sat');
	cmp_ok($oSpot->spot_status ,'==',1,'spot_status');
};






