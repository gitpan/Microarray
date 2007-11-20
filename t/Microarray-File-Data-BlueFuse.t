#! /usr/bin/perl -w

use strict;

use FindBin;
use Test::More tests=>7;
use Test::Group;
use Test::Differences;
use Test::Deep;

BEGIN {
	use_ok('Microarray::File::Data::BlueFuse');
}

my ($oFile,$oSpot);

my $file = $FindBin::Bin.'/../test_files/bluefuse_output.xls';
begin_skipping_tests "The test-file 'bluefuse_output.xls' could not be found" unless (-e $file);  

test "Object creation" => sub {
	ok($oFile = bluefuse_file->new($file),'object creation');
	isa_ok($oFile,'bluefuse_file','bluefuse_file object');
};

test "Header information" => sub {
	is($oFile->analysis_software,'BlueFuse for Microarrays v3.5(6446)','analysis_software');
	is($oFile->build,'6446','build');
	is($oFile->date,'9/26/2007','date');
	is($oFile->experiment,'new scanner Male vs Female');
	is($oFile->channel1_image_file,'C:\Documents and Settings\Tanya Lebi\Desktop\UCL\2007-09-26_12006_532.tif','channel1_image_file');
	is($oFile->channel2_image_file,'C:\Documents and Settings\Tanya Lebi\Desktop\UCL\2007-09-26_12006_635.tif','channel2_image_file');
	is($oFile->frame_ch1,'1','frame_ch1');
	is($oFile->frame_ch2,'1','frame_ch2');
	is($oFile->gal_file,'C:\Documents and Settings\Tanya Lebi\Desktop\12.gal','gal_file');
	is($oFile->clone_file,'C:\Documents and Settings\Tanya Lebi\Desktop\ParaffinArrays.cfg','clone_file');
	is($oFile->clone_text,'Generated for demo purposes GRS 29/01/07','clone_text');
	is($oFile->replicate_field,'Name','replicate_field');
	test "Confidence flag range" => sub {
		my $hFlags = $oFile->confidence_flag_range;
		eq_or_diff $hFlags,{A=>['0.90','1.00'],B=>['0.70','0.90'],C=>['0.30','0.70'],D=>['0.10','0.30'],E=>['0.00','0.10']},'confidence_flag_range hashref';
		my @aFlags = $oFile->confidence_flag_range;
		eq_or_diff \@aFlags,[('0.00','0.10','0.30','0.70','0.90','1.00')],'confidence_flag_range list';
		my @aRangeA = $oFile->confidence_flag_range('A');
		eq_or_diff \@aRangeA,[('0.90','1.00')],'confidence_flag_range list A';
		my @aRangeB = $oFile->confidence_flag_range('B');
		eq_or_diff \@aRangeB,[('0.70','0.90')],'confidence_flag_range list B';
		my @aRangeC = $oFile->confidence_flag_range('C');
		eq_or_diff \@aRangeC,[('0.30','0.70')],'confidence_flag_range list C';
		my @aRangeD = $oFile->confidence_flag_range('D');
		eq_or_diff \@aRangeD,[('0.10','0.30')],'confidence_flag_range list D';
		my @aRangeE = $oFile->confidence_flag_range('E');
		eq_or_diff \@aRangeE,[('0.00','0.10')],'confidence_flag_range list E';
		my $aRangeA = $oFile->confidence_flag_range('A');
		eq_or_diff $aRangeA,['0.90','1.00'],'confidence_flag_range arrayref A';
		my $aRangeB = $oFile->confidence_flag_range('B');
		eq_or_diff $aRangeB,['0.70','0.90'],'confidence_flag_range arrayref B';
		my $aRangeC = $oFile->confidence_flag_range('C');
		eq_or_diff $aRangeC,['0.30','0.70'],'confidence_flag_range arrayref C';
		my $aRangeD = $oFile->confidence_flag_range('D');
		eq_or_diff $aRangeD,['0.10','0.30'],'confidence_flag_range arrayref D';
		my $aRangeE = $oFile->confidence_flag_range('E');
		eq_or_diff $aRangeE,['0.00','0.10'],'confidence_flag_range arrayref E';
	};
};
test "Array QC" => sub {
	test "Confidence flag percen" => sub {
		my $hFlags = $oFile->confidence_flag_percen;
		eq_or_diff $hFlags,{A=>59.28,B=>18.28,C=>8.44,D=>1.25,E=>12.74},'confidence_flag_percen hashref';
		my @aFlags = $oFile->confidence_flag_percen;
		eq_or_diff \@aFlags,[(59.28,18.28,8.44,1.25,12.74)],'confidence_flag_percen list';
		cmp_ok($oFile->confidence_flag_percen('A'),'==',59.28,'confidence_flag_percen A');
		cmp_ok($oFile->confidence_flag_percen('B'),'==',18.28,'confidence_flag_percen B');
		cmp_ok($oFile->confidence_flag_percen('C'),'==',8.44,'confidence_flag_percen C');
		cmp_ok($oFile->confidence_flag_percen('D'),'==',1.25,'confidence_flag_percen D');
		cmp_ok($oFile->confidence_flag_percen('E'),'==',12.74,'confidence_flag_percen E');
	};
	cmp_ok($oFile->log_ratio_sd,'==',0.63,'log_ratio_sd');
	cmp_ok($oFile->rep_median_sd,'==',0.26,'rep_median_sd');
	cmp_ok($oFile->mean_ch1_amp,'==',2044.17,'mean_ch1_amp');
	cmp_ok($oFile->mean_ch2_amp,'==',1654.71,'mean_ch2_amp');
	cmp_ok($oFile->sbr_ch1,'==',8.05,'sbr_ch1');
	cmp_ok($oFile->sbr_ch2,'==',13.52,'sbr_ch2');
	eq_or_diff $oFile->bad_flags,{ 'C'=>'1', 'D'=>'1', 'E'=>'1' },'bad_flags';
};
test "Data" => sub {
	test "Spot 1" => sub {
		cmp_ok($oFile->block_row(0),'==',1,'block_row');
		cmp_ok($oFile->block_col(0),'==',1,'block_col');
		cmp_ok($oFile->spot_row(0),'==',1,'spot_row');
		cmp_ok($oFile->spot_col(0),'==',1,'spot_col');
		cmp_ok($oFile->spot_index(0),'==',1,'spot_index');
		is($oFile->feature_id(0),'RP11-1137J16','feature_id');
		is($oFile->synonym_id(0),'RP11-1137J16','synonym_id');
		cmp_ok($oFile->confidence(0),'==',0.01,'confidence');
		is($oFile->flag_id(0),'E','flag_id');
		is($oFile->man_excl(0),'no','man_excl');
		is($oFile->auto_excl(0),'yes','auto_excl');
		cmp_ok($oFile->channel1_signal(0),'==',436.126,'channel1_signal');
		cmp_ok($oFile->channel2_signal(0),'==',114.35,'channel2_signal');
		cmp_ok($oFile->ratio_ch1ch2(0),'==',3.814,'ratio_ch1ch2');
		cmp_ok($oFile->log2ratio_ch1ch2(0),'==',1.931,'log2ratio_ch1ch2');
		cmp_ok($oFile->log10ratio_ch1ch2(0),'==',0.581,'log10ratio_ch1ch2');
		cmp_ok($oFile->ratio_ch2ch1(0),'==',0.262,'ratio_ch2ch1');
		cmp_ok($oFile->log2ratio_ch2ch1(0),'==',-1.931,'log2ratio_ch2ch1');
		cmp_ok($oFile->log10ratio_ch2ch1(0),'==',-0.581,'log10ratio_ch2ch1');
		cmp_ok($oFile->sumch1ch2(0),'==',550.476,'sumch1ch2');
		cmp_ok($oFile->log2sum(0),'==',9.105,'log2sum');
		cmp_ok($oFile->log10sum(0),'==',2.741,'log10sum');
		cmp_ok($oFile->product(0),'==',49871.043,'product');
		cmp_ok($oFile->log2product(0),'==',15.606,'log2product');
		cmp_ok($oFile->log10product(0),'==',4.698,'log10product');
		cmp_ok($oFile->y_pos(0),'==',138,'y_pos');
		cmp_ok($oFile->x_pos(0),'==',128,'x_pos');
		cmp_ok($oFile->channel1_quality(0),'==',0.15,'channel1_quality');
		cmp_ok($oFile->channel2_quality(0),'==',0.13,'channel2_quality');
		cmp_ok($oFile->spot_diameter(0),'==',6.38,'spot_diameter');
		cmp_ok($oFile->uniformity(0),'==',0.29,'uniformity');
		cmp_ok($oFile->circularity(0),'==',0.54,'circularity');
		cmp_ok($oFile->grid_offset(0),'==',20.1,'grid_offset');
		cmp_ok($oFile->quality(0),'==',0,'quality');
		is($oFile->chromosome(0),'7','chromosome');
		cmp_ok($oFile->position(0),'==',152034734,'position');
		is($oFile->cyto_locn(0),undef,'cyto_locn');
		is($oFile->display(0),undef,'display');
		is($oFile->omim(0),undef,'omim');
		is($oFile->disease(0),undef,'disease');
		is($oFile->gc_content(0),undef,'gc_content');
	};
	test "Spot 14288" => sub {
		cmp_ok($oFile->block_row(14287),'==',6,'block_row');
		cmp_ok($oFile->block_col(14287),'==',1,'block_col');
		cmp_ok($oFile->spot_row(14287),'==',10,'spot_row');
		cmp_ok($oFile->spot_col(14287),'==',14,'spot_col');
		cmp_ok($oFile->spot_index(14287),'==',14288,'spot_index');
		is($oFile->feature_id(14287),'RP11-683G6','feature_id');
		is($oFile->synonym_id(14287),'RP11-683G6','synonym_id');
		cmp_ok($oFile->confidence(14287),'==',0.97,'confidence');
		is($oFile->flag_id(14287),'A','flag_id');
		is($oFile->man_excl(14287),'no','man_excl');
		is($oFile->auto_excl(14287),'no','auto_excl');
		cmp_ok($oFile->channel1_signal(14287),'==',3625.97,'channel1_signal');
		cmp_ok($oFile->channel2_signal(14287),'==',2938.945,'channel2_signal');
		cmp_ok($oFile->ratio_ch1ch2(14287),'==',1.234,'ratio_ch1ch2');
		cmp_ok($oFile->log2ratio_ch1ch2(14287),'==',0.303,'log2ratio_ch1ch2');
		cmp_ok($oFile->log10ratio_ch1ch2(14287),'==',0.091,'log10ratio_ch1ch2');
		cmp_ok($oFile->ratio_ch2ch1(14287),'==',0.811,'ratio_ch2ch1');
		cmp_ok($oFile->log2ratio_ch2ch1(14287),'==',-0.303,'log2ratio_ch2ch1');
		cmp_ok($oFile->log10ratio_ch2ch1(14287),'==',-0.091,'log10ratio_ch2ch1');
		cmp_ok($oFile->sumch1ch2(14287),'==',6564.915,'sumch1ch2');
		cmp_ok($oFile->log2sum(14287),'==',12.681,'log2sum');
		cmp_ok($oFile->log10sum(14287),'==',3.817,'log10sum');
		cmp_ok($oFile->product(14287),'==',10656525,'product');
		cmp_ok($oFile->log2product(14287),'==',23.345,'log2product');
		cmp_ok($oFile->log10product(14287),'==',7.028,'log10product');
		cmp_ok($oFile->y_pos(14287),'==',4930,'y_pos');
		cmp_ok($oFile->x_pos(14287),'==',548,'x_pos');
		cmp_ok($oFile->channel1_quality(14287),'==',1,'channel1_quality');
		cmp_ok($oFile->channel2_quality(14287),'==',1,'channel2_quality');
		cmp_ok($oFile->spot_diameter(14287),'==',22.34,'spot_diameter');
		cmp_ok($oFile->uniformity(14287),'==',0.62,'uniformity');
		cmp_ok($oFile->circularity(14287),'==',0.74,'circularity');
		cmp_ok($oFile->grid_offset(14287),'==',2,'grid_offset');
		cmp_ok($oFile->quality(14287),'==',1,'quality');
		is($oFile->chromosome(14287),'17','chromosome');
		cmp_ok($oFile->position(14287),'==',44271611,'position');
		is($oFile->cyto_locn(14287),undef,'cyto_locn');
		is($oFile->display(14287),undef,'display');
		is($oFile->omim(14287),undef,'omim');
		is($oFile->disease(14287),undef,'disease');
		is($oFile->gc_content(14287),undef,'gc_content');
	};
	test "Spot 33696" => sub {
		cmp_ok($oFile->block_row(33695),'==',12,'block_row');
		cmp_ok($oFile->block_col(33695),'==',4,'block_col');
		cmp_ok($oFile->spot_row(33695),'==',27,'spot_row');
		cmp_ok($oFile->spot_col(33695),'==',26,'spot_col');
		cmp_ok($oFile->spot_index(33695),'==',33696,'spot_index');
		is($oFile->feature_id(33695),'','feature_id');
		is($oFile->synonym_id(33695),'','synonym_id');
		cmp_ok($oFile->confidence(33695),'==',0.01,'confidence');
		is($oFile->flag_id(33695),'E','flag_id');
		is($oFile->man_excl(33695),'no','man_excl');
		is($oFile->auto_excl(33695),'yes','auto_excl');
		cmp_ok($oFile->channel1_signal(33695),'==',13.687,'channel1_signal');
		cmp_ok($oFile->channel2_signal(33695),'==',230.077,'channel2_signal');
		cmp_ok($oFile->ratio_ch1ch2(33695),'==',0.059,'ratio_ch1ch2');
		cmp_ok($oFile->log2ratio_ch1ch2(33695),'==',-4.071,'log2ratio_ch1ch2');
		cmp_ok($oFile->log10ratio_ch1ch2(33695),'==',-1.226,'log10ratio_ch1ch2');
		cmp_ok($oFile->ratio_ch2ch1(33695),'==',16.81,'ratio_ch2ch1');
		cmp_ok($oFile->log2ratio_ch2ch1(33695),'==',4.071,'log2ratio_ch2ch1');
		cmp_ok($oFile->log10ratio_ch2ch1(33695),'==',1.226,'log10ratio_ch2ch1');
		cmp_ok($oFile->sumch1ch2(33695),'==',243.764,'sumch1ch2');
		cmp_ok($oFile->log2sum(33695),'==',7.929,'log2sum');
		cmp_ok($oFile->log10sum(33695),'==',2.387,'log10sum');
		cmp_ok($oFile->product(33695),'==',3149.083,'product');
		cmp_ok($oFile->log2product(33695),'==',11.621,'log2product');
		cmp_ok($oFile->log10product(33695),'==',3.498,'log10product');
		cmp_ok($oFile->y_pos(33695),'==',10862,'y_pos');
		cmp_ok($oFile->x_pos(33695),'==',3656,'x_pos');
		cmp_ok($oFile->channel1_quality(33695),'==',0,'channel1_quality');
		cmp_ok($oFile->channel2_quality(33695),'==',0.15,'channel2_quality');
		cmp_ok($oFile->spot_diameter(33695),'==',5.04,'spot_diameter');
		cmp_ok($oFile->uniformity(33695),'==',0.2,'uniformity');
		cmp_ok($oFile->circularity(33695),'==',0.37,'circularity');
		cmp_ok($oFile->grid_offset(33695),'==',10,'grid_offset');
		cmp_ok($oFile->quality(33695),'==',0,'quality');
		is($oFile->chromosome(33695),undef,'chromosome');
		is($oFile->position(33695),undef,'position');
		is($oFile->cyto_locn(33695),undef,'cyto_locn');
		is($oFile->display(33695),undef,'display');
		is($oFile->omim(33695),undef,'omim');
		is($oFile->disease(33695),undef,'disease');
		is($oFile->gc_content(33695),undef,'gc_content');
	};
	test "ch1_mean_f and ch1_median_b" => sub {
		cmp_ok($oFile->ch1_mean_f(1),'==',$oFile->channel1_signal(1),'ch1_mean_f');
		cmp_ok($oFile->ch2_mean_f(1),'==',$oFile->channel2_signal(1),'ch2_mean_f');
		cmp_ok($oFile->ch1_median_b(1),'==',1,'ch1_mean_b');
		cmp_ok($oFile->ch2_median_b(1),'==',1,'ch2_mean_b');
	};
};
test "Spot test" => sub {
	ok($oSpot = $oFile->spot_object(1),'get spot 1');
	cmp_ok($oSpot->block_row,'==',1,'block_row');
	cmp_ok($oSpot->block_col,'==',1,'block_col');
	cmp_ok($oSpot->spot_row,'==',1,'spot_row');
	cmp_ok($oSpot->spot_col,'==',1,'spot_col');
	cmp_ok($oSpot->spot_index,'==',1,'spot_index');
	is($oSpot->feature_id,'RP11-1137J16','feature_id');
	is($oSpot->synonym_id,'RP11-1137J16','synonym_id');
	cmp_ok($oSpot->channel1_signal,'==',436.126,'channel1_signal');
	cmp_ok($oSpot->channel2_signal,'==',114.35,'channel2_signal');
	cmp_ok($oSpot->ch1_mean_f,'==',436.126,'ch1_mean_f');
	cmp_ok($oSpot->ch2_mean_f,'==',114.35,'ch2_mean_f');
	cmp_ok($oSpot->ch1_median_b,'==',1,'ch1_median_b');
	cmp_ok($oSpot->ch2_median_b,'==',1,'ch2_median_b');
	cmp_ok($oSpot->channel1_quality,'==',0.15,'channel1_quality');
	cmp_ok($oSpot->channel2_quality,'==',0.13,'channel2_quality');
	cmp_ok($oSpot->x_pos,'==',128,'x_pos');
	cmp_ok($oSpot->y_pos,'==',138,'y_pos');
	cmp_ok($oSpot->spot_diameter,'==',6.38,'spot_diameter');
	is($oSpot->flag_id,'E','flag_id');
};
test "Array Layout" => sub {
	cmp_ok($oFile->array_rows,'==',12,'array_rows');
	cmp_ok($oFile->array_columns,'==',4,'array_columns');
	cmp_ok($oFile->spot_rows,'==',27,'spot_rows');
	cmp_ok($oFile->spot_columns,'==',26,'spot_columns');
};

end_skipping_tests;
