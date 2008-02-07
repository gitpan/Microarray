#!/usr/bin/perl -w

use strict;

use FindBin;
use Test::More tests=>21;
use Test::Group;
use Test::Differences;
use Test::Deep;

BEGIN {
	use_ok('Microarray');
	use_ok('Microarray::Spot');
	use_ok('Microarray::Reporter');
	use_ok('Microarray::Analysis');
	use_ok('Microarray::Analysis::CGH');
	use_ok('Microarray::Image');
	use_ok('Microarray::File');
	use_ok('Microarray::File::Image');
	use_ok('Microarray::File::Clone_Locns');
	use_ok('Microarray::File::Data::GenePix');
	use_ok('Microarray::File::Data::Agilent');
	use_ok('Microarray::File::Data::Manor_Output');
}

my ($oArray,$oReporter,$oReporter1,$oReporter2,$oReporter3,$oReporter4,$oReporter5,$aReporter_Objects,$aReporter_Names,$hReporters,$file,$aSpots);

$file = $FindBin::Bin.'/../test_files/quantarray.csv';
begin_skipping_tests "The test-file 'quantarray.csv' could not be found" unless (-e $file);  

test "Object creation" => sub {
	ok($oArray = microarray->new('10001',$file),'object creation');
	isa_ok($oArray,'microarray','microarray object');
};

$oArray->set_param(min_snr=>2,high_signal=>65000,low_signal=>500,signal_quality=>50);

test "getting and setting params for spot quality assessment criteria" => sub {
    is($oArray->low_signal,'500','low_signal');
    is($oArray->high_signal,'65000','high_signal');
    is($oArray->percen_sat,'10','percen_sat');
    is($oArray->min_snr,'2','min_snr');
    is($oArray->signal_quality,'50','signal_quality');
    is($oArray->min_diameter,'80','min_diameter');
    is($oArray->max_diameter,'150','max_diameter');
    is($oArray->target_diameter,'100','target_diameter'); 
    is($oArray->max_diameter_deviation,'10','max_diameter_deviation');
    is($oArray->normalisation,'yes','normalisation');
    is($oArray->signal_normalisation,'yes','signal_normalisation');
    is($oArray->genetic_data_source,'data_file','genetic_data_source');
};

$oArray->set_reporter_data;

test "get reporter data" => sub {
	ok($oReporter = $oArray->get_reporter('RP13-827M24'),'get_reporter'); 
	ok($aReporter_Objects = $oArray->get_reporter_objects,'get_reporter_objects'); 
	ok($aReporter_Names = $oArray->get_reporter_ids,'get_reporter_ids'); 
	ok($hReporters = $oArray->get_all_reporters,'get_all_reporters'); 
	eq_or_diff $hReporters->{ 'RP13-827M24' }, $oReporter, "comparing reporter objects";
};

test "check reporter objects" => sub {
	ok($oReporter = $oArray->get_reporter('RP13-827M24'),'get_reporter');  # returns a single reporter object
	isa_ok($oReporter, 'array_reporter','get_reporter object');
	ok($aReporter_Objects = $oArray->get_reporter_objects,'get_reporter_objects');  # returns a list of reporter objects
	for my $oRep (@$aReporter_Objects){
		isa_ok($oRep, 'array_reporter','get_reporter_objects objects'); 
	}	
	ok($aReporter_Names = $oArray->get_reporter_ids,'get_reporter_ids');  # returns a list of reporter ids
	for my $id (@$aReporter_Names){
		ok(defined($id), 'get_reporter_ids defined'); 
		ok($id ne '', 'get_reporter_ids empty string'); 
	}	
	ok($hReporters = $oArray->get_all_reporters,'get_all_reporters');  # returns a hash of reporters; key=reporter_id, value=reporter object
	while( my($id, $oRep) = each(%$hReporters)) {
		ok(defined($id), 'get_all_reporters defined'); 
		ok($id ne '', 'get_all_reporters empty string');
		isa_ok($oRep, 'array_reporter','get_all_reporters object');  
    }
};

test "CTD-2023C19" => sub {
	my $oCTD_2023C19;
    ok($oCTD_2023C19 = $oArray->get_reporter('CTD-2023C19'),'get_reporter CTD-2023C19');
	cmp_ok($oCTD_2023C19->reporter_id,'eq','CTD-2023C19','reporter_id');
    cmp_ok($oCTD_2023C19->get_reporter_replicates,'==',3,'get_reporter_replicates');
	cmp_ok($oCTD_2023C19->spots_passed_qc,'==',2,'CTD-2023C19 spots passed QC');
    cmp_ok($oCTD_2023C19->mean_ch1,'==',4490.5,'mean_ch1');
    cmp_ok($oCTD_2023C19->mean_ch2,'==',5698.5,'mean_ch2');
    cmp_ok($oCTD_2023C19->mean_ratios,'eq',0.788117500091488,'mean_ratios');
    cmp_ok($oCTD_2023C19->ratio_means,'eq',0.788014389751689,'ratio means');
    
    my $aSpots = $oCTD_2023C19->get_reporter_spots;
    is(@$aSpots,3,'3 spots');
    my ($oSpot1,$oSpot2,$oSpot3) = @$aSpots;
	isa_ok($oSpot1,'array_spot','spot object 1');
	isa_ok($oSpot2,'array_spot','spot object 2');
	isa_ok($oSpot3,'array_spot','spot object 3');
    cmp_ok($oSpot1->spot_index,'==',21438,'spot1 index');
    cmp_ok($oSpot2->spot_index,'==',21621,'spot2 index');
    cmp_ok($oSpot3->spot_index,'==',21805,'spot3 index');
    cmp_ok($oSpot1->spot_status,'==',0,'spot1 failed QC');
    cmp_ok($oSpot2->spot_status,'==',1,'spot2 passed QC');
    cmp_ok($oSpot3->spot_status,'==',1,'spot3 passed QC');
    cmp_ok($oSpot1->channel1_signal,'==',1964,'spot1 ch1 signal');
    cmp_ok($oSpot2->channel1_signal,'==',4533,'spot2 ch1 signal');
    cmp_ok($oSpot3->channel1_signal,'==',4448,'spot3 ch1 signal');
    cmp_ok($oSpot1->channel2_signal,'==',2382,'spot1 ch2 signal');
    cmp_ok($oSpot2->channel2_signal,'==',5796,'spot2 ch2 signal');
    cmp_ok($oSpot3->channel2_signal,'==',5601,'spot3 ch2 signal');
};
test "RP11-40I8" => sub {
	my $oRP11_40I8;
    ok($oRP11_40I8 = $oArray->get_reporter('RP11-40I8'),'get_reporter RP11-40I8');
	cmp_ok($oRP11_40I8->reporter_id,'eq','RP11-40I8','reporter_id');
    cmp_ok($oRP11_40I8->get_reporter_replicates,'==',3,'get_reporter_replicates');
	is($oRP11_40I8->spots_passed_qc,undef,'RP11-40I8 spots passed QC');
    cmp_ok($oRP11_40I8->mean_ch1,'==',0,'mean_ch1');
    cmp_ok($oRP11_40I8->mean_ch2,'==',0,'mean_ch2');
    cmp_ok($oRP11_40I8->mean_ratios,'==',0,'mean_ratios');
    is($oRP11_40I8->ratio_means,undef,'ratio means');

    my $aSpots = $oRP11_40I8->get_reporter_spots;
    is(@$aSpots,3,'3 spots');
    my ($oSpot1,$oSpot2,$oSpot3) = @$aSpots;
	isa_ok($oSpot1,'array_spot','spot object 1');
	isa_ok($oSpot2,'array_spot','spot object 2');
	isa_ok($oSpot3,'array_spot','spot object 3');
    cmp_ok($oSpot1->spot_index,'==',1785,'spot1 index');
    cmp_ok($oSpot2->spot_index,'==',1969,'spot2 index');
    cmp_ok($oSpot3->spot_index,'==',2152,'spot3 index');
    cmp_ok($oSpot1->spot_status,'==',0,'spot1 failed QC');
    cmp_ok($oSpot2->spot_status,'==',0,'spot2 failed QC');
    cmp_ok($oSpot3->spot_status,'==',0,'spot3 failed QC');
    cmp_ok($oSpot1->channel1_signal,'==',1203,'spot1 ch1 signal');
    cmp_ok($oSpot2->channel1_signal,'==',2657,'spot2 ch1 signal');
    cmp_ok($oSpot3->channel1_signal,'==',576,'spot3 ch1 signal');
    cmp_ok($oSpot1->channel2_signal,'==',1448,'spot1 ch2 signal');
    cmp_ok($oSpot2->channel2_signal,'==',2362,'spot2 ch2 signal');
    cmp_ok($oSpot3->channel2_signal,'==',895,'spot3 ch2 signal');
};
test "RP11-558G22" => sub {
	my $oRP11_558G22;
    ok($oRP11_558G22 = $oArray->get_reporter('RP11-558G22'),'get_reporter RP11-558G22');
	cmp_ok($oRP11_558G22->reporter_id,'eq','RP11-558G22','reporter_id');
    cmp_ok($oRP11_558G22->get_reporter_replicates,'==',3,'get_reporter_replicates');
	cmp_ok($oRP11_558G22->spots_passed_qc,'==',3,'RP11-558G22 spots passed QC');
    cmp_ok($oRP11_558G22->mean_ch1,'==',8983,'mean_ch1');
    cmp_ok($oRP11_558G22->mean_ch2,'eq',11071.6666666667,'mean_ch2');
    cmp_ok($oRP11_558G22->mean_ratios,'eq',0.811711465229086,'mean_ratios');
    cmp_ok($oRP11_558G22->ratio_means,'eq',0.811350293542074,'ratio means');

    my $aSpots = $oRP11_558G22->get_reporter_spots;
    is(@$aSpots,3,'3 spots');
    my ($oSpot1,$oSpot2,$oSpot3) = @$aSpots;
	isa_ok($oSpot1,'array_spot','spot object 1');
	isa_ok($oSpot2,'array_spot','spot object 2');
	isa_ok($oSpot3,'array_spot','spot object 3');
    cmp_ok($oSpot1->spot_index,'==',9952,'spot1 index');
    cmp_ok($oSpot2->spot_index,'==',10136,'spot2 index');
    cmp_ok($oSpot3->spot_index,'==',10319,'spot3 index');
    cmp_ok($oSpot1->spot_status,'==',1,'spot1 passed QC');
    cmp_ok($oSpot2->spot_status,'==',1,'spot2 passed QC');
    cmp_ok($oSpot3->spot_status,'==',1,'spot3 passed QC');
    cmp_ok($oSpot1->channel1_signal,'==',8401,'spot1 ch1 signal');
    cmp_ok($oSpot2->channel1_signal,'==',8778,'spot2 ch1 signal');
    cmp_ok($oSpot3->channel1_signal,'==',9770,'spot3 ch1 signal');
    cmp_ok($oSpot1->channel2_signal,'==',11307,'spot1 ch2 signal');
    cmp_ok($oSpot2->channel2_signal,'==',10580,'spot2 ch2 signal');
    cmp_ok($oSpot3->channel2_signal,'==',11328,'spot3 ch2 signal');
};
test "RP11-622I2" => sub {
	my $oRP11_622I2;
    ok($oRP11_622I2 = $oArray->get_reporter('RP11-622I2'),'get_reporter RP11-622I2');
	cmp_ok($oRP11_622I2->reporter_id,'eq','RP11-622I2','reporter_id');
    cmp_ok($oRP11_622I2->get_reporter_replicates,'==',3,'get_reporter_replicates');
	cmp_ok($oRP11_622I2->spots_passed_qc,'==',2,'RP11-622I2 spots passed QC');
    cmp_ok($oRP11_622I2->mean_ch1,'==',10582.5,'mean_ch1');
    cmp_ok($oRP11_622I2->mean_ch2,'==',10189,'mean_ch2');
    cmp_ok($oRP11_622I2->mean_ratios,'eq',0.973928228448155,'mean_ratios');
    cmp_ok($oRP11_622I2->ratio_means,'eq',1.03862008047895,'ratio means');

    my $aSpots = $oRP11_622I2->get_reporter_spots;
    is(@$aSpots,3,'3 spots');
    my ($oSpot1,$oSpot2,$oSpot3) = @$aSpots;
	isa_ok($oSpot1,'array_spot','spot object 1');
	isa_ok($oSpot2,'array_spot','spot object 2');
	isa_ok($oSpot3,'array_spot','spot object 3');
    cmp_ok($oSpot1->spot_index,'==',1741,'spot1 index');
    cmp_ok($oSpot2->spot_index,'==',1925,'spot2 index');
    cmp_ok($oSpot3->spot_index,'==',2108,'spot3 index');
    cmp_ok($oSpot1->spot_status,'==',0,'spot1 failed QC');
    cmp_ok($oSpot2->spot_status,'==',1,'spot2 passed QC');
    cmp_ok($oSpot3->spot_status,'==',1,'spot3 passed QC');
    cmp_ok($oSpot1->channel1_signal,'==',4904,'spot1 ch1 signal');
    cmp_ok($oSpot2->channel1_signal,'==',17251,'spot2 ch1 signal');
    cmp_ok($oSpot3->channel1_signal,'==',3914,'spot3 ch1 signal');
    cmp_ok($oSpot1->channel2_signal,'==',5188,'spot1 ch2 signal');
    cmp_ok($oSpot2->channel2_signal,'==',15809,'spot2 ch2 signal');
    cmp_ok($oSpot3->channel2_signal,'==',4569,'spot3 ch2 signal');
};
test "RP11-753E17" => sub {
	my $oRP11_753E179;
    ok($oRP11_753E179 = $oArray->get_reporter('RP11-753E17'),'get_reporter RP11-753E17');
	cmp_ok($oRP11_753E179->reporter_id,'eq','RP11-753E17','reporter_id');
    cmp_ok($oRP11_753E179->get_reporter_replicates,'==',3,'get_reporter_replicates');
	cmp_ok($oRP11_753E179->spots_passed_qc,'==',3,'RP11-753E17 spots passed QC');
    cmp_ok($oRP11_753E179->mean_ch1,'eq',13229.3333333333,'mean_ch1');
    cmp_ok($oRP11_753E179->mean_ch2,'eq',13865.3333333333,'mean_ch2');
    cmp_ok($oRP11_753E179->mean_ratios,'eq',0.955659050793768,'mean_ratios');
    cmp_ok($oRP11_753E179->ratio_means,'eq',0.954130204827387,'ratio means');

    my $aSpots = $oRP11_753E179->get_reporter_spots;
    is(@$aSpots,3,'3 spots');
    my ($oSpot1,$oSpot2,$oSpot3) = @$aSpots;
	isa_ok($oSpot1,'array_spot','spot object 1');
	isa_ok($oSpot2,'array_spot','spot object 2');
	isa_ok($oSpot3,'array_spot','spot object 3');
    cmp_ok($oSpot1->spot_index,'==',27247,'spot1 index');
    cmp_ok($oSpot2->spot_index,'==',27430,'spot2 index');
    cmp_ok($oSpot3->spot_index,'==',27614,'spot3 index');
    cmp_ok($oSpot1->spot_status,'==',1,'spot1 passed QC');
    cmp_ok($oSpot2->spot_status,'==',1,'spot2 passed QC');
    cmp_ok($oSpot3->spot_status,'==',1,'spot3 passed QC');
    cmp_ok($oSpot1->channel1_signal,'==',15137,'spot1 ch1 signal');
    cmp_ok($oSpot2->channel1_signal,'==',12098,'spot2 ch1 signal');
    cmp_ok($oSpot3->channel1_signal,'==',12453,'spot3 ch1 signal');
    cmp_ok($oSpot1->channel2_signal,'==',16158,'spot1 ch2 signal');
    cmp_ok($oSpot2->channel2_signal,'==',12647,'spot2 ch2 signal');
    cmp_ok($oSpot3->channel2_signal,'==',12791,'spot3 ch2 signal');
};


end_skipping_tests;
