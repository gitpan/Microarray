#!/usr/bin/perl -w

use strict;


use Test::More tests => 12;
BEGIN { use_ok('Microarray') };
BEGIN { use_ok('Microarray::Feature') };
BEGIN { use_ok('Microarray::Analysis') };
BEGIN { use_ok('Microarray::Analysis::CGH') };
BEGIN { use_ok('Microarray::Spot') };
BEGIN { use_ok('Microarray::Image') };
BEGIN { use_ok('Microarray::File') };
BEGIN { use_ok('Microarray::File::Image') };
BEGIN { use_ok('Microarray::File::Clone_Locns') };
BEGIN { use_ok('Microarray::File::Data::GenePix') };
BEGIN { use_ok('Microarray::File::Data::Agilent') };
BEGIN { use_ok('Microarray::File::Data::Manor_Output') };

