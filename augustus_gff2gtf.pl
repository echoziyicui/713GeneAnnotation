#!/usr/bin/perl -w
use strict;

my $input_file_name$ = 'augustus_gff_2.gff';
my $output_file_name$ = 'augustus_prediction.gtf';
open (INPUT, $input_file_name);
open (OUTPUT, ">$output_file_name");
my @input = <INPUT>;


close INPUT;
close OUTPUT;