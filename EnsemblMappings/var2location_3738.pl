#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use List::MoreUtils qw(first_index);


#die "Input file with variant name in 1st column needed as input\n\n" unless scalar(@ARGV) ==1;
die "\n\tUsage: var2location_3738.pl [variant list name (in 'variants/') folder] [GRCh37 || GRCh38]\n\n" unless scalar(@ARGV) ==2;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all("$ARGV[1]/ensembl.registry");

my $v_adaptor  = $registry->get_adaptor('homo_sapiens', 'variation', 'variation');
my $vf_adaptor = $registry->get_adaptor('homo_sapiens', 'variation', 'variationfeature');
## also return variants failing Ensembl QC
$v_adaptor->db()->include_failed_variations(1);

open my $log, ">", "$ARGV[1]/$ARGV[0]\.log" ||die "Failed to open log file  : $!\n";
open my $out, ">", "$ARGV[1]/$ARGV[0]\.out" ||die "Failed to open output file  : $!\n";

open my $infile, "<", "$ARGV[1]/variants/$ARGV[0]\.txt" ||die "Failed to open file of variant names, $ARGV[0] : $!\n";

#my $first_line = <$infile>;
#my @headers = split /\t/, $first_line;
#my $rsid_idx = first_index { $_ eq 'variant_id' } @headers;
my $rsid_idx = 0;

while(<$infile>){
 
  my $varname = (split)[$rsid_idx];
  next unless $varname =~/rs/;

  my $var = $v_adaptor->fetch_by_name($varname);
  unless( $var){
    print $log "No data for $varname\n";
    next;
  }

  my $vfs = $vf_adaptor->fetch_all_by_Variation($var); 

  unless( scalar(@{$vfs}) >0 ){
    print $log "No location data for $varname\n";
    next;
  }

  foreach my $vf (@{$vfs}){
    print $out $varname ."\t". $vf->variation_name() . "\t". $vf->allele_string(). "\t". $vf->seq_region_name() ."\t". $vf->start() . "\t". $vf->end() . "\n";
  }
}
