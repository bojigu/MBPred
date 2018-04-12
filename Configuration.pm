#!/usr/bin/env perl

package Configuration;
require Exporter;

use warnings;

our @ISA    = qw / Exporter /;
our @EXPORT = qw / $path_to_TMBPred
                   $path_to_Rscript
                   $path_to_hhblits
                   $path_to_hhblits_database
                   $path_to_polyphobius
                   $path_to_freecontact 
                   $path_to_Lips 
                   $path_to_cluster
                   $path_to_aaindex /;

$path_to_TMBPred              = "/home/students/zeng/TMBPred1.0";
$path_to_Rscript              = "/home/students/zeng/usr/local/bin/Rscript";
$path_to_hhblits              = "/home/students/zeng/usr/opt/hhsuite/bin/hhblits";
$path_to_hhblits_database     = "/scratch/zeng/uniprot20_2015_06/uniprot20_2015_06";
$path_to_polyphobius          = "/home/students/zeng/TMBPred1.0/Scripts/phobius.pl";
$path_to_freecontact          = "/home/software/bin/freecontact";
$path_to_Lips                 = "/home/students/zeng/TMBPred1.0/Scripts/lips.pl";
##give the particular cluster name that the jobs will be submitted to
$path_to_cluster              = "all.q";
$path_to_aaindex              = "/home/students/zeng/TMBPred1.0/AAindex/AAindex";


1;
