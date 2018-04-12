#!/usr/bin/env perl

##################################################################################################
# Script to calculate the features as inputs of Random Forest for TMPs with fasta sequence input
# 2017/03/02
##################################################################################################
#!/usr/bin/perl -w

use strict;
use warnings;
use TMBPred;
use Common_Util;
use Getopt::Long;
use feature "switch";
use Math::Complex;
use Getopt::Std;
use Time::HiRes;
use Configuration;


print STDERR " \n Run TMBPred as follows: \n";
print STDERR  " perl TMBPredTest.pl -i QuePro.fasta -c AllFit.rf -o PredOut.txt\n";

getopts('i:c:ho:');
our ($opt_i, $opt_c,$opt_h, $opt_o); 


#-------------------------------- Prints help ----------------------------------

if ($opt_h) {
    print STDERR "\nUsage:\n"
               . "perl TMBPredTest.pl [flags]\n\n"
               . "flags:\n"
               . " -i '/some/file' : specify the input the fasta file location.\n" 
               . "      this input is compulsory.\n"
			   . " -c '/some/classifier file' : locate the classifier used.\n"
			   . "      this flag is optional, without the flag the classifier will be TMBPred1.0/RandomForest/AllFit.rf.\n"
               . " -o '/some/file' : writes the output to a user provided output file.\n" 
               . "      this flag is optional, wthout the flag the output will be given to TMBPred1.0/TestData/PredOut.txt \n"            
               . " -h : this help\n\n";
    exit;
}



##################################################################################################
#				create hash to save setting parameters			                                 #
##################################################################################################
my %set;
# Get the basic directory where TMBPred is located, remember first to cd /some path/TMBPred/. 
$set{"base_dir"}=$path_to_TMBPred ;
#get the Rscript path
$set{"R_script"}=$path_to_Rscript ;
#get the hhblits path
$set{"hhblits"}=$path_to_hhblits;
#get the hhblits database path
$set{"hhblits_database"}=$path_to_hhblits_database;
#get the polyphobius script path
$set{"Polyphobius"}=$path_to_polyphobius;
#get the freecontact path
$set{"freecontact"}=$path_to_freecontact ;
#set the LIPS script path
$set{"lips"}=$path_to_Lips  ;
#set the cluster name where jobs will be submitted to 
$set{"cluster_name"}=$path_to_cluster;
#set where the physical parameter file to hash
$set{"aaindex"}=$path_to_aaindex;
#set where the fasta input file path to hash
$set{"fastaf"}=$opt_i;
#set where the Randonm Forest file to hash
if($opt_c){
	$set{"Train_RF"}=$opt_c;
}
else{
	my $current_dir=getcwd();
	$set{"Train_RF"}=$current_dir. "/RandomForest/AllFit.rf";
}
#set where the prediction output file to hash
if ($opt_o) 
{
	$set{"Pred_Out"}=$opt_o;
}
else{
	my $current_dir=getcwd();
	$set{"Pred_Out"}=$current_dir. "/PredOUt.txt";
}



##################################################################################################
#			set conditions to define which functions will be run 			                     #
##################################################################################################
$set{"run_hhblits"}=1;
$set{"oa3m_parser"}=1;
$set{"run_pssm_entropy"}=1;
$set{"run_freecontact"}=1;
$set{"run_phobius"}=1;
$set{"segment_parser"}=1;
$set{"freecontact_parse"}=1;
$set{"run_relative_position"}=1;
$set{"run_lips_score"}=1;
#$set{"bind_convert"}=0;
$set{"merge_file"}=1;
$set{"add_pp"}=1;
$set{"create_testdata"}=1;
$set{"remove_uncomplete_lines"}=1;
$set{"run_test"}=1;

# Gets out the actual protein name from the provided fasta file path
if ($opt_i =~ m/\//) {
    my $inputfile = (split/\//, $opt_i)[-1];
	$set{"protein"}= (split/\./, $inputfile)[-2];
}



##################################################################################################
#					Main Program						 #
##################################################################################################
if($set{"run_hhblits"})
{
	Run_Hhblits(\%set);
	sleep(1) while not -f $set{"base_dir"}. "/TestData/". $set{"protein"}. ".out.oa3m";
}

if($set{"oa3m_parser"})
{
Oa3m_Parser(\%set);
}

if($set{"run_pssm_entropy"})
{
	Pssm_Entropy_Calculate(\%set);
}

if($set{"run_freecontact"})
{
	Run_Freecontact(\%set);
	sleep(1) while not -f $set{"base_dir"}. "/TestData/". $set{"protein"}. ".freecontact";
}


if($set{"run_phobius"})
{
	Run_Phobius(\%set);
}


if($set{"segment_parser"})
{
	my $seq=Segment_Parser(\%set);
	print "$seq\n";
}

if($set{"freecontact_parse"})
{
	Freecontact_Parser(\%set);
}

if($set{"run_relative_position"})
{
	Relative_Position_Calc(\%set);
}

if($set{"run_lips_score"})
{
	Lips_Score_Calc(\%set);
}

if($set{"merge_file"})
{
	MergeFile(\%set);	
}

if($set{"bind_convert"})
{
	Convert_Bind(\%set);	
}

if($set{"add_pp"})
{	
	Add_Physical_Property(\%set);
}

if($set{"create_testdata"})
{	
	Create_TestData(\%set);
}


if($set{"remove_uncomplete_lines"})
{
	Remove_Uncomplete_Lines(\%set);
}


if($set{"run_test"})
{
	Run_Test(\%set);
}

