#!/usr/bin/env perl

package MBPred;
require Exporter;
use warnings;

our @ISA       = qw / Exporter /;

                     
our @EXPORT    = qw / Run_Hhblits
                      Oa3m_Parser
                      Pssm_Entropy_Calculate
                      Run_Freecontact
                      Run_Phobius
                      Read_Fasta
                      Segment_Parser
                      Get_Seg_Start_End
                      Get_Res_Di_Mi_Hash
                      Freecontact_Parser
                      Relative_Position_Calc
                      Lips_Score_Calc
					  COnvert_Bind
                      MergeFile 
					  Add_Physical_Property
					  Create_TestData
					  Remove_Uncomplete_Lines
					  Run_Test/;
               

#-------------------------------------------------------------------------------    
#                               - Subroutines -   
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Purpose   : run Hhblits to extract homologous from uniprot database
# Usage     : Run_Hhblits($Fastaf,$Oa3mf) 
# Arguments : the input fasta file and output oa3m file
# Returns   : the subroutine returns nothing but the homologous file will be saved
# Globals   : none
#******************

sub Run_Hhblits
{
	print STDERR " ... start to run hhblits ... ";
	my $set=$_[0];
	my %set=%$set;
	my $fasta_file = $set{"fastaf"};
	my $output_oa3m_file=$set{"base_dir"}. "/TestData/". $set{"protein"}. ".out.oa3m";
	if($set{"cluster_name"} ne ""){
	print "\n ... using server cluster ". $set{"cluster_name"}. " to run hhblits ... \n";
	system("qsub -b y -q $set{\"cluster_name\"} -N \"HhblitsJob\" $set{\"hhblits\"} -i $fasta_file  -oa3m $output_oa3m_file -d $set{\"hhblits_database\"} -Z 999999999	 -B 999999999 -maxfilt 999999999 -id 99 -diff inf");
	}
	else
	{
	    print "\n ... run hhblits at local computer... \n";
	    system("$set{\"hhblits\"} -i $fasta_file  -oa3m $output_oa3m_file -d $set{\"hhblits_database\"}");
	}
}


#-------------------------------------------------------------------------------
# Purpose   : parse hhblits output oa3m file, the output multiple sequense
#			  alignments will have no lower-case character and no gaps in 
#             the query sequence.  
# Usage     : Oa3m_Parser($Oa3mf,$A3mf)
#             sequences don't contain '>'
# Arguments : hhblits output oa3m file, and the parsed MSAs a3m output filename
# Returns   : nothing
# Globals   : none
#******************

sub Oa3m_Parser
{
	my $set=$_[0];
	my %set=%$set;
	my $oa3m_file = $set{"base_dir"}. "/TestData/". $set{"protein"}. ".out.oa3m";
	my $output_a3m_file=$set{"base_dir"}. "/TestData/". $set{"protein"}. ".a3m";
	system("grep -v \"^>\" $oa3m_file |sed 's/[a-z]//g' >$output_a3m_file");
}



#-------------------------------------------------------------------------------
# Purpose   : calculate residue Position Specific Scoring Matrix (PSSM)
#			  and residue Conservation (here is entropy)
# Usage     : Pssm_Entropy_Calculate($A3mf,$Out_Pssm_file,$Out_Entropy_file)
# Arguments : parsed a3m file as MSAs input, and the pssm and entropy output files
# Returns   : nothing
# Globals   : none
#******************

sub Pssm_Entropy_Calculate
{
	print STDERR " ... start to calculate pssm and entropy ... ";
	my $set=$_[0];
	my %set=%$set;
	my $pssm_output_file = $set{"base_dir"}. "/TestData/". $set{"protein"}. ".pssm";
	my $entropy_output_file = $set{"base_dir"}. "/TestData/". $set{"protein"}. ".entropy";
	open Pssm_Output, ">$pssm_output_file" or die $!;
	open Entropy_Output, ">$entropy_output_file" or die $!;
	my $pssm_head = join(" ","Res_num","Res_name","A","I","L","V","F","W","Y","N","C","Q","M","S","T","D","E","R","H","K","G","P");
	my $entropy_head = join(" ","Res_num","Res_name","Entropy");
	print Pssm_Output "$pssm_head\n";
	print Entropy_Output "$entropy_head\n";
	$set{'a3mf'}=$set{"base_dir"}. "/TestData/". $set{"protein"}. ".a3m";
	open FILE,"$set{'a3mf'}" or die $!;
	my @filedata=<FILE>;
	close (FILE);
	my @pssm_array;
	foreach my $line (@filedata) 
	{
		my @row = split ("", $line);
		push(@pssm_array, [@row]);
	}
	my $collen = scalar @pssm_array;
	my $rowlen = length($filedata[0]) - 1;  #$rowlength include \n, so minus 1
	my @res_seq = split("", $filedata[0]);
	my @aanum= (0) x 20;
	for (my $j = 0; $j < $rowlen; $j++)
	{
		for (my $i = 0; $i < $collen; $i++)
		{
				if($pssm_array[$i][$j] eq 'A')	{$aanum[0] = $aanum[0] + 1;}
				if($pssm_array[$i][$j] eq 'I')	{$aanum[1] = $aanum[1] + 1;}
				if($pssm_array[$i][$j] eq 'L')	{$aanum[2] = $aanum[2] + 1;}
				if($pssm_array[$i][$j] eq 'V')	{$aanum[3] = $aanum[3] + 1;}
				if($pssm_array[$i][$j] eq 'F')	{$aanum[4] = $aanum[4] + 1;}
				if($pssm_array[$i][$j] eq 'W')	{$aanum[5] = $aanum[5] + 1;}
				if($pssm_array[$i][$j] eq 'Y')	{$aanum[6] = $aanum[6] + 1;}
				if($pssm_array[$i][$j] eq 'N')	{$aanum[7] = $aanum[7] + 1;}
				if($pssm_array[$i][$j] eq 'C')	{$aanum[8] = $aanum[8] + 1;}
				if($pssm_array[$i][$j] eq 'Q')	{$aanum[9] = $aanum[9] + 1;}
				if($pssm_array[$i][$j] eq 'M')	{$aanum[10] = $aanum[10] + 1;}
				if($pssm_array[$i][$j] eq 'S')	{$aanum[11] = $aanum[11] + 1;}
				if($pssm_array[$i][$j] eq 'T')	{$aanum[12] = $aanum[12] + 1;}
				if($pssm_array[$i][$j] eq 'D')	{$aanum[13] = $aanum[13] + 1;}
				if($pssm_array[$i][$j] eq 'E')	{$aanum[14] = $aanum[14] + 1;}
				if($pssm_array[$i][$j] eq 'R')	{$aanum[15] = $aanum[15] + 1;}
				if($pssm_array[$i][$j] eq 'H')	{$aanum[16] = $aanum[16] + 1;}
				if($pssm_array[$i][$j] eq 'K')	{$aanum[17] = $aanum[17] + 1;}
				if($pssm_array[$i][$j] eq 'G')	{$aanum[18] = $aanum[18] + 1;}
				if($pssm_array[$i][$j] eq 'P')	{$aanum[19] = $aanum[19] + 1;}
		}
		###pssm calculate
		my $res_num=$j+1;
		my @t = ();
   		foreach my $n (@aanum){
      		if ($n == 0){
         	$n = '0';
      		}
      		my $t = $n/$collen;
     		push(@t,$t);
		}
		print Pssm_Output "$res_num\t$res_seq[$j]\t@t\n";
		###entropy calculate
		my $entropy=0;
		for(my $k=0;$k<20;$k++)
		{
			$aanum[$k]/=$collen;
			if($aanum[$k]!=0)
			{
				$entropy+=-1*$aanum[$k]*(log($aanum[$k])/log(2));
			}
		}
		print Entropy_Output "$res_num\t$res_seq[$j]\t$entropy\n";
		
	}
	close(Pssm_Output);
	close(Entropy_Output);
	print STDERR "  .... pssm and entropy calculation done ... ";
}



#-------------------------------------------------------------------------------
# Purpose   : run freecontact to get residue-residue co-evolution information
#			  Mutual information and Direct Interact score were calculated
# Usage     : Run_Freecontact($A3mf,$Out_Freecontact_file)
# Arguments : input MSAs file and output co-evolution file
# Returns   : nothing
# Globals   : none
#******************

sub Run_Freecontact
{
	print STDERR "  .... start to run freecontact ... ";
	my $set=$_[0];
	my %set=%$set;
	my $a3m_file=$set{"base_dir"}. "/TestData/". $set{"protein"}. ".a3m";
	my $freecontact_file=$set{"base_dir"}. "/TestData/". $set{"protein"}. ".freecontact";
	if($set{"cluster_name"} ne "")
	{
        print "\n ... using server cluster ". $set{"cluster_name"}. " to run freecontact ... \n";
        system("qsub -b y -q $set{\"cluster_name\"} -o $freecontact_file -N \"FreecontactJob\" $set{\"freecontact\"} -f $a3m_file");
	}
	else
	{
	    print "\n ... run freecontact at local computer... \n";
	    system("$set{\"freecontact\"} -f $a3m_file > $freecontact_file");
	}
}



#-------------------------------------------------------------------------------
# Purpose   : read fasta format file and get the protein sequence
# Usage     : Read_Fasta($Fastaf)
# Arguments : input fasta file
# Returns   : protein sequence string 
# Globals   : none
#******************

sub Read_Fasta
{
	my $fasta;
	my $fastafile=$_[0];
	open FILE, "$fastafile" or die $!;
	while(<FILE>)
	{
		chomp;
		if($_=~/^>/)
		{
			next;
		}
		else
		{ 
			$fasta .=$_;
		}
	}
	close(FILE);
	return $fasta;
}


#-------------------------------------------------------------------------------
# Purpose   : run phobius prediction to get transmembrane, cytoplasmic and 
#             extra-cellular segments 
# Usage     : Run_Phobius($Fastaf,$MemSegF)
# Arguments : input fasta file and output phobius prediciton results
# Returns   : nothing
# Globals   : none
#******************

sub Run_Phobius
{
	print STDERR "\n .... start to run phobus prediction ... \n";
	my $set=$_[0];
	my %set=%$set;
	my $fasta_file=Read_Fasta($set{"fastaf"});
	my $seg_file=$set{"base_dir"}. "/TestData/". $set{"protein"}. ".seg";
	system("perl $set{\"phobius\"} $fasta_file >$seg_file");
	print STDERR "\n ...phobius finished running ... \n";
}



#-------------------------------------------------------------------------------
# Purpose   : parse phobius output file, to get transmembrane, cytoplasmic and extracellular
#	          subunits, if protein have multiple pass membranes, then multiple above segments 
#			  sequences will be got and saved seperately  
# Usage     : Segment_Parser($set{"protein"},$MemSegF,$set{"base_dir"},$fasta_seq)
# Arguments : input protein name, phobius output file, segments output file save directory,
#             and the protein sequence string got from subruting Read_Fasta()
# Returns   : nothing
# Globals   : none
#******************

sub Segment_Parser
{
print STDERR "  .... start to parse segments .... \n";
	my $set=$_[0];
	my %set=%$set;
	my $protein=$set{"protein"};
	my $fastaseq=Read_Fasta($set{"fastaf"});
	my $seg_file=$set{"base_dir"}. "/TestData/". $set{"protein"}. ".seg";
	my $i=0;
	my $j=0;
	my $k=0;
    open FILE,"$seg_file" or die $!;
	while(<FILE>)
	{
		chomp;
		my @array=split(" ",$_);
        my $seg;
        my $seg_output_file;
        my $segbeg;
        my $segend;
        if($_=~/NON/)
        {
            $seg = "out";
            $i = $i +1;
            $seg_output_file=$set{"base_dir"}. "/TestData/".  $protein. ".". $seg. $i;
            $segbeg=$array[2];
            $segend=$array[3];
            open SegOut, ">$seg_output_file" or die $!;
            for(my $Seg=$segbeg;$Seg<=$segend;$Seg++)
            {
                my $char=substr($fastaseq,$Seg-1,1);
                print SegOut "$Seg\t$char\t$seg\n";
            }
            close(SegOut);
        }
        if($_=~/TRANSMEM/)
        {
            $seg = "mem";
            $j = $j +1;
            $seg_output_file=$set{"base_dir"}. "/TestData/".  $protein. ".". $seg. $j;
            $segbeg=$array[2];
            $segend=$array[3];
            open SegOut, ">$seg_output_file" or die $!;
            for(my $Seg=$segbeg;$Seg<=$segend;$Seg++)
            {
                my $char=substr($fastaseq,$Seg-1,1);
                print SegOut "$Seg\t$char\t$seg\n";
            }
            close(SegOut);
        }
        if($_=~/CYTOPLASMIC./)
        {
            if($array[4]=~/CYTOPLASMIC./)
            {
                $seg = "in";
                $k = $k +1;
                $seg_output_file=$set{"base_dir"}. "/TestData/".  $protein. ".". $seg. $k;
                $segbeg=$array[2];
                $segend=$array[3];
                open SegOut, ">$seg_output_file" or die $!;
                for(my $Seg=$segbeg;$Seg<=$segend;$Seg++)
                {
                	my $char=substr($fastaseq,$Seg-1,1);
                	print SegOut "$Seg\t$char\t$seg\n";
                }
                close(SegOut);
            }
        }
	}
	close(FILE);

	#combine segment subunits
    my @segs = qw('in' 'mem' 'out');
	foreach (@segs)
	{
		my $Seg_file=$set{"base_dir"}. "/TestData/". $protein. ".". $_;
		my $cat_file=$Seg_file. "[0-30]";  #cat all the sub segments ,here assuming each protein has less than 30 tm helix
		system("cat $cat_file >$Seg_file");
	}
	print STDERR ".... Segments parser finished ... \n";
}


#-------------------------------------------------------------------------------
# Purpose   : parse the segment file and get the segment start and end position. 
# Usage     : Get_Seg_Start_End($segmentfile)
# Arguments : single segment file
# Returns   : segment start and end positions
# Globals   : none
#******************

sub Get_Seg_Start_End
{
	my $segstart;
	my $segend;
	my $i=0;
	open FILE, "$_[0]" or die $!;
	while(<FILE>)
	{
		chomp;
		my @array=split(" ",$_);
		if($i==0)
		{
			$segstart=$array[0];
			$segend=$array[0];
		}
		else
		{
			$segend=$array[0];
		}
		$i++;
	}
	close(FILE);
	return $segstart,$segend;
}



#-------------------------------------------------------------------------------
# Purpose   : parse the freecontact output file, and for each residue, stores the
#             mutual information and direct interact scores with anyother residues along
#			  to the single segment or full TM sequence. 
# Usage     : Get_Res_Di_Mi_Hash($freecontact_outputfile,$segment_file_path)
# Arguments : freecontact output file and segment file directory
# Returns   : four hashes store co-evolutionary information for each residue
# Globals   : none
#******************

sub Get_Res_Di_Mi_Hash
{
	my %Di_Hash_Seq;
	my %Mi_Hash_Seq;
	my %Di_Hash_Seg;
	my %Mi_Hash_Seg;
	my ($segstart,$segend)=Get_Seg_Start_End($_[1]);
	open FreecontactFile, "$_[0]" or die $!;
	while(<FreecontactFile>)
	{
		chomp;
		if($_=~/^Your/)
		{
			next;
		}
		my @array=split(" ",$_);
		#Di hash in the full seq
		if(!exists $Di_Hash_Seq{$array[0]})
		{
			$Di_Hash_Seq{$array[0]}="$array[5]";
		}
		else
		{
			$Di_Hash_Seq{$array[0]}="$Di_Hash_Seq{$array[0]}_$array[5]";
		}
		if(!exists $Di_Hash_Seq{$array[2]})
		{
			$Di_Hash_Seq{$array[2]}="$array[5]";
		}
		else
		{
			$Di_Hash_Seq{$array[2]}="$Di_Hash_Seq{$array[2]}_$array[5]";
		}
		#mi hash in the full seq
		if(!exists $Mi_Hash_Seq{$array[0]})
		{
			$Mi_Hash_Seq{$array[0]}="$array[4]";
		}
		else
		{
			$Mi_Hash_Seq{$array[0]}="$Mi_Hash_Seq{$array[0]}_$array[4]";
		}
		if(!exists $Mi_Hash_Seq{$array[2]})
		{
			$Mi_Hash_Seq{$array[2]}="$array[4]";
		}
		else
		{
			$Mi_Hash_Seq{$array[2]}="$Mi_Hash_Seq{$array[2]}_$array[4]";
		}

		if(($segstart <= $array[0] and $array[0] <= $segend) and ($segstart <= $array[2] and $array[2] <= $segend))
		{
			#Di hash in the segment
			if(!exists $Di_Hash_Seg{$array[0]} )
			{
				$Di_Hash_Seg{$array[0]}="$array[5]";
			}
			else
			{
				$Di_Hash_Seg{$array[0]}="$Di_Hash_Seg{$array[0]}_$array[5]";
			}
			if(!exists $Di_Hash_Seg{$array[2]})
			{
				$Di_Hash_Seg{$array[2]}="$array[5]";
			}
			else
			{
				$Di_Hash_Seg{$array[2]}="$Di_Hash_Seg{$array[2]}_$array[5]";
			}
			#mi hash in the segment
			if(!exists $Mi_Hash_Seg{$array[0]})
			{
				$Mi_Hash_Seg{$array[0]}="$array[4]";
			}
			else
			{
				$Mi_Hash_Seg{$array[0]}="$Mi_Hash_Seg{$array[0]}_$array[4]";
			}
			if(!exists $Mi_Hash_Seg{$array[2]})
			{
				$Mi_Hash_Seg{$array[2]}="$array[4]";
			}
			else
			{
				$Mi_Hash_Seg{$array[2]}="$Mi_Hash_Seg{$array[2]}_$array[4]";
			}
		}
	}
	close(FreecontactFile);
	return (\%Di_Hash_Seq,\%Mi_Hash_Seq,\%Di_Hash_Seg,\%Mi_Hash_Seg);
}



#-------------------------------------------------------------------------------
# Purpose   : parse freecontact output file to calculate 8 types of co-evolutionary strengthes
#             for each residue.8 scores are cumulative MI or DI score along the segment or full
#             sequence. plus maximum MI and DI scores along segment and full sequence. 
# Usage     : Freecontact_Parser($Out_Freecontact_file)
# Arguments : the freecontact outut file as input
# Returns   : nothing, but the cumulative coevolution scores will be saved 
# Globals   : none
#******************

sub Freecontact_Parser
{
	print STDERR "  - Start to parse freecontact output and calculate CCES ... ";
	my $set=$_[0];
	my %set=%$set;
	my $Fastaf=$set{"fastaf"};
	my $protein=$set{"protein"};
	opendir(DIR, $set{"base_dir"}. "/TestData/");
	my @files = grep(/$protein\.out\d{1,2}$|$protein\.mem\d{1,2}$|$protein\.in\d{1,2}$/,readdir(DIR));
	closedir(DIR);
	#my $file=$files[1];
	foreach my $file (@files) {
		my $Cumu_Str_out_file=$set{"base_dir"}. "/TestData/". $file. ".CumStr";
		open CumuOut, ">$Cumu_Str_out_file" or die $!;
		print CumuOut "Res_num\tRes_name\tDiseq\tMiseq\tDiseg\tMiseg\tDiseqmax\tMiseqmax\tDisegmax\tMisegmax\n";
		my $seg_file_path=$set{"base_dir"}. "/TestData/". $file;
		my ($SegStart,$SegEnd)=Get_Seg_Start_End($seg_file_path);
		##remember not available to return two hashes in subrouting but with scalar references like below
		my $freecontact_output_file=$set{"base_dir"}. "/TestData/". $set{"protein"}. ".freecontact";
		my ($DI_HASH_SEQ,$MI_HASH_SEQ,$DI_HASH_SEG,$MI_HASH_SEG) = Get_Res_Di_Mi_Hash($freecontact_output_file,$seg_file_path);
		my %DI_HASH_SEQ = %$DI_HASH_SEQ;
		my %MI_HASH_SEQ = %$MI_HASH_SEQ;
		my %DI_HASH_SEG = %$DI_HASH_SEG;
		my %MI_HASH_SEG = %$MI_HASH_SEG;	
		#return %DI_HASH_SEQ;
		foreach my $key ( sort {$a<=>$b} keys  %DI_HASH_SEQ)
		{
			##only consider segment residues here
			if($key >= $SegStart and $key <= $SegEnd)
			{	
				#di and mi array in full sequence
				my @Di_Seq_arr=split("_",$DI_HASH_SEQ{$key});
				my @Di_Seq_Sort_arr=sort {$b<=>$a} @Di_Seq_arr;
				my @Mi_Seq_arr=split("_",$MI_HASH_SEQ{$key});
				my @Mi_Seq_Sort_arr=sort {$b<=>$a} @Mi_Seq_arr;

				#di and mi array in segment
				my @Di_Seg_arr=split("_",$DI_HASH_SEG{$key});
				my @Di_Seg_Sort_arr=sort {$b<=>$a} @Di_Seg_arr;
				my @Mi_Seg_arr=split("_",$MI_HASH_SEG{$key});
				my @Mi_Seg_Sort_arr=sort {$b<=>$a} @Mi_Seg_arr;

				##di and mi sum in the full seq and segment
				my $disumseq=0;
				my $disumseg=0;
				my $misumseq=0;
				my $misumseg=0;
				for(my $i=0;$i<8;$i++)
				{
					$disumseq+=$Di_Seq_Sort_arr[$i];
					$disumseg+=$Di_Seg_Sort_arr[$i];
					$misumseq+=$Mi_Seq_Sort_arr[$i];
					$misumseg+=$Mi_Seg_Sort_arr[$i];
				}

				##di and mi max in the full seq and segment
				my $dimaxseq=$Di_Seq_Sort_arr[0];
				my $dimaxseg=$Di_Seg_Sort_arr[0];
				my $mimaxseq=$Mi_Seq_Sort_arr[0];
				my $mimaxseg=$Mi_Seg_Sort_arr[0];
			
				my $fasta=Read_Fasta($Fastaf);
				my $str=substr($fasta,$key-1,1);
				my @Cumu_arr=($key,$str,$disumseq,$disumseg,$misumseq,$misumseg,$dimaxseq,$dimaxseg,$mimaxseq,$mimaxseg);
				my $Cumu_arr_str=join("\t",@Cumu_arr). "\n";
				print CumuOut $Cumu_arr_str;
				#return @Cumu_arr;
			}
		}
		#close(CumuOut);
	}
	print STDERR "  - Cumulative Coevolution Strength Calculation fininshed ... ";
}





#-------------------------------------------------------------------------------
# Purpose   : calculate the residue relative position in the full sequence or
#             the correponding segment
# Usage     : Relative_Position_Calc();
# Arguments : none
# Returns   : nothing
# Globals   : none
#******************

sub Relative_Position_Calc
{
	print STDERR "  - start to calculate relative postions ... ";
	my $set=$_[0];
	my %set=%$set;
	my $Fastaf=$set{"fastaf"};
	my $protein=$set{"protein"};
	opendir(DIR, $set{"base_dir"}. "/TestData/");
	my @files = grep(/$protein\.out\d{1,2}$|$protein\.mem\d{1,2}$|$protein\.in\d{1,2}$/,readdir(DIR));
	closedir(DIR);
	#my $file=$files[1];
	foreach my $file (@files) 
	{
		my $Rel_Pos_out_file=$set{"base_dir"}. "/TestData/". $file. ".RelPos";
		open RelPosOut, ">$Rel_Pos_out_file" or die $!;
		print RelPosOut "Res_num\tRes_name\tRelPos1\tRelPos2\n";
		my $seg_file_path=$set{"base_dir"}. "/TestData/". $file;
		my ($SegStart,$SegEnd)=Get_Seg_Start_End($seg_file_path);
		my $fasta=Read_Fasta($Fastaf);
		my $seg_len=$SegEnd-$SegStart+1;
		my $seq_len=length($fasta);
		open SegFile, "$seg_file_path" or die $!;
		my $row=1;
		while(<SegFile>)
		{
			chomp;
			my @seg_arr=split(" ",$_);
			my ($res_num,$res_name)=($seg_arr[0],$seg_arr[1]);
			my $rel_pos1=$row/$seg_len;
			my $rel_pos2=$res_num/$seq_len;
			print RelPosOut "$res_num\t$res_name\t$rel_pos1\t$rel_pos2\n";
			$row=$row+1;
		}
	}
	print STDERR "  - Relative Position calculation finished ... ";
}







#-------------------------------------------------------------------------------
# Purpose   : this subruting includes three seperate methods, first from MSAs 
#			  to remove rows with gaps, second no gaps MSAs was inputted into 
#			  LIPS method function and got 7 predicted helix-oriendted surfaces each
#			  surface with LIPS score, finally parse the LIPS output to get for each 
#			  residue the entropy and lipophilicity scores.   
# Usage     : Lips_Score_Calc();
# Arguments : none
# Returns   : nothing
# Globals   : none
#******************

sub Lips_Score_Calc
{
	print STDERR "  - Running and Parsing LIPS method ... ";
	my $set=$_[0];
	my %set=%$set;
	my $protein=$set{"protein"};
	opendir(DIR, $set{"base_dir"}. "/TestData/");
	my $a3m_file=$set{"base_dir"}. "/TestData/". $set{"protein"}. ".a3m";
	#only for membrane region calculation
	my @files = grep(/$protein\.mem\d{1,2}$/,readdir(DIR));
	closedir(DIR);
	foreach my $file (@files) 
	{
	my $seg_file_path=$set{"base_dir"}. "/TestData/". $file;
	my ($SegStart,$SegEnd)=Get_Seg_Start_End($seg_file_path);
	my $A3m_Nogap_out_file=$set{"base_dir"}. "/TestData/". $file. ".lipinput";
	#parser a3m file to get TM segment homology and nogaps in TM homology
	open A3mNogapOut, ">$A3m_Nogap_out_file" or die $!;
	open A3mFile, "$a3m_file" or die $!;
	while(<A3mFile>)
		{
			chomp;
			my $memseg=substr($_,$SegStart-1,$SegEnd-$SegStart+1);
			if($memseg!~/-/)
			{
				print A3mNogapOut "$memseg\n";
			}
		}
		close(A3mFile);
		close(A3mNogapOut);

		##run lips.pl 
		my $Lips_out_file=$set{"base_dir"}. "/TestData/". $file. ".lipout";
		system("perl $set{\"lips\"} $A3m_Nogap_out_file > $Lips_out_file");
		
		###parse Lips out 
		my %ResHash;
		my $Lips_Score_out_file=$set{"base_dir"}. "/TestData/". $file. ".lipsscore";
		open LipsScoreOut, ">$Lips_Score_out_file" or die $!;
		print LipsScoreOut "Res_num\tRes_name\tLipophilicity\tLip_entropy\n";
		open FILE, "$Lips_out_file" or die $!;
		while(<FILE>)
		{
			chomp;
			if($_=~/^\s+\d+\s+[A-Z]/)
				{
					my @array=split(" ",$_);
					if(!exists $ResHash{$array[0]})
					{
						$ResHash{$array[0]}=1;		
						print LipsScoreOut "$_\n";
					}
				}
		}
		close(FILE);
		close(LipsScoreOut);
		my $Lips_Score_out_file_sort=$set{"base_dir"}. "/TestData/". $file. ".lipsscore.sort";
		system("sort -g $Lips_Score_out_file >$Lips_Score_out_file_sort");
		system("mv $Lips_Score_out_file_sort $Lips_Score_out_file");
	}
	print STDERR " .... LIPS method done ... ";
}


#-------------------------------------------------------------------------------
# Purpose   : Convert bind file as full sequence contact or non-contact list 
# Usage     : COnvert_Bind(/%set);
# Arguments : setting hash
# Returns   : nothing
# Globals   : none
#******************

sub Convert_Bind
{
	my $set=$_[0];
	my %set=%$set;
	my $protein=$set{"protein"};
	my %bindhash;
	my $bindfile=$set{"base_dir"}. "/Bind/". $set{"protein"}. ".". $set{"contact_type"}. "bind";
	open FILE, "$bindfile" or die $!;
	while(<FILE>)
	{
		chomp;
		my @array=split(" ",$_);
		$bindhash{$array[0]}=1;
	}
	close(FILE);
	my $fastaseq=Read_Fasta($set{"fastaf"});
	my @seqarr=split("",$fastaseq);
	my $seqlen=length($fastaseq);
	my $Bind_file_out=$set{"base_dir"}. $set{"protein"}. ".". $set{"contact_type"}. "BIND";
	open BindOut,">$Bind_file_out" or die $!;
	print BindOut "Res_num\tRes_name\tBind\n";
	for(my $i=1;$i<=$seqlen;$i++)
	{
		if(exists $bindhash{$i})
		{
			print BindOut "$i\t$seqarr[$i-1]\t1\n";
		}
		else
		{
			print BindOut "$i\t$seqarr[$i-1]\t0\n";
		}
	}
	close(BindOut);

}


#-------------------------------------------------------------------------------
# Purpose   : merge all features into one file, including merging all segments 
#             into one  
# Usage     : MergeFile();
# Arguments : none
# Returns   : nothing
# Globals   : none
#******************

sub MergeFile
{
	print STDERR "  - start to merge features ... ";
	my $set=$_[0];
	my %set=%$set;
	my $protein=$set{"protein"};
	opendir(DIR, $set{"base_dir"}. "/TestData/");
	my @files = grep(/$protein\.out\d{1,2}$|$protein\.mem\d{1,2}$|$protein\.in\d{1,2}$/,readdir(DIR));
	closedir(DIR);
	#my $file=$files[0];
	#print "$file\n";
	foreach my $file (@files) 
	{
		opendir(DIR, $set{"base_dir"}. "/TestData/");
		#my $bind_type=$set{"contact_type"}."BIND";
		my @files1=grep(/^$protein\.pssm$|^$protein\.entropy$|^$file\.CumStr$|^$file\.RelPos$/,readdir(DIR));
		my @sortfiles1=sort { lc($a) cmp lc($b) } @files1;
		closedir(DIR);
		my $SegDataOutFile=$set{"base_dir"}. "/TestData/". $file. ".".  "data";
		open SegDataOut,">$SegDataOutFile" or die$!;
		my %res_num;
		my $i=0;
		foreach my $file1 (@sortfiles1)
		{
			my %hash_kk;
			my $file1_path=$set{"base_dir"}. "/TestData/". $file1;
			open my $infh, '<', $file1_path or die $!;
			while (<$infh>) 
			{
				chomp;
				my @array = split(" ",$_);
				my @arr_val=splice(@array,2,$#array-1);
				my $kk=$array[0]. "\t". $array[1];
				$hash_kk{$kk}=1;
				if($i==0)
				{
					$res_num{$kk}=[@arr_val];
				}
				else
				{
					if(exists $res_num{$kk})
					{
						my @arr =(@{$res_num{$kk}},@arr_val);
						$res_num{$kk}=[@arr];
					}
				}
	
	    	}
			$i++;
			while( my( $k, $value ) = each %res_num )
			{
				if(!exists $hash_kk{$k})
				{
					delete $res_num{$k};
				}
			}
	
		}
	
	
		foreach my $key (sort {$a<=>$b} keys  %res_num)
		{
			print SegDataOut "$key\t@{$res_num{$key}}\n";
		}
	}
}



#-------------------------------------------------------------------------------
# Purpose   : adding physical parameters into the traindata, the physical parameters are got form 
#             the AAindex database
# Usage     : Add_Physical_Property(\%set);
# Arguments : setting parameters
# Returns   : nothing
# Globals   : none
#******************

sub Add_Physical_Property
{
	print STDERR "  - start to add physical property ... \n";
	my $set=$_[0];
	my %set=%$set;
	my $protein=$set{"protein"};
	my $aaindexfile= $set{"aaindex"};
	my %aaindexhash;
	open FILE, "$aaindexfile" or die $!;
	while(<FILE>)
	{
		chomp;
		my @array=split(" ",$_);
		if($_=~/^PP/)
		{
			next;
		}
		my @arr=splice(@array,1,6);
		$aaindexhash{$array[0]}=[@arr];
	}
	close(FILE); 

	opendir(DIR, $set{"base_dir"}. "/TestData/");
	my $bind_type= "data";
	my @files = grep(/^$protein\.out\d{1,2}\.$bind_type$|^$protein\.mem\d{1,2}\.$bind_type$|^$protein\.in\d{1,2}\.$bind_type$/,readdir(DIR));
	closedir(DIR);

	foreach my $file (@files) 
	{
		my ($pro,$seg,$contacttype)=split(/\./,$file);
		my $outputfile=$set{"base_dir"}. "/TestData/". $pro. ".". $seg. ".". "phypardata";
		open OUTPUT, ">$outputfile" or die $!;
		my $binddatafile=$set{"base_dir"}. "/TestData/". $file;
		open FILE1, "$binddatafile" or die $!;
		while(<FILE1>)
		{
			chomp;
			if($_=~/Res_num/)
			{
				print OUTPUT "$_\thydrophobicity\tpolarity\tcharge\tresidue_volume\tisoelectric_point\tsteric_parameter\n";
			}
			my @array1=split(" ",$_);
			if(exists $aaindexhash{$array1[1]})
			{	
				print OUTPUT "$_\t@{$aaindexhash{$array1[1]}}\n";
			}
		}
		close(OUTPUT);
		close(FILE1); 
	}


	print STDERR "  - adding physical property done ... \n";
}

#-------------------------------------------------------------------------------
# Purpose   : join all segment phypardata together to form the testdata
# Usage     : Create_TestData(\%set);
# Arguments : setting parameters
# Returns   : nothing
# Globals   : none
#******************

sub Create_TestData
{
	my $set=$_[0];
	my %set=%$set;
	my $protein=$set{"protein"};
	my $i=0;
	my $TestFile=$set{"base_dir"}. "/TestData/". $protein. ".".  "testdata";
	open OUTPUT, ">$TestFile" or die $!;
	opendir(DIR, $set{"base_dir"}. "/TestData/");
	my $bind_type= "phypardata";
	my @files = grep(/$protein.*\.$bind_type$/,readdir(DIR));
	closedir(DIR);

	foreach my $file (@files) 
	{
		my $ppdatafile=$set{"base_dir"}. "/TestData/". $file;
		open FILE1, "$ppdatafile" or die $!;
		while(<FILE1>)
		{
			chomp;
			if($_=~/Res_num/)
			{
				if($i==0)
				{
					#my @array=split(" ",$_);
					#push (@array, splice(@array, 31, 1));
					#my $row=join(" ",@array);
					print OUTPUT "$_\tSeg\n";
				}
				else{

					next;
				}
				$i+=1;
			}
			else
			{	
				# my @array1=split(" ",$_);
				# push (@array1, splice(@array1, 31, 1));
				# my $row1=join(" ",@array1);
				if($file=~/.mem/)
				{
					print OUTPUT "$_\tM\n";
				}
				if($file=~/.in/)
				{
					print OUTPUT "$_\tI\n";
				}
				if($file=~/.out/)
				{
					print OUTPUT "$_\tO\n";
				}
			}
		}
		close(FILE1); 
	}
	close(OUTPUT);
}
	

#-------------------------------------------------------------------------------
# Purpose   : remove uncomplete liens in the traindata (more or less 12 uncomplete lines found in seg TrainData)
# Usage     : Remove_Uncomplete_Lines(\%set);
# Arguments : setting parameters
# Returns   : nothing
# Globals   : none
#******************

sub Remove_Uncomplete_Lines
{
	print STDERR "  - start to remove uncomplete lines ... \n";
	my $set=$_[0];
	my %set=%$set;
	opendir(DIR, $set{"base_dir"}. "/TestData/");
	#my $bind_type=$set{"contact_type"}. "phypardata";
	my @files = grep(/.*testdata$/,readdir(DIR));
	my @sortfiles=sort { lc($a) cmp lc($b) } @files;
	closedir(DIR);

	foreach my $file (@sortfiles) 
	{
		my $collen;
		my $ppdatafile=$set{"base_dir"}. "/TestData/". $file;
		my $newTestData=$set{"base_dir"}. "/TestData/". $file. "new";
		open OUTPUT ,">$newTestData" or die $!;
		open FILE1, "$ppdatafile" or die $!;
		while(<FILE1>)
		{
			chomp;
			if($_=~/Res_num/)
			{
				my @array=split(" ",$_);
				$collen=$#array;
				print OUTPUT "$_\n";
			}
			else
			{	
				my @array1=split(" ",$_);
				my $collen1=$#array1;
				if($collen eq $collen1)
				{
					print OUTPUT "$_\n";
				}
			}
		}
		close(FILE1); 
		close(OUTPUT);
	}
	
	print STDERR "  - Removing uncomplete lines finished  ... \n";
}



sub Run_Test
{
	my $set=$_[0];
	my %set=%$set;
	my $protein=$set{"protein"};

	my $Rscript_loc = $set{"R_script"};
    my $Random_Forest_R_code_file=$set{"base_dir"}. "/RandomForest/RF_test.R";
    #my $train_data_file=$set{"base_dir"}. "/FASTA/RandomForest/AllFit.rf";
	my $AllFit=$set{"Train_RF"};
	my $test_data_file=$set{"base_dir"}. "/TestData/". $set{"protein"}. ".testdatanew";
	#my $pred_out=$set{"base_dir"}. "/FASTA/TestData/". $set{"protein"}. ".predout";
	my $pred_out=$set{"Pred_Out"};
	system("$Rscript_loc $Random_Forest_R_code_file $AllFit $test_data_file $pred_out");
	
}
1;
