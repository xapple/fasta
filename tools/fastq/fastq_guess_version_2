#!/usr/bin/perl
use strict;
use warnings;
my $MAX_LIM = 50000;  ## these many reads will be sampled


	if( !(scalar(@ARGV)) || ($ARGV[0] eq "-h") || ($ARGV[0] =~ m/-help/) ){
	print_help();
	exit 0;
	}

open FH1,$ARGV[0] or die "\n can not open file $ARGV[0] \n";
##open FH2,">$ARGV[1]";  ## output file
my $counter = 0;
my (%qualHash,$qualValue,@quals,$min,$max);
print "\nsampling..\n";
	while(<FH1>){
	$_ = <FH1>;
	$_ = <FH1>;
	$_ = <FH1>;

	s/\n//;
	s/\r//;
		if($counter >= $MAX_LIM){
		last;
		}

	my @a1 = split('',$_);
		for(my $i = 0;$i<scalar(@a1);++$i){
		$qualValue = ord($a1[$i]); 
		$qualHash{$qualValue} = 1;
		##print FH2 ord($a1[$i]),"\n";
		}
	$counter  += 1;
	}  ## while(<FH1>) ends
close FH1;
##close FH2;


	for my $key (keys %qualHash){
	push(@quals,$key);
	}

@quals = sort { $a <=> $b  }@quals;

##print "\n",join(" ",@quals),"\n";

$min = $quals[1];  ## second
$max = $quals[scalar(@quals)-2];  ## second from last
print "observed range: ",$quals[0],"-",$quals[scalar(@quals)-2],"\n";

print "format: ",which_qual_format($min,$max),"\n\n";

	sub which_qual_format{
	my($min,$max)  = @_;
	my $qualFormat;
		if( ($min >= 33) && ($max <= 76) ){
		$qualFormat = "Sanger or Illumina 1.9+ (offset by 33)";
		}
		elsif( ($min >= 67) && ($max <= 104) ){
		$qualFormat = "Illumina 1.5+ (offset by 64)";
		}
		
		elsif( ($min >= 54) && ($max <= 104) ){
		$qualFormat = "Illumina 1.3+ (offset by 64)";
		}
		
		elsif( ($min >= 59) && ($max <= 104) ){
		$qualFormat = "Solexa (offset by 64)";
		}
		else{
		$qualFormat = "Unknown";
		}
	return($qualFormat);
	}  ## function ends
	###########################
	sub print_help{
	print "\nUsage: perl GetQualityFormat.pl <input fastq>\n";
	print "possible formats: Sanger( or Illumina 1.9+) , Solexa, Illumina 1.3+ and Illumina 1.5+\n\n";
	}  ## function ends
	###########################
	
