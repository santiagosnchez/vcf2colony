#!/usr/bin/perl

my $gzip=0;
my $vcf;
my $name;
my $par_offs_file;
my $al;
my %par_offs=();
my @exclude=();
my @indnms=();
my %data_ind=();
my %data_loc=();
my %loc_tran=();
my @loci=();
my @indiv;
my $locn=0;
my $warn=0;
my $recode=0;
my $usage = "Usage:
perl vcf2colony.pl -v <vcf-file> -o <name> -f <parent-offspring-map> [ --recode ]

example:
perl vcf2colony.pl -v snps.vcf -o mydataset -f parent-offspring.txt --recode

The --recode option will generate a condensed vcf based on the samples found.

Parent-offspring file example:
ind1<tab>O1
ind2<tab>O2
ind3<tab>O3
ind4<tab>M1
ind5<tab>M2
ind6<tab>F1
ind7<tab>F2\n";

# arguments

if ((scalar(@ARGV) == 0) or (grep { /-+he{0,1}l{0,1}p{0,1}/ } @ARGV)){
	die $usage;
}

if (my ($indV) = grep { $ARGV[$_] =~ /^-v$/ } 0 .. $#ARGV){
	$vcf = $ARGV[$indV+1];
	if ($vcf =~ m/\.gz$/){
		print "VCF is compressed\n";
		open(VCF, "gunzip -c $vcf |") or die "Can't open gunzip pipe to $vcf\n";
	} else {
		open(VCF, "<", $vcf) or die "Can't open $vcf\n";
	}
} else {
	die $usage;
}

if (my ($indP) = grep { $ARGV[$_] =~ /^-f$/ } 0 .. $#ARGV){
	$par_offs_file = $ARGV[$indP+1];
	open(PO, "<", $par_offs_file) or die "Can't open $par_offs_file\n";
} else {
	die $usage;
}

if (my ($indN) = grep { $ARGV[$_] =~ /^-o$/ } 0 .. $#ARGV){
	$name = $ARGV[$indN+1];
} else {
	$name = "colony_output";
	print "No name provided with -o ... Using default \"colony_output\"\n";
}

if (my ($indN) = grep { $ARGV[$_] =~ /^--recode$/ } 0 .. $#ARGV){
	$recode = 1;
	print "A new and condensed VCF file will be generated\n";
}

# read and process pop file

while (<PO>){
	chomp $_;
	@line = split /\t/, $_;
	$par_offs{$line[0]} = $line[1];
}
close PO;

# read and process VCF file

open WAR, "> warnings.txt";
local $| = 1;
while(<VCF>){
	if (/^##/){
		$head .= $_;
		if ($recode == 1){
			$ovcf = $vcf;
			$ovcf =~ s/\.vcf.*//;
			open OVCF, "> $ovcf.recode.vcf";
			print OVCF $head;
		}
	}
	elsif (/^#CHROM/){
		chomp($_);
		@HEAD = split /\t/, $_;
		($format_ind) = grep { @HEAD[$_] =~ m/FORMAT/ } 0 .. $#HEAD;
		@indiv = keys %par_offs;
		@myind = getind(\@HEAD,\@indiv);
		if ($recode == 1){
			print OVCF join("\t", (@HEAD[ 0 .. $format_ind ],@HEAD[ @myind ])) . "\n";
		}
		$cline = $. + 1;
		$al = 0;
	}
	else {
		++$al;
		chomp($_);
		@DAT = split /\t/, $_;
		if (scalar(split(/,/, $DAT[4])) > 1){
			print WAR "Warning: allele $DAT[0]:$DAT[1] is not bialellic. Skipping...\n";
			$warn += 1;
			next;
		}
		@site = @DAT[ @myind ];
		map { s/:.*// } @site;
		map { s/\/|\|//g } @site;
		$gentmp = join("", @site);
		$gentmp =~ tr/[01\.]/[120]/;
		if (($gentmp =~ tr/1//) == 0 or ($gentmp =~ tr/2//) == 0){
			print WAR "Warning: allele $DAT[0]:$DAT[1] is monomorphic. Skipping...\n";
			$warn += 1;
			next;
		}
		$locn += 1;
		push @loci, "LOCUS.$locn";
		$loc_tran{"LOCUS.$locn"} = join(":",@DAT[0..4]);
		@genind = ($gentmp =~ m/(\d\d)/g);
		foreach (0 .. $#indiv){
			$data_ind{$indiv[$_]} .= $genind[$_];
		}
		$data_loc{"LOCUS.$locn"} = $gentmp;
		if ($recode == 1){
			my @H = @DAT[ 0 .. 8 ];
			splice(@H,7,1,reformatInfo($gentmp));
			print OVCF join("\t", (@H,@DAT[ @myind ])) . "\n";
		}
		push @nalleles, scalar(split(/,/, $DAT[4])) + 1;
		print "Processing SNP $al\r";
	}
}
close VCF;
close WAR;
close OVCF;
local $| = 0;

print "\nData read. Saving to file: $name.txt\n";

my $rn = 999 + int(rand(9999 - 999));;
my $nl = scalar(keys %data_loc);

# Get Offspring, Mothers, and Fathers from hash table

my @offspring = sort {$a <=> $b} grep { $par_offs{$_} =~ m/^O\d+/ } keys %par_offs;
@offspring = grep { $data_ind{$_} } @offspring;
my @mothers = sort {$a <=> $b} grep { $par_offs{$_} =~ m/^F\d+/ } keys %par_offs;
@mothers = grep { $data_ind{$_} } @mothers;
my @fathers = sort {$a <=> $b} grep { $par_offs{$_} =~ m/^M\d+/ } keys %par_offs;
@fathers = grep { $data_ind{$_} } @fathers;
my $no = scalar(@offspring);

print "Final number of loci: $nl\n";
print scalar(@offspring) . " offspring; " . scalar(@mothers) . " mothers; " . scalar(@fathers) . " fathers\n";
if ($warn == 0){
	print "No warnings recorded\n";
	`rm warnings.txt`;
} else {
	print "$warn warnings recorded in warnings.txt\n";
}


# Prepare output

my $output = 
sprintf("%-12s", $name) . " ! Dataset name\n" .
sprintf("%-12s", $name) . " ! Output file name\n" .
sprintf("%-12s", $no) .   " ! Number of offspring in the sample\n" .
sprintf("%-12s", $nl) .   " ! Number of loci\n" .
sprintf("%-12s", $rn) .   " ! Seed for random number generator\n" .
"0           ! 0/1=Not updating/updating allele frequency
2           ! 2/1=Dioecious/Monoecious species
0           ! 0/1=No inbreeding/inbreeding
0           ! 0/1=Diploid species/HaploDiploid species
0  0        ! 0/1=Polygamy/Monogamy for males & females
0           ! 0/1=Clone inference =No/Yes
1           ! 0/1=Full sibship size scaling =No/Yes
1 1.0 1.0   ! 0,1,2,3=No,weak,medium,strong sibship size prior; mean paternal & maternal sibship size
1           ! 0/1=Unknown/Known population allele frequency\n" . 
join(' ', map { $_ = 2 } 0 .. ($nl-1)) . "  !Number of alleles per locus\n";

open TRAN, "> loci_translation.txt";

foreach (@loci){
	@fr = frequency($data_loc{$_});
	$output .= sprintf("% 7s% 7s\n", (1,2)) . sprintf("% 7s% 7s\n", @fr);
	print TRAN $_ . "\t" . $loc_tran{$_} . "\n";
}
close TRAN;

$output .= 
"\n\n1         ! Number of runs
2         ! Length of run
0         ! 0/1=Monitor method by Iterate#/Time in second
100000    ! Monitor interval in Iterate# / in seconds
0         ! non-Windows version
1         ! Fulllikelihood
3         ! 1/2/3=low/medium/high Precision for Fulllikelihood\n\n";

map { $_ = sprintf("% 15s", $_) } @loci;

$output .= join("", @loci) . "   !Marker names\n";

my @codominant = map { sprintf("% 4s", 0) } 0 .. ($nl-1);

$output .= join("", @codominant) . "   !Marker types, 0/1 = codominant/dominant\n";

my @dropout = map { sprintf("% 8s", sprintf("%.4f",0)) } 0 .. ($nl-1);

$output .= join("", @dropout) . "   !Allelic dropout rate\n";

my @false_rate = map { sprintf("% 8s", 0.0001) } 0 .. ($nl-1);

$output .= join("", @false_rate) . "   !false allele rate\n\n";

if (scalar(@offspring) != 0){
	for $offs (@offspring){
		if ($data_ind{$offs}){
			my @line = split //, $data_ind{$offs};
			unshift @line, $par_offs{$offs};
			$output .= join("", map { sprintf("% 5s", $_) } @line) . "\n";
		}
	}
}

if (scalar(@mothers) == 0 and scalar(@fathers) == 0){
	$output .= "\n0.0  0.0     ! prob. of dad/mum included in the candidates\n";
	$output .= "0   0         ! numbers of candiadte males & females\n\n";
}
elsif (scalar(@mothers) != 0 and scalar(@fathers) == 0){
	$output .= "\n0.0  0.5     !prob. of dad/mum included in the candidates\n";
	$output .= 0 . "  " . scalar(@mothers) . "     ! numbers of candiadte males & females\n\n";
}
elsif (scalar(@mothers) == 0 and scalar(@fathers) != 0){
	$output .= "\n0.5  0.0     !prob. of dad/mum included in the candidates\n";
	$output .= scalar(@fathers) . "  " . 0 . "     ! numbers of candiadte males & females\n\n";
} else {
	$output .= "\n0.5  0.5     !prob. of dad/mum included in the candidates\n";
	$output .= scalar(@fathers) . "  " . scalar(@mothers) . "     ! numbers of candiadte males & females\n\n";
}

if (scalar(@fathers) != 0){
	for $fath (@fathers){
		if ($data_ind{$fath}){
			my @line = split //, $data_ind{$fath};
			unshift @line, $par_offs{$fath};
			$output .= join("", map { sprintf("% 5s", $_) } @line) . "\n";
		}
	}
}

if (scalar(@mothers) != 0){
	for $mums (@mothers){
		if ($data_ind{$mums}){
			my @line = split //, $data_ind{$mums};
			unshift @line, $par_offs{$mums};
			$output .= join("", map { sprintf("% 5s", $_) } @line) . "\n";
		}
	}
	$output .= "\n";
}

$output .= "


0  0          ! #known fater-offspring dyads, paternity exclusion threshold


0  0          ! #known moter-offspring dyads, maternity exclusion threshold


0             ! #known paternal sibship with unknown fathers 


0             ! #known maternal sibship with unknown mothers 


0             ! #known paternity exclusions


0             ! #known maternity exclusions 


0             ! #known paternal sibship exclusions


0             ! #known maternal sibship exclusions\n\n";

# save output to file

open(OUT, "> $name.txt");
print OUT $output;
close OUT;
print "done. ;‹)\n";

# subroutines

sub getind {
	my ($line,$names) = @_;
	my @ind = map { $el = $_; grep { $$line[$_] =~ m/^$el$/ } 0 .. $#$line } @$names;
	return @ind;
}

sub frequency {
	my ($str) = @_;
	my $tot = $str =~ tr/[12]//;
	my $one = $str =~ tr/1//;
	my $two = $str =~ tr/2//;
	my $freq1 = sprintf("%.3f", $one/$tot);
	my $freq2 = sprintf("%.3f", $two/$tot);
	return ($freq1,$freq2);
}

sub reformatInfo {
	my ($str) = @_;
	my $tot = $str =~ tr/[12]//;
	my $in = $tot/2;
	my $two = $str =~ tr/2//;;
	my $freq2 = sprintf("%.3f", $two/$tot);
	my $info = "NS=$in;AF=$freq2";
	return $info;
}

