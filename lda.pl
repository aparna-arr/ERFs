#!/usr/bin/env perl
use warnings;
use strict;

my ($mapwigDir, $genes_file, @classes) = @ARGV;
die "usage: $0 <mapwig directory> <mm9 genes> <bed filenames of classes [ex. top_by_ChIPseq middle_by_ChIPseq bottom_by_ChIPseq]>\n" unless @ARGV;

opendir(DIR, $mapwigDir);
my @files = readdir(DIR);
closedir(DIR);

my %hash;

my $header = "";
foreach my $file (sort @files)
{
	next if (-d $file || $file =~ /rep[12]/ || $file !~/\.mapwig$/);

	my ($wig, $bed) = $file =~ /(.+?(?:_GSE\d+|_merge)*)_(.+?)(?:_combined_reps|\.)/;

	my $next = 1;
	foreach my $class (@classes)
	{
		if ($bed eq $class)
		{
			$next = 0;
			last;
		}
	}	

	next if ($next == 1);

	#  make sure it is sorted
	`sort -k 1,1 -k 2,2n $mapwigDir/$file > tmp ; mv tmp $file.sort`;

	push(@{$hash{$bed}{files}}, "$file.sort");

	$hash{$bed}{header} .= "$wig\t";
}
$header .= $hash{$classes[0]}{header} . "length\tgene_num\tsample";

my @beds;

foreach my $bed (sort keys %hash)
{	
	my %lines;

	`cut -f 4 $hash{$bed}{files}[0] > tmp`;
	for (my $i = 1; $i < @{$hash{$bed}{files}}; $i++)
	{
		`cut -f 4 $hash{$bed}{files}[$i] | paste tmp - > tmp2 ; mv tmp2 tmp`;
	}

	# FIXME this is wrong, those cols are no longer start and end
#	`awk '{print \$0 "\t" \$3 - \$2 "\t$bed"}' tmp > tmp2 ; mv tmp2 tmp_$bed`;

	warn "file is [$hash{$bed}{files}[0]]\n";
#	`awk '{print \$3 - \$2}' $hash{$bed}{files}[0] > test`;

# length of peaks
	`awk '{print \$3 - \$2}' $hash{$bed}{files}[0] | paste tmp - > tmp2; mv tmp2 tmp_$bed`;

# number of genes
	`bedtools intersect -c -a $hash{$bed}{files}[0] -b $genes_file | awk '{print \$NF "\t$bed"}' | paste tmp_$bed - > tmp2; mv tmp2 tmp_$bed`;


	push (@beds, $bed);
}

`cp tmp_$beds[0] tmp`;
for (my $j = 1; $j < @beds; $j++)
{
	`cat tmp tmp_$beds[$j] > tmp2; mv tmp2 tmp`
}

`echo \"$header\" | cat - tmp > lda_data.txt`;

`Rscript pca.R`;
