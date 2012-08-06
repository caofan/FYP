#!/usr/bin/perl
print "create exon record frm refgene: <CHR>\n";
unless(open (REF,"refGene.txt")) {die};
unless(open (ANN,$ARGV[0].".ann.location")) {die};

#unless (open (OUTP,">$ARGV[0].+")) {die};
#unless (open (OUTPB,">$ARGV[0].+.bin")) {die};
#binmode OUTPB;
#unless (open (OUTM,">$ARGV[0].-")) {die};
#unless (open (OUTMB,">$ARGV[0].-.bin")) {die};
#binmode OUTMB;
#$Chr=$ARGV[0];

my $Gen_Count=0;my %Chromosomes;my @ChrPlus;my @ChrMinus;my @ChrPlusR;my @ChrMinusR;$i=0;
while ($Chr=<ANN>)
{
	my $OUT;
	$Chr=<ANN>; chomp $Chr;
	$Chromosomes{$Chr}=$i;
	open ($ChrPlus[$i],">$Chr.+.bin") || die;
	binmode $ChrPlus[$i];#$OUT;
	#push (@ChrPlus,$OUT);
	open ($ChrPlusR[$i],">$Chr.+") || die;
	#push (@ChrPlusR,$OUT);
	open ($ChrMinus[$i],">$Chr.-.bin") || die;
	binmode $ChrMinus[$i];#$OUT;
	#push (@ChrMinus,$OUT);
	open ($ChrMinusR[$i],">$Chr.-") || die;
	#push (@ChrMinusR,$OUT);
	$i++;
}
while( $Record=<REF>)
{
	chomp $Record;
	split ('\t',$Record);
	$geneName=$_[0];#           "Name of gene as it appears in Genome Browser."
	$name=$_[1]; #              "Name of gene"
	$chrom=$_[2]; #             "Chromosome name"
	$Chr=$Chromosomes{$chrom};

	#if ($Chr eq $chrom)
	{
		$strand=$_[3]; #            "+ or - for strand"
		$txStart=$_[4]; #           "Transcription start position"
		$txEnd=$_[5]; #             "Transcription end position"
		$cdsStart=$_[6]; #          "Coding region start"
		$cdsEnd=$_[7]; #            "Coding region end"
		$exonCount=$_[8]; #         "Number of exons"
		@exonStarts=split(',',$_[9]); #"Exon start positions"
		@exonEnds=split(',',$_[10]); #  "Exon end positions"
		if ($strand eq '+')
		{
			$PlusR= $ChrPlusR[$Chr];
			$Plus= $ChrPlus[$Chr];
			print  $PlusR "@\n";
			for ($i=0;$i<$exonCount;$i++)
			{
				print $PlusR "$exonStarts[$i]\n";
				print $Plus pack("V",$exonStarts[$i]);
				print $PlusR "$exonEnds[$i]\n";
				print $Plus pack("V",$exonEnds[$i]);
			}
		}
		else
		{
			$MinusR=$ChrMinusR[$Chr];
			$Minus=$ChrMinus[$Chr];
			print  $MinusR "@\n";
			for ($i=0;$i<$exonCount;$i++)
			{
				print $MinusR "$exonStarts[$i]\n";
				print $Minus pack("V",$exonStarts[$i]);
				print $MinusR "$exonEnds[$i]\n";
				print $Minus pack("V",$exonEnds[$i]);
			}
		}
	}
}
