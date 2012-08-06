#!/usr/bin/perl
print "Decode Formatter.....\n";
if (($ARGV[0] eq "-h") || ($ARGV[0] eq "-?"))
{
	print "Format output of decode \n";
	print "Command line :\n";
	print "fmt.pl <input file> <format file> <output file>\n";
	exit;
}

open FILE, $ARGV[0] or die;
open FORMATFILE, $ARGV[1] or die;
open(OUTPUT,'>'.$ARGV[2]) or die;

$Line=<FILE>;
@Line=split(/\t/,$Line);
if ($Line[0] != '@') {die;}
$StrlenH= $Line[1];
$StrlenT= $Line[2];
chomp $StrlenT;
chomp $StrlenH;
$Strlen=$StrlenH;

while ($FLine=<FORMATFILE>) 
{
	chop ($FLine);
	$FLine =~ s/^\s*//;     
	if ( ($FLine !~ /^#/) && ($FLine ne "") )
	{    
		($Name, $Value) = split (/=/, $FLine);    
		$Name =~ s/\s*$//;
		$Hash{$Name} = $Value;                  
	}
}
foreach $Key (keys(%Hash))
{
$Hash{$Key}=~ s/\\t/\t/g;
$Hash{$Key}=~ s/\\n/\n/g;
}

$Head=$Hash{"HEAD"};
$Match=$Hash{"BODY"};
$Hit=$Hash{"HIT"};
print OUTPUT $Hash{"TITLE"};

if(!$StrlenT) {Print_Normal();} else {Print_Pair();}

sub Print_Pair()
{
	while ($Line =<FILE>)
	{
		@Header=split('\t',$Line,2);
		$Description=$Header[0];#DES
		$Tag=$Header[1];#TAG
		chomp $Tag;

		$Formatted_Head=$Head;
		$Formatted_Head=~s/%DES%/$Description/g;
		$Formatted_Head=~s/%TAG%/$Tag/g;
		$Formatted_Head=~s/%LENH%/$StrlenH/g;
		$Formatted_Head=~s/%LENT%/$StrlenT/g;
		print OUTPUT "$Formatted_Head";

		$Line=<FILE>;
		while($Line ne "@\n")
		{
			@Data_Line=split("\t",$Line,5);
			$Tag_Number= $Data_Line[0];#TAGNUMBER
			$Mismatch_Count= $Data_Line[1];#MISMATCHNO
			$Head_Tail= $Data_Line[2];#HT
			if($Head_Tail=='H') {$Strlen=$StrlenH;} else {$Strlen=$StrlenT;}
			$Orientation= $Data_Line[3];#SIGN
			@Rest=split /\|/,$Data_Line[4];
			$Mismatch_Description= $Rest[0];#MISINFO
			@Last=split("\t",$Rest[1],2);
			$Hit_Count=$Last[0];#HITCOUNT
			$Locations= $Last[1]; chomp $Locations;#LOC

			$Formatted_Match=$Match;
			$Formatted_Match=~s/%HT%/$Head_Tail/g;
			$Formatted_Match=~s/%DES%/$Description/g;
			$Formatted_Match=~s/%SIGN%/$Orientation/g;
			$Formatted_Match=~s/%TAG%/$Tag/g;
			$Formatted_Match=~s/%TAGNUMBER%/$Tag_Number/g;
			$Formatted_Match=~s/%MISMATCHNO%/$Mismatch_Count/g;
			$Formatted_Match=~s/%MISINFO%/$Mismatch_Description/g;
			$Formatted_Match=~s/%HITCOUNT%/$Hit_Count/g;
			$Formatted_Match=~s/%LOC%/$Locations/g;
			$Formatted_Match=~s/%LENH%/$StrlenH/g;
			$Formatted_Match=~s/%LENT%/$StrlenT/g;
			$Formatted_Match=~s/%LEN%/$Strlen/g;


			print OUTPUT "$Formatted_Match";
			@Hits=split /;/,$Locations;
			foreach (@Hits)
			{
				@Location=split /:/,$_;
				$Chromosome=$Location[0];
				if(!$Chromosome) {$Offset=$Location[0];} else {$Offset=$Location[1];}
				$Formatted_Hit=$Hit;
				$Formatted_Hit=~s/%HT%/$Head_Tail/g;
				$Formatted_Hit=~s/%OFFSET%/$Offset/g;
				$Formatted_Hit=~s/%CHR%/$Chromosome/g;
				$Formatted_Hit=~s/%DES%/$Description/g;
				$Formatted_Hit=~s/%SIGN%/$Orientation/g;
				$Formatted_Hit=~s/%TAG%/$Tag/g;
				$Formatted_Hit=~s/%TAGNUMBER%/$Tag_Number/g;
				$Formatted_Hit=~s/%MISMATCHNO%/$Mismatch_Count/g;
				$Formatted_Hit=~s/%MISINFO%/$Mismatch_Description/g;
				$Formatted_Hit=~s/%HITCOUNT%/$Hit_Count/g;
				$Formatted_Hit=~s/%LOC%/$Locations/g;
				$Formatted_Hit=~s/%HIT%/$_/g;
				$Formatted_Hit=~s/%LENH%/$StrlenH/g;
				$Formatted_Hit=~s/%LENT%/$StrlenT/g;
				$Formatted_Hit=~s/%LEN%/$Strlen/g;

				print OUTPUT "$Formatted_Hit";
			}
			if (!($Line=<FILE>)) {exit};
		}	
	}
}

sub Print_Normal()
{
	while ($Line =<FILE>)
	{
		@Header=split('\t',$Line);
		$Description=$Header[0];#DES
		chomp $Header[1];
		$Orientation=chop($Header[1]);#SIGN
		$Tag=$Header[1];#TAG
		chomp $Tag;

		$Formatted_Head=$Head;
		$Formatted_Head=~s/%DES%/$Description/g;
		$Formatted_Head=~s/%SIGN%/$Orientation/g;
		$Formatted_Head=~s/%TAG%/$Tag/g;
		$Formatted_Head=~s/%LEN%/$Strlen/g;
		print OUTPUT "$Formatted_Head";

		$Line=<FILE>;
		while($Line ne "@\n")
		{
			@Data_Line=split("\t",$Line,3);
			$Tag_Number= $Data_Line[0];#TAGNUMBER
			$Mismatch_Count= $Data_Line[1];#MISMATCHNO
			@Rest=split /\|/,$Data_Line[2];
			$Mismatch_Description= $Rest[0];#MISINFO
			@Last=split("\t",$Rest[1],2);
			$Hit_Count=$Last[0];#HITCOUNT
			$Locations= $Last[1]; chomp $Locations;#LOC

			$Formatted_Match=$Match;
			$Formatted_Match=~s/%DES%/$Description/g;
			$Formatted_Match=~s/%SIGN%/$Orientation/g;
			$Formatted_Match=~s/%TAG%/$Tag/g;
			$Formatted_Match=~s/%TAGNUMBER%/$Tag_Number/g;
			$Formatted_Match=~s/%MISMATCHNO%/$Mismatch_Count/g;
			$Formatted_Match=~s/%MISINFO%/$Mismatch_Description/g;
			$Formatted_Match=~s/%HITCOUNT%/$Hit_Count/g;
			$Formatted_Match=~s/%LOC%/$Locations/g;
			$Formatted_Match=~s/%LEN%/$Strlen/g;

			print OUTPUT "$Formatted_Match";
			@Hits=split /;/,$Locations;
			foreach (@Hits)
			{
				$Formatted_Hit=$Hit;
				@Location=split /:/,$_;
				$Chromosome=$Location[0];
				if(!$Chromosome) {$Offset=$Location[0];} else {$Offset=$Location[1];}

				$Formatted_Hit=~s/%OFFSET%/$Offset/g;
				$Formatted_Hit=~s/%CHR%/$Chromosome/g;
				$Formatted_Hit=~s/%DES%/$Description/g;
				$Formatted_Hit=~s/%SIGN%/$Orientation/g;
				$Formatted_Hit=~s/%TAG%/$Tag/g;
				$Formatted_Hit=~s/%TAGNUMBER%/$Tag_Number/g;
				$Formatted_Hit=~s/%MISMATCHNO%/$Mismatch_Count/g;
				$Formatted_Hit=~s/%MISINFO%/$Mismatch_Description/g;
				$Formatted_Hit=~s/%HITCOUNT%/$Hit_Count/g;
				$Formatted_Hit=~s/%LOC%/$Locations/g;
				$Formatted_Hit=~s/%HIT%/$_/g;
				$Formatted_Hit=~s/%LEN%/$Strlen/g;

				print OUTPUT "$Formatted_Hit";
			}
			if (!($Line=<FILE>)) {exit};
		}	
	}
}

