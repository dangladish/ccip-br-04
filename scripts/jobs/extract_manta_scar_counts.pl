#!/usr/bin/perl
# Extracts the per-tow feeding-scar count (ScarsCount) from the Manta sheet of the COTS manta/cull/RHIS
# workbook and writes data/raw/manta_tows/manta_scar_counts.csv (tow_id, scar_code, scars_count).
# The count is recorded for tows from 14 June 2025; only tows with a count are written.
# The sheet XML is streamed directly because reading it with readxl is very slow.
#
# Usage (needs perl and unzip; run from the repository root):
#   perl scripts/jobs/extract_manta_scar_counts.pl data/raw/manta_tows/260529-COTS-Manta-Cull-RHIS.xlsx \
#        > data/raw/manta_tows/manta_scar_counts.csv
#
# Manta sheet columns used: A CrownOfThornsStarfishSurveillanceId, M FeedingScarCountRangeCode, N ScarsCount.
use strict; use warnings;
my $xlsx = shift or die "usage: $0 <workbook.xlsx>\n";
sub unz { my $part = shift; open(my $fh, "-|", "unzip", "-p", $xlsx, $part) or die "cannot read $part: $!"; local $/; my $s = <$fh>; close $fh; return $s; }
# locate the worksheet part of the sheet named "Manta"
my $wbxml = unz("xl/workbook.xml");
my ($rid) = $wbxml =~ /<sheet\b[^>]*\bname="Manta"[^>]*\br:id="([^"]+)"/ or die "no sheet named Manta\n";
my $target;
for my $rel (unz("xl/_rels/workbook.xml.rels") =~ /<Relationship\b[^>]*>/g) {
  if ($rel =~ /\bId="\Q$rid\E"/ && $rel =~ /\bTarget="([^"]+)"/) { $target = $1; last; }
}
die "no worksheet part for sheet Manta\n" unless $target;
$target =~ s{^/?(xl/)?}{xl/};
my @ss;
{
  open(my $fh, "-|", "unzip", "-p", $xlsx, "xl/sharedStrings.xml") or die "cannot read sharedStrings: $!";
  local $/ = "</si>";
  while (<$fh>) { my $s = ""; while (/<t[^>]*>([^<]*)<\/t>/g) { $s .= $1 } push @ss, $s if /<si>/; }
}
open(my $in, "-|", "unzip", "-p", $xlsx, $target) or die "cannot read $target: $!";
local $/ = "</row>";
print "tow_id,scar_code,scars_count\n";
my $n = 0;
while (<$in>) {
  my %v;
  while (/<c r="([A-Z]+)(\d+)"([^>]*?)(?:\/>|>(?:<v>([^<]*)<\/v>)?<\/c>)/g) {
    my ($col, $r, $attr, $val) = ($1, $2, $3, $4);
    next unless defined $val;
    $val = $ss[$val] if $attr =~ /t="s"/;
    $v{$col} = $val; $v{row} = $r;
  }
  next unless defined $v{row} && $v{row} > 1;
  next unless defined $v{N} && $v{N} =~ /^\d+(\.\d+)?$/ && defined $v{A} && defined $v{M};   # text such as "NA" = not recorded
  print join(",", $v{A}, $v{M}, $v{N}), "\n"; $n++;
}
print STDERR "wrote $n tows with a ScarsCount\n";
