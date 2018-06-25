#!/usr/bin/perl
use strict;
use warnings;
use Statistics::Basic qw/mean median/;

use PerlIO::gzip;
use Getopt::Long;

my $usage="
	Options:
		-case <STR>		The tags file of case;
		-control <STR>		The tags file of control;
		-o <STR>		The output file;
		-h|?			Help!
";
my ($case,$control,$out,$help);
GetOptions(
	'case=s' => \$case,
	'control=s' => \$control,
	'o=s' => \$out,
	'w=s' => \$win,
	'h|?' => \$help
);

if($help or !$out){die "$usage\n";}
