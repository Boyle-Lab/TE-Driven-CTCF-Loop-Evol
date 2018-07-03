#!/usr/bin/perl

##  File:
##      @(#) beds_to_bedpe.pl
##  Author:
##      Adam G. Diehl  adadiehl@umich.edu
##  Description:
##      Given two bed files corresponding to individual paired end reads,
##      producea a bedpe file with paried reads on one line.
##
#******************************************************************************
#* Copyright (C) Alan Boyle Lab, Prepared by Adam Diehl
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
#
# ChangeLog
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

beds_to_bedpe.pl - Prepare a bedpe file from two bed files.

=head1 SYNOPSIS

  beds_to_bedpe.pl [-version] <pair_1.bed> <pair_2.bed2 [OPTIONS]

=head1 DESCRIPTION

  Given two bed files corresponding to individual paired end reads,
  produce a a bedpe file describing the read pairs.

The options are:

=over 4

=item -version

Displays the version of the program

=item -unpaired

Write unpaired reads to the given file name. By default these are
discarded.

=item -help

Display this message

=back

=head1 COPYRIGHT

Copyright 2016 Adam Diehl

=head1 AUTHOR

Adam Diehl <adadiehl@umich.edu>

=cut

use strict;
use warnings;
use Getopt::Long;

my $Version = 0.01;

my ($bed1_f, $bed2_f) = @ARGV;

if ($#ARGV+1 < 2) {
    usage();
}

my $unpaired_f;

my @getopt_args = (
    '-version',                    # Print out the version and exit
    '-unpaired=s' => \$unpaired_f, # File name for unpaired reads
    '-help'                        # Print usage message
    );

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
    usage();
}

if ( $options{'version'} ) {
    print "$Version\n";
    exit 1;
}

if ( $options{'help'} ) {
    usage();
}

# A hash to keep track of (as yet) unpaired reads:
my %unpaired;

open my $bed_1, '<', $bed1_f || die "Cannot read $bed1_f: $!\n";
open my $bed_2, '<', $bed2_f || die "Cannot read $bed2_f: $!\n";

while (!eof($bed_1) and !eof($bed_2)) {
    my $l1 = <$bed_1>;
    my $l2 = <$bed_2>;

    chomp $l1;
    chomp $l2;
    
    my @tmp1 = split /\s+/, $l1;
    my @tmp2 = split /\s+/, $l2;

    if ($tmp1[3] eq $tmp2[3]) {
	# We have a valid read pair
	my @out = ($tmp1[0], $tmp1[1], $tmp1[2], $tmp2[0], $tmp2[1], $tmp2[2], $tmp1[3], '.', $tmp1[5], $tmp2[5]);
	print_array(\@out);
    } else {
	check_unpaired(\%unpaired, undef, \@tmp1);
	check_unpaired(\%unpaired, undef, \@tmp2);
    }
}

if (!eof($bed_1)) {
    while (<$bed_1>) {
	check_unpaired(\%unpaired, $_, undef);
    }
}

if (!eof($bed_2)) {
    while (<$bed_2>) {
	check_unpaired(\%unpaired, $_, undef);
    }
}

close $bed_1;
close $bed_2;

if (defined($unpaired_f)) {
    open my $out, '>', $unpaired_f || die "Cannot write to $unpaired_f: $!\n";
    foreach my $key (keys(%unpaired)) {
	print_array($unpaired{$key}, undef, $out);
    }
}

exit 0;

sub usage {
    print "$0 - $Version\n";
    exec "pod2text $0";
    exit( 1 );
}

sub check_unpaired {
    my ($unpaired, $l, $lr) = @_;
    my $tmp;
    if (defined($lr)) {
	$tmp = $lr;
    } else {
	chomp $l;
	$tmp = [split /\s+/, $l];
    }
    if ( exists($unpaired{$$tmp[3]}) ) {
	my @out = ($$tmp[0], $$tmp[1], $$tmp[2], $unpaired{$$tmp[3]}->[0], $unpaired{$$tmp[3]}->[1], $unpaired{$$tmp[3]}->[2], $$tmp[3], '.', $$tmp[5], $unpaired{$$tmp[3]}->[5]);
	print_array(\@out);
	delete($unpaired{$$tmp[3]});
    } else {
	$unpaired{$$tmp[3]} = $tmp;
    }
}

sub print_array {
    my ($array, $delim, $fh) = @_;

    if (!defined($delim)) {
	$delim = "\t";
    }

    if (defined($fh)) {
	select $fh;
    }
    
    print join $delim, @{$array}, "\n";

    if (defined($fh)) {
	select STDOUT;
    }
    return 0;
}
