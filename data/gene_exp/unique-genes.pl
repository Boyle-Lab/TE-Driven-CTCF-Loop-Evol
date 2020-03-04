#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $usage = "
unique_genes.pl

Given a bed file with multiple gene isoforms, produce a single, 
maximal-length bed entry for each unique gene.

USAGE:

unique_genes.pl <input.bed>

<input.bed>
    Standard tab-delimited bed file. Only the name, chrom,
    chromStart and chromEnd fields are considered. The
    remaining fields will be taken from the first record for
    each gene.

OPTIONS:

--help
    Show this message.
--convert-ids <conversion_table.txt>
    Convert IDs in input BED based on a conversion table. The first
    column should contain input IDs with output IDs in the second
    column.
\n";

my $infile = $ARGV[0];

my $help = 0;
my $conv_f;
GetOptions(
    "help" => \$help,
    "convert-ids=s" => \$conv_f
    );


if ($help || $#ARGV+1 < 1) {
    die "$usage\n";
}

my %conv_table;
if (defined($conv_f)) {
    %conv_table = read_into_hash($conv_f, 0);
}

my %data = read_into_hash_array($infile, 3, \%conv_table);

my %res;
foreach my $gene (sort(keys(%data))) {
    my @row = @{$data{$gene}->[0]};
    my $start = $row[1];
    my $end = $row[2];    
    for (my $i = 0; $i <= $#{$data{$gene}}; $i++) {
	if ($data{$gene}->[$i]->[1] < $start) {
	    $start = $data{$gene}->[$i]->[1];
	}
	if ($data{$gene}->[$i]->[2] > $end) {
	    $end = $data{$gene}->[$i]->[2];
	}
    }
    $row[1] = $start;
    $row[2] = $end;
    $res{$gene} = \@row;
}

print_hash_of_arrays(\%res);

sub read_into_hash {
    my ($file, $key_idx, $delim) = @_;

    if (!defined($delim)) {
        $delim = "\t";
    }
    my %hash;
    open my $INFILE, '<', $file || die "Cannot read $file: $!\n";
    while (<$INFILE>) {
        if ($_ =~ /^#/ || $_ !~ m/\S+/g) {
            # Skip blank and commented lines
            next;
        }
        chomp;
        my @tmp = split /$delim/, $_;
        my $key = $tmp[$key_idx];
        $hash{$key} = \@tmp;
    }
    close $INFILE;
    return %hash;
}

sub read_into_hash_array {
    # Read a set of records into a hash where duplicate key rows are put into
    # an array.
    my ($file, $key_idx, $conv_table, $delim) = @_;
    if (!defined($delim)) {
        $delim = "\t";
    }
    my %hash;
    open my $INFILE, '<', $file || die "Cannot read $file: $!\n";
    while (<$INFILE>) {
        if ($_ =~ /^#/ || $_ !~ m/\S+/g) {
            # Skip blank and commented lines      
            next;
        }
        chomp;
        my @tmp = split /$delim/, $_;
        my $key = $tmp[$key_idx];
	if (defined($conv_table)) {
	    $key = $conv_table->{$key}->[1];
	    $tmp[$key_idx] = $key;
	}
        if (exists($hash{$key})) {
            push @{$hash{$key}}, \@tmp;
        } else {
            my @t = (\@tmp);
            $hash{$key} = \@t;
        }
    }
    close $INFILE;
    return %hash;
}

sub print_array {
    my ($array, $delim, $fh, $term, $nd_char) = @_;
    # input array, delimiter char, filehandle, terminator char, not-defined char
    if (!defined($delim)) {
        $delim = "\t";
    }
    if (!defined($term)) {
        $term = "\n";
    }
    if (defined($fh)) {
        select $fh;
    } else {
        select STDOUT;
    }
    my $i;
    for ($i = 0; $i < $#{$array}; $i++) {
        if (defined(${$array}[$i])) {
            print "${$array}[$i]$delim";
        } else {
            if (defined($nd_char)) {
                print "$nd_char$delim";
            } else {
                print STDERR "WARNING: Some fields in output array not defined with no default. Skipping!\n";
                next;
            }
        }
    }
    if (defined(${$array}[$#{$array}])) {
        print "${$array}[$#{$array}]$term";
    } else {
        if (defined($nd_char)) {
            print "$nd_char$term";
        } else {
            print STDERR "WARNING: Some fields in output array not defined with no default. Skipping!\n";
            print "$term";
        }
    }
    if (defined($fh)) {
        select STDOUT;
    }
    return 0;
}

sub print_hash_of_arrays {
    # Print a hash of arrays to given filehandle with given delimiting
    my ($hash, $delim, $fh) = @_;
    if (!defined($delim)) {
        $delim = "\t";
    }
    if (!defined($fh)) {
        $fh = \*STDOUT;
    }

    if (exists(${$hash}{header})) {
        print_array(${$hash}{header}, $delim, $fh);
    }
    foreach my $key (sort(keys(%{$hash}))) {
        if ($key eq "header") {
            next;
        }
        print_array(${$hash}{$key}, $delim, $fh);
    }
}
