#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

my $usage = "usage: $0 A.bed B.bed C.bed ...\n\n";

unless (@ARGV) {
    die $usage;
}

my @bed_files = @ARGV;

if (scalar @bed_files == 1) {

    if (-s $bed_files[0]) {
        # allow for a file listing the various files.
        @bed_files = `cat $bed_files[0]`;
        chomp @bed_files;
    }
    else {
        die $usage;
    }
}

=header_format

0       chromosome_start_end_strand
1       score


=cut


main: {


    my %data;

    foreach my $file (@bed_files) {

        open (my $fh, $file) or die "Error, cannot open file $file";
        my $header = <$fh>; # ignore it
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $acc = $x[0];
            my $count = $x[1];
            $data{$acc}->{$file} = $count;
        }
        close $fh;
    }

    my @filenames = @bed_files;
    foreach my $file (@filenames) {
        $file = basename($file);
    }


    print join("\t", "", @filenames) . "\n";
    foreach my $acc (keys %data) {

        print "$acc";

        foreach my $file (@bed_files) {

            my $count = $data{$acc}->{$file};
            unless (defined $count) {
                $count = "0";
            }

            print "\t$count";

        }

        print "\n";

    }


    exit(0);
}


