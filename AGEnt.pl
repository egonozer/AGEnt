#!/usr/bin/perl

my $version = "0.1.3";

##Changes from v0.1.2 -> 0.1.3
# Upgraded nucmer_difference.pl to v0.5.1

##Changes from v0.1 -> 0.1.2
# Accept files with all types of line endings (Unix, Mac, PC)

my $license = "
    AGEnt.pl
    Copyright (C) 2014 Egon A. Ozer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see [http://www.gnu.org/licenses/].
";

use strict;
use warnings;
use Cwd 'abs_path';
use File::Which;

$|++;

my $usage = "
AGEnt.pl

This is a wrapper script to run the AGEnt algorithm

PREREQUISITES:
- Perl 5.10 or above
- MUMmer version 3.22 or above
- Mac OSX or Linux. No guarantees that this will work on Windows or other
  operating systems.

REQUIRED:
  -q    file of query sequence(s) in Fasta or Genbank format. If an annotated
        Genbank formatted file is used, AGEnt will try to extract CDS
        coordinates to separate genes into core and accessory groups.
            AGEnt will try to guess what type of file you have entered based on
            the suffix (Fasta if suffix is .fasta or .fa, Genbank if suffix is
            .gbk or .gb). If you would like to set this manually, use the -Q
            option (see below).
  -r    file of reference sequence(s) in Fasta or Genbank format
            AGEnt will try to guess what type of file you have entered based on
            the suffix (Fasta if suffix is .fasta or .fa, Genbank if suffix is
            .gbk or .gb). If you would like to set this manually, use the -R
            option (see below).
        
OPTIONS:
  -a    minimum overlap with accessory coordinates, in percent, for a query
        gene to be called accessory
        (default: 50).
  -c    path to file containing names and coordinates of genes in the query
        genome. This will output a file separating genes into core or
        accessory categories (\"_orfs.txt\").
        File should be in \"Glimmer\" format, i.e.
            >contig_name_1 
            orf_ID_1<space(s)>start_coordinate<space(s)>stop_coordinate
            orf_ID_2<space(s)>start_coordinate<space(s)>stop_coordinate
            >contig_name_2
            orf_ID_3<space(s)>start_coordinate<space(s)>stop_coordinate
            etc...
            
        Contig names should match those in file given by option -q.
        All coordinates are 1-based.
        Coordinates assuming a circular contig that cross the origin will
            give incorrect results.
        Best if all ORF IDs are unique (i.e. don't restart count
            every contig).
        If an annotated Genbank query file is given (-q), the gene coordinates
            entered here will override the Genbank file annotations
  -l    print license information and quit
  -m    minimum alignment identity between query and reference to be called
        core, in percent
        (default: 85)
  -n    full path to folder containing MUMmer scripts and executables,
        i.e. /home/applications/MUMmer/bin
        (default: tries to find MUMmer in your PATH)
  -o    prefix for output files
        (default: \"output\")
  -p    prefix for output sequences
        (default: same as given by option -o)
  -Q    Manual override of query file type. Enter \"F\" for Fasta or \"G\" for
        Genbank
  -R    Manual override of reference file type. Enter \"F\" for Fasta or \"G\" for
        Genbank
  -s    minimum size of fragment to output
        (default: 10)
  -v    print version information and quit

";

# command line processing
use Getopt::Std;
our ($opt_q, $opt_r, $opt_a, $opt_c, $opt_l, $opt_m, $opt_n, $opt_o, $opt_p, $opt_Q, $opt_R, $opt_s, $opt_v);
getopts('q:r:a:c:lm:n:o:p:Q:R:s:v');
die "version $version\n" if $opt_v;
die "$license\n" if $opt_l;
die $usage unless ($opt_q and $opt_r);

my $qfile   = $opt_q if $opt_q;
my $rfile   = $opt_r if $opt_r;
my $minover = $opt_a ? $opt_a : 50;
my $coords  = $opt_c if $opt_c;
my $minalgn = $opt_m ? $opt_m : 85;
my $nucpath = $opt_n if $opt_n;
my $pref    = $opt_o ? $opt_o : "output";
my $spref   = $opt_p ? $opt_p : $pref;
my $man_q   = $opt_Q if $opt_Q;
my $man_r   = $opt_R if $opt_R;
my $minsize = $opt_s ? $opt_s : 10;

die ("ERROR: Minimum overlap (-a) must be a non-negative number\n") if ($minover =~ m/[^0123456789.]/);
die ("ERROR: Minimum overlap (-a) must be between 0 and 100\n") if ($minover > 100);
die ("ERROR: Minimum aligment (-m) must be a non-negative number\n") if ($minalgn =~ m/[^0123456789.]/);
die ("ERROR: Minimum aligment (-m) must be between 1 and 100\n") if ($minalgn > 100 or $minalgn < 1);
die ("ERROR: Minimum size (-s) must be a non-negative integer\n") if ($minsize =~ m/[^0123456789]/);
die ("ERROR: Minimum size (-s) must be a non-negative integer\n") if ($minsize == 0);
if ($man_q){
    $man_q = uc(substr($man_q, 0, 1)); #take first letter of input only and make it uppercase
    die ("ERROR: Query file type (-Q) must be either \"F\" (Fasta) or \"G\" (Genbank)\n") if ($man_q =~ m/[^FG]/);
}
if ($man_r){
    $man_r = uc(substr($man_r, 0, 1)); #take first letter of input only and make it uppercase
    die ("ERROR: Reference file type (-R) must be either \"F\" (Fasta) or \"G\" (Genbank)\n") if ($man_r =~ m/[^FG]/);
}

#if path to nucmer was given, check that the executables are present
my $nuc_loc = which("nucmer");
if ($nucpath){
    if (-e "$nucpath/nucmer"){
        $nuc_loc = "$nucpath/nucmer";
    } else {
        print STDERR "WARNING: Could not find nucmer at $nucpath. Searching PATH...\n";
    }
}
die "ERROR: Could not find nucmer in PATH. Make sure MUMmer is installed and executable.\n" if !$nuc_loc;
print STDERR "nuc_loc: $nuc_loc\n";

#check that nucmer_difference is present and accessible
die "ERROR: Perl must be installed and in your PATH.\n" unless (which("perl"));
my $home_dir = abs_path($0); #get the absolute path to AGEnt.pl
$home_dir =~ s/\/[^\/]*$//; #strip off AGEnt.pl
die "ERROR: Can't find the required file \"nucmer_difference.pl\". Make sure it is in the \"scripts\" directory with AGEnt.pl ($home_dir/scripts) and has not been renamed.\n" unless (-e "$home_dir/scripts/nucmer_difference.pl");
print STDERR "home_dir = $home_dir\n";

#Make sure we're not going to delete anything unintentionally
my @temp_list = qw (temp_ref.fasta temp_qry.fasta temp_qry.coords.txt temp_align.delta);
for my $i (0 .. $#temp_list){
    my $file = $temp_list[$i];
    die "ERROR: The temporary file $file already exists in this directory. Please delete, rename, or move the file before running AGEnt\n" if -e $file;
}

#check input sequence file types (trying to be flexible)
my %types = (
    FASTA   => 1,
    FA      => 1,
    FNA     => 1,
    FSA     => 1,
    FFN     => 1,
    GBK     => 1,
    GB      => 1,
    GENBANK => 1
);
if (!$man_q){
    (my $suffix) = $qfile =~ m/\.([^.]+)$/;
    die ("ERROR: Cannot determine query file type. No suffix found. Please set -Q with sequence file type\n") if (!$suffix);
    my $usuffix = uc($suffix);
    die ("ERROR: Query file suffix ($suffix) not recognized, cannot determine file type. Please set -Q with sequence file type\n") if (!$types{$usuffix});
    $man_q = substr($usuffix, 0, 1);
}
if (!$man_r){
    (my $suffix) = $rfile =~ m/\.([^.]+)$/;
    die ("ERROR: Cannot determine reference file type. No suffix found. Please set -R with appropriate sequence file type\n") if (!$suffix);
    my $usuffix = uc($suffix);
    die ("ERROR: Reference file suffix ($suffix) not recognized, cannot determine file type. Please set -R with approrpriate sequence file type\n") if (!$types{$usuffix});
    $man_r = substr($usuffix, 0, 1);
}

#read in reference file
if ($man_r eq "G"){
    print STDERR "Processing reference genbank file ...\n";
    my $rpref = "temp_ref";
    my $status = gbk_convert($rfile, $rpref, "R");
    die "ERROR: Reference Genbank file does not contain DNA sequence. Please check file.\n" if ($status == 2);
} else {
    print STDERR "Processing reference fasta file ...\n";
    open (my $r_in, "<$rfile") or die "ERROR: Can't open reference sequence file $rfile: $!\n";
    open (my $r_out, ">temp_ref.fasta") or die "ERROR: Can't open temporary output file: $!\n";
    while (my $fline = <$r_in>){
        $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
        my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
        while (@lines){
            my $line = shift @lines;
            print $r_out "$line\n";
        }
    }
    close ($r_in);
    close ($r_out);
}
my $r_source = "temp_ref.fasta";

#read in query file
open (my $q_in, "<$qfile") or die "ERROR: Can't open query sequence file $qfile: $!\n";
if ($man_q eq "G"){
    print STDERR "Processing query genbank file ...\n";
    my $qpref = "temp_qry";
    my $status = gbk_convert($qfile, $qpref, "Q");
    die "ERROR: Query Genbank file does not contain DNA sequence. Please check file.\n" if ($status == 2);
} else {
    print STDERR "Processing query fasta file ...\n";
    open (my $q_in, "<$qfile") or die "ERROR: Can't open query sequence file $rfile: $!\n";
    open (my $q_out, ">temp_qry.fasta") or die "ERROR: Can't open temporary output file: $!\n";
    while (my $fline = <$q_in>){
        $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
        my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
        while (@lines){
            my $line = shift @lines;
            print $q_out "$line\n";
        }
    }
    close ($q_in);
    close ($q_out);
}
my $q_source = "temp_qry.fasta";

#Align with nucmer
my $n_command = "$nuc_loc --maxmatch -p temp_align $r_source $q_source";
print STDERR "\nRunning $n_command ...\n";
my @n_result = `$n_command > /dev/null 2>&1`;
my $error = $?;
die "ERROR: Nucmer failure (", join(",",@n_result), ")\n" if $error;

#Run nucmer_difference
print STDERR "Running nucmer_difference ...\n";
my $return;
{
    our ($opt_d);
    local $opt_d = "temp_align.delta";
    local $opt_m = $minalgn;
    local $opt_s = $minsize;
    local $opt_q = $q_source;
    local $opt_c = $coords if $coords;
    local $opt_a = $minover;
    local $opt_o = $pref;
    local $opt_p = $spref;
    local $opt_v = "";
    $return = do "$home_dir/scripts/nucmer_difference.pl";
}
unless ($return){
    die "ERROR: Couldn't run nucmer_difference.pl: $@\n" if $@;
    die "ERROR: Couldn't run nucmer_difference.pl: $!\n" unless defined $return;
}

print STDERR "Cleaning up ...\n";
for my $i (0 .. $#temp_list){
    my $file = $temp_list[$i];
    unlink ($file) if -e $file;
}
print STDERR "\nFinished!\n";

#------------------------
sub gbk_convert{
    my $file = shift;
    my $filepref = shift;
    my $filetype = shift;
    open (my $seqout, "> $filepref.fasta") or die "Can't open temporary file: $!\n";
    my $crdstring;
    open (my $gbkin, "<", $file) or die "ERROR: Can't open $file: $!\n";
    my (@list);
    my @seqlist;
    my $loccount = 0;
    my ($c_id, $c_seq);
    my ($start, $stop);
    my $is_origin;
    while (my $fline = <$gbkin>){
        $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
        my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
        while (@lines){
            my $line = shift @lines;
            if ($line =~ m/^LOCUS\s+([^\s]+)/){
                $is_origin = "";
                if ($c_id){ #error if there is no sequence in the genbank file
                    if (!$c_seq) {
                        close ($seqout);
                        unlink ("$filepref.fasta");
                        return (2);
                    } 
                    print $seqout ">$c_id\n$c_seq\n";
                    $c_seq = "";
                    ($start, $stop) = ("") x 2;
                }
                $c_id = $1;
                $crdstring .= ">$c_id\n";
                next;
            }
            if ($line =~ m/^\s+CDS\s+/){
                $line =~ m/(\d+)\.\.(\d+)/;
                ($start, $stop) = ($1, $2);
                ($start, $stop) = ($stop, $start) if $line =~ m/complement\(/; #reverse coordinates if the record is complemented
                next;
            }
            if ($start and $line =~ m/^\s+\/locus_tag="/){
                $loccount++;
                $line =~ m/"([^"]+)"/;
                my $lid = $1;
                $crdstring .= "$lid\t$start\t$stop\n";
                ($start, $stop) = ("") x 2;
                next;
            }
            if ($line =~ m/^ORIGIN/){
                $is_origin = 1;
                next;
            }
            if ($is_origin and $line =~ m/^\s*\d+\s(.*)/){
                my $seqline = $1;
                $seqline =~ s/\s//g;
                if ($seqline =~ m/[^ACTGNactgn]/){
                    ###since sequences may contain some IUPAC ambiguity codes other than "N" I'll ignore this for now
                    #return (4);
                }
                $c_seq .= $seqline;
                next;
            }
        }
    }
    if ($c_id){
        if (!$c_seq){
            close ($seqout);
            unlink ("$filepref.fasta");
            return (2);
        }
        print $seqout ">$c_id\n$c_seq\n";
        $c_seq = "";
        ($start, $stop) = ("") x 2;
    }
    close ($gbkin);
    close ($seqout);
    if ($filetype eq "Q" and !$coords){
        if ($loccount > 0){
            open (my $crdout, "> $filepref.coords.txt") or die "Can't open temporary file: $!\n";
            print $crdout "$crdstring";
            close ($crdout);
            $coords = "$filepref.coords.txt";
        } else { #if no CDS records and/or no locus_ids were found in the genbank file
            print STDERR
            "
            ****
            FYI: No CDS annotations or no locus_id identifiers were found in
            the query genbank file and no gene coordinate file was given with
            option -c. No core/accessory gene information will be output.
            If you need to add locus_id tags to an annotated genbank file
            (for example, if RAST was used to annotate the sequence), please
            visit http://vfsmspineagent.fsm.northwestern.edu/gbk_reformat.cgi
            for online or download-able file conversion tool.
            ****
            ";
        }
    }
    return (0);
}

