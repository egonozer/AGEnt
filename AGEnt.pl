#!/usr/bin/perl

my $version = "0.2";

##Changes from v0.1.3 to v0.2
# Updated output file formatting to match Spine v0.2
# Accept genemark, prodigal, and gff ORF files (in addition to glimmer)
# Web status updates directly from AGEnt
# Add gff coordinate input option

##Changes from v0.1.2 -> 0.1.3
# Upgraded nucmer_difference.pl to v0.5.1

##Changes from v0.1 -> 0.1.2
# Accept files with all types of line endings (Unix, Mac, PC)

my $license = "
    AGEnt.pl
    Copyright (C) 2016 Egon A. Ozer

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
  -q    file of query genome sequence(s) in Fasta or Genbank format. If an
        annotated Genbank formatted file is used, AGEnt will try to extract CDS
        coordinates to separate genes into core and accessory groups.
            AGEnt will try to guess what type of file you have entered based on
            the suffix (Fasta if suffix is .fasta or .fa, Genbank if suffix is
            .gbk or .gb). If you would like to set this manually, use the -Q
            option (see below).
  -r    file of core genome sequence(s) in Fasta or Genbank format. Example of
            input file would be the \"backbone.fasta\" file output by Spine.
            AGEnt will try to guess what type of file you have entered based on
            the suffix (Fasta if suffix is .fasta or .fa, Genbank if suffix is
            .gbk or .gb). If you would like to set this manually, use the -R
            option (see below).
        
OPTIONS:
  -c    path to file containing names and coordinates of ORFs in the query
        genome. This will output a file separating genes into core or
        accessory categories (\"_loci.txt\").
        Default file format is \"Glimmer\" format, i.e.
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
        
        If entering ORF-calling results from GeneMark or Prodigal instead,
            please indicate this with option -f
            
        If an annotated Genbank query file is given (-q), the gene coordinates
            entered here will override the Genbank file annotations
  -f    format of ORF coordinate file given to -c. Options are:
            'glimmer'
            'genemark'
            'prodigal' ('gbk' or web format) 
            'gff'
        (default: glimmer)
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
our ($opt_q, $opt_r, $opt_a, $opt_c, $opt_l, $opt_m, $opt_n, $opt_o, $opt_p, $opt_Q, $opt_R, $opt_s, $opt_v, $opt_f, $opt_w);
getopts('q:r:a:c:lm:n:o:p:Q:R:s:vf:w');
die "version $version\n" if $opt_v;
die "$license\n" if $opt_l;
die $usage unless ($opt_q and $opt_r);

my $qfile   = $opt_q if $opt_q;
my $rfile   = $opt_r if $opt_r;
#my $minover = $opt_a ? $opt_a : 50;
#  -a    minimum overlap with accessory coordinates, in percent, for a query
#        gene to be called accessory
#        (default: 50).
my $coords  = $opt_c if $opt_c;
my $coordf  = $opt_f ? $opt_f : "glimmer";
my $minalgn = $opt_m ? $opt_m : 85;
my $nucpath = $opt_n if $opt_n;
my $pref    = $opt_o ? $opt_o : "output";
my $spref   = $opt_p ? $opt_p : $pref;
my $man_q   = $opt_Q if $opt_Q;
my $man_r   = $opt_R if $opt_R;
my $minsize = $opt_s ? $opt_s : 10;

#die ("ERROR: Minimum overlap (-a) must be a non-negative number\n") if ($minover =~ m/[^0123456789.]/);
#die ("ERROR: Minimum overlap (-a) must be between 0 and 100\n") if ($minover > 100);
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
die ("ERROR: ORF file format (-f) must be either 'glimmer', 'genemark', 'prodigal', or 'gff'\n") unless ($coordf eq "glimmer" or $coordf eq "genemark" or $coordf eq "prodigal" or $coordf eq "gff");
$spref =~ s/^.*\///;

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
print STDERR "nuc_loc: $nuc_loc\n" unless $opt_w;

#check that nucmer_difference is present and accessible
unless ($opt_w){
    die "ERROR: Perl must be installed and in your PATH.\n" unless (which("perl"));
}
my $home_dir = abs_path($0); #get the absolute path to AGEnt.pl
$home_dir =~ s/\/[^\/]*$//; #strip off AGEnt.pl
die "ERROR: Can't find the required file \"nucmer_difference.pl\". Make sure it is in the \"scripts\" directory with AGEnt.pl ($home_dir/scripts) and has not been renamed.\n" unless (-e "$home_dir/scripts/nucmer_difference.pl");
print STDERR "home_dir = $home_dir\n" unless $opt_w;

#Make sure we're not going to delete anything unintentionally
my @temp_list = qw (temp_ref.fasta temp_qry.fasta temp_qry.coords.txt temp_align.delta);
for my $i (0 .. $#temp_list){
    my $file = $temp_list[$i];
    die "ERROR: The temporary file $file already exists in this directory. Please delete, rename, or move the file before running AGEnt\n" if -e $file;
}

#Read in coordinates file (if given)
if ($coords){
    open (my $cin, "<$coords") or die "ERROR: Can't open $coords: $!\n";
    open (my $crdout, ">temp_qry.coords.txt") or die "Can't open temporary file: $!\n";
    my $defline; #for Prodigal
	my ($lid, $lstart, $lstop, $ldir); #for Prodigal
	my $count = 0; #for Prodigal
    while (my $fline = <$cin>){
        $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
        my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
        while (@lines){
            my $line = shift @lines;
            if ($coordf eq "prodigal"){
                if ($line =~ m/seqhdr="([^"]*)/ or $line =~ m/DEFINITION\s+>(.*)/){
                    if ($lstart){
                        my $a_count = sprintf("%05d", $count);
                        my $orfid = "$spref-$a_count";
                        $orfid = $lid if $lid;
                        print $crdout "$orfid\t$defline\t$lstart\t$lstop\t$ldir\t\n";
                        ($lid, $lstart, $lstop, $ldir) = ("") x 4;
                    }
                    my $id = $1;
                    $id =~ s/\s.*$//g;
                    $defline = $id;
                    next;
                }
                if ($line =~ m/^\s+(\S+)\s+(complement\()*<*(\d+)\.\.(\d+)>*\)*\s*$/){
                    my ($type, $start, $stop) = ($1, $3, $4);
                    my $dir = "+";
                    $dir = "-" if $2;
                    $count++;
                    if ($lstart){
                        my $a_count = sprintf("%05d", $count);
                        my $orfid = "$spref-$a_count";
                        $orfid = $lid if $lid;
                        print $crdout "$orfid\t$defline\t$lstart\t$lstop\t$ldir\t\n";
                        $lid = ("");
                    }
                    ($lstart, $lstop, $ldir) = ($start, $stop, $dir);
                    next;
                }
                if ($line =~ m/^\s*\/note="ID=([^;]+)/){
                    $lid = $1;
                    next;
                }
            } elsif ($coordf eq "genemark"){
                $line =~ s/^\s*//; #removes any spaces at the beginning of the line
                $line =~ s/\s*$//; #removes any spaces at the end of the line
                if ($line =~ m/FASTA definition line: (.*)/){
                    $defline = $1;
                    $defline =~ s/\s.*$//;
                    next;
                }
                if ($line =~ m/^\s*\d/){
                    my ($id, $strand, $start, $stop, $x1, $x2) = split(' ', $line);
                    $id = sprintf("%05d", $id);
                    $start =~ s/[<>]//;
                    $stop =~ s/[<>]//;
                    my $out1 = "$spref-$id";
                    print $crdout "$out1\t$defline\t$start\t$stop\t$strand\t\n";
                }
            } elsif ($coordf eq "gff"){
                next if $line =~ m/^#/;
                next if $line =~ m/^\s*$/;
                my ($contig, $x1, $type, $start, $stop, $x2, $dir, $x3, $rest) = split("\t", $line);
                next unless $type eq "CDS";
                $count++;
                my $a_count = sprintf("%05d", $count);
                $contig =~ s/^.*\|//;
                my $orfid = "$spref-$a_count";
                if ($rest =~ m/ID="*([^;"]+)/){
                    $orfid = $1;
                }
                print $crdout "$orfid\t$contig\t$start\t$stop\t$dir\t";
                if ($rest =~ m/product="*([^;"]+)/){
                    print $crdout "$1";
                }
                print $crdout "\n";
            } else { #Glimmer
                if ($line =~ m/^>\s*(\S+)/){
                    $defline = $1;
                    next;
                }
                next if $line =~ m/^\s*$/; #skip blank lines
                my @tmp = split(' ', $line);
                my ($orf, $start, $stop) = @tmp;
                my $out1 = "$spref-$orf";
                my $dir = "+";
                if ($start > $stop){
                    $dir = "-";
                    ($start, $stop) = ($stop, $start);
                }
                print $crdout "$out1\t$defline\t$start\t$stop\t$dir\t\n";
            }
        }
    }
    if ($coordf eq "prodigal"){
        if ($lstart){
            my $a_count = sprintf("%05d", $count);
            my $orfid = "$spref-$a_count";
            $orfid = $lid if $lid;
            print $crdout "$orfid\t$defline\t$lstart\t$lstop\t$ldir\t\n";
        }
    }
    close $coords;
    $coords = "temp_qry.coords.txt";
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
    print STDERR "<br>\n" if $opt_w;
    my $rpref = "temp_ref";
    my $status = gbk_convert($rfile, $rpref, "R");
    die "ERROR: Reference Genbank file does not contain DNA sequence. Please check file.\n" if ($status == 2);
} else {
    print STDERR "Processing reference fasta file ...\n";
    print STDERR "<br>\n" if $opt_w;
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
    print STDERR "Processing query genome genbank file ...\n";
    print STDERR "<br>\n" if $opt_w;
    my $qpref = "temp_qry";
    my $status = gbk_convert($qfile, $qpref, "Q");
    die "ERROR: Query Genbank file does not contain DNA sequence. Please check file.\n" if ($status == 2);
    die "ERROR: CDS records in query file missing \"locus_tag\" tags. Please check file and visit http://vfsmspineagent.fsm.northwestern.edu/gbk_reformat.cgi for conversion tool.\n" if ($status == 3);
} else {
    print STDERR "Processing query genome fasta file ...\n";
    print STDERR "<br>\n" if $opt_w;
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
print STDERR "<br>\n<h1>" if $opt_w;
print STDERR "\nRunning Nucmer ...\n";
print STDERR "</h1><br>\n" if $opt_w;
my @n_result = `$n_command > /dev/null 2>&1`;
my $error = $?;
die "ERROR: Nucmer failure (", join(",",@n_result), ")\n" if $error;

#Run nucmer_difference
print STDERR "<h1>" if $opt_w;
print STDERR "Running nucmer_difference ...\n";
print STDERR "</h1><br>\n" if $opt_w;
my $return;
{
    our ($opt_d);
    local $opt_d = "temp_align.delta";
    local $opt_m = $minalgn;
    local $opt_s = $minsize;
    local $opt_q = $q_source;
    local $opt_c = $coords if $coords;
    #local $opt_a = $minover;
    local $opt_o = $pref;
    local $opt_p = $spref;
    local $opt_v = "";
    local $opt_w = 1 if $opt_w;
    $return = do "$home_dir/scripts/nucmer_difference.pl";
}
unless ($return){
    die "ERROR: Couldn't run nucmer_difference.pl: $@\n" if $@;
    die "ERROR: Couldn't run nucmer_difference.pl: $!\n" unless defined $return;
}

print STDERR "Cleaning up ...\n";
print STDERR "<br>\n" if $opt_w;
for my $i (0 .. $#temp_list){
    my $file = $temp_list[$i];
    unlink ($file) if -e $file;
}
print STDERR "<br>\n<strong>" if $opt_w;
print STDERR "\nFinished!\n";
print STDERR "</strong><br>\n" if $opt_w;

#------------------------
sub gbk_convert{
    my $file = shift;
    my $filepref = shift;
    my $filetype = shift;
    open (my $seqout, ">$filepref.fasta") or die "Can't open temporary file: $!\n";
    my $crdout;
    if ($filetype eq "Q"){
        open ($crdout, ">$filepref.coords.txt") or die "Can't open temporary file: $!\n";
    }
    my $crdstring;
    open (my $gbkin, "<", $file) or die "ERROR: Can't open $file: $!\n";
    my (@list);
    #my @seqlist;
    my $loccount = 0;
    my $seqcount = 0;
    my ($c_id, $c_seq);
    my ($start, $stop);
    my $is_origin;
    my $is_feature;
    my ($is_cds, $is_prod);
    my %tags;
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
                }
                $seqcount++;
                if ($line =~ m/^LOCUS\s+(\S+)\s+\d+ bp/){
                    $c_id = $1;
                } else {
                    $c_id = "rec$seqcount";
                }
                next;
            }
            if ($line =~ m/^FEATURES\s+Location\/Qualifiers/){
                $is_feature = 1;
                next;
            }
            next unless $is_feature;
            
            if ($line =~ m/^\s+(\S+)\s+(complement\()*<*(\d+)\.\.(\d+)>*\)*\s*$/){
                my ($type, $start, $stop) = ($1, $3, $4);
                my $dir = "+";
                $dir = "-" if $2;
                if (%tags and $filetype eq "Q"){
                    return(3) unless $tags{'locus_tag'}; #no locus_tag was present on the last record
                    my ($o_id, $o_start, $o_stop, $o_dir) = ($tags{'locus_tag'}, $tags{'start'}, $tags{'stop'}, $tags{'dir'});
                    my $o_prod = $tags{'product'} if $tags{'product'};
                    print $crdout "$o_id\t$c_id\t$o_start\t$o_stop\t$o_dir\t";
                    print $crdout "$o_prod" if $o_prod;
                    print $crdout "\n";
                    $loccount++;
                }
                undef %tags;
                $is_prod = "";
                $is_cds = "";
                if ($type eq "CDS"){
                    ($start, $stop) = ($stop, $start) if $start > $stop;
                    $tags{'start'} = $start;
                    $tags{'stop'} = $stop;
                    $tags{'dir'} = $dir;
                    $is_cds = 1;
                }
                next;
            }
            if ($is_cds and $line !~ m/^ORIGIN/){
                if ($line =~ m/^\s+\/(\S+)=\"*([^"]*)\"*/){
                    my ($key, $val) = ($1, $2);
                    if ($key eq "locus_tag" or $key eq "product"){
                        $val =~ s/\s*$//;
                        $tags{$key} = $val;
                        $is_prod = 1 if $key eq "product";
                        next;
                    } else {
                        $is_prod = "";
                    }
                } else {
                    if ($is_prod){
                        $line =~ s/^\s*//;
                        $line =~ s/\"//g;
                        $tags{'product'} .= " $line";
                    }
                }
                next;
            }
            if ($line =~ m/^ORIGIN/){
                if (%tags  and $filetype eq "Q"){
                    return(3) unless $tags{'locus_tag'}; #no locus_tag was present on the last record
                    my ($o_id, $o_start, $o_stop, $o_dir) = ($tags{'locus_tag'}, $tags{'start'}, $tags{'stop'}, $tags{'dir'});
                    my $o_prod = $tags{'product'} if $tags{'product'};
                    print $crdout "$o_id\t$c_id\t$o_start\t$o_stop\t$o_dir";
                    print $crdout "\t$o_prod" if $o_prod;
                    print $crdout "\n";
                    $loccount++;
                }
                undef %tags;
                $is_prod = "";
                $is_cds = "";
                #if ($start){ #no locus_tag was present on the last record
                #    return (3);
                #}
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
        #($start, $stop) = ("") x 2;
    }
    close ($gbkin);
    close ($seqout);
    if ($filetype eq "Q" and !$coords){
        if ($loccount > 0){
            $coords = "$filepref.coords.txt";
        } else { #if no CDS records and/or no locus_ids were found in the genbank file
            if ($opt_w){
                print STDERR "<strong>****</strong><br>\n";
                print STDERR "<p>FYI: No CDS annotations or no locus_id identifiers were found in the query genbank file and no gene coordinate file was given. No core/accessory gene information will be output. If you need to add locus_id tags to an annotated genbank file (for example, if RAST was used to annotate the sequence), please visit http://vfsmspineagent.fsm.northwestern.edu/gbk_reformat.cgi for online or download-able file conversion tool.</p><\n>";
                print STDERR "<strong>****</strong><br>\n";
            } else {
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
        close ($crdout);
    }
    return (0);
}

