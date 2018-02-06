#!/usr/bin/perl

my $version = "0.3.1";

##Changes from v0.3 -> v0.3.1
# Improved gff/gff3 file processing
# Improved genbank file processing. Now can accept genbank files where locus_ids are in the gene records, but not in the CDS records.

##Changes from v0.2.4 -> v0.3
# Can accept multiple files for the query genome and query genome coordinates. This is useful for multi-chromosome organisms or inclusion of plasmids.
# Updated gff compatibility to read ensembl-formatted gff3 files
# Made core output optional. This step can be slow and perhaps of limited value for most analyses.
# statistics.txt file will now show the AGEnt.pl wrapper script version, not the nucmer_difference.pl version
# make temporary file names unique so that multiple instances of AGEnt can be run in the same directory

##Changes from v0.2.3 -> v0.2.4
# Removed File::Which dependency. Added subroutine to test for whether executable is in PATH that uses only core Perl modules

##Changes from v0.2.2 to v0.2.3
# Fixed bug in genbank file parsing where some genes that span the end of a contig might not appear in results

##Changes from v0.2.1 to v0.2.2
# Fixed bug in nucmer_difference where error message could be produced if all-N region was encountered

##Changes from v0.2 to v0.2.1
# Fixed bug in nucmer_difference.pl where first CDS on each contig was not being output

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
    Copyright (C) 2016-2018 Egon A. Ozer

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
use File::Basename;
use File::Spec::Functions qw ( catfile path );

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
  -q    file(s) of query genome sequence(s) in Fasta or Genbank format. If an
        annotated Genbank formatted file is used, AGEnt will try to extract CDS
        coordinates to separate genes into core and accessory groups.
            AGEnt will try to guess what type of file you have entered based on
            the suffix (Fasta if suffix is .fasta or .fa, Genbank if suffix is
            .gbk or .gb). If you would like to set this manually, use the -Q
            option (see below).
        ** If your genome is split into multiple files (i.e. multiple
        chromosomes or plasmids), you can include all files here separated by
        commas (no spaces).
        Example:
        -q /path/to/chrom_I.fasta,/path/to/chrom_II.fasta,/path/to/plasmid.fasta
            
  -r    file of core genome sequence(s) in Fasta or Genbank format. Example of
            input file would be the \"backbone.fasta\" file output by Spine.
            AGEnt will try to guess what type of file you have entered based on
            the suffix (Fasta if suffix is .fasta or .fa, Genbank if suffix is
            .gbk or .gb). If you would like to set this manually, use the -R
            option (see below).
        
OPTIONS:
  -b    also output core (i.e. non-accessory) sequences and coordinates
        (default: only accessory sequences and coordinates are output)
  -c    path(s) to file(s) containing names and coordinates of ORFs in the
        query genome. This will output a file separating genes into core or
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
        
        If entering ORF-calling results from GeneMark, Prodigal, or GFF/GFF3
            instead, please indicate this with option -f
            
        If an annotated Genbank query file is given (-q), the gene coordinates
            entered here will override the Genbank file annotations
            
        ** If you provided multiple sequence files as query files (-q), you can
            also provide multiple coordinate files here, again separated by
            commas with no spaces
            Example:
            -c /path/to/chrom_I.gff3,/path/to/chrom_II.gff3,/path/to/plasmid.gff3
        
  -f    format of ORF coordinate file given to -c. Options are:
            'glimmer'
            'genemark'
            'prodigal' ('gbk' or web format) 
            'gff' (reads gff or gff3 formats)
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
our ($opt_q, $opt_r, $opt_b, $opt_c, $opt_l, $opt_m, $opt_n, $opt_o, $opt_p, $opt_Q, $opt_R, $opt_s, $opt_v, $opt_f, $opt_w);
getopts('q:r:a:bc:lm:n:o:p:Q:R:s:vf:w');
die "version $version\n" if $opt_v;
die "$license\n" if $opt_l;
die $usage unless ($opt_q and $opt_r);

my $qfile   = $opt_q if $opt_q;
my $rfile   = $opt_r if $opt_r;
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
$coordf = lc($coordf);
die ("ERROR: ORF file format (-f) must be either 'glimmer', 'genemark', 'prodigal', or 'gff'\n") unless ($coordf eq "glimmer" or $coordf eq "genemark" or $coordf eq "prodigal" or $coordf eq "gff");
$spref =~ s/^.*\///;

#if path to nucmer was given, check that the executables are present
my $nuc_loc = is_path("nucmer");
if ($nucpath){
    if (-x "$nucpath/nucmer"){
        $nuc_loc = "$nucpath/nucmer";
    } else {
        print STDERR "WARNING: Could not find nucmer at $nucpath. Searching PATH...\n";
    }
}
die "ERROR: Could not find nucmer in PATH. Make sure MUMmer is installed and executable.\n" if !$nuc_loc;
print STDERR "nuc_loc: $nuc_loc\n" unless $opt_w;

#check that nucmer_difference is present and accessible
unless ($opt_w){
    die "ERROR: Perl must be installed and in your PATH.\n" unless (is_path("perl"));
}
my $home_dir = abs_path($0); #get the absolute path to AGEnt.pl
$home_dir =~ s/\/[^\/]*$//; #strip off AGEnt.pl
die "ERROR: Can't find the required file \"nucmer_difference.pl\". Make sure it is in the \"scripts\" directory with AGEnt.pl ($home_dir/scripts) and has not been renamed.\n" unless (-e "$home_dir/scripts/nucmer_difference.pl");
print STDERR "home_dir = $home_dir\n" unless $opt_w;

#Make sure we're not going to delete anything unintentionally
#my @temp_list = qw (temp_ref.fasta temp_qry.fasta temp_qry.coords.txt temp_align.delta);
#for my $i (0 .. $#temp_list){
#    my $file = $temp_list[$i];
#    die "ERROR: The temporary file $file already exists in this directory. Please delete, rename, or move the file before running AGEnt\n" if -e $file;
#}
my @temp_list;

#Read in coordinates file (if given)
if ($coords){
    my @cfiles = split(",", $coords);
    open (my $crdout, ">$pref\_$spref\_temp_qry.coords.txt") or die "Can't open temporary file: $!\n";
    push @temp_list, "$pref\_$spref\_temp_qry.coords.txt";
    my $count = 0;
    foreach my $cfile (@cfiles){
        open (my $cin, "<$cfile") or die "ERROR: Can't open $cfile: $!\n";
        my $defline; #for Prodigal
        my %crecs; #crecs{contig}{start}{stop}{dir}: [0]is_cds,[1]locus_id,[2]product
        my %ctg_order;
        my $ctg_num = 0;
        my ($lstart, $lstop, $ldir, $lid);
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
                            @{$crecs{$defline}{$lstart}{$lstop}{$ldir}} = (1, $orfid);
                            ($lid, $lstart, $lstop, $ldir) = ("") x 4;
                        }
                        my $id = $1;
                        $id =~ s/\s.*$//g;
                        die "ERROR: Sequence ID '$defline' is used more than once in the annotation file\n" if $ctg_order{$defline};
                        $ctg_num++;
                        $ctg_order{$defline} = $ctg_num;
                        $defline = $id;
                        next;
                    }
                    if ($line =~ m/^\s+(\S+)\s+(complement\()*<*(\d+)<*\.\.>*(\d+)>*\)*\s*$/){
                        my ($type, $start, $stop) = ($1, $3, $4);
                        my $dir = "+";
                        $dir = "-" if $2;
                        $count++;
                        if ($lstart){
                            my $a_count = sprintf("%05d", $count);
                            my $orfid = "$spref-$a_count";
                            $orfid = $lid if $lid;
                            @{$crecs{$defline}{$lstart}{$lstop}{$ldir}} = (1, $orfid);
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
                        die "ERROR: Sequence ID '$defline' is used more than once in the annotation file\n" if $ctg_order{$defline};
                        $ctg_num++;
                        $ctg_order{$defline} = $ctg_num;
                        next;
                    }
                    if ($line =~ m/^\s*\d/){
                        my ($id, $strand, $start, $stop, $x1, $x2) = split(' ', $line);
                        $id = sprintf("%05d", $id);
                        $start =~ s/[<>]//;
                        $stop =~ s/[<>]//;
                        my $out1 = "$spref-$id";
                        $count++;
                        @{$crecs{$defline}{$start}{$stop}{$strand}} = (1, $out1);
                    }
                } elsif ($coordf eq "gff"){
                    ## gff is tough because there seems to be very little standardization of tags for locus IDs and gene products
                    ## Need to try to set some priorities. 
                    next if $line =~ m/^#/;
                    next if $line =~ m/^\s*$/;
                    my ($contig, $x1, $type, $start, $stop, $x2, $dir, $x3, $rest) = split("\t", $line);
                    unless ($ctg_order{$contig}){
                        $ctg_num++;
                        $ctg_order{$contig} = $ctg_num;
                    }
                    unless ($crecs{$contig}{$start}{$stop}{$dir}){
                        @{$crecs{$contig}{$start}{$stop}{$dir}} = (0); #initialize the record as non-cds
                    }
                    if ($type eq "gene"){ # in Ensembl gff3 files, records corresponding to locus_id and product in gbk files are found in the gene record, not the CDS record
                        if ($rest =~ m/(?:\A|;)gene_id="*([^;"]+)/){
                            ${$crecs{$contig}{$start}{$stop}{$dir}}[1] = $1;
                        }
                        if ($rest =~ m/(?:\A|;)description="*([^;"]+)/){
                            ${$crecs{$contig}{$start}{$stop}{$dir}}[2] = $1;
                        }
                    }
                    next unless $type eq "CDS";
                    $count++;
                    ${$crecs{$contig}{$start}{$stop}{$dir}}[0] = 1;
                    #if the CDS record has "locus" or "name" records, will use these as the locus id and product names
                    if ($rest =~ m/(?:\A|;)locus="*([^;"]+)/){
                        ${$crecs{$contig}{$start}{$stop}{$dir}}[1] = $1;
                    }
                    if ($rest =~ m/(?:\A|;)name="*([^;"]+)/){
                        ${$crecs{$contig}{$start}{$stop}{$dir}}[2] = $1;
                    }
                    #assign a locus ID if none was found
                    unless (${$crecs{$contig}{$start}{$stop}{$dir}}[1]){
                        my $a_count = sprintf("%05d", $count);
                        $contig =~ s/^.*\|//;
                        my $orfid = "$spref-$a_count";
                        if ($rest =~ m/(?:\A|;)ID="*([^;"]+)/){
                            $orfid = $1;
                        }
                        ${$crecs{$contig}{$start}{$stop}{$dir}}[1] = $orfid;
                    }
                } else { #Glimmer
                    if ($line =~ m/^>\s*(\S+)/){
                        $defline = $1;
                        die "ERROR: Sequence ID '$defline' is used more than once in the annotation file\n" if $ctg_order{$defline};
                        $ctg_num++;
                        $ctg_order{$defline} = $ctg_num;
                        next;
                    }
                    next if $line =~ m/^\s*$/; #skip blank lines
                    $count++;
                    my @tmp = split(' ', $line);
                    my ($orf, $start, $stop) = @tmp;
                    my $out1 = "$spref-$orf";
                    my $dir = "+";
                    if ($start > $stop){
                        $dir = "-";
                        ($start, $stop) = ($stop, $start);
                    }
                    @{$crecs{$defline}{$start}{$stop}{$dir}} = (1, $out1);
                }
            }
        }
        if ($coordf eq "prodigal"){
            if ($lstart){
                my $a_count = sprintf("%05d", $count);
                my $orfid = "$spref-$a_count";
                $orfid = $lid if $lid;
                @{$crecs{$defline}{$lstart}{$lstop}{$ldir}} = (1, $orfid);
            }
        }
        close ($cfile);
        #print to crdout
        foreach my $cid (sort{$ctg_order{$a} <=> $ctg_order{$b}} keys %ctg_order){
            foreach my $start (sort{$a <=> $b} keys %{$crecs{$cid}}){
                foreach my $stop (sort{$a <=> $b} keys %{$crecs{$cid}{$start}}){
                    foreach my $dir (sort keys %{$crecs{$cid}{$start}{$stop}}){
                        my ($is_cds, $lid, $prod) = @{$crecs{$cid}{$start}{$stop}{$dir}};
                        if ($is_cds){
                            print $crdout "$lid\t$cid\t$start\t$stop\t$dir\t";
                            print $crdout "$prod" if $prod;
                            print $crdout "\n";
                        }
                    }
                }
            }
        }
        
    }
    close ($crdout);
    print STDERR "Read in $count total annotation records from coordinate file(s)\n";
    print STDERR "<br>\n" if $opt_w;
    $coords = "$pref\_$spref\_temp_qry.coords.txt";
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
    my @qfiles = split(",",$qfile);
    foreach my $qf (@qfiles){
        (my $suffix) = $qf =~ m/\.([^.]+)$/;
        die ("ERROR: Cannot determine query file type for '$qf'. No suffix found. Please set -Q with sequence file type\n") if (!$suffix);
        my $usuffix = uc($suffix);
        die ("ERROR: Query file '$qf' suffix ($suffix) not recognized, cannot determine file type. Please set -Q with sequence file type\n") if (!$types{$usuffix});
        if ($man_q){
            die("ERROR: All query files must be of the same type (fasta or genbank)\n") if $man_q ne substr($usuffix, 0, 1);
        }
        $man_q = substr($usuffix, 0, 1);
    }
}
if (!$man_r){
    (my $suffix) = $rfile =~ m/\.([^.]+)$/;
    die ("ERROR: Cannot determine reference file type. No suffix found. Please set -R with appropriate sequence file type\n") if (!$suffix);
    my $usuffix = uc($suffix);
    die ("ERROR: Reference file suffix ($suffix) not recognized, cannot determine file type. Please set -R with approrpriate sequence file type\n") if (!$types{$usuffix});
    $man_r = substr($usuffix, 0, 1);
}

#read in reference file
push @temp_list, "$pref\_$spref\_temp_ref.fasta";
if ($man_r eq "G"){
    print STDERR "Processing reference genbank file ...\n";
    print STDERR "<br>\n" if $opt_w;
    my $rpref = "$pref\_$spref\_temp_ref";
    my $status = gbk_convert($rfile, $rpref, "R");
    die "ERROR: Reference Genbank file does not contain DNA sequence. Please check file.\n" if ($status == 2);
} else {
    print STDERR "Processing reference fasta file ...\n";
    print STDERR "<br>\n" if $opt_w;
    open (my $r_in, "<$rfile") or die "ERROR: Can't open reference sequence file $rfile: $!\n";
    open (my $r_out, ">$pref\_$spref\_temp_ref.fasta") or die "ERROR: Can't open temporary output file: $!\n";
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
my $r_source = "$pref\_$spref\_temp_ref.fasta";

#read in query file(s)
push @temp_list, "$pref\_$spref\_temp_qry.fasta";
my @qfiles = split(",", $qfile);
foreach my $qf (@qfiles){
    if ($man_q eq "G"){
        print STDERR "Processing query genome genbank file ...\n";
        print STDERR "<br>\n" if $opt_w;
        my $status = gbk_convert($qf, "$pref\_$spref\_temp_qry", "Q");
        die "ERROR: Query Genbank file does not contain DNA sequence. Please check file.\n" if ($status == 2);
        die "ERROR: CDS records in query file missing \"locus_tag\" tags. Please check file and visit http://vfsmspineagent.fsm.northwestern.edu/gbk_reformat.cgi for conversion tool.\n" if ($status == 3);
    } else {
        print STDERR "Processing query genome fasta file ...\n";
        print STDERR "<br>\n" if $opt_w;
        open (my $q_in, "<$qf") or die "ERROR: Can't open query sequence file $qf: $!\n";
        open (my $q_out, ">>$pref\_$spref\_temp_qry.fasta") or die "ERROR: Can't open temporary output file: $!\n";
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
}
my $q_source = "$pref\_$spref\_temp_qry.fasta";

#Align with nucmer
push @temp_list, "$pref\_$spref\_temp_align.delta";
my $n_command = "$nuc_loc --maxmatch -p $pref\_$spref\_temp_align $r_source $q_source";
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
    our ($opt_d, $opt_V);
    local $opt_d = "$pref\_$spref\_temp_align.delta";
    local $opt_m = $minalgn;
    local $opt_s = $minsize;
    local $opt_q = $q_source;
    local $opt_c = $coords if $coords;
    local $opt_o = $pref;
    local $opt_p = $spref;
    local $opt_v = "";
    local $opt_w = 1 if $opt_w;
    $opt_r = "";
    local $opt_r = 1 if $opt_b;
    local $opt_V = $version;
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
sub is_path {
    ## Subroutine based on StackOverflow post by Sinan Unur (https://stackoverflow.com/a/8243770)
    my $exe = shift;
    my @path = path;
    my @pathext = ( q{} );
    if ($^O eq 'MSWin32'){
        push @pathext, map { lc } split /;/, $ENV{PATHEXT};
    }
    for my $dir (@path){
        for my $ext (@pathext){
            my $f = catfile $dir, "$exe$ext";
            return ($f) if -x $f;
        }
    }
    return();
}

sub gbk_convert{
    my $file = shift;
    my $filepref = shift;
    my $filetype = shift;
    open (my $gbkin, "<", $file) or die "ERROR: Can't open $file: $!\n";
    open (my $seqout, ">>$filepref.fasta") or die "Can't open temporary file: $!\n";
    my $loccount = 0;
    my $seqcount = 0;
    my ($c_id, $c_seq);
    my $is_prod;
    my @tags;
    my %crecs;
    my @ctg_order;
    my $reading = 1; # 1 = front material, 2 = annotations, 3 = sequence
    
    while (my $fline = <$gbkin>){
        $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
        my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
        while (@lines){
            my $line = shift @lines;
            next if $line =~ m/^\s*$/;
            if ($line =~ m/^LOCUS\s+\S*\s+\d+\sbp/){
                if ($reading == 2){ #no ORIGIN sequence record was found between LOCUS records
                    close ($seqout);
                    unlink ("$filepref.fasta");
                    return (2);
                }
                if ($reading == 3){
                    if ($c_seq and $c_id){
                        print $seqout ">$c_id\n$c_seq\n";
                        $c_seq = "";
                        $reading = 1;
                    } else {
                        close ($seqout);
                        unlink ("$filepref.fasta");
                        return (2);
                    }
                }
            }
            if ($line =~ m/^\/\//){ #reached the end of the file (or record)
                if ($c_seq and $c_id){
                    print $seqout ">$c_id\n$c_seq\n";
                    $c_seq = "";
                    $reading = 1;
                } else {
                    close ($seqout);
                    unlink ("$filepref.fasta");
                    return (2);
                }
            }
            if ($reading == 1){
                if ($line =~ m/^LOCUS\s+([^\s]+)/){
                    $seqcount++;
                    if ($line =~ m/^LOCUS\s+(\S+)\s+\d+ bp/){
                        $c_id = $1;
                    } else {
                        $c_id = "rec$seqcount";
                    }
                    push @ctg_order, $c_id;
                    next;
                }
                if ($line =~ m/^FEATURES\s+Location\/Qualifiers/){
                    $reading = 2;
                    next;
                }
            } elsif ($reading == 2){
                if ($line =~ m/^\s+(\S+)\s+(complement\()*[<>]*(\d+)<*\.\.[<>]*(\d+)>*\)*\s*$/){
                    $is_prod = "";
                    my ($type, $start, $stop) = ($1, $3, $4);
                    my $dir = "+";
                    $dir = "-" if $2;
                    unless ($crecs{$c_id}{$start}{$stop}{$dir}){
                        @{$crecs{$c_id}{$start}{$stop}{$dir}} = (0);
                    }
                    if ($type eq "CDS"){
                        ${$crecs{$c_id}{$start}{$stop}{$dir}}[0] = 1;
                        $loccount++;
                    }
                    if (@tags){
                        my ($o_start, $o_stop, $o_dir) = @tags;
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[3] if $tags[3];
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[2] = $tags[4] if $tags[4];
                        $loccount++;
                    }
                    @tags = ($start, $stop, $dir);
                    next;
                }
                if ($line =~ m/^ORIGIN\s*$/){
                    $is_prod = "";
                    if (@tags){
                        my ($o_start, $o_stop, $o_dir) = @tags;
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[3] if $tags[3];
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[2] = $tags[4] if $tags[4];
                        $loccount++;
                    }
                    undef @tags;
                    $reading = 3;
                    next
                }
                if ($line =~ m/^\s+\/(\S+)=\"*([^"]*)\"*/){
                    $is_prod = "";
                    my ($key, $val) = ($1, $2);
                    if ($key eq "locus_tag"){
                        $tags[3] = $val;
                    }
                    if ($key eq "product"){
                        $tags[4] = $val;
                        $is_prod = 1;
                    }
                    next;
                }
                if ($is_prod){
                    $line =~ s/^\s*//;
                    $line =~ s/"*\s*$//;
                    $tags[4] .= " $line";
                }
            } elsif ($reading == 3){
                $line =~ s/\d//g;
                $line =~ s/\s//g;
                $c_seq .= $line;
                next;
            }
        }
    }
    if ($c_seq and $c_id){
        print $seqout ">$c_id\n$c_seq\n";
        $c_seq = "";
        $reading = 1;
    }
    close ($gbkin);
    close ($seqout);
    
    if ($filetype eq "Q" and !$opt_c){
        if ($loccount > 0){
            $coords = "$filepref.coords.txt";
            push @temp_list, "$filepref.coords.txt";
            open (my $crdout, ">>$coords") or die "Can't open temporary file: $!\n";
            foreach my $cid (@ctg_order){
                foreach my $start (sort{$a <=> $b} keys %{$crecs{$cid}}){
                    foreach my $stop (sort{$a <=> $b} keys %{$crecs{$cid}{$start}}){
                        foreach my $dir (sort keys %{$crecs{$cid}{$start}{$stop}}){
                            my ($is_cds, $lid, $prod) = @{$crecs{$cid}{$start}{$stop}{$dir}};
                            if ($is_cds){
                                unless ($lid){
                                    close ($crdout);
                                    unlink ($coords);
                                    print STDERR "ERROR: CDS at $start..$stop on contig $cid has no locus_id\n";
                                    print STDERR "<br>\n" if $opt_w;
                                    return(3);
                                }
                                print $crdout "$lid\t$cid\t$start\t$stop\t$dir\t";
                                print $crdout "$prod" if $prod;
                                print $crdout "\n";
                            }
                        }
                    }
                }
            }
            close ($crdout);
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
    }
    return (0);
}

