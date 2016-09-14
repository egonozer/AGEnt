#!/usr/bin/perl

my $license = "
    nucmer_difference.pl
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

my $version = "0.5.1";

## Changes from version 0.5
## Fixed a bug where contigs in the query that were aligned to the reference by nucmer, but where none of the alignment identities were above the given threshold (-m) were not included in the accessory genome.

## Changes from version 0.4
## Overhauled the procedure for parsing nucmer alignment data to identify aligned and non-aligned regions. The old version was missing duplicate sequence in the query, i.e. if two query seqs aligned to the same reference seq, the sequence would still be output as "non-core"
## No longer outputs segments that only consist of N's
## -q option is no longer required. Program will use query file listed in the delta file header unless another file is specified in the options.
## Outputs the percentage of the reference sequence that was found in the query
## Improved speed with large, single-contig genomes by keeping %crds in memory and not re-loading between accessory segments if the contig id is the same as the last

## Changes from version 0.3
## Added option to identify genes contained in the core or accessory segments based on an input file of genome coordinates of orfs in the query genome. This had previously been perfomred by a separate program "nucmer_difference_to_orfs.pl"
## Changed default minimum fragment size to 10 from 500

## Changes from version 0.2
## Ability to filter out sequences with low percentage alignment (but I don't recommend it)
## If a reference sequence or portion of reference sequence aligns to multiple query sequences, prioritizes longest alignments

## Changes from version 0.1
##  output key file with query sequence borders 

use strict;
use Cwd;
use warnings;

## Usage
my $usage = "
nucmer_difference.pl (version: $version)

Description: takes nucmer .delta file as input and outputs coordinates and 
             sequences from query sequences not present in the reference

Copyright (C) 2014  Egon A. Ozer
This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you
are welcome to redistribute it under certain conditions;
see LICENSE.txt for details

Usage: nucmer_difference.pl <options>
required:
    -d      nucmer .delta file

optional:
    -m      minimum alignment identity, in % (default: 1)
    -s      minimum size of fragment to output (default: 10)
    -q      fasta file of query sequences (if different than what was input to
            nucmer)    
    -c      path to file containing names and coordinates of genes in the query
            genome. This will output a file separating genes into core or
            accessory categories.
            File should be in \"Glimmer\" format, i.e.
                >contig_name_1 
                orf_ID_1<tab>start_coordinate<tab>stop_coordinate
                orf_ID_2<tab>start_coordinate<tab>stop_coordinate
                >contig_name_2
                orf_ID_3<tab>start_coordinate<tab>stop_coordinate
                etc...
                
                contig names should match those in file given by option -q
                all coordinates are 1-based
                coordinates assuming a circular contig that cross the origin will
                    give incorrect results
    -a      minimum overlap with accessory coordinates, in percent, for a gene
            in file given by -c to be called accessory (default: 50)
    -o      prefix for output files (default: \"out\")
    -p      sequence prefix (default: same as given by option -o)
    -v      verbose output

";

## command line processing.
use Getopt::Std;
use vars qw( $opt_d $opt_q $opt_o $opt_m $opt_c $opt_a $opt_s $opt_p $opt_v );
getopts('d:q:o:m:s:p:c:a:v');
die $usage unless ($opt_d);

my ($infile, $inquery, $minalign, $minsize, $pref, $seqpref, $coord, $minacc);

$infile     = $opt_d if $opt_d;
$inquery    = $opt_q if $opt_q;
$minalign   = $opt_m ? $opt_m : 1 ;
$minsize    = $opt_s ? $opt_s : 10 ;
$pref       = $opt_o ? $opt_o : "out" ;
$seqpref    = $opt_p ? $opt_p : $pref ;
$coord      = $opt_c if $opt_c;
$minacc     = $opt_a ? $opt_a : 50 ;

## Read in delta file
my @matches;
my %seen;
open (my $in, "<", $infile);
my ($files, $head_id, $head_leng, $ref_id, $ref_leng);
my $inref;
print STDERR "Reading in delta file\n";
while (my $line = <$in>){
    chomp $line;
    if (!$files){ #takes the file paths of the reference and query sequences from the delta file header line;
        my ($ref, $qry) = split(' ', $line);
        $inquery = $qry if (!$opt_q);
        $inref = $ref;
        $files = 1;
        next;
    }
    if ($line =~ /^>/) {
        my @head = split(' ', $line);
        ($head_id, $head_leng) = ($head[1], $head[3]);
        ($ref_id, $ref_leng) = (substr($head[0], 1), $head[2]);
        next;
    }
    if ($head_id){
        my @slice = split(' ', $line);
        if (@slice == 1){
            next;
        }
        my ($rstart, $rend, $qstart, $qend, $err, $t2, $t3) = @slice;
        my ($start, $end);
        if ($qend > $qstart){    
            ($start, $end) = ($qstart, $qend);
        } else {
            ($start, $end) = ($qend, $qstart);
        }
        my $dist = ($rend - $rstart) + 1;
        my $ident = 100*(($dist - $err)/$dist);
        print STDERR "id: $head_id, length: $head_leng, start: $start, end: $end, ident: $ident\n" if ($opt_v);
        if ($ident < $minalign){
            print STDERR "\tAlignment identity below $minalign%, skipping\n" if ($opt_v);
            next;
        }
        push @matches, ([$head_id, $head_leng, $start, $end, $ident, $ref_id, $rstart, $rend, $dist]);
	$seen{$head_id}++;
    }
}
close ($in);

#identify non-aligning regions in the query sequence
@matches = sort{$a->[0] cmp $b->[0] || $a->[2] <=> $b->[2] || $b->[8] <=> $a->[8]} @matches; #sort by query name, query start, and overlap length

my ($ostart, $oend, $last_id, $last_leng);
my @to_output;
for my $i (0 .. $#matches){
    my ($head_id, $head_leng, $start, $end, $ident, $ref_id, $rstart, $rend, $dist) = @{$matches[$i]};
    print "sorted matches: id: $head_id, length: $head_leng, start: $start, end: $end, ident: $ident\n" if ($opt_v);
    if (!$last_id or $head_id ne $last_id){ #if this is a new query contig...
        if ($oend and $oend != $last_leng){ #if the last alignment on the last contig did not extend to the end of the contig...
            my $outstart = ($oend + 1);
            my $oleng = ($last_leng - $outstart) + 1; #length of sequence from the last overlap end to the end of the contig
            push @to_output, ([$last_id, $last_leng, $outstart, $last_leng, $oleng]) if $oleng >= $minsize; #add to array if fragment to output is at least the minimum size
            print "\tto_output: $last_id, $last_leng, $outstart, $last_leng, $oleng\n" if ($opt_v and $oleng >= $minsize);
        }
        ($ostart, $oend) = ($start, $end);
        if ($ostart > 1){ #if the alignment does not start at the beginning of the contig...
            my $outend = ($ostart - 1);
            my $oleng = $outend;
            push @to_output, ([$head_id, $head_leng, 1, $outend, $oleng]) if $oleng >= $minsize; #add to array if fragment to output is at least the minimum size
            print "\tto_output: $head_id, $head_leng, 1, $outend, $oleng\n" if ($opt_v and $oleng >= $minsize);
        }
    } else {
        if ($end > $oend){
            if ($start <= $oend){
                $oend = $end;
            } else { #if the next alignment is outside the boundaries of the previous alignment(s)...
                my $outstart = $oend + 1;
                my $outend = $start - 1;
                my $oleng = ($outend - $outstart) + 1;
                push @to_output, ([$head_id, $head_leng, $outstart, $outend, $oleng]) if $oleng >= $minsize; #add to array if fragment to output is at least the minimum size
                print "\tto_output: $head_id, $head_leng, $outstart, $outend, $oleng\n" if ($opt_v and $oleng >= $minsize);
                ($ostart, $oend) = ($start, $end);
            }
        }
    }
    ($last_id, $last_leng) = ($head_id, $head_leng);
}
if ($oend != $last_leng){ #if the last alignment on the last contig did not extend to the end of the contig...
    my $outstart = ($oend + 1);
    my $oleng = ($last_leng - $outstart) + 1; #length of sequence from the last overlap end to the end of the contig
    push @to_output, ([$last_id, $last_leng, $outstart, $last_leng, $oleng]) if $oleng >= $minsize; #add to array if fragment to output is at least the minimum size
    print "\tto_output: $last_id, $last_leng, $outstart, $last_leng, $oleng\n" if ($opt_v and $oleng >= $minsize);
}

my ($id, $seq);
my %query_seqs;
open (my $inq, "<", $inquery);
while (my $line = <$inq>){
    chomp $line;
    if ($line =~ /^>/){
        if ($id){
            $query_seqs{$id} = $seq;
            my $seq_leng = length($seq);
            if (!$seen{$id}){
                push @to_output, ([$id, $seq_leng, 1, $seq_leng, $seq_leng]) if ($seq_leng >= $minsize);
                print "\tto_output: $id, $seq_leng, 1, $seq_leng, $seq_leng\n" if ($opt_v and $seq_leng >= $minsize);
            }
            ($id, $seq) = ("") x 2;
        }
        $id = substr($line, 1);
        $id =~ s/\s.*//; #remove everything after first whitespace
        next;
    }
    $line =~ s/\s//g;
    $seq .= $line; 
}
$query_seqs{$id} = $seq;
my $seq_leng = length($seq);
if (!$seen{$id}){
    push @to_output, ([$id, $seq_leng, 1, $seq_leng, $seq_leng]) if ($seq_leng >= $minsize);
    print "\tto_output: $id, $seq_leng, 1, $seq_leng, $seq_leng\n" if ($opt_v and $seq_leng >= $minsize);
}
close ($inq);

my @to_output_sorted = sort{$b->[4] <=> $a->[4]} @to_output;

#input genes, if requested (option -c)
my $cout;
my %coords;
my @orf_list;
if ($opt_c){
    open (my $cin, "<", $coord) or warn "** Couldn't open $coord: $!\n** No gene output file will be produced\n";
    if ($cin){
        my %c_hash;
        my $ctg_id;
        my $orf_count = 0;
        while (my $line = <$cin>){
            chomp $line;
            if ($line =~ m/^>/){
                if (%c_hash){
                    $coords{$ctg_id} = [%c_hash];
                    %c_hash = ();
                }
                $ctg_id = substr($line, 1);
                $ctg_id =~ s/\s.*//; #removes everything from the first space onwards
                next;
            } else {
                next if (!$ctg_id); #skips any header lines, i.e. as might be output by Glimmer3
                my ($id, $start, $stop) = split(' ', $line);
                my ($temp_id) = "$ctg_id\_$id";
                if ($start > $stop){
                    my ($tstart, $tstop) = ($start, $stop);
                    ($start, $stop) = ($tstop, $tstart);
                }
                my $dist = ($stop - $start) + 1;
                push @orf_list, ([$ctg_id, $id, $temp_id, $dist]);
                for my $i ($start .. $stop){
                    my @tmp;
                    if ($c_hash{$i}){
                        @tmp = @{$c_hash{$i}};
                    }
                    push @tmp, $temp_id;
                    $c_hash{$i} = [@tmp];
                }
                $orf_count++;
                print STDERR "\rInputting gene# $orf_count" if ($orf_count % 100 == 0);
            }
        }
        if (%c_hash){
            $coords{$ctg_id} = [%c_hash];
            %c_hash = ();
        }
        print STDERR "\rInputting gene# $orf_count\n";
    }
    close ($cin);
}

#output accessory sequences and key
my (%acc_orfs, %crds);
my ($seq_num, $tot_size) = (0) x 2;
$last_id = "";
open (my $out, "> $pref\_non-core.fasta");
open (my $key, "> $pref\_non-core.key.txt");
print $key "id\tlength\tgc%\tscaffold_parent\tscaffold_parent_leng\tl_border_coord\tr_border_coord\n";
for my $i (0 .. $#to_output_sorted){
    my ($id, $qleng, $start, $stop, $leng) = @{$to_output_sorted[$i]};
    my $seq = $query_seqs{$id};
    my $subseq = substr($seq, ($start - 1), $leng);
    next if ($subseq !~ m/[ACTGactg]/); #skip if the sequence only contains ambiguous bases (N's)
    print "outputting_sorted: $id, $qleng, $start, $stop, $leng\n" if ($opt_v);
    $seq_num++;
    my $lz_sn = sprintf("%05d", $seq_num);
    my $gc = gc_content($subseq);
    my $new_id = "$seqpref\_non.core$lz_sn\_leng.$leng";
    print $out ">$new_id\n$subseq\n";
    $tot_size += $leng;
    my ($l_bc, $r_bc) = ("?") x 2;
    if ($start > 1){
        $l_bc = $start - 1;
    } 
    if ($stop < $qleng){
        $r_bc = $stop + 1;
    }
    print $key "$new_id\t$leng\t$gc\t$id\t$qleng\t$l_bc\t$r_bc\n";
    if ($l_bc eq "?"){
        $l_bc = 1;
    } else {
        $l_bc++;
    }
    if ($r_bc eq "?"){
        $r_bc = $qleng;
    } else {
        $r_bc--;
    }
    if ($coords{$id}){
        print STDERR "\tSearching for accessory genes... \n" if ($opt_v);
        if (!$last_id or $last_id ne $id){
            %crds = @{$coords{$id}};
        }
        for my $j ($l_bc .. $r_bc){
            print STDERR "\r\t$j" if ($opt_v);
            if ($crds{$j}){
                my @tmp = @{$crds{$j}};
                for my $k (0 .. $#tmp){
                    $acc_orfs{$tmp[$k]}++;
                }
            }
        }
        print STDERR "\n" if ($opt_v);
    }
    $last_id = $id;
}

print STDERR "Output $seq_num non-core sequences with total length of $tot_size\n";

# Determine how much of the reference sequence aligned to the query (mostly for the sake of curiosity?)
open (my $rin, "<", $inref) or warn "Can't open reference sequence file $inref:$!\nNo reference coverage statistics can be output\n";
if ($rin){
    @matches = sort{$a->[5] cmp $b->[5] || $a->[6] <=> $b->[6] || $b->[8] <=> $a->[8]} @matches; #sort by reference name, reference start, and overlap length
    
    my ($ostart, $oend, $last_id);
    my $ref_cov = 0;
    for my $i (0 .. $#matches){
        my ($head_id, $head_leng, $start, $end, $ident, $ref_id, $rstart, $rend, $dist) = @{$matches[$i]};
        print "sorted reference matches: id: $ref_id, start: $rstart, end: $rend\n" if ($opt_v);
        if (!$last_id or $ref_id ne $last_id){ #if this is a new reference contig...
            if ($oend){
                my $odist = ($oend - $ostart) + 1;
                $ref_cov += $odist;
                print "\t$last_id: $ostart - $oend ($odist) total: $ref_cov\n" if ($opt_v);
            }
            ($ostart, $oend) = ($rstart, $rend);
        } else {
            if ($rend > $oend){
                if ($rstart <= $oend){
                    $oend = $rend;
                } else { #if the next alignment is outside the boundaries of the previous alignment(s)...
                    my $odist = ($oend - $ostart) + 1;
                    $ref_cov += $odist;
                    print "\t$ref_id: $ostart - $oend ($odist) total: $ref_cov\n" if ($opt_v);
                    ($ostart, $oend) = ($rstart, $rend);
                }
            }
        }
        $last_id = $ref_id;
    }
    my $odist = ($oend - $ostart) + 1;
    $ref_cov += $odist;
    print "\t$last_id: $ostart - $oend ($odist) total: $ref_cov\n" if ($opt_v);
    
    my $ref_size = 0;
    while (my $line = <$rin>){
        chomp $line;
        next if ($line =~ m/^>/);
        $line =~ s/\s*$//;
        my $line_leng = length($line);
        $ref_size += $line_leng;
    }
    print "Total reference size = $ref_size\n" if ($opt_v);
    close ($rin);
    my $pct_ref_cov = sprintf("%.2f", (100*($ref_cov / $ref_size)));
    
    print STDERR "$ref_cov / $ref_size ($pct_ref_cov\%) of the reference sequence was found in the query\n";
}

close ($out);
close ($key);

if (%coords){
    %coords = ();
    my (@accy, @core);
    my ($core_count, $acc_count) = (0) x 2;
    for my $i (0 .. $#orf_list){
        my ($ctg_id, $id, $temp_id, $dist) = @{$orf_list[$i]};
        if ($acc_orfs{$temp_id}){
            my $cov = $acc_orfs{$temp_id};
            my $pct_cov = sprintf("%.2f", 100*($cov/$dist));
            if ($pct_cov >= $minacc){
                $acc_count++;
                push @accy, ([$ctg_id, $id, $pct_cov]);
            } else {
                $core_count++;
                push @core, ([$ctg_id, $id, $pct_cov]);
                #print "core: ";
            }
            #print "$ctg_id $id $pct_cov\n";
            next;
        }
        $core_count++;
        push @core, ([$ctg_id, $id, 0.00]);
        #print "core: $ctg_id $id 0.00\n";
    }
    open (my $out, "> $pref\_orfs.txt");
    print $out "## Accessory ORFS ##\n";
    print $out "## Number of accessory ORFs: $acc_count ##\n";
    for my $i (0 .. $#accy){
        my ($ctg_id, $id, $pct_cov) = @{$accy[$i]};
        print $out "$ctg_id\t$id\t$pct_cov\n";
    }
    print $out "## Core ORFS ##\n";
    print $out "## Number of core ORFs: $core_count ##\n";
    for my $i (0 .. $#core){
        my ($ctg_id, $id, $pct_cov) = @{$core[$i]};
        print $out "$ctg_id\t$id\t$pct_cov\n";
    }
    close ($out);
    print "Number of accessory ORFs: $acc_count\n";
    print "Number of core ORFs: $core_count\n";
}

#----------------------------------------------
sub gc_content {
    my $seq = $_[0];
    my $count = 0;
    while ($seq =~/G|C|g|c/g) {
        $count++;
    }
    my $n_count = 0;
    while ($seq =~/[^ACTGactg]/g) {
        $n_count++;
    }
    my $length = length($seq);
    my $leng_no_ns = $length - $n_count;
    $leng_no_ns = 1 if ($leng_no_ns == 0);
    my $gc_long = 100*($count/$leng_no_ns);
    my $gc = sprintf("%.2f", $gc_long);
    return ($gc);
}
