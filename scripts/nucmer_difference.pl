#!/usr/bin/perl

my $license = "
    nucmer_difference.pl
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

my $version = "0.6";

## Changes from verion 0.5.1
## Updated output to match spine v0.2 style output including core sequence and annotation files that include gene products
## Web-style formatting

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

Copyright (C) 2016  Egon A. Ozer
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
            File should be in the following format, one gene per line:
                gene_ID<tab>contig_ID<tab>start_coordinate<tab>stop_coordinate<tab>strand<tab>(optional)gene product
                
                contig IDs should match those in file given by option -q
                all coordinates are 1-based
                coordinates assuming a circular contig that cross the origin will
                    give incorrect results
    -o      prefix for output files (default: \"out\")
    -p      sequence prefix (default: same as given by option -o)
    -v      verbose output

";

## command line processing.
use Getopt::Std;
use vars qw( $opt_d $opt_q $opt_o $opt_m $opt_c $opt_a $opt_s $opt_p $opt_v $opt_w );
getopts('d:q:o:m:s:p:c:a:vw');
die $usage unless ($opt_d);

my ($infile, $inquery, $minalign, $minsize, $pref, $seqpref, $coord, $minacc);

$infile     = $opt_d if $opt_d;
$inquery    = $opt_q if $opt_q;
$minalign   = $opt_m ? $opt_m : 1 ;
$minsize    = $opt_s ? $opt_s : 10 ;
$pref       = $opt_o ? $opt_o : "out" ;
$seqpref    = $opt_p ? $opt_p : $pref ;
$coord      = $opt_c if $opt_c;
#$minacc     = $opt_a ? $opt_a : 50 ;
#    -a      minimum overlap with accessory coordinates, in percent, for a gene
#            in file given by -c to be called accessory (default: 50)

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
my @accessory;
my @core;
my @to_output;
for my $i (0 .. $#matches){
    my ($head_id, $head_leng, $start, $end, $ident, $ref_id, $rstart, $rend, $dist) = @{$matches[$i]};
    push @core, ([$start, $end, $head_id]);
    print "sorted matches: id: $head_id, length: $head_leng, start: $start, end: $end, ident: $ident\n" if ($opt_v);
    if (!$last_id or $head_id ne $last_id){ #if this is a new query contig...
        if ($oend and $oend != $last_leng){ #if the last alignment on the last contig did not extend to the end of the contig...
            my $outstart = ($oend + 1);
            push @accessory, ([$outstart, $last_leng, $last_id]);
            #my $oleng = ($last_leng - $outstart) + 1; #length of sequence from the last overlap end to the end of the contig
            #push @to_output, ([$last_id, $last_leng, $outstart, $last_leng, $oleng]) if $oleng >= $minsize; #add to array if fragment to output is at least the minimum size
            #print "\tto_output: $last_id, $last_leng, $outstart, $last_leng, $oleng\n" if ($opt_v and $oleng >= $minsize);
        }
        ($ostart, $oend) = ($start, $end);
        if ($ostart > 1){ #if the alignment does not start at the beginning of the contig...
            my $outend = ($ostart - 1);
            push @accessory, ([1, $outend, $head_id]);
            #my $oleng = $outend;
            #push @to_output, ([$head_id, $head_leng, 1, $outend, $oleng]) if $oleng >= $minsize; #add to array if fragment to output is at least the minimum size
            #print "\tto_output: $head_id, $head_leng, 1, $outend, $oleng\n" if ($opt_v and $oleng >= $minsize);
        }
    } else {
        if ($end > $oend){
            if ($start <= $oend){
                $oend = $end;
            } else { #if the next alignment is outside the boundaries of the previous alignment(s)...
                my $outstart = $oend + 1;
                my $outend = $start - 1;
                push @accessory, ([$outstart, $outend, $head_id]);
                #my $oleng = ($outend - $outstart) + 1;
                #push @to_output, ([$head_id, $head_leng, $outstart, $outend, $oleng]) if $oleng >= $minsize; #add to array if fragment to output is at least the minimum size
                #print "\tto_output: $head_id, $head_leng, $outstart, $outend, $oleng\n" if ($opt_v and $oleng >= $minsize);
                ($ostart, $oend) = ($start, $end);
            }
        }
    }
    ($last_id, $last_leng) = ($head_id, $head_leng);
}
if ($oend != $last_leng){ #if the last alignment on the last contig did not extend to the end of the contig...
    my $outstart = ($oend + 1);
    push @accessory, ([$outstart, $last_leng, $last_id]);
    #my $oleng = ($last_leng - $outstart) + 1; #length of sequence from the last overlap end to the end of the contig
    #push @to_output, ([$last_id, $last_leng, $outstart, $last_leng, $oleng]) if $oleng >= $minsize; #add to array if fragment to output is at least the minimum size
    #print "\tto_output: $last_id, $last_leng, $outstart, $last_leng, $oleng\n" if ($opt_v and $oleng >= $minsize);
}

my ($id, $seq);
my %query_seqs;
my %contig_lengs;
open (my $inq, "<", $inquery);
while (my $line = <$inq>){
    chomp $line;
    if ($line =~ /^>/){
        if ($id){
            $query_seqs{$id} = $seq;
            my $seq_leng = length($seq);
            $contig_lengs{$id} = $seq_leng;
            if (!$seen{$id}){
                push @accessory, ([1, $seq_leng, $id]);
                #push @to_output, ([$id, $seq_leng, 1, $seq_leng, $seq_leng]) if ($seq_leng >= $minsize);
                #print "\tto_output: $id, $seq_leng, 1, $seq_leng, $seq_leng\n" if ($opt_v and $seq_leng >= $minsize);
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
$contig_lengs{$id} = $seq_leng;
if (!$seen{$id}){
    push @accessory, ([1, $seq_leng, $id]);
    #push @to_output, ([$id, $seq_leng, 1, $seq_leng, $seq_leng]) if ($seq_leng >= $minsize);
    #print "\tto_output: $id, $seq_leng, 1, $seq_leng, $seq_leng\n" if ($opt_v and $seq_leng >= $minsize);
}
close ($inq);

@core = sort{$a->[2] cmp $b->[2] || $a->[0] <=> $b->[0]} @core;
@accessory = sort{$a->[2] cmp $b->[2] || $a->[0] <=> $b->[0]} @accessory;

#my @to_output_sorted = sort{$b->[4] <=> $a->[4]} @to_output;

#input genes, if requested (option -c)
my @loci_order;
my %locuslengs;
my %loci_starts;
my %loci_stops;
if ($opt_c){
    open (my $cin, "<", $coord) or warn "** Couldn't open $coord: $!\n** No gene output file will be produced\n";
    if ($cin){
        my %c_hash;
        my $ctg_id;
        my $orf_count = 0;
        while (my $line = <$cin>){
            chomp $line;
            next if $line =~ m/^\s*$/;
            my @tmp = split("\t", $line);
            my ($lid, $contig, $start, $stop) = @tmp;
            my $leng = ($stop - $start + 1);
            $locuslengs{$lid}{$contig} = $leng;
            push @loci_order, ([@tmp]);
        }
        close ($cin);
        @loci_order = sort{$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]} @loci_order; #loci should be in order by start coordinates, but will sort just in case
        for my $i (0 .. $#loci_order){
            my ($x1, $contig, $start, $stop) = @{$loci_order[$i]};
            $loci_starts{$contig}{$start} = $i unless $loci_starts{$start};
            $loci_stops{$contig}{$stop} = $i;
        }
    }
}

#Output core and accessory sequences, coordinates, and genes (if given)
my $accessory_stats = post_process("accessory", \@accessory);
my $core_stats = post_process("core", \@core);

open (my $statout, ">$pref.statistics.txt") or die "ERROR: Can't open $pref.statistics.txt: $!\n";
if ($opt_w){
    print $statout "AGEnt";
} else {
    print $statout "$0";
}
print $statout " version $version\n";
print $statout "inputs -m $minalign -s $minsize\n\n";
print $statout "source\ttotal_bp\tgc_%\tnum_segs\tmin_seg\tmax_seg\tavg_leng\tmedian_leng";
print $statout "\tnum_cds" if @loci_order;
print $statout "\n";
print $statout "$accessory_stats\n$core_stats\n";
close $statout;

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
    print STDERR "<br>\n" if $opt_w;
}

#----------------------------------------------
sub stats{
	my @lengths = @{$_[0]};
	my ($sum, $num, $min, $maxi, $mean, $median, $mode, $mode_freq);
        return (0,0,0,0,0,0,0,0) if (scalar @lengths == 0);
	my %seen;
	my @sorted_leng = sort {$a <=> $b} @lengths;
	for my $i (0 .. $#sorted_leng){
		$sum += $sorted_leng[$i];
		$seen{$sorted_leng[$i]}++;
	}
	$num = $#sorted_leng + 1;
	$min = $sorted_leng[0];
	$maxi = $sorted_leng[$#sorted_leng];
	$mean = $sum/$num;
	my $rounded_mean = sprintf("%.2f", $mean);
	my @modes;
	foreach my $leng (sort {$seen{$b} <=> $seen{$a}} keys %seen) {
		push @modes, ([$leng, $seen{$leng}]);
	}
	$mode = $modes[0][0];
	$mode_freq = $modes[0][1];
	my $mid = int @sorted_leng/2;
	if (@sorted_leng % 2){
		$median = $sorted_leng[$mid];
	} else {
		$median = ($sorted_leng[$mid-1] + $sorted_leng[$mid])/2;
	}
	return ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq);
}

sub gc_content {
    my $seq = $_[0];
    return ("0.0000") unless $seq;
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

sub post_process {
    my $type = shift;
    my @array = @_;
    
    my $fileid = "$pref.$seqpref.$type";
    my $seqid = "$seqpref\_$type\_";
    
    open (my $out_cor, ">$fileid\_coords.txt") or return "ERROR: Can't open $fileid\_coords.txt: $!\n";
    print $out_cor "contig_id\tcontig_length\tstart\tstop\tout_seq_id\n" unless ($type eq "out" or $type eq "pan");
    open (my $out_seq, ">$fileid.fasta") or return "ERROR: Can't open $fileid.fasta: $!\n";

    my %loci;
    my @lengs;
    my $count = 0;
    my $tot_seq = "";
    my $last_contig = " ";
    my $seq;
    for my $i (0 .. $#array){
        my ($o_start, $o_stop, $c_id) = @{$array[$i]};
        my $cleng = $contig_lengs{$c_id};
        my $b_leng = ($o_stop - $o_start) + 1;
        #output sequences
        my $out_id;
        if ($b_leng >= $minsize){
            if ($c_id ne $last_contig){
                $seq = $query_seqs{$c_id};
            }
            $last_contig = $c_id;
            my $b_seq = substr($seq, ($o_start - 1), $b_leng);
            #trim ambiguous bases from front and back of sequence. If the size of the remaining sequence is smaller than the minimum length, drop it.
            my ($lead, $tail) = (0) x 2;
            (my $lead_seq) = $b_seq =~ m/^([NnXx-]*)/;
            if ($lead_seq){
                $lead = length($lead_seq);
            }
            if ($lead == $b_leng){ #skip if the entire sequence was N's
                #output coordinates
                print $out_cor "$c_id\t$cleng\t$o_start\t$o_stop\t$out_id\n";
                next;
            }
            (my $tail_seq) = $b_seq =~ m/([NnXx-]*$)/;
            if ($tail_seq){
                $tail = length($tail_seq);
            }
            if ($lead > 0 or $tail > 0){
                my $dist = length($b_seq) - ($lead + $tail);
                if ($dist < $minsize){
                    #output coordinates
                    print $out_cor "$c_id\t$cleng\t$o_start\t$o_stop\t$out_id\n";
                    next;
                }
                $b_seq = substr($b_seq, $lead, $dist);
                $o_start += $lead;
                $o_stop -= $tail;
            }
            push @lengs, length($b_seq);
            $count++;
            $out_id = $seqid.sprintf("%04d", $count)."_length\_$b_leng";
            print $out_seq ">$out_id\n$b_seq\n";
            $tot_seq .= $b_seq;

            #output coordinates
            print $out_cor "$c_id\t$cleng\t$o_start\t$o_stop\t$out_id\n";

            #collect annotation information, if present
            my @locusids;
            if (@loci_order){ #first, will set a range of loci to screen so that we don't have to search the entire array each time
                my ($start_rec, $stop_rec);
                for my $i (reverse 1 .. $o_start){
                    if ($loci_stops{$c_id}{$i}){
                        $start_rec = $loci_stops{$c_id}{$i};
                        last;
                    }
                }
                if (!$start_rec){
                    for my $i (reverse 1 .. $o_start){
                        if ($loci_starts{$c_id}{$i}){
                            $start_rec = $loci_starts{$c_id}{$i};
                            last;
                        }
                    }
                }
                if (!$start_rec){
                    for my $i ($o_start .. $o_stop){
                        if ($loci_starts{$c_id}{$i}){
                            $start_rec = $loci_starts{$c_id}{$i};
                            last;
                        }
                    }
                }
                if ($start_rec){ #no point in looking the other way if there's no start_rec
                    for my $i ($o_stop .. $cleng){
                        if ($loci_starts{$c_id}{$i}){
                            $stop_rec = $loci_starts{$c_id}{$i};
                            last;
                        }
                    }
                    if (!$stop_rec){
                        for my $i (reverse $o_start .. $o_stop){
                            if ($loci_starts{$c_id}{$i}){
                                $stop_rec = $loci_starts{$c_id}{$i};
                                last;
                            }
                        }
                    }
                    $stop_rec = $start_rec if !$stop_rec;
                    @locusids = @loci_order[$start_rec .. $stop_rec];
                }
            }
            if (@locusids){
                foreach my $slice (@locusids){
                    my ($lid, $c_id, $l_start, $l_stop) = @{$slice};
                    next if $l_stop < $o_start;
                    next if $l_start > $o_stop;
                    my ($over_front, $over_back) = (0) x 2; #keep track of how much a gene overhangs the end of the segment
                    if ($l_start < $o_start){
                        $over_front = $o_start - $l_start;
                        if ($l_stop <= $o_stop){
                            push @{$loci{$lid}{$c_id}}, ([$out_id, $o_start, $o_start, $l_stop, $over_front, $over_back]);
                        } else {
                            $over_back = $l_stop - $o_stop;
                            push @{$loci{$lid}{$c_id}}, ([$out_id, $o_start, $o_start, $o_stop, $over_front, $over_back]);
                        }
                        next;
                    }
                    if ($l_stop > $o_stop){
                        $over_back = $l_stop - $o_stop;
                        push @{$loci{$lid}{$c_id}}, ([$out_id, $o_start, $l_start, $o_stop, $over_front, $over_back]);
                        next;
                    }
                    push @{$loci{$lid}{$c_id}}, ([$out_id, $o_start, $l_start, $l_stop, $over_front, $over_back]); #don't need an if statement. Everything left after the above ifs are loci with borders within o_star and o_stop, inclusive
                }
            }
        } else {
            #output coordinates
            print $out_cor "$c_id\t$cleng\t$o_start\t$o_stop\n";
        }
    }
    close ($out_cor);
    close ($out_seq);
    my $loci_count = 0;
    if (%loci){
        #sort and output the annotation information to file
        open (my $out_gen, ">$fileid\_loci.txt") or die "ERROR: Can't open $fileid\_loci.txt: $!\n";
        print $out_gen "locus_id\tgen_contig_id\tgen_contig_start\tgen_contig_stop\tstrand\tout_seq_id\tout_seq_start\tout_seq_stop\tpct_locus\toverhangs\tproduct\n" unless ($type eq "out" or $type eq "pan");
        foreach my $slice (@loci_order){
            my @slice_arr = @{$slice};
            #$lid, $contig, $start, $stop
            my ($locus, $contig) = @slice_arr;
            next unless $loci{$locus}{$contig};
            my $dir = $slice_arr[4];
            my $prod = "NoAnnotation";
            $prod = $slice_arr[5] if $slice_arr[5];
            my $locusleng = $locuslengs{$locus}{$contig};
            my @hits = @{$loci{$locus}{$contig}};
            @hits = sort{$a->[0] cmp $b->[0]}@hits;
            #calculate total percentage of hit
            my $totleng = 0;
            foreach (@hits){
                $totleng += (${$_}[3] - ${$_}[2] + 1);
            }
            my $tot_pct = 0;
            foreach my $hit (@hits){
                my ($id, $hoffset, $hstart, $hstop, $over_front, $over_back) = @{$hit};
                my $pct = sprintf("%.2f", 100 * (($hstop - $hstart + 1) / $locusleng));
                $tot_pct += $pct;
                my ($cstart, $cstop) = ($hstart - $hoffset + 1, $hstop - $hoffset + 1);
                print $out_gen "$locus\t$contig\t$hstart\t$hstop\t$dir\t$id\t$cstart\t$cstop\t$pct\t$over_front,$over_back\t$prod\n";
            }
            $loci_count++ if $tot_pct >= 50; #$minacc; #for the purposes of counting, will use a 50% cutoff (so that only if the gene is in the majority in either core of accessory it will be counted as belonging to that group)
        }
        close ($out_gen);
    }
    #get stats
    my $stat_string = "x\n";
    my ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq) = stats(\@lengs);
    my $gc = gc_content($tot_seq);
    if ($type eq "accessory"){
        if ($opt_w){
            print "<br><strong>Accessory Genome Statistics:</strong><br>\n";
        } else {
            print "\nAccessory Genome Statistics:/n";
        }
    }
    print "Number of segments >= $minsize bp: $num\n";
    print "<br>\n" if $opt_w;
    print "Total bp: $sum\n";
    print "<br>\n" if $opt_w;
    print "GC content: $gc%\n";
    print "<br>\n" if $opt_w;
    print "Shortest segment length: $min\n";
    print "<br>\n" if $opt_w;
    print "Longest segment length: $maxi\n";
    print "<br>\n" if $opt_w;
    print "Average segment length: $rounded_mean\n";
    print "<br>\n" if $opt_w;
    print "Median segment length: $median\n";
    print "<br>\n" if $opt_w;
    print "Mode segment length: $mode, Frequency: $mode_freq\n";
    print "<br><br>\n" if $opt_w;
    
    $stat_string = "$type\t$sum\t$gc\t$num\t$min\t$maxi\t$rounded_mean\t$median";
    if (%loci){
        $stat_string .= "\t$loci_count";
    }
    $stat_string .= "\n";
    return($stat_string);
}
