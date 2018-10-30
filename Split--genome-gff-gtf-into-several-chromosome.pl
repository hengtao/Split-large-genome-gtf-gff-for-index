use strict;
use warnings;
use lib "/libPath";
use File::Path;
use Data::Dumper;
use Getopt::Long;
use File::Path;
use PerlIO::gzip;
use File::Basename;
#use Bio::SeqIO;
#use Bio::Seq;

my ($gtf, $outdir, $dict, $help);
GetOptions(
    'gtf=s'         => \$gtf,
    'outdir=s'      => \$outdir,
    'dict=s'        => \$dict,
    'help|h!'       => \$help
);
my $genomefa = $ARGV[0];

if ($help)
{
    usage();
    exit(1);
}

#my $in = Bio::SeqIO->new(-file => "$genomefa", -format = "Fasta");
#my $annot = Bio::SeqIO->new(
my $threshold = 512000000;
my %part = (
    "i1" => "A",
    "i2" => "B",
    "i3" => "C",
    "i4" => "D",
    "i5" => "E"
);
my %dict = ();
my %chromdict = ();
open DICT, "$dict" or die "cannot open $dict:$!\n";
open GTF, "$gtf" or die "cannot open $gtf:$!\n";
open GENOME, "$genomefa" or die "cannot open $genomefa:$!\n";
while(<DICT>){
    next if($_!~/SN\:/);
    my @array = split(/\t/, $_);
    my $chrom = $array[1] =~ s/SN://r;
    my $length = $array[2] =~ s/LN://r;
    $length = int($length);
    $chromdict{$chrom} = $length;

    my $nsplit = int($length/$threshold) + 1;

    $dict{$chrom} = {};
    if($nsplit == 1){
        my $newchrom = $chrom.$part{"i1"};
        $dict{$chrom}{$newchrom} = $length;
    }else{
        my ($k,$newchrom);
        my $perlength = int($length/$nsplit);
        for(my $i = 1; $i<=$nsplit-1; $i++){
            $k = "i"."$i";
            $newchrom = $chrom.$part{$k};
            my $binlength = $perlength*$i;
            $dict{$chrom}{$newchrom} = $binlength;
        }
        $k = "i"."$nsplit";
        $newchrom = $chrom.$part{$k};
        $dict{$chrom}{$newchrom} = $length;
    }
}
close DICT;


foreach my $chrom (keys %dict){
    foreach my $key (keys $dict{$chrom}){
        print $key."\t".$chrom."\t".$dict{$chrom}{$key}."\n";
    }
}


my @suffixlist=qw(.gtf .gff .gff3);
my ($name, $path, $suffix)=fileparse($gtf, @suffixlist);
my $newgtf = $outdir."/$name.new$suffix";
open OUTGTF, ">$newgtf" or die "cannot write into $newgtf:$!\n";
while(<GTF>){
    chomp;
    if(/^#|\tchromosome\t/){
        print OUTGTF $_."\n";
        next;
    }
    my ($chrom, $gene, $start, $end) = (split(/\t/, $_))[0,2,3,4];
    my $line = $_;
    $start = int($start);
    $end = int($end);
    foreach my $key (sort keys $dict{$chrom}){
        if($end < $dict{$chrom}{$key}){
            $line =~s/$chrom/$key/;
            print OUTGTF $line."\n";
            last;
        }elsif($end > $dict{$chrom}{$key} && $start < $dict{$chrom}{$key}){
            $line =~s/$chrom/$key/;
            print OUTGTF $line."\n";
            $dict{$chrom}{$key} = $end;
            last;
        }
        #}elsif($start >$dict{$chrom}{$key}){
            #$line =~ s/$chrom/$key/;
            #print OUTGTF $line."\n";
            #last;
        #}
    }
}
close GTF;
close OUTGTF;

print "------\n";
foreach my $chrom (keys %dict){
    foreach my $key (keys $dict{$chrom}){
        print $key."\t".$chrom."\t".$dict{$chrom}{$key}."\n";
    }
}

#foreach my $key (keys %chromdict){
    #my @chrompart = keys $dict{$key};
    #my $nums = @chrompart;
    #my $tmplength = 0;
    #if($nums > 1){
        #for(my $j = 0;$j<$nums-1;$j++){
            #$tmplength += $dict{$key}{$chrompart[$j]};
        #}
    #}
    #my $lastlength = $chromdict{$key} - $tmplength;
    #$dict{$key}{$chrompart[-1]} = $lastlength;
#}

@suffixlist=qw(.fa .fasta .fna);
($name, $path, $suffix)=fileparse($genomefa, @suffixlist);
my $newfa = $outdir."/$name.new$suffix";
open OUTGENOME, ">$newfa" or die "cannot write into $newfa:$!\n";
my $dnaseq = '';
my @chromlist = ();
while(<GENOME>){
    chomp;
    if(/^>(.*?)\s/){
        my $chrom = $1;
        push(@chromlist, $chrom);
        my $num = @chromlist;
        if($num > 1) {
            my $tmppos = 0;
            my $tmplen = 0;
            my $seq = "";
            for my $key(sort keys $dict{$chromlist[-2]}){
                my $length = $dict{$chromlist[-2]}{$key} - $tmplen;
                print "$key\t$dict{$chromlist[-2]}{$key}\t$length\n";
                $seq = substr($dnaseq, $tmppos, $length);
                print OUTGENOME ">$key\n";
                print OUTGENOME "$seq\n";
                $tmppos = $dict{$chromlist[-2]}{$key};
                $tmplen = $dict{$chromlist[-2]}{$key};
            }
        }
        $dnaseq = "";
    }else{
        $dnaseq .= $_;
    }
}
close GENOME;
my $tmppos = 0;
my $tmplen = 0;
my $seq = "";
for my $key(sort keys $dict{$chromlist[-1]}){
    my $length = $dict{$chromlist[-1]}{$key} - $tmplen;
    $seq = substr($dnaseq, $tmppos, $length);
    print OUTGENOME ">$key\n";
    print OUTGENOME "$seq\n";
    $tmppos = $dict{$chromlist[-1]}{$key};
    $tmplen = $dict{$chromlist[-1]}{$key};
}
close OUTGENOME;

sub usage{
    print <<"EOF";
This script was used to GO enrichment for DE genes.\n
Usage:\n
perl  Split-genome-gff-gtf-into-several-chromosome.pl <args1> <args2> <...> genome.fa

Args:
    --gtf           <genome.gtf(gff3)>
    --dict          <genome.dict>
    --outdir        <outdir>
EOF
}


