use strict;
use warnings;
use lib "/libPath";
use File::Path;
use Data::Dumper;
use Getopt::Long;
use File::Path;
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
        my @range = (1, $length);
        $dict{$chrom}{$newchrom} = \@range;
    }else{
        my ($k,$newchrom);
        my $perlength = int($length/$nsplit);
        my $binstartlength;
        my $binendlength;
        for(my $i = 1; $i<=$nsplit-1; $i++){
            $k = "i"."$i";
            $newchrom = $chrom.$part{$k};
            $binstartlength = $perlength*($i-1)+1;
            $binendlength = $perlength*$i;
            my @range = ($binstartlength,$binendlength);
            $dict{$chrom}{$newchrom} = \@range;
        }
        $k = "i"."$nsplit";
        $newchrom = $chrom.$part{$k};
        my @range = ($binendlength+1,$length);
        $dict{$chrom}{$newchrom} = \@range;
    }
}
close DICT;


foreach my $chrom (keys %dict){
    foreach my $key (keys $dict{$chrom}){
        print $key."\t".$chrom."\t".join("\t", @{$dict{$chrom}{$key}})."\n";
    }
}

my @suffixlist=qw(.gtf .gff .gff3);
my ($name, $path, $suffix)=fileparse($gtf, @suffixlist);
my $newgtf = $outdir."/$name.parts$suffix";
my %hash;
my %transcript;
my $newchr;
open OUTGTF, ">$newgtf" or die "cannot write into $newgtf:$!\n";
while(<GTF>){
    chomp;
    if(/^#|\tchromosome\t/){
        print OUTGTF $_."\n"; 
        next;
    }
    my ($chrom, $gene, $start, $end) = (split(/\t/, $_))[0,2,3,4];
    next if ($gene =~/biological/i);
    my $newstart;
    my $newend;
    if($gene =~/gene|lnc_RNA|ncRNA_gene|pre_miRNA|pseudogenic_transcript|RNase_MRP_RNA|rRNA|snoRNA|snRNA|SRP_RNA|tRNA/i){
        my $line = $_;
        my $geneid = (split(/\t/, $_))[8];
        if($geneid =~/ID=(?:gene:|transcript:)?(.*?);/){
            $geneid = $1;
            $transcript{$geneid} = ();
        }else{
            exit("Wrong gff/gtf format, no gene marker detected, please check your annotation file...");
        }
        #print($geneid."\n");
        $start = int($start);
        $end = int($end);
        my @tmp = (sort keys $dict{$chrom});
        my $length = scalar(@tmp);
        #print($length);
        #print($tmp[$#tmp]);
        for(my $i=0;$i<=$#tmp;$i++){
        #foreach my $key (sort keys $dict{$chrom}){
            if($end <= @{$dict{$chrom}{$tmp[$i]}}[1]){
                #print("$tmp[$i]\n");
                $newstart = $start - @{$dict{$chrom}{$tmp[$i]}}[0]  + 1;
                $newend = $end - @{$dict{$chrom}{$tmp[$i]}}[0]  + 1;
                my @range = (@{$dict{$chrom}{$tmp[$i]}}[0], $tmp[$i]);
                $hash{$geneid} = \@range;
                $line =~s/^$chrom\t/$tmp[$i]\t/;
                $line =~s/\t$start\t/\t$newstart\t/;
                $line =~s/\t$end\t/\t$newend\t/;
                print OUTGTF "$line\n";
                $newchr = $tmp[$i];
                last;
            }elsif($end > @{$dict{$chrom}{$tmp[$i]}}[1] && $start <= @{$dict{$chrom}{$tmp[$i]}}[1]){
                #print("$tmp[$i+1]\n");
                @{$dict{$chrom}{$tmp[$i]}}[1] = $start - 1;
                @{$dict{$chrom}{$tmp[$i+1]}}[0] = $start;
                $newstart = $start - @{$dict{$chrom}{$tmp[$i+1]}}[0]  + 1;
                $newend = $end - @{$dict{$chrom}{$tmp[$i+1]}}[0]  + 1;
                my @range = (@{$dict{$chrom}{$tmp[$i+1]}}[0], $tmp[$i+1]);
                $hash{$geneid} = \@range;
                $line =~s/^$chrom\t/$tmp[$i+1]\t/;
                $line =~s/\t$start\t/\t$newstart\t/;
                $line =~s/\t$end\t/\t$newend\t/;
                print OUTGTF "$line\n";
                $newchr = $tmp[$i+1];
                last;
            }else{
                next;
            }
        }
=pod
    #}elsif($gene =~ /mRNA/i){
    	my $line = $_;
	my $geneid;
        my $transcriptid = (split(/\t/, $_))[8];
        if($transcriptid =~/Parent=(?:gene:|transcript:)?(.*?);/){
            $geneid = $1;
        }else{
            exit("Wrong gff/gtf format, no gene marker detected, please check your annotation file...");
        }
        if($transcriptid =~/ID=(?:gene:|transcript:)?(.*?);/){
            $transcriptid = $1;
        }else{
            exit("Wrong gff/gtf format, no transcript marker detected, please check your annotation file...");
        }
        if(exists $transcript{$geneid}){
	    push(@{$transcript{$geneid}}, $transcriptid);   
        }else{
            exit("$transcriptid has no parent geneid,please check your annotation file..."); 
        }
        my $newstart = $start - @{$hash{$geneid}}[0] + 1;
        my $newend = $end - @{$hash{$geneid}}[0] + 1;
        $line =~s/^$chrom\t/@{$hash{$geneid}}[1]\t/;
        $line =~s/\t$start\t/\t$newstart\t/;
        $line =~s/\t$end\t/\t$newend\t/;
        print OUTGTF "$line\n";
=cut
    }else{
        my $line = $_;
        #my $cdsid = (split(/\t/, $_))[8];
        #if($cdsid=~/Parent=(.*?);/){
            #$cdsid = $1;
        #}
        $start = int($start);
        $end = int($end);
        #foreach my $geneid (keys %hash){
            #foreach my $transcript (@{$transcript{$geneid}}){
                #if($cdsid =~ /$geneid|$transcript/){
        #my $newstart = $start - @{$hash{$geneid}}[0] + 1;
        #my $newend = $end - @{$hash{$geneid}}[0] + 1;
        my $newstart = $start - @{$dict{$chrom}{$newchr}}[0] + 1;
        my $newend = $end - @{$dict{$chrom}{$newchr}}[0] + 1;
        $line =~s/^$chrom\t/$newchr\t/;
        $line =~s/\t$start\t/\t$newstart\t/;
        $line =~s/\t$end\t/\t$newend\t/;
        print OUTGTF "$line\n";
                    #last;
                
            
        
    }
}
close GTF;
close OUTGTF;

print "------\n";


foreach my $chrom (keys %dict){
    foreach my $key (keys $dict{$chrom}){
        print $key."\t".$chrom."\t".join("\t", @{$dict{$chrom}{$key}})."\n";
    }
}



my $bed = $outdir."/$name.bed";
open(BED, ">$bed") or die $!;
foreach my $chrom (sort keys %dict){
    foreach my $key (sort keys $dict{$chrom}){
        print BED $key."\t".$chrom."\t".join("\t",@{$dict{$chrom}{$key}})."\n";
    }
}
close BED;




@suffixlist=qw(.fa .fasta .fna);
($name, $path, $suffix)=fileparse($genomefa, @suffixlist);
my $newfa = $outdir."/$name.parts$suffix";
open OUTGENOME, ">$newfa" or die "cannot write into $newfa:$!\n";
my $dnaseq = '';
my @chromlist = ();
my $length;
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
                
                my $start = @{$dict{$chromlist[-2]}{$key}}[0];
                my $end = @{$dict{$chromlist[-2]}{$key}}[1];
                #print "$key\t$dict{$chromlist[-2]}{$key}\t$length\n";
                $length = $end - $start +1;
                $seq = substr($dnaseq, $start-1, $length);
                print OUTGENOME ">$key\n";
                print OUTGENOME "$seq\n";
                #$tmppos = $dict{$chromlist[-2]}{$key};
                #$tmplen = $dict{$chromlist[-2]}{$key};
            }
        }
        $dnaseq = "";
    }else{
        $dnaseq .= $_;
    }
}
close GENOME;
#my $tmppos = 0;
#my $tmplen = 0;
my $seq = "";
for my $key(sort keys $dict{$chromlist[-1]}){
    my $start = @{$dict{$chromlist[-1]}{$key}}[0];
    my $end = @{$dict{$chromlist[-1]}{$key}}[1];
    $length = $end - $start +1;
    $seq = substr($dnaseq, $start-1, $length);
    print OUTGENOME ">$key\n";
    print OUTGENOME "$seq\n";
    #$tmppos = $dict{$chromlist[-1]}{$key};
    #$tmplen = $dict{$chromlist[-1]}{$key};
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

