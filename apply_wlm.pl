#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
my $cm=0.5;
my $cp=0.5;
my $mp=0;
my $child=undef;
my $mother=undef;
my $father=undef;
my $outfile=undef;
GetOptions(
  "mp=f"  => \$mp,
  "cm=f"  => \$cm,
  "cp=f"  => \$cp,
  "child_file=s" => \$child,
  "mother_file=s" => \$mother,
  "father_file=s" => \$father,
  "out_file=s" => \$outfile,
);
if($child && $mother && $father && $outfile){
open(my $in,"zcat ".$child." | ") or die $!;
my %fetal=();
<$in>;
while(<$in>){
  chomp;
  my @F=split(' ');
  $fetal{$F[0]}=$F[1]."\t".$F[2]."\t".$F[3]."\t".$F[4]."\t".$F[5];
}
close($in);
open($in,"zcat ".$mother." | ") or die $!;
my %maternal=();
<$in>;
while(<$in>){
  chomp;
  my @F=split(' ');
  $maternal{$F[0]}=$F[1]."\t".$F[2]."\t".$F[3]."\t".$F[4]."\t".$F[5];
}
close($in);
open($in,"zcat ".$father." | ") or die $!;
my %paternal=();
<$in>;
while(<$in>){
  chomp;
  my @F=split(' ');
  $paternal{$F[0]}=$F[1]."\t".$F[2]."\t".$F[3]."\t".$F[4]."\t".$F[5];
}
close($in);
open(my $out," | gzip -c > ".$outfile);
print $out "MarkerName\tEffect_allele\tOther_allele\tbeta_fetal\tse_fetal\tp_fetal\tbeta_maternal\tse_maternal\tp_maternal\tbeta_paternal\tse_paternal\tp_paternal\twlm_beta_fetal\twlm_se_fetal\twlm_beta_maternal\twlm_se_maternal\twlm_beta_paternal\twlm_se_paternal";
print $out "\n";
foreach my $key (keys %fetal){
  if(exists $maternal{$key} && exists $paternal{$key}){
    my @F=split(' ',$fetal{$key});
    my @M=split(' ',$maternal{$key});
    my @P=split(' ',$paternal{$key});
    print $out $key."\t".$fetal{$key};
    if($F[0] eq $M[0] && $F[1] eq $M[1]){
      print $out "\t".join("\t",$M[2],$M[3],$M[4]);
    }elsif($F[0] eq $M[1] && $F[1] eq $M[0]){
      $M[2]=-$M[2];
      print $out "\t".join("\t",$M[2],$M[3],$M[4]);
    }else{
      print "maternal\t".$key."\n";
    }
    if($F[0] eq $P[0] && $F[1] eq $P[1]){
      print $out "\t".join("\t",$P[2],$P[3],$P[4]);
    }elsif($F[0] eq $P[1] && $F[1] eq $P[0]){
      $P[2]=-$P[2];
      print $out "\t".join("\t",$P[2],$P[3],$P[4]);
    }else{
      print "paternal\t".$key."\n";
    }
    my $beta=2*$F[2]-$M[2]-$P[2];
    my $se=sqrt(4*$F[3]**2+$M[3]**2+$P[3]**2+2*$mp*sqrt(($M[3]**2)*($P[3]**2))-4*$cm*sqrt(($M[3]**2)*($F[3]**2))-4*$cp*sqrt(($P[3]**2)*($F[3]**2)));
    print $out "\t".$beta."\t".$se;
    $beta=(3*$M[2]-2*$F[2]+$P[2])/2;
    $se=sqrt((9*$M[3]**2)/4+$F[3]**2+($P[3]**2)/4-$cp*sqrt(($F[3]**2)*($P[3]**2))-3*$cm*sqrt(($M[3]**2)*($F[3]**2))+3*$mp*sqrt(($P[3]**2)*($M[3]**2))/2);
    print $out "\t".$beta."\t".$se;
    $beta=(3*$P[2]-2*$F[2]+$M[2])/2;
    $se=sqrt((9*$P[3]**2)/4+$F[3]**2+($M[3]**2)/4-$cm*sqrt(($F[3]**2)*($M[3]**2))-3*$cp*sqrt(($P[3]**2)*($F[3]**2))+3*$mp*sqrt(($P[3]**2)*($M[3]**2))/2);
    print $out "\t".$beta."\t".$se;
    print $out "\n";
  }
}
close($in);
}else{
  print "\n\tUsage:\n\t--cm\toverlap of child/maternal\n\t--cp\toverlap of child/paternal\n\t--mp\toverlap of maternal/paternal\n\t--child_file\t<path_to_child_file>\n\t--mother_file\t<path_to_mother_file>\n\t--father_file\t<path_to_father_file>\n\t--out_file\toutput filename\n\n";
}
