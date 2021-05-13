#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
my $cm=2;	# fetal-maternal overlap multiplied by the phenotypic correlation
my $cp=2;	# fetal-paternal overlap multiplied by the phenotypic correlation
my $mp=2;	# maternal-paternal overlap multiplied by the phenotypic correlation
my $child=undef;	# location of child file
my $mother=undef;	# location of maternal file
my $father=undef;	# location of paternal file
my $outfile=undef;	# destination filename
# Read in options
GetOptions(
  "mp=f"  => \$mp,
  "cm=f"  => \$cm,
  "cp=f"  => \$cp,
  "child_file=s" => \$child,
  "mother_file=s" => \$mother,
  "father_file=s" => \$father,
  "out_file=s" => \$outfile,
);

sub print_usage{	# subroutine to print usage information
  print "\n\tUsage:\n";
  print "\t--child_file\t<path_to_child_file>\n";
  print "\t--mother_file\t<path_to_mother_file>\n";
  print "\t--father_file\t<path_to_father_file> (only required if allowing non-zero paternal effects)\n";
  print "\t--out_file\toutput filename\n";
  print "\t--cm\toverlap of child/maternal\n";
  print "\t--cp\toverlap of child/paternal (only required if allowing non-zero paternal effects)\n";
  print "\t--mp\toverlap of maternal/paternal (only required if allowing non-zero paternal effects)\n\n";
}
# check that minimal options are supplied
if($child && $mother && $outfile && $cm<=1 && $cm>=0){
  # check if any of the father options are supplied, if so unless all are supplied then error
  if(!(!$father && $cp==2 && $mp==2)){
    # we have at least one paternal effect, so check that they're all supplied
    if(!($father && $cp<=1 && $cp>=0 && $mp<=1 && $mp>=0)){
      print "ERROR: at least one paternal argument provided. Please either supply only fetal-maternal arguments, or supply all fetal, maternal and paternal options\n";
      print_usage();
      exit(2);
    }
  }

  # read in the child file and store in a hash
  open(my $in,"zcat ".$child." | ") or die $!;
  my %fetal=();
  <$in>;
  while(<$in>){
    chomp;
    my @F=split(' ');
    my @a=sort(uc($F[1]), uc($F[2]));
    $fetal{$F[0].":".$a[0].":".$a[1]}=uc($F[1])."\t".uc($F[2])."\t".$F[3]."\t".$F[4]."\t".$F[5];
  }
  close($in);
  # read in the mother file and store in a hash
  open($in,"zcat ".$mother." | ") or die $!;
  my %maternal=();
  <$in>;
  while(<$in>){
    chomp;
    my @F=split(' ');
    my @a=sort(uc($F[1]), uc($F[2]));
    $maternal{$F[0].":".$a[0].":".$a[1]}=uc($F[1])."\t".uc($F[2])."\t".$F[3]."\t".$F[4]."\t".$F[5];
  }
  close($in);
  # if the father options are provided the read in and store the father estimates
  my %paternal=();
  if($father){
    open($in,"zcat ".$father." | ") or die $!;
    <$in>;
    while(<$in>){
      chomp;
      my @F=split(' ');
      my @a=sort(uc($F[1]), uc($F[2]));
      $paternal{$F[0].":".$a[0].":".$a[1]}=uc($F[1])."\t".uc($F[2])."\t".$F[3]."\t".$F[4]."\t".$F[5];
    }
  }
  close($in);
  # now open the output file
  open(my $out," | gzip -c > ".$outfile);
  # print the headers
  print $out "MarkerName\tEffect_allele\tOther_allele\t";
  print $out "beta_fetal\tse_fetal\tp_fetal\t";
  print $out "beta_maternal\tse_maternal\tp_maternal\t";
  if($father){
    print $out "beta_paternal\tse_paternal\tp_paternal\t";
  }
  print $out "wlm_beta_fetal\twlm_se_fetal\t";
  print $out "wlm_beta_maternal\twlm_se_maternal\t";
  if($father){
    print $out "wlm_beta_paternal\twlm_se_paternal";
  }
  print $out "\n";
  # now loop over the contents of the fetal file
  foreach my $key (keys %fetal){
    # check if it's in the maternal and paternal (if present) files as we can only run on SNPs in both files
    my $maternal_present=0;
    if(exists $maternal{$key}){
      $maternal_present=1;
    }
    my $paternal_present=0;
    if($father){
      if(exists $paternal{$key}){
        $paternal_present=1;
      }
    }else{
      $paternal_present=1;
    }
    if($maternal_present && $paternal_present){
      # separate the elements of the fetal, maternal (and paternal) estimates
      my @F=split(' ',$fetal{$key});
      my @M=split(' ',$maternal{$key});
      my @P=();
      if($father){
        @P=split(' ',$paternal{$key});
      }
      # print out the fetal raw estimates for this SNP
      print $out $key."\t".$fetal{$key};
      # align the maternal raw estimates for this SNP to the fetal effect allele and print it
      if($F[0] eq $M[0] && $F[1] eq $M[1]){
        print $out "\t".join("\t",$M[2],$M[3],$M[4]);
      }elsif($F[0] eq $M[1] && $F[1] eq $M[0]){
        $M[2]=-$M[2];
        print $out "\t".join("\t",$M[2],$M[3],$M[4]);
      }
      if($father){
        if($F[0] eq $P[0] && $F[1] eq $P[1]){
          print $out "\t".join("\t",$P[2],$P[3],$P[4]);
        }elsif($F[0] eq $P[1] && $F[1] eq $P[0]){
          $P[2]=-$P[2];
          print $out "\t".join("\t",$P[2],$P[3],$P[4]);
        }
      }
      if($father){
        my $beta=2*$F[2]-$M[2]-$P[2];
        my $se=sqrt(4*$F[3]**2+$M[3]**2+$P[3]**2+2*$mp*sqrt(($M[3]**2)*($P[3]**2))-4*$cm*sqrt(($M[3]**2)*($F[3]**2))-4*$cp*sqrt(($P[3]**2)*($F[3]**2)));
        print $out "\t".$beta."\t".$se;
        $beta=(3*$M[2]-2*$F[2]+$P[2])/2;
        $se=sqrt((9*$M[3]**2)/4+$F[3]**2+($P[3]**2)/4-$cp*sqrt(($F[3]**2)*($P[3]**2))-3*$cm*sqrt(($M[3]**2)*($F[3]**2))+3*$mp*sqrt(($P[3]**2)*($M[3]**2))/2);
        print $out "\t".$beta."\t".$se;
        $beta=(3*$P[2]-2*$F[2]+$M[2])/2;
        $se=sqrt((9*$P[3]**2)/4+$F[3]**2+($M[3]**2)/4-$cm*sqrt(($F[3]**2)*($M[3]**2))-3*$cp*sqrt(($P[3]**2)*($F[3]**2))+3*$mp*sqrt(($P[3]**2)*($M[3]**2))/2);
        print $out "\t".$beta."\t".$se;
      }else{
        my $beta=2*(2*$F[2]-$M[2])/3;
        my $se=sqrt(4*(4*$F[3]**2+$M[3]**2-4*$cm*sqrt(($F[3]**2)*($M[3]**2)))/9);
        print $out "\t".$beta."\t".$se;
        $beta=2*(2*$M[2]-$F[2])/3;
        $se=sqrt(4*(4*$M[3]**2+$F[3]**2-4*$cm*sqrt(($F[3]**2)*($M[3]**2)))/9);
        print $out "\t".$beta."\t".$se;
      }
      print $out "\n";
    }
  }
  close($in);
}else{
  print_usage();
}
