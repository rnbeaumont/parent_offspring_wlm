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
my $poe;	# estimate POE term rather than direct paternal term in trios
# Read in options
GetOptions(
	"mp=f"	=> \$mp,
	"cm=f"	=> \$cm,
	"cp=f"	=> \$cp,
	"child_file=s" => \$child,
	"mother_file=s" => \$mother,
	"father_file=s" => \$father,
	"out_file=s" => \$outfile,
	"poe"	=> \$poe,
);

sub print_usage{	# subroutine to print usage information
	print "\n\tUsage:\n";
	print "\t--child_file\t<path_to_child_file>\n";
	print "\t--mother_file\t<path_to_mother_file>\n";
	print "\t--father_file\t<path_to_father_file> (only required if allowing non-zero paternal effects)\n";
	print "\t--out_file\toutput filename\n";
	print "\t--cm\toverlap of child/maternal\n";
	print "\t--cp\toverlap of child/paternal (only required if allowing non-zero paternal effects)\n";
	print "\t--mp\toverlap of maternal/paternal (only required if allowing non-zero paternal effects)\n";
	print "\t--poe\testimate parent of origin term rather than direct paternal term (optional, requires all paternal parameters to run)\n\n";
}

sub read_file{	# subroutine to read in each file
	my $fname=$_[0];
	my $hash=$_[1];
	open(my $in,"zcat ".$fname." | ") or die $!;
	<$in>;
	while(<$in>){
		chomp;
		my @F=split(' ');
		my @a=sort(uc($F[1]), uc($F[2]));
		$$hash{$F[0].":".$a[0].":".$a[1]}=uc($F[1])."\t".uc($F[2])."\t".$F[3]."\t".$F[4]."\t".$F[5];
	}
	close($in);
}

sub print_headers{	# subroutine to print file header
	my $father=$_[0];
	my $out=$_[1];
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
}

sub wlm_pairs{	# subroutine to calculate WLM in pars
	my $fbeta=$_[0];
	my $mbeta=$_[1];
	my $fse=$_[2];
	my $mse=$_[3];
	my $cm=$_[4];
	my $out=$_[5];
	my $beta=2*(2*$fbeta-$mbeta)/3;
	my $se=sqrt(4*(4*$fse**2+$mse**2-4*$cm*sqrt(($fse**2)*($mse**2)))/9);
	print $out "\t".$beta."\t".$se;
	$beta=2*(2*$mbeta-$fbeta)/3;
	$se=sqrt(4*(4*$mse**2+$fse**2-4*$cm*sqrt(($fse**2)*($mse**2)))/9);
	print $out "\t".$beta."\t".$se;
}

sub wlm_trio{	# subrountine to calculate WLM in trios
	my $fbeta=$_[0];
	my $mbeta=$_[1];
	my $pbeta=$_[2];
	my $fse=$_[3];
	my $mse=$_[4];
	my $pse=$_[5];
	my $cm=$_[6];
	my $cp=$_[7];
	my $mp=$_[8];
	my $out=$_[9];
	my $beta=2*$fbeta-$mbeta-$pbeta;
	my $se=sqrt(4*$fse**2+$mse**2+$pse**2+2*$mp*sqrt(($mse**2)*($pse**2))-4*$cm*sqrt(($mse**2)*($fse**2))-4*$cp*sqrt(($pse**2)*($fse**2)));
	print $out "\t".$beta."\t".$se;
	$beta=(3*$mbeta-2*$fbeta+$pbeta)/2;
	$se=sqrt((9*$mse**2)/4+$fse**2+($pse**2)/4-$cp*sqrt(($fse**2)*($pse**2))-3*$cm*sqrt(($mse**2)*($fse**2))+3*$mp*sqrt(($pse**2)*($mse**2))/2);
	print $out "\t".$beta."\t".$se;
	$beta=(3*$pbeta-2*$fbeta+$mbeta)/2;
	$se=sqrt((9*$pse**2)/4+$fse**2+($mse**2)/4-$cm*sqrt(($fse**2)*($mse**2))-3*$cp*sqrt(($pse**2)*($fse**2))+3*$mp*sqrt(($pse**2)*($mse**2))/2);
	print $out "\t".$beta."\t".$se;
}

sub wlm_poe{	# subrountine to calculate WLM in trios
	my $fbeta=$_[0];
	my $mbeta=$_[1];
	my $pbeta=$_[2];
	my $fse=$_[3];
	my $mse=$_[4];
	my $pse=$_[5];
	my $cm=$_[6];
	my $cp=$_[7];
	my $mp=$_[8];
	my $out=$_[9];
	my $beta=4*$fbeta-2*$mbeta-4*$pbeta;
	my $se=sqrt(16*$fse**2+4*$mse**2+16*$pse**2+16*$mp*sqrt(($mse**2)*($pse**2))-16*$cm*sqrt(($mse**2)*($fse**2))-32*$cp*sqrt(($pse**2)*($fse**2)));
	print $out "\t".$beta."\t".$se;
	$beta=2*$mbeta-2*$fbeta+2*$pbeta;
	$se=sqrt(4*$mse**2+4*$fse**2+4*$pse**2-8*$cp*sqrt(($fse**2)*($pse**2))-8*$cm*sqrt(($mse**2)*($fse**2))+8*$mp*sqrt(($pse**2)*($mse**2)));
	print $out "\t".$beta."\t".$se;
	$beta=6*$pbeta-4*$fbeta+2*$mbeta;
	$se=sqrt(36*$pse**2+16*$fse**2+4*$mse**2-16*$cm*sqrt(($fse**2)*($mse**2))-48**$cp*sqrt(($pse**2)*($fse**2))+24*$mp*sqrt(($pse**2)*($mse**2)));
	print $out "\t".$beta."\t".$se;
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
	# if POE is specified make sure all father options are too
	if($poe){
		if(!($father && $cp<=1 && $cp>=0 && $mp<=1 && $mp>=0)){
			print "ERROR: when supplying --poe argument, trios must be supplied\n";
			print_usage();
			exit(3);
		}
	}

	# read in the child file and store in a hash
	my %fetal=();
	read_file($child,\%fetal);
	# read in the mother file and store in a hash
	my %maternal=();
	read_file($mother,\%maternal);
	# if the father options are provided the read in and store the father estimates
	my %paternal=();
	if($father){
		read_file($father,\%paternal);
	}
	# now open the output file
	open(my $out," | gzip -c > ".$outfile);
	# print the headers
	print_headers($father,$out);
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
				if($poe){
					wlm_poe($F[2],$M[2],$P[2],$F[3],$M[3],$P[3],$cm,$cp,$mp,$out);
				}else{
					wlm_trio($F[2],$M[2],$P[2],$F[3],$M[3],$P[3],$cm,$cp,$mp,$out);
				}
			}else{
				wlm_pairs($F[2],$M[2],$F[3],$M[3],$cm,$out);
			}
			print $out "\n";
		}
	}
	close($out);
}else{
	print_usage();
}
