#!/usr/bin/perl -w
use strict;
my $f=shift;



#my @tab_pb=qw(A B C D E F G H I J K L M N O P Z);
my @tab_pb=qw(A B C D E F G H I J K L M N O P);

my %hash_tp;
my %hash_tn;
my %hash_fp;
my %hash_fn;

foreach my $pb (@tab_pb)
{
	$hash_tp{$pb}=0;
	$hash_tn{$pb}=0;
	$hash_fp{$pb}=0;
	$hash_fn{$pb}=0;
}

open(F,"$f");
my @tab_f1=<F>;
close F;

my @tab_f;
foreach my $l (@tab_f1)
{
	if ($l =~/^TRUE:/ or $l =~/^PRED:/)
	{
		push @tab_f,$l;
	}
}

my $num=0;
for (my $i=0 ; $i<=$#tab_f-1; $i=$i+2)
{
	my ($id,$true)=split(/:/,$tab_f[$i]);
	my ($id2,$pred)=split(/:/,$tab_f[$i+1]);
	$true=~s/\s//g;
	$pred=~s/\s//g;

	my @tab_true=split('',$true);
        my @tab_pred=split('',$pred);	
	$num++;
	FOR2:for (my $j =0; $j <= $#tab_true ; $j++)
	{
		my $true_pb=$tab_true[$j];
	        my $pred_pb=$tab_pred[$j];
		#if ($true_pb eq 'Z') {print "$true_pb\n";}
		next FOR2 if ($true_pb eq 'Z');
		if ($true_pb eq $pred_pb)
		{
			$hash_tp{$true_pb}++;
			FOREACH:foreach my $pb (@tab_pb)
			{
				next FOREACH if ($pb eq $true_pb);
				$hash_tn{$pb}++;
			}	
		}
		else
		{
			$hash_fn{$true_pb}++;
			$hash_fp{$pred_pb}++;
			FOREACH:foreach my $pb (@tab_pb)
			{
				next FOREACH if ($pb eq $true_pb or $pb eq $pred_pb );
				$hash_tn{$pb}++;
			}	

		}


	}
}
print("N:$num\n");
printf("%2s %6s %6s %6s %6s %6s   %5s  %5s  %5s  %5s %5s \n"," ","TP","TN","FP","FN","TOT","TPR","TNR","PPV","BER","F1");

my @tab_value_weighted;
my @tab_value;

for (my $i=0 ; $i <= 5 ; $i++)
{
	$tab_value_weighted[$i]=0;
	$tab_value[$i]=0;

}


my $tot_all=0;
foreach my $pb (@tab_pb)
{
		my $tpr=$hash_tp{$pb}/($hash_tp{$pb}+$hash_fn{$pb});
		my $tnr=$hash_tn{$pb}/($hash_tn{$pb}+$hash_fp{$pb});
		my $ppv=$hash_tp{$pb}/($hash_tp{$pb}+$hash_fp{$pb});
		my $ber=($tpr+$tnr)/2;
		my $f1=2*($ppv*$tpr)/($ppv+$tpr); 
	
	
	printf("%2s %6d %6d %6d %6d %6d   %5.2f  %5.2f  %5.2f  %5.2f  %5.2f \n",
	      	$pb,	
		$hash_tp{$pb},  $hash_tn{$pb}, $hash_fp{$pb}, $hash_fn{$pb}, $hash_tp{$pb}+$hash_fn{$pb},
		$tpr,
		$tnr,
		$ppv,
		$ber,
		$f1);

		$tab_value_weighted[0]=$tab_value_weighted[0]+$tpr;
		$tab_value_weighted[1]=$tab_value_weighted[1]+$tnr;
		$tab_value_weighted[2]=$tab_value_weighted[2]+$ppv;
		$tab_value_weighted[3]=$tab_value_weighted[3]+$ber;
		$tab_value_weighted[4]=$tab_value_weighted[4]+$f1;
		$tab_value_weighted[5]=$tab_value_weighted[5]+($hash_tp{$pb}+$hash_fn{$pb});

		my $tot=($hash_tp{$pb}+$hash_fn{$pb});

		$tab_value[5]=$tab_value[5]+($hash_tp{$pb}+$hash_fn{$pb});
		
		$tot_all=$tot_all+$tot;

		$tab_value[0]=$tab_value[0]+$tpr*$tot;
		$tab_value[1]=$tab_value[1]+$tnr*$tot;
		$tab_value[2]=$tab_value[2]+$ppv*$tot;
		$tab_value[3]=$tab_value[3]+$ber*$tot;
		$tab_value[4]=$tab_value[4]+$f1*$tot;
		#$tab_value[5]=$tab_value[5]+($hash_tp{$pb}+$hash_fn{$pb});
}


print "\n";
printf("%2s %6d %6d %6d %6d %6d   ","WM","0", "0", "0", "0",$tab_value[5]);

for (my $i=0 ; $i < 5 ; $i++)
{
	printf("%5.2f  ", $tab_value_weighted[$i]/16);

}
print"\n";
printf("%2s %6d %6d %6d %6d %6d   ","M","0", "0", "0", "0",$tab_value[5]);
for (my $i=0 ; $i < 5 ; $i++)
{
	printf("%5.2f  ", $tab_value[$i]/$tab_value[5]);

}
print "\n";

