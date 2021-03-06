#!/usr/bin/perl -w 
use strict;


sub helper
{
    print STDERR "$0 : a programm to perfom vector for learning\n";
    print STDERR "Usage:   $0 <file profile input> <windows size>";
    print STDERR "Example: $0 1aqv.aammtx 15";
}

helper if ($#ARGV != 1);

my $file=shift;
my $window=shift;
my $limit=($window-1)/2;


my @tab_file;
open(F,"$file") or die "Cannot open file profile :\"$file\" :$!\n";
while(my $line=<F>)
{
    if ($line !~/^>/)
    { 
        push @tab_file,$line;
    }
}
close F;
for (my $i=0 ; $i <= $#tab_file ; $i++)
{
    my $line=$tab_file[$i];
    chomp $line;
    my @tab_line=split(/\s+/,$line);
    #    my @tab_line=split(/,/,$line);
    my $size_vector=$#tab_line;
    #print "size vecotr = $size_vector\n";
    # NTER
    my $output="";
    if ($i < $limit)
    {
        for(my $j=($i-$limit) ; $j < $i ; $j++)
        {
            for(my $k=0; $k <= $size_vector ; $k++)
            {
                $output.=sprintf("0.0000 ");
            }
            $output.=sprintf("1.0000 ");
        }

        for(my $j=$i ; $j <= ($i+$limit) ; $j++)
        {
		if ($j <= $#tab_file)
	    	{
			chomp $tab_file[$j];
	        	$tab_file[$j]=~s/\s+$//;
        		$output.=sprintf("$tab_file[$j] 0.0000 ");
        	}
		else
		{
			# if size of protein is less than size of windows
			#for(my $j=($#tab_file+1) ; $j <= ($i+$limit) ; $j++)
			#{
            			for(my $k=0; $k <= $size_vector ; $k++)
            			{
                			$output.=sprintf("0.0000 ");
            			}
            			$output.=sprintf("1.0000 ");
				#}
		}

	}

        chop $output;
	$output=~s/,/ /g;
        print "$output\n";
    }
    elsif ($i > (($#tab_file)-$limit))
    {
	# if position is infrerio ar limit ant 
        for(my $j=($i-$limit) ; $j <= $#tab_file ; $j++)
        {
            chomp $tab_file[$j];
            $tab_file[$j]=~s/\s+$//;
            $output.=sprintf("$tab_file[$j] 0.0000 ");
        }

        for(my $j=($#tab_file+1) ; $j <= ($i+$limit) ; $j++)
        {
            for(my $k=0; $k <= $size_vector ; $k++)
            {
                $output.=sprintf("0.0000 ");
            }
            $output.=sprintf("1.0000 ");
        }
        chop $output;
	$output=~s/,/ /g;
        print "$output\n";
    }
    else
    {
        for(my $j=($i-$limit) ; $j <= ($i+$limit) ; $j++)
        {
            exit if (not $tab_file[$j]);
            chomp $tab_file[$j];
            $tab_file[$j]=~s/\s+$//;
            $output.=sprintf("$tab_file[$j] 0.0000 ");
        }
         chop $output;
	 $output=~s/,/ /g;
         print "$output\n";
    }

}
