#!/usr/local/bin/perl
use strict;
{

    # use output from plotmany0203loop.pl
    my $head="new0129";
    my $ticks=5;
    open(my $SCRIPT,">$head.m") || die;
     my %results=();
    my %cslabels=();
    my %psi2labels=();
    while(my $line=<STDIN>){
	next if $line=~ /#/;
	chomp($line);

	my($cs2,$psi2,$taustar)=split(' ',$line);
	$results{$cs2}{$psi2}=$taustar;
	$psi2labels{$psi2}=$psi2;
	$cslabels{$cs2}=$cs2;
    }



    my $ncs=(keys %cslabels);
    my $npsi2=(keys %psi2labels);
print "$ncs $npsi2\n";
        for my $i (0..0){
	print_matlab_file($i);



	open(my $OUT,">$head$i.dat") || die;
	foreach my $cs (sort {$a<=>$b} keys %results){
	    my $count=0;
	    foreach my $w (sort {$a<=>$b} keys %{$results{$cs}}){
		print $OUT "$results{$cs}{$w} ";
		$count++;

	    }
	    print $OUT "\n";
	}
	close($OUT);
    }

    sub print_matlab_file{
	my ($num)=@_;
	my $num1=$num+1;
	print $SCRIPT  "
figure($num1);clf;
load('$head$num.dat','-ascii');
thiss=$head$num;
ss=surf(thiss);
ss.EdgeColor='none';
title('new$num');
xlabel('psi2');
ylabel('cs2');
";
	    print $SCRIPT  "xticks([";
	    for my $i (0..$ticks){
		my $nx=int($i*$ncs/$ticks);
		print $SCRIPT  "$nx ";
	    }
	    print $SCRIPT  "]);\n";

	    print $SCRIPT  "yticks([";
	    for my $i (0..$ticks){
		my $ny=int($i*$npsi2/$ticks);
		print $SCRIPT  "$ny ";
	    }
	    print $SCRIPT  "]);\n";

	    print $SCRIPT  "xticklabels({".string_from_labels($ticks,%psi2labels);

	    print $SCRIPT  "});\n";

	print $SCRIPT  "yticklabels({".string_from_labels($ticks,%cslabels);
	    print $SCRIPT  "});\n";
	print $SCRIPT  "colorbar;
savefig('new$num1.fig');
	";

	}




	sub string_from_labels{
	    my($nn,%l)=@_;
	    $nn;
	    my $data=(keys %l);
	    my $n=int($data/$nn);
	    my $count=0;
	    print "$n $data $nn\n";
	    my $str="";
	    foreach my $l (sort {$a<=>$b} keys %l){
#		print "$count $l\n";
		if(int($count/$n)*$n==$count){
		    my $ll=sprintf("%.2e",$l);
		    $str.=" ".$ll;
		}
		$count++;
	    }
	    return "[ $str ]";
	}

    }
