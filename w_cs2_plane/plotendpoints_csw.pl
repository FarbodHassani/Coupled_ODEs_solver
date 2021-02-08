#!/usr/bin/perl
use strict;


## CAREFUL use either ifort_larger or ifort_norpar
 ## usage: plotendpoints.pl >! file.txt;
  ## full_prepare_for_matlab0203ampl1minuscs2.pl
  ## or: full_prepare_for_matlab0203ampl.pl*
  ##depending whether you are interested in 1-cs^2 or cs^2
{

    my $psi2= 1;
    my $nc=200;
    my $nw=200;
    my $ncmin=5e-15;
    my $ncmax=1;


    my $facc=exp(-log($ncmin/$ncmax)/$nc);




    
    print "# zfinal=0\n";
    print "# cs_2 near 0\n";
    print "# Odarke not 0\n";
    print "# psi2 exponentially from 10\n"; 
    for my $cc (0..$nc){
#			my $cs2=$ncmin*$facc**$cc*0.999999999999e0;
	my $cs2=$cc/$nc;
#		my $cs2=1e0-$ncmin*$facc**$cc*0.999999999999e0;
#	print "$cs2\n";next;

 

		for my $ww (0..$nw){
		    my $w=-1.5e0+$ww/$nw*1e0;
		    
#		    print "$psi2 $cs2 $w\n";
		    my $text=`endpoint   $psi2 $cs2 $w`;
		    my($cs22,$psi2,$x,$idid)=split(' ',$text);
		    print "$cs22 $w $x $idid\n";
						 

						 
		    my %data=();
		my $maxx=0;
		my @line=split('\n',$text);
	    for my $line (@line){
		next if $line=~ /=/;
		my($x)=split(' ',$line);

			if($x>$maxx){$maxx=$x;}
	    }

#	    print "$cs2 $psi2 $maxx\n";
	}			# end psi2
    }  # end cs2

}

