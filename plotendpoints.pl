#!/usr/local/bin/perl

use strict;


## CAREFUL use either ifort_larger or ifort_norpar
 ## usage: plotendpoints.pl >! file.txt;
  ## full_prepare_for_matlab0203ampl1minuscs2.pl
  ## or: full_prepare_for_matlab0203ampl.pl*
  ##depending whether you are interested in 1-cs^2 or cs^2
{

    my $w=-0.7E00;
    my $nc=30;
    my $np=30;
    my $ncmin=5e-15;
    my $ncmax=1;


    my $npmin=5e-8;
    my $npmax=1000e0;

    my $facc=exp(-log($ncmin/$ncmax)/$nc);
    my $facp=exp(-log($npmin/($npmax))/$np);



    for my $cc (0..$nc){
		my $cs2=$ncmin*$facc**$cc*0.999999999999e0;
#		my $cs2=1e0-$ncmin*$facc**$cc*0.999999999999e0;
#	print "$cs2\n";next;



		for my $pp (0..$np){

		    my $psi2=$npmin*$facp**$pp*0.999999999999e0;


		    my $text=`endpoint   $psi2 $cs2 $w`;
		    print "$text";
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
