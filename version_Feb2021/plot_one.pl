#!/usr/bin/perl

use strict;
{
    my $psi2;
	my $cs2;
    my $w;
    my $zfinal;
    my $Odarke;
    #    my $text=`ifort_norpar $psi2 $cs2 $w`;
    #        my $text=`a.out $psi2 $cs2 $w`;
    my $ymax=10;
    my $maxx=0;
    my %data=();
    my $maxx=0;
#    my @line=split('\n',$text);
    while( my $line=<STDIN>){
#	print "$line\n";
	if ($line=~ /=/){
	    chomp $line;
	    my($txt,$data)=split('=',$line);
	    $data=sprintf("%.4e",$data);
	    if($txt=~ /psi2/){$psi2=$data;}
	    if($txt=~ /cs2/){$cs2=$data;}
	    if($txt=~ /w/){$w=$data;}
	    if($txt=~ /zfinal/){$zfinal=$data;}
	    if($txt=~ /Odarke/){$Odarke=$data;}

	    next;
	}
	my($x,@y)=split(' ',$line);
	@{$data{$x}}=@y;
	if($x>$maxx){$maxx=$x;}

    }


    print_header();
    print "\@ with g0\n";
    print "\@ title \"psi2=$psi2  cs2=$cs2 w=$w\"\n";
    print "\@ subtitle \"zfinal=$zfinal  Odarke=$Odarke\"\n";
    print "\@ world 0, -$ymax,$maxx+1,$ymax\n";
    print "\@ xaxis tick major 5\n";
    print "\@ xaxis label \"tau\"\n";
    print "\@ yaxis tick major ",$ymax/2,"\n";
    print "\@    legend 1.2, 0.8\n";
    my @c=();
    @{$c[0]}=(5,"psi2");
    @{$c[1]}=(10,"psip2");
    @{$c[2]}=(13,"pi0");
    @{$c[3]}=(14,"pi1");
    @{$c[4]}=(15,"pi2");
    @{$c[5]}=(16,"pi3");
    @{$c[6]}=(17,"pi4");
    @{$c[7]}=(18,"pip0");
    @{$c[8]}=(19,"pip1");
    @{$c[9]}=(20,"pip2");
    @{$c[10]}=(21,"pip3");
    @{$c[11]}=(22,"pip4");
    @{$c[12]}=(23,"taustar");
    @{$c[13]}=(24,"hprime");
      @{$c[14]}=(25,"l2norm");
#       for my $i (4,9){
     for my $i (2,4,7,9){
#	    for my $i (0..14){
	print "\@target G0.S$i\n";
	print "\@ s$i linewidth 2.5\n";
	print "\@  s$i legend \"$c[$i][1]\"\n";
	foreach my $z (sort {$a<=>$b} keys %data){
	    print "$z ",$data{$z}[$c[$i][0]-2]+0,"\n";
	}

    }

    sub print_header(){
	print
'# Grace project file
#
@version 50122
@page size 837, 594
@page scroll 5%
@page inout 5%
@link page off
@map font 0 to "Times-Roman", "Times-Roman"
@map font 1 to "Times-Italic", "Times-Italic"
@map font 2 to "Times-Bold", "Times-Bold"
@map font 3 to "Times-BoldItalic", "Times-BoldItalic"
@map font 4 to "Helvetica", "Helvetica"
@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"
@map font 6 to "Helvetica-Bold", "Helvetica-Bold"
@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"
@map font 8 to "Courier", "Courier"
@map font 9 to "Courier-Oblique", "Courier-Oblique"
@map font 10 to "Courier-Bold", "Courier-Bold"
@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"
@map font 12 to "Symbol", "Symbol"
@map font 13 to "ZapfDingbats", "ZapfDingbats"
@map color 0 to (255, 255, 255), "white"
@map color 1 to (0, 0, 0), "black"
@map color 2 to (255, 0, 0), "red"
@map color 3 to (0, 255, 0), "green"
@map color 4 to (0, 0, 255), "blue"
@map color 5 to (255, 255, 0), "yellow"
@map color 6 to (188, 143, 143), "brown"
@map color 7 to (220, 220, 220), "grey"
@map color 8 to (148, 0, 211), "violet"
@map color 9 to (0, 255, 255), "cyan"
@map color 10 to (255, 0, 255), "magenta"
@map color 11 to (255, 165, 0), "orange"
@map color 12 to (114, 33, 188), "indigo"
@map color 13 to (103, 7, 72), "maroon"
@map color 14 to (64, 224, 208), "turquoise"
@map color 15 to (0, 139, 0), "green4"
@reference date 0
@date wrap off
@date wrap year 1950
@default linewidth 1.0
@default linestyle 1
@default color 1
@default pattern 1
@default font 0
@default char size 1.000000
@default symbol size 1.000000
@default sformat "%.8g"
@background color 0
@page background fill on
@timestamp off
@timestamp 0.03, 0.03
@timestamp color 1
@timestamp rot 0
@timestamp font 0
@timestamp char size 1.000000
@timestamp def "Tue Feb  2 13:47:02 2021"
@r0 off
@link r0 to g0
@r0 type above
@r0 linestyle 1
@r0 linewidth 1.0
@r0 color 1
@r0 line 0, 0, 0, 0
@r1 off
@link r1 to g0
@r1 type above
@r1 linestyle 1
@r1 linewidth 1.0
@r1 color 1
@r1 line 0, 0, 0, 0
@r2 off
@link r2 to g0
@r2 type above
@r2 linestyle 1
@r2 linewidth 1.0
@r2 color 1
@r2 line 0, 0, 0, 0
@r3 off
@link r3 to g0
@r3 type above
@r3 linestyle 1
@r3 linewidth 1.0
@r3 color 1
@r3 line 0, 0, 0, 0
@r4 off
@link r4 to g0
@r4 type above
@r4 linestyle 1
@r4 linewidth 1.0
@r4 color 1
@r4 line 0, 0, 0, 0
@g0 on
@g0 hidden false
@g0 type XY
@g0 stacked false
@g0 bar hgap 0.000000
@g0 fixedpoint off
@g0 fixedpoint type 0
@g0 fixedpoint xy 0.000000, 0.000000
@g0 fixedpoint format general general
	@g0 fixedpoint prec 6, 6
	';
    }
}
