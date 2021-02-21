#!/usr/bin/perl

use strict;
{
    my $psi2= shift || die 'psi2';
    my $cs2= shift || die 'cs2  0 is not allowed. Type 0.0';
    my $w= shift || die 'w';
    open(my $OUT,">temp.file") || die 'cannot open temp.file';
    my $text=`one_orbit_no_taustar $psi2 $cs2 $w`;
    print $OUT "$text";
    close $OUT;
    open(my $OUT,">temp.xgr") || die 'cannot open temp.xgr';
    my $text=`perl plot_one.pl <  temp.file`;
    print $OUT "$text";
    close $OUT;
    exec '/usr/local/bin/xmgrace','-param','/home/eckmann/.grace/parameters.par', '-geometry',  '1400x1100','temp.xgr';

}
