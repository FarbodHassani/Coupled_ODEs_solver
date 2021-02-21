
figure(1);clf;
load('new01290.dat','-ascii');
thiss=new01290;
ss=surf(thiss);
ss.EdgeColor='none';
title('new0');
xlabel('psi2');
ylabel('cs2');
xticks([0 6 12 18 24 31 ]);
yticks([0 6 12 18 24 31 ]);
xticklabels({[  5.00e-08 5.74e-06 6.60e-04 7.58e-02 8.71e+00 1.00e+03 ]});
yticklabels({[  5.00e-15 3.62e-12 2.63e-09 1.90e-06 1.38e-03 1.00e+00 ]});
colorbar;
savefig('new1.fig');
	