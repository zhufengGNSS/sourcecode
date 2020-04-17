%calculate R0 - average of contact number over 5 days
close all
clear all
res(1) = calcR0(@getDataPortugal);
res(2) = calcR0(@getDataNorte);
res(3) = calcR0(@getDataCentro);
res(4) = calcR0(@getDataSul);
res(5) = calcR0(@getDataAlentejo);
res(6) = calcR0(@getDataAlgarve);
res(7) = calcR0(@getDataMAdeira);
res(8) = calcR0(@getDataAcores);


fprintf('%12s %7s  %7s  %12s  %12s  %12s\n',...
        'Country','R0','stdR0','N','stdN','nday');
for n = 1:length(res)
    rr = res(n);
    fprintf('%12s %7.3f  %7.3f  %12d  %12d  %4d\n',...
        rr.country,rr.R0,rr.stdR0,fix(rr.N), fix(rr.stdN),fix(rr.nday));
end
