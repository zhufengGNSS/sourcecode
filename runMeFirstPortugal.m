close all
res = fitVirusCV19Portugal(@getDataPortugal,'prn','on','nsp',3);
n = 12; % from day onward
out = analyseCV19Portugal(@getDataPortugal,n);