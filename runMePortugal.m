close all
resportugal = fitVirusCV19Portugal(@getDataPortugal,'prn','on','nsp',3);
nportugal = 12; % from day onward
outportugal = analyseCV19Portugal(@getDataPortugal,nportugal);

