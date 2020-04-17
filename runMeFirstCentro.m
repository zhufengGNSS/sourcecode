close all
rescentro = fitVirusCV19Centro(@getDataCentro,'prn','on','nsp',3);
ncentro = 12; % from day onward
outcentro = analyseCV19Centro(@getDataCentro,ncentro);