close all
resnorte = fitVirusCV19Norte(@getDataNorte,'prn','on','nsp',3);
nnorte = 20; % from day onward
outnorte = analyseCV19Norte(@getDataNorte,nnorte);