close all
resmadeira = fitVirusCV19Madeira(@getDataMadeira,'prn','on','nsp',3);
nmadeira = 18; % from day onward
outmadeira = analyseCV19Madeira(@getDataMadeira,nmadeira);