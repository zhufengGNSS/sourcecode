close all
resacores = fitVirusCV19Acores(@getDataAcores,'prn','on','nsp',3);
nacores = 18; % from day onward
outacores = analyseCV19Acores(@getDataAcores,nacores);