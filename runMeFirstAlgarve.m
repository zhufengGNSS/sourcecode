close all
resalgarve = fitVirusCV19Algarve(@getDataAlgarve,'prn','on','nsp',3);
nalgarve = 18; % from day onward
outalgarve = analyseCV19Algarve(@getDataAlgarve,nalgarve);