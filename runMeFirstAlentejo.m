close all
resalentejo = fitVirusCV19Alentejo(@getDataAlentejo,'prn','on','nsp',3);
naletenjo = 18; % from day onward
outalentejo = analyseCV19Alentejo(@getDataAlentejo,naletenjo);