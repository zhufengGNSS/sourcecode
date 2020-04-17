close all
ressul = fitVirusCV19Sul(@getDataSul,'prn','on','nsp',3);
nsul = 18; % from day onward
outsul = analyseCV19Sul(@getDataSul,nsul);