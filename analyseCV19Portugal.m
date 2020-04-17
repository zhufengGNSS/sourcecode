function out = analyseCV19Portugal(getDataPortugal,day)
%ANALISE Plot evaluation of R0, N, Cend. If start day is not giben then it
%is taken ad the middle of data.
%
% Optional output:
%   out -- structure
%         out.R0 -  reproduction number
%         out.Rn -  reproduction number
%         out.N  -  polulation size
%         out.Cend - epidemy size
%         out.date0 - start day
%         out.nday - number of days
%
%Example: data set from 10 days till current set
%   analsys(@getDataItaly,10)
%
    narginchk(1,2)
    nargoutchk(0,1)

    % get data
    res = fitVirusCV19Portugal(getDataPortugal,'plt','off');
    date0 = res.date0;
    nday = res.day;
    if nargin < 2
        day = ceil(0.5*nday);
    end

    ndat = nday - day ;
    if ndat < 1
        if nargout > 0
            out = [];
        end
        error('Número Inválido de dias.')
    end
    R0 = NaN(ndat,1);
    N  = NaN(ndat,1);
        Rn  = NaN(ndat,1);
    Cend = NaN(ndat,1);
    Sc = NaN(ndat,1);    
    k = 0;
    date0 = date0 + day - 1;
    for n = day:nday
        res = fitVirusCV19Portugal(getDataPortugal,'day',n,'plt','off');
        k = k + 1;
        if ~isempty(res)
            R0(k) = res.R0;
            N(k) = res.N; %????????????????????????????
             Rn(k) = res.Rn; %????????????????????????????
            Cend(k) = res.Clim;
            Sc(k)  = res.Sc;
        else
            fprintf('Fail day %d\n',n)
        end
    end
    t = 0:k-1;

    % plot results
    figure
    set(gcf,'Position',[0 0 832 832])  % 642
    hold on

    % plot R0 ---------------------
    subplot(4,1,1)
    hold on

   % RR = log10(R0);
    RR = R0;
    plot(t+date0,RR,'k','LineWidth',2)
    scatter(t+date0,RR,50,'k','filled')
    scatter(t+date0,RR,30,'w','filled')

    %... limits
    xlim([t(1),t(end)]+date0);

    %... what kind of thicks?
    datetick('x',19,'keeplimits')

    %... add title
    title({sprintf('Modulação Epidemiológica do Vírus SARS-CoV2 pelo modelo SIR em %s. Estimativas Diárias.',res.country),...
        'R_0 - N.º Básico de Reprodução'})

    %... add axis labels
    ylabel('R_0 - N.º Básico de Reprodução')
    xlabel('Data')

    %... add grid
    grid on

    % plot R0 ---------------------
    subplot(4,1,2)
    hold on

   % RR = log10(R0);
    RR = Rn;
    plot(t+date0,RR,'k','LineWidth',2)
    scatter(t+date0,RR,50,'k','filled')
    scatter(t+date0,RR,30,'w','filled')

    %... limits
    xlim([t(1),t(end)]+date0);

    %... what kind of thicks?
    datetick('x',19,'keeplimits')

    %... add title
    title('Taxa Diária de R-->R_O')

    %... add axis labels
    ylabel('R')
    xlabel('Data')

    %... add grid
    grid on
    
    % plot Clim ---------------------
    subplot(4,1,3)
    hold on
    %...set scale
    if max(Cend) > 1000
        sf = 1000;
    else
        sf = 1;
    end
    plot(t+date0,Cend/sf,'r','LineWidth',2)
    scatter(t+date0,Cend/sf,50,'k','filled')
    scatter(t+date0,Cend/sf,30,'w','filled')

    %... limits
    xlim([t(1),t(end)]+date0);

    %... what kind of thicks?
    datetick('x',19,'keeplimits')

    % add title
    title({sprintf('Modulação Epidemiológica do Vírus SARS-CoV2 pelo modelo SIR em %s. Estimativas Diárias.',res.country),...
        'Número de Casos Confirmados'})

    % add axis labels
    xlabel('Data')
    if sf == 1
        ylabel('N.º de Casos')
    else
        ylabel('N.º de Casos (x1000)')
    end

    % add grid
    grid on

    % plot Clim ---------------------
    subplot(4,1,4)
    hold on
    %...set scale
    if max(N) > 1000
        sf = 1000;
    else
        sf = 1;
    end
    plot(t+date0,N/sf,'b','LineWidth',2)
    scatter(t+date0,N/sf,50,'k','filled')
    scatter(t+date0,N/sf,30,'w','filled')
    
%     plot(t+date0,Sc/sf,'g','LineWidth',2)
%     scatter(t+date0,Sc/sf,50,'k','filled')
%     scatter(t+date0,Sc/sf,30,'w','filled')

    %... limits
    xlim([t(1),t(end)]+date0);

    %... what kind of thicks?
    datetick('x',19,'keeplimits')

    % add title
    title({sprintf('Modulação Epidemiológica do Vírus SARS-CoV2 pelo modelo SIR em %s. Estimativas Diárias.',res.country),...
        'N.º Inicial de Indivíduos Susceptíveis '})

    % add axis labels
    xlabel('Date')
    if sf == 1
        ylabel('N.º de Indivíduos Susceptíveis')
    else
        ylabel('N.º Ind. Susceptíveis (x1000)')
    end

    % add grid
    grid on

    hold off
    
    if nargout > 0
        out.R0 = R0;
        out.N  = N;
        out.Cend = Cend;
        out.date0 = date0;
        out.nday = nday - day ;
    end
end