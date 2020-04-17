function plotGrowthRatePortugal(res)
%PLOTGROWTHRATE Plot growth rate of epidemic from data in structure res.
    
    %check data
    narginchk(1,1)
    
    if isempty(res)
        fprintf('Empty data.\n');
        return
    end
    
    % collect data
    % ... predictions
    Ce = res.Ce;
    dCe = diff(res.Ce);
    %... actual
    dC = diff(res.C);
    C = res.C;
    %... time span
    t = 0:length(Ce)- 1; %res.t;
    
    % start new figure
    figure
    hold on
    
    set(gca, 'YScale', 'log')
    
    %... plot prediction
    plot(t(1:end-1)+res.date0,(dCe./Ce(1:end-1))*100,'LineWidth',2);
    
    %... plot actual
    scatter(t(1:end-1)+res.date0,(dC./C(1:end-1))*100,50,'k','filled')
    h = scatter(t(1:end-1)+res.date0,(dC./C(1:end-1))*100,30,'w','filled');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    %... plot limits
    plot(t(1:end-1)+res.date0,(5*ones(length(Ce)-1,1)),...
        'r','LineWidth',2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
     text(t(1)+res.date0,5,'5%','VerticalAlignment','bottom');
    
    %... label axes
    datetick('x',20,'keeplimits')
    xlabel('Data')
    ylabel('\DeltaC/C (x100%)')
    
    %... add legend
   legend('Previsão','Actual','Posição','Aproximação')
    
    %... add title
    title(sprintf('Taxa de Crescimento em %s (%s)',...
        res.country,datestr(res.date0+length(t))));
    
    %... add grid
    grid on
    
    %... end plot
    hold off
end