function COVID19Modelingv2(country)
%Coronavirus Tracker - Country Tracker - Joshua McGee
%Created to track the simulate the spread of Coronavirus (COVID-19)
%Data is stored online and is provided via JHU CSSE from various sources including:
%"the World Health Organization (WHO), DXY.cn. Pneumonia. 2020, BNO News,
%National Health Commission of the People? Republic of China (NHC),
%China CDC (CCDC), Hong Kong Department of Health, Macau Government, Taiwan CDC, US CDC,
%Government of Canada, Australia Government Department of Health,
%European Centre for Disease Prevention and Control (ECDC) and Ministry of
%Health Singapore (MOH)"
countryall = country;
for i = 1:length(countryall)
    country = countryall(i); %specify country to model
    %Obtaining and formating data - courtesy of Toshi Takeuchi - https://www.mathworks.com/matlabcentral/profile/authors/951521
    resultadospt=readtable('confirmados.txt');
    obitospt = readtable('obitos.txt');
    writetable(resultadospt,'result.txt','WriteVariableNames',false);
    writetable(obitospt,'deathresult.txt','WriteVariableNames',false);
    opts = detectImportOptions('result.txt', "TextType","string");
    opts1 = detectImportOptions('deathresult.txt', "TextType","string");
    
    first_day = datetime(2020,1,22);
    day_add = size(resultadospt);
    last_day = first_day+days(day_add(2)-5);
    time = first_day:last_day;
    
    C = cell(1,day_add(2));
    for i = 1:day_add(2)
        if i == 1
            C(i) = {'Province_State'};
        elseif i == 2
            C(i) = {'Country_Region'};
        elseif i == 3
            C(i) = {'Lat'};
        elseif i == 4
            C(i) = {'Long'};
        else
            formatOut = 'xmm_dd_yy';
            C(i) = {sprintf('%s',datestr(datenum(time(i-4)),formatOut))};
        end
    end
    times_conf = readtable('result.txt',opts);
    times_conf1 = readtable('deathresult.txt',opts1);
    matlab.lang.makeValidName(C);
    times_conf.Properties.VariableNames = C;
    times_conf1.Properties.VariableNames = C;
    times_conf.("Country_Region")(times_conf.("Country_Region") == "China") = "Mainland China";
    times_conf.("Country_Region")(times_conf.("Country_Region") == "Czechia") = "Czech Republic";
    times_conf.("Country_Region")(times_conf.("Country_Region") == "Iran (Islamic Republic of)") = "Iran";
    times_conf.("Country_Region")(times_conf.("Country_Region") == "Republic of Korea") = "Korea, South";
    times_conf.("Country_Region")(times_conf.("Country_Region") == "Republic of Moldova") = "Moldova";
    times_conf.("Country_Region")(times_conf.("Country_Region") == "Russian Federation") = "Russia";
    times_conf.("Country_Region")(times_conf.("Country_Region") == "Taipei and environs") = "Taiwan";
    times_conf.("Country_Region")(times_conf.("Country_Region") == "Taiwan*") = "Taiwan";
    times_conf.("Country_Region")(times_conf.("Country_Region") == "United Kingdom") = "UK";
    times_conf.("Country_Region")(times_conf.("Country_Region") == "Viet Nam") = "Vietnam";
    times_conf.("Country_Region")(times_conf.("Province_State") == "St Martin") = "St Martin";
    times_conf.("Country_Region")(times_conf.("Province_State") == "Saint Barthelemy") = "Saint Barthelemy";
    times_conf1.("Country_Region")(times_conf1.("Country_Region") == "China") = "Mainland China";
    times_conf1.("Country_Region")(times_conf1.("Country_Region") == "Czechia") = "Czech Republic";
    times_conf1.("Country_Region")(times_conf1.("Country_Region") == "Iran (Islamic Republic of)") = "Iran";
    times_conf1.("Country_Region")(times_conf1.("Country_Region") == "Republic of Korea") = "Korea, South";
    times_conf1.("Country_Region")(times_conf1.("Country_Region") == "Republic of Moldova") = "Moldova";
    times_conf1.("Country_Region")(times_conf1.("Country_Region") == "Russian Federation") = "Russia";
    times_conf1.("Country_Region")(times_conf1.("Country_Region") == "Taipei and environs") = "Taiwan";
    times_conf1.("Country_Region")(times_conf1.("Country_Region") == "Taiwan*") = "Taiwan";
    times_conf1.("Country_Region")(times_conf1.("Country_Region") == "United Kingdom") = "UK";
    times_conf1.("Country_Region")(times_conf1.("Country_Region") == "Viet Nam") = "Vietnam";
    times_conf1.("Country_Region")(times_conf1.("Province_State") == "St Martin") = "St Martin";
    times_conf1.("Country_Region")(times_conf1.("Province_State") == "Saint Barthelemy") = "Saint Barthelemy";
    vars = times_conf.Properties.VariableNames;
    vars1 = times_conf1.Properties.VariableNames;
    times_conf_country = groupsummary(times_conf,"Country_Region",{'sum'},vars(3:end));
    times_conf_country1 = groupsummary(times_conf1,"Country_Region",{'sum'},vars1(3:end));
    vars = times_conf_country.Properties.VariableNames;
    vars = regexprep(vars,"^(sum_)(?=L(a|o))","remove_");
    vars = erase(vars,{'sum_'});
    times_conf_country.Properties.VariableNames = vars;
    vars1 = times_conf_country1.Properties.VariableNames;
    vars1 = regexprep(vars1,"^(sum_)(?=L(a|o))","remove_");
    vars1 = erase(vars1,{'sum_'});
    times_conf_country1.Properties.VariableNames = vars1;
    infectedtable = removevars(times_conf_country,[{'GroupCount'},vars(contains(vars,"remove_"))]);
    countrytable = infectedtable(strcmp(infectedtable.("Country_Region"),country), :);
    deathtable = removevars(times_conf_country1,[{'GroupCount'},vars1(contains(vars1,"remove_"))]);
    countrytable1 = deathtable(strcmp(deathtable.("Country_Region"),country), :);
    countrytable = countrytable(:,2:end);
    countrytable1 = countrytable1(:,2:end);
    cols1 = size(countrytable);
    cols2 = size(countrytable1);
    Countrytotaldead = zeros(1,cols2(2));
    Countrytotalinfected = zeros(1,cols1(2));
    for i = 1:cols1(2)
        Countrytotalinfected(i) = table2array(countrytable(1,i));
        Countrytotaldead(i) = table2array(countrytable1(1,i));
    end
    infected = [Countrytotalinfected];
    deathrate = Countrytotaldead(end)/Countrytotalinfected(end);
    warning('off')
    %algorithm to determine how many slow-growth phase points to include
    startidx = zeros(length(infected));
    for i = 1:length(infected)
        if i == length(infected)
            break
        end
        if max(infected) < 10000
            if infected(i+1) > 1.5*infected(i) && infected(i) > 10
                startidx(i) = i;
            end
        else
            if infected(i+1) > 1.5*infected(i) && infected(i) > 50
                startidx(i) = i;
            end
        end
    end
    startidx = find(startidx~=0, 1, 'first');
    first_day=datenum('2020/01/22'); % spread start date - China (do not change)
    first_day = first_day+startidx;
    fprintf('Virus Spread Prediction for: %s\n',country);
    sampledTime = 0:1:length(infected(startidx:end))-1;
    sampledDate = first_day + sampledTime;
    [b0] = iniGuess(infected(startidx:end));
    if isempty(b0) || b0(2) == Inf || b0(3) == Inf
        fprintf('***Warning: Fail to calculate initial quess. Use default.\n');
        b0 = [max(infected) 0.5 max(infected)]';
    end
    K0 = b0(1);
    r0 = b0(2);
    A0 = b0(3);
    fprintf('  Initial guess K = %g  r = %g  A = %g\n',K0,r0,A0);
    infected = infected(startidx:end);
    infectedidx = length(infected);
    K      = NaN(length(infected),1);
    r      = NaN(length(infected),1);
    A      = NaN(length(infected),1);
    C0     = NaN(length(infected),1);
    tpeak  = NaN(length(infected),1);
    dpeak  = NaN(length(infected),1);
    dend   = NaN(length(infected),1);
    dCpeak = NaN(length(infected),1);
    tau    = NaN(length(infected),1);
    err    = NaN(length(infected),1);
    R2     = NaN(length(infected),1);
    n0 = ceil(0.25*length(infected));
    opts = optimoptions('lsqcurvefit','Display','off',...
        'SpecifyObjectiveGradient',true);
    for n = n0:length(infected)
        [b,resnorm,~,~,~] = lsqcurvefit(@fun,b0,...
            sampledTime(1:n),infected(1:n),[0 0 0],[],opts);
        err(n) = resnorm;
        R2(n) = calcR2(sampledTime(1:n),infected(1:n),b);
        K(n)   = fix(b(1));
        r(n)   = b(2);
        A(n)   = b(3);
        C0(n)  = fix(K(n)/(A(n) + 1));
        tpeak(n)  = fix(log(A(n))/r(n));
        dpeak(n)  = tpeak(n) + first_day;
        dend(n)   = 2*tpeak(n) + first_day;
        dCpeak(n) = fix(r(n)*K(n)/4);
        tau(n)    = fix(4/r(n));
    end
    fprintf('Regression parameters for complete data set\n')
    mdl = fitnlm(sampledTime(1:n),infected(1:n),@fun,b,...
        'CoefficientNames',{'K','r','A'})
    coef = mdl.Coefficients.Estimate;
    if abs(coef(1)/K(infectedidx)) > 2
        fprintf('***Warning: results of lsqcurvefit and fitnlm differs significantly. \n');
        fprintf('   Knlm/Klsq = %g\n',coef(1)/K(infectedidx));
        fprintf('   rnlm/rlsq = %g\n',coef(2)/r(infectedidx));
        fprintf('   Anlm/Alsq = %g\n',coef(3)/A(infectedidx));
    end
    fprintf('\nEvaluation of model parameters for %s\n',country)
    fprintf('%4s %10s %8s %8s %7s %5s %5s %10s %5s %5s %10s %5s\n',...
        'day','date','C','K','r','C0','Tau','end','dCdt','tpeak','peak','R2')
    fprintf('%4s %10s %8s %8s %7s %5s %5s %10s %5s %5s %10s\n',...
        ' ',' ','(cases)','(cases)','(1/day)','(cases)','(day)',' ','(c/day)','(day)',' ')
    for n = n0:length(infected)
        fprintf('%4d %10s %8d %8d %7.3f %5d %5d %10s %5d %5d %10s %5.3f\n',...
            n,datestr(sampledDate(n-1)),infected(n),K(n),r(n),C0(n),...
            tau(n),datestr(dend(n)),dCpeak(n),tpeak(n),datestr(dpeak(n)),R2(n));
    end
    time = 0:1:max([ceil(tpeak(n)+3*fix(2/r(n))),sampledTime(n)]);
    date = first_day + time;
    Ce = fun(b,time);
    [R2,~,RMSE,~,~] = calcR2a(infected,Ce(1:length(infected)));
    if R2 < 0.9
        fprintf('***Warning: %s R2 = %g\n',country, R2)
    end
    
    figure
    subplot(2,1,1)
    hold on
    % ...plot prediction
    plot(date,fun(b,time)/1000,'k','LineWidth',2)
    %... plot +/-SDE
    nsd = 3;
    sf = 1000;
    h = plot(date,(Ce + nsd*RMSE)/sf,'r','LineWidth',1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    Ct = (Ce - nsd*RMSE)/sf;
    Ct(Ct<0) = 0;
    h = plot(date,Ct,'r','LineWidth',1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % for i = 1:length(time)
    %     fprintf('%s  %g  %g\n',datestr(date(i)),time(i),fun(b,time(i)));
    % end
    %ax = gca;
    %ax.XTick = linspace(date(1),date(end),7);
    %------------------------------
    % ...plot limits epidemy phases
    ylm = get(gca,'Ylim');  % get y-axes limits
    xlm = get(gca,'Xlim');  % get x-axes limits
    www = xlm(2);
    hhh = ylm(2);
    tp2 = tpeak(n) - fix(2/r(n)) + first_day;
    tp3 = tpeak(n) + fix(2/r(n)) + first_day;
    %tp4 = 2*tpeak(n) + first_day;  % end date
    tp4 = tpeak(n) + 2*fix(2/r(n)) + first_day;  % end date
    tp = tpeak(n) + first_day;
    h = plot([tp,tp],[0,hhh],'r','LineWidth',1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % tp = tpeak(n) - fix(2/r(n)) + first_day;
    % h = plot([tp,tp],[0,hhh],'r','LineWidth',1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h = fill([tp2,tp3,tp3,tp2],[0 0 hhh hhh],'r','FaceAlpha',0.15,'EdgeColor','none');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %h = plot([tp3,tp3],[0,hhh],'r','LineWidth',1);
    %h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h = fill([tp3,tp4,tp4,tp3],[0 0 hhh hhh],'y','FaceAlpha',0.15,'EdgeColor','none');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %h = plot([tp4,tp4],[0,hhh],'r','LineWidth',1);
    %h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h = fill([tp4,www,www,tp4],[0 0 hhh hhh],'g','FaceAlpha',0.15,'EdgeColor','none');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %------------
    %... plot cases limits
    h = plot(date,K(n)*ones(length(date),1)/1000,'g--','LineWidth',1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %... plot data
    scatter(sampledDate,infected/1000,50,'k','filled')
    h = scatter(sampledDate,infected/1000,30,'w','filled');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % for i = 1:length(sampledDate)
    %     fprintf('%s  %g\n',datestr(sampledDate(i)),sampleC(i));
    % end
    %-------------------------
    % ...add axes labels
    %xtk = get(gca, 'XTick');
    datetick('x',2,'keepticks')
    %xlim([date(1) date(end)])
    %set(gca, 'XTick',date(1):7:date(end))
    legend('Prediction','Actual','Location','best')
    ylabel('Infected (x1000)')
    xlabel('Date')
    txt1 = sprintf('Coronavirus Epidemic in %s on: %s (Logistic Model)',(country),datestr(sampledDate(infectedidx-1)));
    txt2 = sprintf('K = %d  r =  %g   C_0 =  %d   RMSE = %d',...
        fix(K(n)),r(n),fix(C0(n)),fix(RMSE));
    title({txt1,txt2},'FontWeight','normal')
    grid on
    % add growth rate
    subplot(2,1,2)
    hold on
    plot(date,dfun(b,time),'k','LineWidth',2)
    ddC=diff(infected);
    ddC(ddC<0) = 0;
    bar(sampledDate(1:(end-1)),ddC,'b'); %,30,'b','filled')
    plot(date,dfun(b,time),'r','LineWidth',1)
    % dd = 7;
    % xtk = get(gca, 'XTick')
    % set(gca, 'XTick', [date(1):dd:date(end)])
    ylabel('Infected/day')
    datetick('x',2,'keepticks')
    %... add legend
    legend('Predicted','Actual','Location','best')
    % ...add title
    % ...add grid
    grid on
    %... add tick marks
    % % xtk = get(gca, 'XTick')
    % dd = 7;
    % set(gca, 'XTick', [date(1):dd:date(end)])
    tx1 = sprintf('Infection Rate in: %s on: %s',country,datestr(sampledDate(infectedidx-1)));
    tx2 = sprintf('Total Infected: %0.0f, Total Dead: %0.0f',max(Countrytotalinfected),max(Countrytotaldead));
    title({tx1,tx2},'FontWeight','normal')
    xlabel('Date')
    hold off
    
    % plot evaluation of final size
    figure
    hold on
    bar(sampledDate(n0:n),K(n0:n)/1000); %,'LineWidth',2)
    datetick('x',20,'keeplimits')
    ylabel('Infected (x1000)')
    xlabel('Date')
    title(sprintf('Daily estimated final size of epidemic for %s',country))
    mxx = max(K(n0:n)/1000);
    if mxx > 2*K(n)/1000
        ylim([0 2*K(n)/1000]);
    end
    grid on
    hold off
    
    predict()
    warning('on')
    fprintf('Summary\n')
    fprintf('-------\n')
    fprintf(' date: %10s  day: %3d\n',datestr(sampledDate(infectedidx-1)),infectedidx);
    fprintf(' start date: %s \n',datestr(first_day));
    fprintf(' number of cases: %d\n',infected(infectedidx));
    fprintf(' number of deaths: %d\n',max(Countrytotaldead));
    fprintf(' estimated epidemic size (cases): %d\n',K(infectedidx));
    fprintf(' estimated epidemic rate (1/day): %d\n',r(infectedidx));
    fprintf(' estimated initial state (cases): %d\n',C0(infectedidx));
    fprintf(' estimated initial doubling time (day): %3.1f\n',round(log(2)/r(infectedidx),1));
    fprintf(' estimated duration of fast growth phase (day): %d\n',tau(infectedidx));
    fprintf(' estimated peak date: %s  day: %d \n',datestr(dpeak(infectedidx)),tpeak(infectedidx));
    fprintf(' estimated peak rate (cases/day): %d\n',dCpeak(infectedidx));
    fprintf(' estimated end of transition phase: %s  day: %d\n',datestr(dend(infectedidx)),2*tpeak(infectedidx));
    tp2 = tpeak(infectedidx) - fix(2/r(infectedidx));
    tp3 = tpeak(infectedidx) + ceil(2/r(infectedidx)) ;
    tp4 = 2*tpeak(infectedidx);  % end date
    tp = tpeak(infectedidx) ;
    if infectedidx < tp2
        phase = 1;
        txt = 'slow growth';
    elseif infectedidx < tp
        phase = 2;
        txt = 'fast growth acceleration phase';
    elseif infectedidx < tp3
        phase = 3;
        txt = 'fast growth deceleration phase';
    elseif infectedidx < tp4
        phase = 4;
        txt = 'slow growth transition phase';
    else
        phase = 5;
        txt = 'ending phase';
    end
    fprintf(' epidemic phase: %g/5 %s\n',phase, txt);
end
    function [y,J] = fun(b,t)
        % Logistic model
        y = b(1)./(1 + b(3)*exp(-b(2)*t));
        if nargout > 1
            nt = length(t);
            J = zeros(3,nt);
            J(1,1:nt) = 1./(1 + b(3)*exp(-b(2)*t));
            J(3,1:nt) = -b(1)*exp(-b(2)*t).*J(1,1:nt).^2;
            J(2,1:nt) = -b(3)*t.*J(3,1:nt);
            J = J';
        end
    end
    function y = dfun(b,t)
        y = b(1)*b(2)*b(3)*exp(-b(2)*t)./(1 + b(3)*exp(-b(2)*t)).^2;
    end
    function predict
        %PREDICT - short time forcasting.
        nc = infectedidx;
        K0 = K(nc);
        r0 = r(nc);
        A0 = A(nc);
        
        fprintf('\nShort-term forecasting for %s\n',country)
        fprintf('%4s %12s %10s %10s %10s %10s %10s  %10s\n',...
            'day','date','actual','predict','error %','c./day act.','c./day pred.','error %')
        for n = nc-5:nc+6
            
            if n <= nc
                c0 = infected(n) - infected(n-1);
                c1 = round(funC(n-1),0)-round(funC(n-2),0);
                fprintf('%4d %12s %10d %10d %10.2f %10d %10d  %10.2f\n',n,datestr(n-2+first_day),infected(n),round(funC(n-1),0),...
                    abs(round(funC(n-1),0)/infected(n)-1)*100,...
                    c0,c1,abs(c1/c0-1)*100) ;
            else
                if round(funC(n-1),0) > infected(length(infected))
                    fprintf('%4d %12s %10s %10d %10s %10s %10d\n',n,datestr(n-2+first_day),'-',...
                        round(funC(n-1),0),'-','-',round(funC(n-1),0)-round(funC(n-2),0));
                else
                    fprintf('%4d %12s %10s\n',n,datestr(n-2+first_day),'**** Fail: actual > predict');
                    break
                end
            end
        end
        
        %=======================================
        
        function C = funC(t)
            C = K0/(1+A0*exp(-r0*t));
        end
    end
    function R2 = calcR2(t,C,b)
        K00 = b(1);
        r00 = b(2);
        A00 = b(3);
        zbar = sum(C)/length(C);
        SStot = sum((C - zbar).^2);
        SSres = sum((C - funC(t)).^2);
        R2 = 1 - SSres/SStot;
        function Ca = funC(t)
            Ca = K00./(1+A00*exp(-r00*t));
        end
    end
    function [R2, AdjR2, RMSE, Fval,pval] = calcR2a(y,ye)
        %CALCR2 Calculate the coefficient of determination
        %
        % Input:
        %   y  -- actual values
        %   ye -- estimated values
        %
        % Output:
        %   R2 -- the coefficient of determination
        %   AdjR2 -- adjusted R2
        %
        % References:
        %   https://en.wikipedia.org/wiki/Coefficient_of_determination
        %
        n = length(y);  % number of data points
        p = 4;          % number of explanatory terms in a model
        
        ybar  = sum(y)/n;
        SStot = sum((y - ybar).^2);
        SSres = sum((y - ye).^2);    %
        R2    = 1 - SSres/SStot;
        
        % calculate adjusted R2
        if nargout > 1
            AdjR2 = 1 - (1 - R2)*(n - 1)/(n - p - 1);
        end
        
        if nargout > 2
            % http://facweb.cs.depaul.edu/sjost/csc423/documents/f-test-reg.htm
            SSM =  sum((ye - ybar).^2); % sum of squares for regression
            SSE = SSres;                % sum of squares for residuals
            SST = SStot;                % sample variance x (n-1)
            dfm = p - 1;            % Corrected Degrees of Freedom for Model
            dfe = n - p;            %Degrees of Freedom for Error
            dft = n - 1;            %Corrected Degrees of Freedom Total
            MSM = SSM/dfm;          % Mean of Squares for Model
            MSE = SSE/dfe;          % Mean of Squares for Error (variance of the residuals)
            %MST = SST/DFT;         % Mean of Squares Total (sample variance)
            RMSE = sqrt(MSE); %sqrt(SStot/(n - p));  % standard error of estimate
            % calculate F statistics
            Fval = MSM/MSE;
            pval = fcdf(1/max(0,Fval),dfe,dfm); %????
        end
        
    end

    function [b0] = iniGuess(C)
        b0 = [];
        n = length(C);
        if mod(n,2) == 0
            k1 = 1;
            k3  = n-1;
        else
            k1 = 1;
            k3 = n;
        end
        k2 = (k1 + k3)/2;
        m = k2 - k1 -1;
        q = C(k2)^2 - C(k3)*C(k1);
        if q <= 0
            return
        end
        p = C(k1)*C(k2) - 2*C(k1)*C(k3) + C(k2)*C(k3);
        if p <= 0
            return
        end
        K = C(k2)*p/q;
        r = log(C(k3)*(C(k2) - C(k1))/C(k1)/(C(k3) - C(k2)))/m;
        if r < 0
            return
        end
        A = (C(k3) - C(k2))*(C(k2) - C(k1))/q*...
            (C(k3)*(C(k2) - C(k1))/C(k1)/(C(k3) - C(k2)))^((k3-m)/m);
        if A <= 0
            return
        end
        b0 = [K r A]';
    end
end