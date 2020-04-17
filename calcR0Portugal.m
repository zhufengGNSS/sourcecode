function out = calcR0Portugal(getDataPortugal,nday)
%CALCR0 Calculate basic reproduction number 
%
% Optional input:
%       nday -- sample length (def 14)
%
% Output: struucture
%       out.R0 -- basic reproduction number;
%       out.stdR0 -- std of R0
%       out.N  -- initial population size
%       out.stdN -- stdev of N
%
%Example: 
%   calcR0(@getDataItaly)
%
    narginchk(1,2)
    nargoutchk(0,1)

    % get data
    res = fitVirusCV19Portugal(getDataPortugal,'plt','off');
    nnday = res.day;
    
    % set sample length
    if nargin < 2
        nday = 14;
    end
    nday = min(nnday,nday);

    R0 = NaN(nday,1);
    N  = NaN(nday,1);
 
    nn = 0;
    for n = 1:nnday
        res = fitVirusCV19Portugal(getDataPortugal,'day',n,'plt','off');
        if ~isempty(res)
            if abs(res.R0) < 1e-3 || res.R0 > 6
                continue
            end
            nn = nn + 1;
            R0(nn) = res.R0;
            N(nn) = res.N;
            if nn > 2
            mR0 = mean(R0(1:nn));
            stdR0 = std(R0(1:nn));
            if stdR0/mR0 < 0.4 && nn > 4
                break
            end
            end
            if nn >= nday
                break
            end
        else
            fprintf('Erro - Dia %d\n',n)
        end
    end
    if nn == 0
        if nargout > 0
            out = [];
        end
        fprintf('Erro no Cálculo.\n')
        return
    end
    
    % clean data
    R0 = R0(1:nn);
    N = N(1:nn);
    
    % pack output
    if nargout > 0
        out.country = res.country;
        out.R0 = mean(R0);
        out.stdR0 = std(R0);
        out.N  = mean(N);
        out.stdN = std(N);
        out.nday = length(R0);
    end
end