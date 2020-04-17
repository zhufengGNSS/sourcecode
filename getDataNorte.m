function [country,C,date0] = getDataNorte()
%GETDATA Coronavirus data
%  https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_Portugal
country = 'Norte';
%Norte
date0=datenum('2020/03/18'); % start date
C = [
    289
381
506
644
825
1007
1130
1517
1858
2443
3035
3550
3801
4452
4910
5338
5899
6280
6530
6706
7052
7386
8102 %20/04/09
8897
9264
9747
9984
10302
10751
11237
11324
%<-------------- add new data here
]';
end

