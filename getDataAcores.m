function [country,C,date0] = getDataAcores()
%GETDATA Coronavirus data
%  https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_Portugal
country = 'Açores';
%Acores
date0=datenum('2020/03/23'); % start date
C = [
    11
12
17
24
24
30
33
41
48
52
57
63
63
67
68
68
70
91
94
94
94
94
100
100
102
102
%<-------------- add new data here
]';
end

