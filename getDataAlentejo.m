function [country,C,date0] = getDataAlentejo()
%GETDATA Coronavirus data
%  https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_Portugal
country = 'Alentejo';
%Alentejo
date0=datenum('2020/03/25'); % start date
C = [
    12
20
30
34
41
45
50
54
59
62
63
82
84
85
93
94
125
130
139
140
155
155
156
158
%<-------------- add new data here
]';
end

