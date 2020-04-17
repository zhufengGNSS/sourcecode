function [country,C,date0] = getDataCentro()
%GETDATA Coronavirus data
%  https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_Portugal
country = 'Centro';
%Centro
date0=datenum('2020/03/10'); % start date
C = [
    1
3
5
6
8
10
31
51
74
86
106
137
180
238
293
365
435
520
647
709
784
911
1043
1161
1286
1372
1442
1521
1766
1865
1905
2197
2327
2426
2477
2549
2629
2756
2778
%<-------------- add new data here
]';
end

