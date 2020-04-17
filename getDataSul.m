function [country,C,date0] = getDataSul()
%GETDATA Coronavirus data
%  https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_Portugal
country = 'Sul';
%Sul
date0=datenum('2020/03/10'); % start date
C = [
    10
17
23
46
73
116
142
180
243
278
364
448
534
737
852
992
1082
1110
1287
1478
1577
1799
1998
2207
2347
2513
2904
3070
3185
3424
3451
3821
3834
3841
3896
3994
4102
4237
4302
%<-------------- add new data here
]';
end

