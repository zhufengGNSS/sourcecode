function [country,C,date0] = getDataAlgarve()
%GETDATA Coronavirus data
%  https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_Portugal
country = 'Algarve';
%Algarve
date0=datenum('2020/03/25'); % start date
C = [
    62
89
99
106
108
116
137
146
164
179
182
201
229
234
251
260
279
279
279
284
289
295
305
%<-------------- add new data here
]';
end

