function h = histfit_HXB(data,linespec,is_bar)
% revised histfit by Huixin Bi
% is_bar = n.a. : plot only approximated normal distribution
%        = 1    : plot approximated normal and actual bar distributions
%        = 2    : plot only actual density curve 
% date: 05/18/2010


[row,col] = size(data);
if min(row,col) > 1
   error('stats:histfit:VectorRequired','DATA must be a vector.');
end

if row == 1
  data = data(:);
end
row = sum(~isnan(data));

inter = 0.1;    % width of the interval

if nargin ==2
    mr = nanmean(data); % Estimates the parameter, MU, of the normal distribution.
    sr = nanstd(data);  % Estimates the parameter, SIGMA, of the normal distribution.
    x=(-3*sr+mr:inter*sr:3*sr+mr)';% Evenly spaced samples of the expected data range.
    y = normpdf(x,mr,sr); 
    h = plot(x,y,linespec,'LineWidth',2);
    
    % axis_tmp = axis;
    % text(mr+2*sr,(axis_tmp(3)+axis_tmp(4))/2,textspec)

elseif nargin ==3
    
    % show the bar graph of simulated b*
    if is_bar == 1
        mr = nanmean(data); % Estimates the parameter, MU, of the normal distribution.
        sr = nanstd(data);  % Estimates the parameter, SIGMA, of the normal distribution.
        x=(-3*sr+mr:inter*sr:3*sr+mr)';% Evenly spaced samples of the expected data range.

        [count xbin]=histc(data,x);
        h1 = bar(x,count/row/(0.1*sr),'histc');
        hold on

        y = normpdf(x,mr,sr); 
        h2 = plot(x,y,'r','LineWidth',2);
        hold off
        h = [h1 h2];
        
    elseif is_bar ==2
        mr = nanmean(data); % Estimates the parameter, MU, of the normal distribution.
        sr = nanstd(data);  % Estimates the parameter, SIGMA, of the normal distribution.
        x=(-3*sr+mr:inter*sr:3*sr+mr)';% Evenly spaced samples of the expected data range.

        [count xbin]=histc(data,x);
        h3 = plot(x, count/row/(0.1*sr),linespec,'LineWidth',2);
        h = h3;
    end
end



% % test
% load bstar_Benchmark
% data = bstar_vec_Benchmark/0.25;
% numbins = 100;
% n = length(data);
%   binwidth = range(data)/numbins;
%   edg= min(data):binwidth:max(data);
%   [count, bin] = histc(data, edg);
%   h = bar(edg, count./(n*binwidth), 'histc');
%   set(h, 'facecolor', [0.8 0.8 1]); % change the color of the bins
%   set(h, 'edgecolor', [0.8 0.8 1]);
% p = count./(n);
% figure
% plot(edg,p);    % plot(p,edg) is the smooth curve representing the probability density function you are looking for.
% 
%  