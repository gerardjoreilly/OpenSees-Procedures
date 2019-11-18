% This is a testing script for the computation of AvgSa

%% Clean up 
addpath('/Users/gor/GitHub/Matlab-Functions')
% run('startup.m');
close all
clear 
clc

%% Compute the spectrum
[~,~,Sa,~,~,~,T] = RsSpNewmark(load('/Users/gor/GitHub/Matlab-Functions/eq_1g_0.01s.txt'),0.01,0.05,0.01,0.01,3,0.5,0.25);

%% Compute the spectral values at set periods
Tlist = [0.27 0.46 0.65 0.83 1.02 1.21 1.40 1.58 1.77 1.96];

for i = 1:length(Tlist)
    [~,~,Sai(i),~,~,~]=RsSpNewmarkT(load('/Users/gor/GitHub/Matlab-Functions/eq_1g_0.01s.txt'),0.01,0.05,Tlist(i),0.5,0.25);
end

Sai
n = length(Sai)
AvgSa = prod(Sai)^(1/n)
AvgSa = exp(sum(log(Sai))/n)
AvgSa = geomean(Sai)

%% Plot the spectrum
figure
grid on; hold on;
plot(T,Sa,'-k')
plot(Tlist,Sai,'or');
xlabel('T')
ylabel('Sa(T)')
