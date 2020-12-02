% This is a replication of the 2017 Levin-Reisman model

close all; clear all;clc

%% define parameters 
K = 1;
MIC = 1;
c = 10;
n = 2;
uMin = 2;
uMax = 1;
Nmin = 1e-3;
u = .003;
T = 60;
t = T-1;
tau = 10;

%% assumed distributions to create a negative relationship btwn survival and u
u = exprnd(1e-8,1,1000); % mutation rate (low --> target, high --> metabolic)
KHALF = 1e-14;
n=2;
sigman = .05*n;      %noise factor for Michaelis-Menten n (5%)
sigmak = .05*KHALF;  %noise factor for Michaelis-Menten K (5%)
nuse = sigman*randn(1,1000)+n;
KHALFuse = sigmak*randn(1,1000)+KHALF;
Sr = KHALFuse.^nuse./(KHALFuse.^nuse+u.^nuse);

%% calculate probability of establishing a mutation
P = @(u,Nmin,t,T,Sr)u * Nmin * 2 .^ t .* (1 - (1-Sr) .^ (2 .^ (T-t)));
prob1 = P(u,Nmin,t,T,Sr);

%% plot
figure; scatter(u,prob1)
set(gca,'xscale','log','yscale','log')
xlabel('mutation frequency')
ylabel('P')


%% Now over multiple generations
uMean = logspace(-10,-8,5);
figure, hold on
for ii = 1:length(uMean)
    
    u = exprnd(uMean(ii),1,1000);
    nuse = sigman*randn(1,1000)+n;
    KHALFuse = sigmak*randn(1,1000)+KHALF;
    Sr = KHALFuse.^nuse./(KHALFuse.^nuse+u.^nuse);
    
    % Calculate and sum probabilities over all T generations
    prob1 = zeros(size(u));
    for tt = 0:(T-1)
        P = @(u,Nmin,t,T,Sr)u * Nmin * 2 .^ t .* (1 - (1-Sr) .^ (2 .^ (T-t)));
        prob1 = prob1 + P(u,Nmin,tt,T,Sr);
        prob1(prob1>1) = 1;
    end
    
    scatter(u,prob1,'filled')
    
end
set(gca,'xscale','log','yscale','log','linewidth',4.0,'fontsize',40)
xlim([10^-15 10^-5])
ylim([10^-6 10^0])
set(gca,'ytick',[10^-6 10^-3 10^0])
xlabel('mutation frequency')
ylabel('\mu')