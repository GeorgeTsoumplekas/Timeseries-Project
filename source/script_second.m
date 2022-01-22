clc
clear

%Load the second timeseries
file=importdata('VideoViews.xlsx');
Y2=file(:,15);
size = length(Y2);

figure_num=1;

%Part 0: Visualize the original timeseries
figure(figure_num)
plot(Y2,'.-');
title('History diagram of the second time series');
xlabel('Days');
ylabel('View count');
figure_num=figure_num+1;
%Finished with the history diagramms. By looking at it we could say that
%the timeseries is not stationary.

%End of Part 0


%Part 1: We create an autocorrelation diagram, to get information about the
%stationarity of the timeseries

%We choose maxtau_r=100 as an appropriate tau. It is large enough to give
%us the info we want(whether there is stationarity or not)
maxtau_r=100;                                            
alpha=0.05;
zalpha = norminv(1-alpha/2);
autlim = zalpha/sqrt(size);

rt2 = autocorrelation(Y2, maxtau_r);
figure(figure_num)
figure_num=figure_num+1;
hold on
for ii=1:maxtau_r
    plot(rt2(ii+1,1)*[1 1],[0 rt2(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau_r+1],[0 0],'k','linewidth',1.5)
%We form the horizontal lines, to have a visual perspective on a more
%subjective way to tell if the rt is large or not (by comparing with
%statistical significance)
plot([0 maxtau_r+1],autlim*[1 1],'--c','linewidth',1.5)                 
plot([0 maxtau_r+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title('The autocorrelation diagram of the second timeseries');
%We observe that the autocorrelation values are very large and they
%dissipate slowly as tau gets larger. We also observe the existence of
%spikes and periodicity in this plot.
%This makes me believe that the timeseries has both a trend and some periodicity


%Remove the trend using a moving average filter. The order of the filter
%(2q+1) is set to 5 since it looks to work better for this value than
%others
maorder_trend=5;

%Trend contains the mean values 
trend= movingaveragesmooth2(Y2,maorder_trend);              
Y2_detrend=Y2-trend;
figure(figure_num)
figure_num=figure_num+1;
plot(Y2_detrend,'.-');
title('History diagram of the second time series, detrended');
xlabel('Days');
ylabel('View change');

%Autocorrelation diagram of the detrended timeseries
rt2_det = autocorrelation(Y2_detrend, maxtau_r);
figure(figure_num)
figure_num=figure_num+1;
hold on
for ii=1:maxtau_r
    plot(rt2_det(ii+1,1)*[1 1],[0 rt2_det(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau_r+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau_r+1],autlim*[1 1],'--c','linewidth',1.5)                  
plot([0 maxtau_r+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title('The autocorrelation diagram of the detrended timeseries');
%By looking at this plot, it is clear that we have some seasonality.
%The period seems to be 7 days.

%Remove seasonality using a seasonal components filter with period = 7
maorder_seasonality=7;
s2 = seasonalcomponents(Y2_detrend,maorder_seasonality);

%X2 is the equivalent stationary timeseries we are going to work with from
%now on.
X2=Y2_detrend-s2;             

%Plot of X2
figure(figure_num)
figure_num=figure_num+1;
plot(X2,'.-');
title('History diagram of the stationary timeseries');
xlabel('Days');
ylabel('View change');

%Autocorrelation diagram of the stationary timeseries
rt2_WSS = autocorrelation(X2, maxtau_r);
figure(figure_num)
figure_num=figure_num+1;
hold on
for ii=1:maxtau_r
    plot(rt2_WSS(ii+1,1)*[1 1],[0 rt2_WSS(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau_r+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau_r+1],autlim*[1 1],'--c','linewidth',1.5)                
plot([0 maxtau_r+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title('The autocorrelation diagram of the detrended and deseasoned timeseries');

%From the autocorrelation diagram, we can see that ,without large doubt,
%the timeseries looks stationary since both trend and periodicity have been
%eliminated.We also see that the autocorrelation is not one that a white
%noise would have had. So, we conclude that the timeseries is not white noise.

%End of Part 1.


%Part 2: Fit an ARMA model to the timeseries

%Plot of the partial autocorrelation, to help to decide on the (p,q) parameters.
phi2 = parautocor(X2,maxtau_r);
figure(figure_num)
figure_num=figure_num+1;
hold on
for ii=1:maxtau_r
    plot(rt2_WSS(ii+1,1)*[1 1],[0 phi2(ii)],'b','linewidth',1.5)
end
plot([0 maxtau_r+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau_r+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau_r+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('\phi_{\tau,\tau}')
title('The partial autocorrelation diagram of the stationary timeseries');


%Judging from these plots, we say that a proper q would have to be at least
%4 and the proper p should be around 16. We decided to make some
%brute force computations, to find the (p,q) which give a model with the
%minimum aic criterion value.We see that the returned value is (10,6) which
%is not a very simple model but apparently is the most appropriate
%(considering that the aic takes into consideration the complexity of the
%model too). Moreover, we see that
%for some values of p and q a non-stationary or non-reversible model may
%occur. However, this does not impose a problem and comparing the aic
%values of the feasible models we can choose which one of them is the best.

%We have put a limit, so that p,q won't be any larger than 10.
p=1;
aic=inf;
for i=1:10
    q=1;
    for j=1:10
        [~,~,~,~,aicS,~]=fitARMA(X2,p,q,10);
        if (aicS<aic)
            aic=aicS;
            p_final=p;
            q_final=q;
        end
        q=q+1;
    end
    p=p+1;
end

fprintf('p = %f\n',p_final);
fprintf('q = %f\n',q_final);
fprintf('The most appropriate ARMA model is ARMA(%d,%d)\n',p_final,q_final);



%Fitting the ARMA model
p=10;
q=6;
Tmax=5;
[nrmseV,phiallV,thetaallV,SDz,aicS,fpeS,sarmamodel2]=fitARMA(X2,p,q,Tmax);
sarmamodel2;
fprintf('===== ARMA model ===== \n');
fprintf('Estimated coefficients of phi(B):\n');
for ip=1:p+1
    fprintf('phi%d = %f\n',ip-1,phiallV(ip));
end
fprintf('\nEstimated coefficients of theta(B):\n');
for iq=1:q
    fprintf('theta%d = %f\n',iq,thetaallV(iq));
end
fprintf('\nSD of noise: %f \n',SDz);
fprintf('AIC: %f \n',aicS);
fprintf('FPE: %f \n',fpeS);
fprintf('\n\t T \t\t NRMSE \n');
disp([[1:Tmax]' nrmseV])
if Tmax>3
    figure(figure_num)
    figure_num=figure_num+1;
    plot([1:Tmax]',nrmseV,'.-k')
    hold on
    plot([1 Tmax],[1 1],'r')
    xlabel('T')
    ylabel('NRMSE')
    title(sprintf('ARMA(%d,%d), fitting error',p,q))
    ylim([0 1.1]);
end    

%Apply Portmanteau test to test on whether the residuals are white noise
%or not. If they are not, then there is information that our linear model
%could not describe.
x2preV=predict(sarmamodel2, X2, 1);
errV2=x2preV -X2;
figure(figure_num)
figure_num=figure_num+1;
[hV,pV,qV,xautV]=portmanteauLB(errV2, 15, 0.05, 1);
%Due to the big p-values we cannot reject the hypothesis that the errors'
%time series is white noise. So, the ARMA model seems to fit well.

%Finished with part 2.


%Part 3: Finding the change points. We use the initial fitted model for all
%predictions throughout

T=5;
training_size=400;

%The initial value of n.
%WARNING! THIS VARIABLE IS NOT THE VARIABLE FOR THE SIZE OF THE TRAINING
%SET, just the point after which it starts predecting
n=400;     
change_count=0;

while (n<size-5)
    [s,a]=compute_stat(X2, p, q, n, T, training_size,1.5);
    if(s>a)
        change_count=change_count+1;
        %Vector containing the days in which we have a change
        changes(change_count)=n; 
      
        fprintf('A change in demand has occured!, with s=%f and a=%f for n=%d\n',s,a,n)
        n=n+T;
    else
        n=n+1;
    end
end

%Repeat for the last 5 values in case it skipped them in the loop
n=size-5;
[s,a]=compute_stat(X2, p, q, n, T, training_size,1.5);
if(s>a)
        change_count=change_count+1;
        changes(change_count)=n;
end

%Plot original timeseries with change points marked
figure(figure_num)
figure_num=figure_num+1;
plot(Y2,'.-')
hold on
scatter(changes, Y2(changes), 'k', 'filled')
title('Original timeseries with change points marked');
xlabel('Days');
ylabel('Daily views');

%Plot stationary timeseries with change points marked
figure(figure_num)
figure_num=figure_num+1;
plot(X2,'.-')
hold on
scatter(changes, X2(changes), 'k', 'filled')
title('Detrended timeseries with change points marked');
xlabel('Days');
ylabel('Daily views');

%End of Part 3


%Part 4: Fitting a non-linear chaotic model for the timeseries

%Create the 2-d and 3-d scatterplots
figure(figure_num)
scatter(X2(1:end-1),X2(2:end))
xlabel('X(t-1)');
ylabel('X(t)');
title('2-d scatterplot of the stationary timeseries')

figure_num = figure_num+1;
figure(figure_num)
scatter3(X2(1:end-2),X2(2:end-1),X2(3:end))
xlabel('X(t-2)');
ylabel('X(t-1)');
zlabel('X(t)');
title('3-d scatterplot of the stationary timeseries')

%Set delay tau equal to 1 because we have a discrete system (1 value each
%day)
tau = 1;

%Calculate embedding dimension m
mmax = 10;
escape = 10;
theiler = 0;

figure_num = figure_num+1;
figure(figure_num)
fnnM = falsenearest(X2,tau,mmax,escape,theiler,'2nd timeseries');

%Calculate correlation dimension
fac = 4;
figure_num = figure_num+1;
[~,~,~,~,nuM] = correlationdimension(X2,tau,mmax,'1st timeseries',fac);
figure_num = figure_num+4;

fprintf('\nEstimation of correlation dimension v\n');
for i=1:mmax
    fprintf('For m=%d, v=%f\n',nuM(i,1),nuM(i,4));
end

%We will choose m to be equal to 5 since every value of correlation
%dimension is < 5. Moreover, from the calculation of the false neighbors
%percentage we see that it is already small for m=4 so it will be even
%smaller for m=5 which is something we want.
m=5;

%Local linear prediction model
T=5;

%The initial value of n.
%WARNING! THIS VARIABLE IS NOT THE VARIABLE FOR THE SIZE OF THE TRAINING
%SET, just the point after which it starts predecting
n=training_size+m;

%Find the appropriate number of neighbors K by calculating the nrmse we
%have for different values of K for the 5 first values predicted
iMax = floor(log2(training_size));

nrmse_values = zeros(iMax);

q = 0;
for i=1:iMax
    K = 2^i;
    x_pred = localpredictmultistep(X2,training_size,tau,m,T,K,q);
    nrmse_values(i) = nrmse(X2(n+1:n+T),x_pred);
end

figure_num = figure_num+1;
figure(figure_num);
plot(nrmse_values);
hold on
plot([1 iMax],[1 1],'--r')
xlabel('log_2(K)')
ylabel('NRMSE')
title('NRMSE for different number of neighbors for the first timeseries');

%For m=2 best K is 4
K=4;

%Calculate the change points
changes2 = zeros(1,2);
x_pred = zeros(size-n,1);
change_count = 0;
while (n<size-T)
    [s,a,x_pred_temp] = compute_stat2(X2,K,n,m,T,training_size,1.5);
    x_pred(n-(training_size-1)-m:n+T-training_size-m) = x_pred_temp(:);
    if s>a
        change_count = change_count+1;
        changes2(change_count) = n;
        fprintf('A change in demand has occured!, with s=%f and a=%f for n=%d\n',s,a,n)
        n = n+T;
    else
        n = n+1;
    end
end

%Repeat for the last 5 values in case it skipped them in the loop
n=size-5;
[s,a,x_pred_temp] = compute_stat2(X2,K,n,m,T,training_size,1.5);
x_pred(n-(training_size-1)-m:n+T-training_size-m) = x_pred_temp(:);
if s>a
    change_count = change_count+1;
    changes2(i,change_count) = n;
end

%Plot original timeseries with change points marked
figure_num=figure_num+1;
figure(figure_num)
plot(Y2,'.-')
hold on
scatter(changes2, Y2(changes2), 'k', 'filled')
title('Original timeseries with change points marked');
xlabel('Days');
ylabel('Daily views');

%Plot stationary timeseries with change points marked
figure_num=figure_num+1;
figure(figure_num)
plot(X2,'.-')
hold on
scatter(changes2, X2(changes2), 'k', 'filled')
title('Detrended timeseries with change points marked');
xlabel('Days');
ylabel('Daily views');

%Plot the predicted values
n=training_size+m;
figure_num=figure_num+1;
plotrealpred(X2(n+1:end),x_pred,'predicted');