clc;clear;close all;warning off;

load('Project_3_Data.mat');

%% Question 1. Stimulus Nature

stim=Stimulus;
auto_0_50=autocorr(stim,50);
auto_minus_50_0=autocorr(fliplr(stim),50);
concat=cat(2,fliplr(auto_minus_50_0(1:50)),auto_0_50);

figure;
plot([-50:50],concat);     
xlabel("\tau(ms)");
ylabel("R(\tau)");
title("Autocorrelation Function of Stimulus");

%% Question 2. PSTH and Mean Firing Rate

PSTH = zeros(4, 20000);
for i=1:4
    for j=1:50
        PSTH(i,:) = PSTH(i,:) + histcounts(All_Spike_Times{i,j}*1000,0:20000)*1000/50;
    end
end

figure;
neuron_1=subplot(4,1,1);
plot(PSTH(1,1:20000));
xlabel(neuron_1,'Time(ms)');
ylabel(neuron_1,'Rate(Spikes/s)');
title(['PSTH for the First Neuron']);

neuron_2=subplot(4,1,2);
plot(PSTH(2,1:20000));
xlabel(neuron_2,'Time(ms)');
ylabel(neuron_2,'Rate(Spikes/s)');
title(['PSTH for the Second Neuron']);

neuron_3=subplot(4,1,3);
plot(PSTH(3,1:20000));
xlabel(neuron_3,'Time(ms)');
ylabel(neuron_3,'Rate(Spikes/s)');
title(['PSTH for the Third Neuron']);

neuron_4=subplot(4,1,4);
plot(PSTH(4,1:20000));
xlabel(neuron_4,'Time(ms)');
ylabel(neuron_4,'Rate(Spikes/s)');
title(['PSTH for the Fourth Neuron']);

%% Question 3. Poisson or Non-Poisson

time_bins = [10, 20, 50, 100, 200, 500];

for i = 1:6 
    bin_size = time_bins(i);
    figure
    
    for n=1:4   
        mean_neurons = [];
        variance_neurons = [];
        for j=1:50  
            spike_count_values = histcounts(All_Spike_Times{n,j}*1000,0:bin_size:20000); 
            M = mean(spike_count_values);
            V = var(spike_count_values);
            mean_neurons = [mean_neurons [M]];    
            variance_neurons = [variance_neurons [V]];
        end
 
        max_val = max([max(mean_neurons), max(variance_neurons)]);
        min_val = min([min(mean_neurons), min(variance_neurons)]);
        st_line = [0,max_val];
        p = subplot(2,2,n);
        scatter(variance_neurons, mean_neurons,'bx');
        hold on;
        plot(st_line, st_line,'k');
        hold off;
        title(['Mean vs Variance in Specified Bin Size ' num2str(bin_size) ' ms for Neuron ' num2str(n)]);
        xlabel('Variance');
        ylabel('Mean');
        xlim([min_val max_val]);
        ylim([min_val max_val]);
    end
end

%% Question 4. Spike Triggered Average (STA)

spike_trig_avg = zeros(4, 100);
h = zeros(4, 100);

Rxx = autocorr(Stimulus, 99);  
Css = zeros(100,100);
Css = toeplitz(Rxx);


for i=1:4
    spike_count_total = 0;
    for j=1:50
        spike_count_total = spike_count_total + nnz(All_Spike_Times{i,j} <= 15);
        for k=1:nnz(All_Spike_Times{i,j} <= 15)
            interval_t = round(All_Spike_Times{i,j}(k)*1000);
            stim_values = Stimulus(max([interval_t-99 1]):interval_t);
            partial_STA = [zeros(1,100-length(stim_values)) stim_values];
            spike_trig_avg(i,:) = spike_trig_avg(i,:) + partial_STA;  
        end
    end
    spike_trig_avg(i,:)=spike_trig_avg(i,:)/spike_count_total;
    figure(9);
    subplot(2,2,i);
    plot(0:99,fliplr(spike_trig_avg(i,:)))
    xlabel('Time(ms)');
    ylabel('STA');
    ylim([-0.5 0.5]);
    title(['h(t) Without Correction for Neuron ' num2str(i)]);
      
    h(i,100:-1:1)=(Css\spike_trig_avg(i,:)')';
    figure(10);
    subplot(2,2,i);
    plot(0:99,h(i,:)); 
    xlabel('Time (ms)');
    ylabel('h(t)');
    ylim([-1 1]);
    title(['h(t) With Correction for Neuron ' num2str(i)]);
end

%% Question 5 Determining the output Non-Linearity 
    
    stim = Stimulus;
    psth=PSTH/1000;
    
    pred(1, :) = conv(stim(1:15000), h(1, :)); x(1,1:15000) = pred(1, 1:15000); 
    pred(2, :) = conv(stim(1:15000), h(2, :)); x(2,1:15000) = pred(2, 1:15000);
    pred(3, :) = conv(stim(1:15000), h(3, :)); x(3,1:15000) = pred(3, 1:15000);
    pred(4, :) = conv(stim(1:15000), h(4, :)); x(4,1:15000) = pred(4, 1:15000);
    
 
    y(1, 1:15000) = 1000*psth(1, 1:15000); 
    y(2, 1:15000) = 1000*psth(2, 1:15000); 
    y(3, 1:15000) = 1000*psth(3, 1:15000); 
    y(4, 1:15000) = 1000*psth(4, 1:15000); 
    
    bin_size = 30;
    for i = 1:ceil(15000/bin_size)
        e = i*bin_size;
        if e>15000
            e = 15000;
        end
        x1(i) = mean( x(1, (1+(i-1)*bin_size):e) );
        x2(i) = mean( x(2, (1+(i-1)*bin_size):e) );
        x3(i) = mean( x(3, (1+(i-1)*bin_size):e) );
        x4(i) = mean( x(4, (1+(i-1)*bin_size):e) );
    
        y1(i) = mean( y(1, (1+(i-1)*bin_size):e) );
        y2(i) = mean( y(2, (1+(i-1)*bin_size):e) );
        y3(i) = mean( y(3, (1+(i-1)*bin_size):e) );
        y4(i) = mean( y(4, (1+(i-1)*bin_size):e) );
    end
    
    figure()
    subplot(2,2,1)
    scatter(x1, y1,'kx')
    title('PSTH vs y(t) for Neuron 1')
    xlabel('y(t)= s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    subplot(2,2,2)
    scatter(x2, y2,'k*')
    title('PSTH vs y(t) for Neuron 2')
    xlabel('y(t) = s(t)*h(t)')  
    ylabel('PSTH ~ \lambda(t)')
    subplot(2,2,3)
    scatter(x3, y3,'kx')
    title('PSTH vs y(t) for Neuron 3')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    subplot(2,2,4)
    scatter(x4, y4,'k*')
    title('PSTH vs y(t) for Neuron 4')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
        
    [fit, gof] = Fitting_curve(x1, y1, x2, y2, x3, y3, x4, y4);
    
    save('fits', 'fit', 'gof')
    

    pred1 = conv(stim(15001:20000), h(1, :)); pred1 = pred1(1:5000);
    pred2 = conv(stim(15001:20000), h(2, :)); pred2 = pred2(1:5000);
    pred3 = conv(stim(15001:20000), h(3, :)); pred3 = pred3(1:5000);
    pred4 = conv(stim(15001:20000), h(4, :)); pred4 = pred4(1:5000);

    for i = 1:5000
        pred1(i) = (fit{1,1}.a)/(1+exp(-fit{1,1}.b*(pred1(i)-fit{1,1}.c)));
        pred2(i) = (fit{2,1}.a)/(1+exp(-fit{2,1}.b*(pred2(i)-fit{2,1}.c)));
        pred3(i) = (fit{3,1}.a)/(1+exp(-fit{3,1}.b*(pred3(i)-fit{3,1}.c)));
        pred4(i) = (fit{4,1}.a)/(1+exp(-fit{4,1}.b*(pred4(i)-fit{4,1}.c)));
    end 
 
    gt1 = 1000*psth(1, 15001:20000);
    gt2 = 1000*psth(2, 15001:20000);
    gt3 = 1000*psth(3, 15001:20000);
    gt4 = 1000*psth(4, 15001:20000);
    

%% Question 6. Prediction Performance and Pruning Filter Parameters

R1=corrcoef(gt1,pred1); R_squared1=R1(2)^2 ;
    R1;
    R_squared1;
    
    R2=corrcoef(gt2,pred2); R_squared2=R2(2)^2 ;
    R2;
    R_squared2;
    
    R3=corrcoef(gt3,pred3); R_squared3=R3(2)^2 ;
    R3;
    R_squared3;
    
    R4=corrcoef(gt4,pred4); R_squared4=R4(2)^2;
    R4;
    R_squared4;
    
    figure()
    subplot(2,2,1)
    scatter(gt1,pred1)
    title('PSTH vs Prediction for Neuron 1')
    
 
    subplot(2,2,2)
    scatter(gt2,pred2)
    title('PSTH vs Prediction for Neuron 2')
    
    
    subplot(2,2,3)
    scatter(gt3,pred3)
    title('PSTH vs Prediction for Neuron 3')
    
    
    subplot(2,2,4)
    scatter(gt4,pred4 )
    title('PSTH vs Prediction for Neuron 4')
    
    
    A1=zeros(100,1);
    B1=zeros(100,1);
    
    initial_r_sq=R_squared2;
    old_r_sq=R_squared2;
    new_r_sq=old_r_sq;
    pos=zeros(1,1);
    count1=0;
    while (old_r_sq - new_r_sq) < 0.01
        old_r_sq=new_r_sq;
        min_=1000000000;
        for i = 1:100 
            if abs(h(2,i))<min_ && abs(h(2,i))~=0;
                min_=abs(h(2,i));
                pos=i;
            end
        end
        if pos==0
            break
        end
        h(2,pos)=0;
        pos(1,1)=0;
        pred2 = conv(stim(15001:20000), h(2, :)); pred2 = pred2(1:5000);
        for i = 1:5000
            pred2(i) = (fit{2,1}.a)/(1+exp(-fit{2,1}.b*(pred2(i)-fit{2,1}.c)));
        end
        new_r=corrcoef(gt2,pred2); new_r_sq=new_r(2)^2;
        count1=count1+1;
        A1(count1,1)=count1;
        B1(count1,1)=new_r_sq;
        
    end
    figure()
    scatter(A1,B1,'k*')
    title('Plot of Predicted performance for Neuron 2 (with Iterations)')
    
    max_corr_coef_2=0;
    for i =1:100
        if B1(i,1)>max_corr_coef_2
            max_corr_coef_2=B1(i,1);
        end
    end

    A=zeros(100,1);
    B=zeros(100,1);
    initial_r_sq=R_squared3;
    old_r_sq=R_squared3;
    new_r_sq=old_r_sq;
    pos=zeros(1,1);
    count=0;
    while (old_r_sq - new_r_sq) < 0.01
        old_r_sq=new_r_sq;
        min_=1000000000;
        for i = 1:100 
            if abs(h(3,i))<min_ && abs(h(3,i))~=0;
                min_=abs(h(3,i));
                pos=i;
            end
        end
        h(3,pos)=0;
        pos(1,1)=0;
        pred3 = conv(stim(15001:20000), h(3, :)); pred3 = pred3(1:5000);
        for i = 1:5000
            pred3(i) = (fit{3,1}.a)/(1+exp(-fit{3,1}.b*(pred3(i)-fit{3,1}.c)));
        end
        new_r=corrcoef(gt3,pred3); new_r_sq=new_r(2)^2;
        count=count+1;
        A(count,1)=count;
        B(count,1)=new_r_sq;
    end
    
    figure()
    scatter(A,B,'k*')
    title('Plot of Predicted performance for Neuron 3 (with Iterations)')
    max_corr_coef_3=0;
    for i =1:100
        if B(i,1)>max_corr_coef_3;
            max_corr_coef_3=B(i,1);
        end
    end
 
    figure()
    subplot(2,2,1)
    plot(h(1,:))
    title('Linear filter for Neuron 1')
    subplot(2,2,2)
    plot(h(2,:))
    title('Linear filter for Neuron 2')
    subplot(2,2,3)
    plot(h(3,:))
    title('Linear filter for Neuron 3')
    subplot(2,2,4)
    plot(h(4,:))
    title('Linear filter for Neuron 4')
    
    
    f_t4=fft(h(4,:));
    f_t3=fft(h(3,:));
    f_t2=fft(h(2,:));
    f_t1=fft(h(1,:));
    
    vect=[-50:49];
    vect=vect';
    
    figure()
    subplot(2,2,1)
    plot(vect,f_t1)
    title('FFT of Filter for Neuron 1')
    
    subplot(2,2,2)
    plot(vect,f_t2)
    title('FFT of Filter for Neuron 2')
    
    subplot(2,2,3)
    plot(vect,f_t3)
    title('FFT of Filter for Neuron 3')
    
    subplot(2,2,4)
    plot(vect,f_t4)
    title('FFT of Filter for Neuron 4')
    
    
    fprintf('Maximum predicted performance for the 2nd Neuron is %d.\n',max_corr_coef_2);
    fprintf('Maximum predicted performance for the 3rd Neuron is %d.\n',max_corr_coef_3);
    
    %% (B) Discrimination based on Victor and Purpura (VP) Spike Distance Metric (SDM)
    
q = [0 0.001 0.01 0.1 1 10 100];

MI = zeros(4,100,length(q));

for num_trials = 1:100
    
    t = randperm(19901,8);
    for n = 1:4
        spike_seg = cell(8,50);
        for rept = 1:30
            spikes = All_Spike_Times{n,rept};
            for i = 1:8
                spike_seg{i,rept} = spikes(spikes>=t(i)/1000&spikes<(t(i)+100)/1000);
            end
        end
        
        for m = 1:length(q)
            conf_mat = zeros(8,8);
            for i = 1:8
                for rept = 1:30
                    mean_dist = meandist(i,rept,spike_seg,q(m));
                    mean_dist(i) = mean_dist(i)*50/49;
                    [~,k] = min(mean_dist);
                    conf_mat(i,k) = conf_mat(i,k)+1;
                end
            end
            conf_mat = conf_mat/50;
            MI(n,num_trials,m) = MI(n,num_trials,m) + MutualInfo(conf_mat)/100;
        end
    end
end

ci90 = c(0.9,MI);
a = MI(:,1,:);
b = ci90(:,1,:);
figure()
for n = 1:4
    subplot(2,2,n);
    plot(log10(q), a(n,:),'g',"LineWidth",1);
    hold on
    plot(log10(q), a(n,:)-b(n,:), 'k--',"LineWidth",0.1);
    plot(log10(q), a(n,:)+b(n,:), 'k--',"LineWidth",0.1);
    [~,p] = max(a(n,:));
    p = plot(log10(q(p)), a(n,1,p), 'mo');
    set(p, 'linewidth', 2)
    hold off
    title(['Discrimination - Neuron ' num2str(n)])
    xlabel('log_1_0(q)')
    ylabel('Mean MI(q)')
end

    %% Functions
    function [fitresult, gof] = Fitting_curve(x1, y1, x2, y2, x3, y3, x4, y4)

    fitresult = cell( 4, 1 );
    gof = struct( 'sse', cell( 4, 1 ), ...
        'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

    [xData, yData] = prepareCurveData( x1, y1 );

    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.StartPoint = [0.709754432349746 0.910192467553578 0.978691004156862];

 
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );


    figure( 'Name', 'Fit 1' );
    h = plot( fitresult{1}, xData, yData );
    title('Fit for Neuron 1')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    grid on

    [xData, yData] = prepareCurveData( x2, y2 );

    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [0.515433736542118 0.72193376784407 0.955510096008974];

    [fitresult{2}, gof(2)] = fit( xData, yData, ft, opts );

    figure( 'Name', 'Fit 2' );
    h = plot( fitresult{2}, xData, yData );
    title('Fit for Neuron 2')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    grid on

    [xData, yData] = prepareCurveData( x3, y3 );

    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [0.355142740203084 0.49975884517448 0.624609065987624];

    [fitresult{3}, gof(3)] = fit( xData, yData, ft, opts );

    figure( 'Name', 'Fit 3' );
    h = plot( fitresult{3}, xData, yData );
    title('Fit for Neuron 3')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    grid on

    [xData, yData] = prepareCurveData( x4, y4 );

    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [0.45087821170349 0.723648199208609 0.253875154342024];

    [fitresult{4}, gof(4)] = fit( xData, yData, ft, opts );

    figure( 'Name', 'Fit 4' );
    h = plot( fitresult{4}, xData, yData );
    title('Fit for Neuron 4')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    grid on
    end

 function mean_dist = meandist(i,rep,spike_segments,q)
mean_dist = zeros(1,8);
for i1 = 1:8
    for rep1 = 1:50
        if (rep1 == rep && i1 == i)
            continue
        end
        mean_dist(i1) = mean_dist(i1) + VPSDM(spike_segments{i,rep},spike_segments{i1,rep1},q);
    end
end
mean_dist = mean_dist/50;
end

function ci90 = c(ci,MI)
MIbar = mean(MI,2);
MIstd = std(abs(MI),0,2);
alpha = 1 - ci;
T_multiplier = tinv(1-alpha/2, 99);
ci90 = T_multiplier*MIstd/sqrt(99);
end

function d=VPSDM(tli,tlj,q)
nspi=length(tli);
nspj=length(tlj);

if q==0
   d=abs(nspi-nspj);
   return
elseif q==Inf
   d=nspi+nspj;
   return
end

scr=zeros(nspi+1,nspj+1);
scr(:,1)=(0:nspi)';
scr(1,:)=(0:nspj);
if(nspi && nspj)
   for i=2:nspi+1
      for j=2:nspj+1
         scr(i,j)=min([scr(i-1,j)+1 scr(i,j-1)+1 scr(i-1,j-1)+q*abs(tli(i-1)-tlj(j-1))]);
      end
   end
end
d=scr(nspi+1,nspj+1);
end

function MI = MutualInfo(confusion)
MI=0;
for i=1:size(confusion,1)
    for j=1:size(confusion,2)
        if(confusion(i,j)~=0)
            MI=MI+confusion(i,j)/size(confusion,1)*log2(confusion(i,j)/sum(confusion(:,j)));
        end
    end
end
end
