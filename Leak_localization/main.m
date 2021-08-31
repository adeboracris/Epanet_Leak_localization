clear;close all;clc

%load the data to use in the code
addpath('data')
load('HanoiData2_Du_20.mat', 'TrainingHead','TestingHead','S')
load('coordenades.mat', 'node_coordenades','D','P','Resistence_PipeDistance')

PipeDistance=Resistence_PipeDistance;

% selected sensor with the different scenarios
elec_sensors{1}= [12,18,23,29];
elec_sensors{2} = [6,12,17,23,29,21];
elec_sensors{3} = [6,12,15,17,23,21,27,30];
elec_sensors{4} = [6,9,12,15,17,24,21,22,28,29,31];
% elec_sensors{5} = [1:31];


%the training data
TrainingFlow=TrainingHead(:,33);
TrainingPressure=TrainingHead(:,[1:31]);

%Curent data : 4 days
Flow_leak=TestingHead(:,33,:);
Pressure_leak=TestingHead(:,[1:31],:);

%model calibration
n=1;
for j=1:length(elec_sensors)
    for i=1:length(elec_sensors{1,j})
        alpha(i,:,j) = polyfit(TrainingFlow.^2,TrainingPressure(:,elec_sensors{1,j}(i)),n);
    end
end

%% Clustering
% Will be ns(number of sensors) clusters
Tts=2;%days

AccuracyTime = zeros(Tts*24,length(elec_sensors));
for scenario_s = 1:length(elec_sensors)
    %-------------lambdaausibility vector  distance[lambda]
    [lambda,di]=Plausibility(elec_sensors{1,scenario_s});
    %end----------
    
    ConfusionMatrix = zeros(31);
    Pred = zeros(Tts*24,31,31);
    prob=[]; head_ref_interpolation=[];residual=[];
    for n = 1:31
        for time = 1:Tts*24
            for i=1:length(elec_sensors{1,scenario_s})
                head_ref(i)= polyval(alpha(i,:,scenario_s),Flow_leak(time,:,n)^2);
            end
            %-------------Probability residual
            head_leak=Pressure_leak(time,elec_sensors{1,scenario_s},n); % taking only the sensor measurement
            residual = head_ref-head_leak;
            r_min = min(residual(1:length(elec_sensors{1,scenario_s})));%(eq 9)
            residual(1:length(elec_sensors{1,scenario_s})) = residual(1:length(elec_sensors{1,scenario_s}))-r_min; %(eq 9)
            
            prob(time,:,n) = residual(1:length(elec_sensors{1,scenario_s}))/sum(residual(1:length(elec_sensors{1,scenario_s})));
            
            %---------- pressure residual method
            theta_aux=sum((lambda'.*prob(time,:,n))');%
            theta(time,:,n)=theta_aux;
            [~,Class] = max(theta(time,:,n));
            %end ----------
            % confusionMatrix
            ConfusionMatrix(n,(Class)) = ConfusionMatrix(n,(Class))+1;
        end
    end
    
    Accuracy = 0;
    for i = 1:31
        Accuracy = Accuracy+ConfusionMatrix(i,i);
    end
    Accuracy = 100*Accuracy/(31*24*Tts);
    
    ATD = sum(sum(ConfusionMatrix.*D(1:end-1,1:end-1)))/sum(sum(ConfusionMatrix));
    
    ConfusionMatrixTime = zeros(31,31,Tts*24);
    for i = 1:31
        for j = 1:Tts*24
            Prob(1,1:31) = 1;
            for k = 0:Tts*24-1
                if j+k <= Tts*24
                    Prob = Prob.*theta(j+k,:,i);
                else
                    Prob = Prob.*theta(j+k-Tts*24,:,i);
                end
                Prob = Prob/sum(Prob);
                [~,Class] = max(Prob);
                ConfusionMatrixTime(i,Class,k+1) = ConfusionMatrixTime(i,Class,k+1)+1;
            end
        end
    end
    
    for j = 1:Tts*24
        for i = 1:31
            AccuracyTime(j,scenario_s) = AccuracyTime(j,scenario_s)+ConfusionMatrixTime(i,i,j);
        end
        ATDTime(j,scenario_s) = sum(sum(ConfusionMatrixTime(:,:,j).*D(1:end-1,1:end-1)))/sum(sum(ConfusionMatrixTime(:,:,j)));
    end
    AccuracyTime(:,scenario_s) = 100*AccuracyTime(:,scenario_s)/(31*24*Tts);
    
end


%plot
figure
a1 = plot(ATDTime(:,1),'color',[0 0 0],'linestyle','-','linewidth',2);
hold on
a3 = plot(ATDTime(:,2),'r','linestyle',':','linewidth',2);
a5 = plot(ATDTime(:,3),'b','linestyle','--','linewidth',2);
a6 = plot(ATDTime(:,4),'m','linestyle','-.','linewidth',2);
% a7 = lambdaot(ATDTime(:,5),'k','linestyle','--','linewidth',2);
hold off
h = legend([a1 a3 a5 a6 ],'4 sensors','6 sensors','8 sensors','10 sensors','location','northeast');
set(h,'interpreter','latex','fontsize',12);
xlim([1 24]);
xlabel('Time horizon $\mathrm{[h]}$','interpreter','latex','fontsize',14);
ylabel('ATD $\mathrm{[nodes]}$','interpreter','latex','fontsize',14);
grid
title('Evolution of the ATD')
ylim([0 3.5]);

return