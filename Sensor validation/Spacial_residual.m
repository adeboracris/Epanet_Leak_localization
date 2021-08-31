%UNTITLED Summary 
%   Detailed explanation goes here
clear; close all; clc;
%preparing dataset
% load
load('data_sensor_validation.mat')
elec_sensors{1}= [9,18,23,29]; %selected sensors

%model calibration
n=1;
for j=1:length(elec_sensors)
    for i=1:length(elec_sensors{1,j})
        pol_ref(i,:,j) = polyfit(flow_n.^2,pressure_no_leak(elec_sensors{1,j}(i),:),n);
        p_ref(i,:)=polyval(pol_ref(i,:),flow_n.^2);
    end
end


for n = 1:1
    for time = 1:length(flow_leak)
        for i=1:length(elec_sensors{1,1})
            head_ref(i)= polyval(pol_ref(i,:,1),flow_n_2(time)^2);
        end
        %--Probability residual
        head=pressure_no_leak_2(elec_sensors{1,1},time)';
        residual(time,:) = head_ref-head;
    end
end
for i=1:336-23
   residual_filter(i,:)= mean(residual(i:i+23,:)) ;
end

residual_spacial= [];lable_residual_spacial=[];
for i=1:4
    for j=1:4
        if i~=j && j>i
            lable_residual_spacial= [lable_residual_spacial  sprintf('r_{S%d,S%d};',i,j)];
            residual_spacial=[residual_spacial  residual_filter(:,i)-residual_filter(:,j)];
        end
    end
end
lable_residual_spacial=split(lable_residual_spacial,';');
th_min=1.1*[min(residual_filter)  min(residual_spacial)];
th_max=1.1*[max(residual_filter)  max(residual_spacial)];

 xti=0:24:313;
%% Creating a fault in the 1st sensor
t_fault=100; % time instant that will start the fault
pressure_fault= pressure_no_leak_2;
pressure_fault(elec_sensors{1}(1),t_fault:end)= pressure_no_leak_2(elec_sensors{1}(1),t_fault:end)+0.1;

for n = 1:1
    for time = 1:length(flow_leak)
        for i=1:length(elec_sensors{1,1})
            head_ref(i)= polyval(pol_ref(i,:,1),flow_n_2(time)^2);
        end
        %-------------Probability residual
        head=pressure_fault(elec_sensors{1,1},time)';
        residual_f(time,:) = head_ref-head;
    end
end

for i=1:336-23
   residual_filter_f(i,:)= mean(residual_f(i:i+23,:)) ;
end
residual_spacial_f= [];
for i=1:4
    for j=1:4
        if i~=j && j>i
            residual_spacial_f=[residual_spacial_f  residual_filter_f(:,i)-residual_filter_f(:,j)];
        end
    end
end
residual_spacial_f=[zeros(24,6) ;residual_spacial_f];
residual_filter_f=[zeros(24,4) ;residual_filter_f];

% plot
figure(2)
figure(1)
leak_start=(ones(1,2)*t_fault);
for i=1:10   
    if i<5
        figure(1)
%                 subplot(2,2,i)
        pos=[0.08 0.58 0.08 0.58 ;
            0.65 0.65 0.10 0.1 ]';
        subplot('Position',[pos(i,1) pos(i,2) 0.40 .30])
        plot(residual_filter_f(:,i),'LineWidth',[1]);hold on
        plot(leak_start,[min(residual_filter_f(:,i))-0.1 max(residual_filter_f(:,i))+.1],'g','LineWidth',[1])
        plot(ones(1,length(residual_filter_f))*th_min(i),'r','LineWidth',[1]);plot(ones(1,length(residual_filter_f))*th_max(i),'r','LineWidth',[1]);
        str = sprintf('Sensor %d', i); title(str)
        
        xticks(xti)
        xticklabels({'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '13'  })
        xlabel('[Days]')
        ylabel('[mcw]')
        axis([0 length(residual_filter_f(:,i)),min(residual_filter_f(:,i))-0.1 max(residual_filter_f(:,i))+.1 ])
        if i ==1
            legend('Filtered residual','Start fault','Residual bounds')
        end
    else
        figure(2)
        pos=[0.061 0.52 0.061 0.52 0.061 0.525; 
             0.70 0.70 0.38 0.38 0.080 0.080]';
        subplot('Position',[pos(i-4,1) pos(i-4,2) 0.40 .20])
%         subplot(3,2,i-4)
        plot(residual_spacial_f(:,i-4),'LineWidth',[1]) ;hold on
        plot(leak_start,[min(residual_spacial_f(:,i-4))-0.1 max(residual_spacial_f(:,i-4))+.1],'g','LineWidth',[1])
        plot(ones(1,length(residual_spacial_f))*th_min(i),'r','LineWidth',[1]); plot(ones(1,length(residual_spacial_f))*th_max(i),'r','LineWidth',[1]);
        
        str = sprintf('S%s\n', lable_residual_spacial{i-4}); title(str,'VerticalAlignment', 'middle')
        xticks(xti)
        xticklabels({'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '13'  })
        xlabel('[Days]')
        ylabel('[mcw]')
        axis([0 length(residual_spacial_f(:,i-4)),min(residual_spacial_f(:,i-4))-0.1 max(residual_spacial_f(:,i-4))+.1 ])
        if i ==6
            legend('Spacial residual','Start fault','Residual bounds')
        end
    end
end

return

% this part of the code is to generate a plot without any fault, it not in
%the article, but is good to analyze 

figure(3)
figure(4)
leak_start=(ones(1,2)*t_fault);
for i=1:10   
    if i<5
        figure(3)
%                 subplot(2,2,i)
        pos=[0.08 0.58 0.08 0.58 ;
            0.65 0.65 0.10 0.1 ]';
        subplot('Position',[pos(i,1) pos(i,2) 0.40 .30])
        plot(residual_filter(:,i),'LineWidth',[1]);hold on
%         plot(leak_start,[min(residual_filter(:,i))-0.1 max(residual_filter(:,i))+.1],'g','LineWidth',[1])
        plot(ones(1,length(residual_filter))*th_min(i),'r','LineWidth',[1]);plot(ones(1,length(residual_filter))*th_max(i),'r','LineWidth',[1]);
        str = sprintf('Sensor %d', i); title(str)
        
        xticks(xti)
        xticklabels({'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '13'  })
        xlabel('[Days]')
        ylabel('[mcw]')
        axis([0 length(residual_filter(:,i)),min(residual_filter(:,i))-0.1 max(residual_filter(:,i))+.1 ])
        if i ==1
            legend('Filtered residual','Residual bounds')
        end
    else
        figure(4)
        pos=[0.061 0.54 0.061 0.54 0.061 0.525; 
             0.70 0.70 0.38 0.38 0.080 0.080]';
        subplot('Position',[pos(i-4,1) pos(i-4,2) 0.40 .20])
%         subplot(3,2,i-4)
        plot(residual_spacial(:,i-4),'LineWidth',[1]) ;hold on
        plot(leak_start,[min(residual_spacial(:,i-4))-0.1 max(residual_spacial(:,i-4))+.1],'g','LineWidth',[1])
        plot(ones(1,length(residual_spacial))*th_min(i),'r','LineWidth',[1]); plot(ones(1,length(residual_spacial))*th_max(i),'r','LineWidth',[1]);
        
        str = sprintf('S%s\n', lable_residual_spacial{i-4}); title(str,'VerticalAlignment', 'middle')
        xticks(xti)
        xticklabels({'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '13'  })
        xlabel('[Days]')
        ylabel('[mcw]')
        axis([0 length(residual_spacial(:,i-4)),min(residual_spacial(:,i-4))-0.1 max(residual_spacial(:,i-4))+.1 ])
        if i ==6
            legend('Spacial residual','Start fault','Residual bounds')
        end
    end
end

