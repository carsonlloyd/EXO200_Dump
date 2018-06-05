
%% plot red points (all + cut by detector half)
clear;clc;
load('all-for-alpha-cuts-pos2-noVOLcut-driftTime.mat')
pairs = {};
t0 = time(1);
for i = 2:length(time)
    i
    t = time(i);
    dt = t-t0;%abs(t - t0);
    dz = position(i,3)-position(i-1,3);
    v = dz/dt;
    dU = abs(U(i) - U(i-1));
    R = norm([position(i,1)-position(i-1,1), position(i,2)-position(i-1,2)]);%, position(i,3)-position(i-1,3)]);
    if (dt <= 180) && (dU <= 9.0) && (dt >= 0.0) && (R <= 30.0) %&& (detector_half(i) == detector_half(i-1))
        pairs = [pairs; {i-1,i,dt,dz}]; %#ok<AGROW>
    end
    t0 = t;
end
%plot(pairs{:,8},pairs{:,5})
scatter([pairs{:,4}],[pairs{:,3}],'r')

%% overlay red points (with final cut)
hold on;
load('all-for-alpha-cuts-pos2-noVOLcut-driftTime.mat')
pairs = {};
t0 = time(1);
for i = 2:length(time)
    i
    t = time(i);
    dt = t-t0;%abs(t - t0);
    dz = position(i,3)-position(i-1,3);
    v = dz/dt;
    dU = abs(U(i) - U(i-1));
    R = norm([position(i,1)-position(i-1,1), position(i,2)-position(i-1,2)]);%, position(i,3)-position(i-1,3)]);
    %if (dt <= 180) && (dU <= 9.0) && (dt >= 0.0) && (R <= 30.0) && ((detector_half(i) == detector_half(i-1)) || (detector_half(i) + detector_half(i-1) >= 2))
    if (dt <= 180) && (dU <= 9.0) && (dt >= 0.0) && (R <= 30.0) && (detector_half(i) == detector_half(i-1))
        pairs = [pairs; {i-1,i,dt,dz}]; %#ok<AGROW>
    end
    t0 = t;
end
scatter([pairs{:,4}],[pairs{:,3}],'b')

%% only points that were cut, z1 vs z2
clear; clc; close all;
load('all-for-alpha-cuts-pos2-noVOLcut-driftTime.mat')
first = {};
t0 = time(1);
for i = 2:length(time)
    i
    t = time(i);
    dt = t-t0;%abs(t - t0);
    dz = position(i,3)-position(i-1,3);
    v = dz/dt;
    dU = abs(U(i) - U(i-1));
    R = norm([position(i,1)-position(i-1,1), position(i,2)-position(i-1,2)]);%, position(i,3)-position(i-1,3)]);
    %if (dt <= 180) && (dU <= 9.0) && (dt >= 0.0) && (R <= 30.0) && ((detector_half(i) == detector_half(i-1)) || (detector_half(i) + detector_half(i-1) >= 2))
    if (dt <= 180) && (dU <= 9.0) && (dt >= 0.0) && (R <= 30.0) && (detector_half(i) ~= detector_half(i-1))
        first = [first; {i-1,i,position(i,3),position(i-1,3)}]; %#ok<AGROW>
    end
    t0 = t;
end
%scatter3([first{:,3}],[first{:,4}],[first{:,5}],'b')
%hold on;
%scatter3([second{:,3}],[second{:,4}],[second{:,5}],'r')

%plot3([first{1,3},second{1,3}],[first{1,4},second{1,4}],[first{1,5},second{1,5}],'g')
%plot3([first{2,3},second{2,3}],[first{2,4},second{2,4}],[first{2,5},second{2,5}],'g')
%plot3([first{3,3},second{3,3}],[first{3,4},second{3,4}],[first{3,5},second{3,5}],'g')

scatter([first{:,3}],[first{:,4}])

%% histogram of velocities
clear; clc; close all;
load('all-for-alpha-cuts-pos2-noVOLcut-driftTime.mat')
first = {};
t0 = time(1);
for i = 2:length(time)
    i
    t = time(i);
    dt = t-t0;%abs(t - t0);
    dz = position(i,3)-position(i-1,3);
    v = dz/dt;
    dU = abs(U(i) - U(i-1));
    R = norm([position(i,1)-position(i-1,1), position(i,2)-position(i-1,2)]);%, position(i,3)-position(i-1,3)]);
    %if (dt <= 180) && (dU <= 9.0) && (dt >= 0.0) && (R <= 30.0) && ((detector_half(i) == detector_half(i-1)) || (detector_half(i) + detector_half(i-1) >= 2))
    if (dt <= 180) && (dU <= 9.0) && (dt >= 0.0) && (R <= 30.0) && (detector_half(i) == detector_half(i-1))
        first = [first; {i-1,i,position(i,3),position(i-1,3),v}]; %#ok<AGROW>
    end
    t0 = t;
end
hplot = histfit([first{:,5}],50,'kernel');
axes('Position',[.7 .7 .2 .2])
points = [first{:,5}];
points =  points(points > 0.5 & points < 1.5);
hplot = histfit(points,6,'normal');
[mu,sigma] = normfit(points);
%% drift time vs drift velocity
clear; clc; close all;
load('all-for-alpha-cuts-pos2-noVOLcut-driftTime.mat')
first = {};
second = {};
t0 = time(1);
for i = 2:length(time)
    t = time(i);
    dt = t-t0;%abs(t - t0);
    dz = position(i,3)-position(i-1,3);
    v = dz/dt;
    dU = abs(U(i) - U(i-1));
    R = norm([position(i,1)-position(i-1,1), position(i,2)-position(i-1,2)]);%, position(i,3)-position(i-1,3)]);
    if (dt <= 180) && (dU <= 9.0) && (dt >= 0.0) && (R <= 30.0) && (detector_half(i) == detector_half(i-1))
        if (detector_half(i) == 0) %TPC1
            first = [first; {i-1,i,dt,v,detector_half(i)}]; %#ok<AGROW>
        elseif (detector_half(i) == 1) %TPC2
            second = [second; {i-1,i,dt,v,detector_half(i)}]; %#ok<AGROW>
        end
    end
    t0 = t;
end
hold on;
h1=scatter([first{:,4}],[first{:,3}],'b','DisplayName','TPC1');
h1=scatter([second{:,4}],[second{:,3}],'o','DisplayName','TPC2');
title('Drift Time vs. Drift Velocity');
xlabel('v (mm/s)');
ylabel('\Deltat (s)');
%set(h1,'FontSize',16);
legend(h1,'show');
