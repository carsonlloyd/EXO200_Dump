
%%
clear;clc;close all;
%load('/Users/carsonlloyd/Documents/Drexel/Senior/Spring/SeniorResearch/GoldenMaskedData/Apr17/golden-phase2/surfaceStudies/1CC/30train/30train_70test_30train_ttip.mat')
load('/Users/carsonlloyd/Documents/Drexel/Senior/Spring/SeniorResearch/GoldenMaskedData/Apr17/golden-phase2/surfaceStudies/0CC/30train/30train_70test_30train_ttip.mat')

fprintf('Done loading .pkl file...')
xs = position(:,1);
ys = position(:,2);
zs = position(:,3);

abs_zs = abs(position(:,3));

Rs = zeros(1, length(zs));
for i = 1:length(zs)
    R = norm([position(i,1), position(i,2)]);
    Rs(i) = R;
    fprintf('%s%% done\n',num2str(100*i/length(zs)))
end

%%
clc;
fprintf('R > 172mm: %s\n',num2str(sum(Rs > 172.00)))
fprintf('R < 172mm: %s\n',num2str(sum(Rs < 172.00)))
fprintf('Outside fiducial z:??? %s\n',num2str( length( abs_zs( abs_zs>172.0 | abs_zs<20.0))))
fprintf('20 < |z| < 172mm:??? %s\n',num2str( length( abs_zs( abs_zs>20.0 & abs_zs<172.0 ))))
madeit=0;
for i = 1:length(abs_zs)
    if (abs_zs(i) > 20.0) && (abs_zs(i) < 172.0) && (Rs(i) > 172.0)
        madeit = madeit+1;
    end
end
fprintf('20 < |z| < 172mm AND R>172mm: %s\n',num2str(madeit))
clear madeit;

%% ALL
close all;
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter3(xs,ys,zs,[],zs);
xlabel({'x (mm)'});
zlabel('z (mm)	');
title('title');
ylabel('y (mm)');
view(axes1,[-37.5 30]);
grid(axes1,'on');
set(axes1,'FontSize',16);
legend(axes1,'show');

%% OUTSIDE FIDUCIAL Z
close all;
madeit = zeros(length(abs_zs),3);
for i = 1:length(abs_zs)
    if (abs_zs(i)<20.0) || (abs_zs(i)>172.0)
        madeit(i,:) = [xs(i),ys(i),zs(i)];
        fprintf('%s%% done\n',num2str(100*i/length(abs_zs)))
    end
end
new = madeit;
madeit(~any(new,2),:) = [];

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter3(madeit(:,1),madeit(:,2),madeit(:,3),[],madeit(:,3));
xlabel('x (mm)');
zlabel('z (mm)');
title('title');
ylabel('y (mm)');
view(axes1,[-37.5 30]);
grid(axes1,'on');
set(axes1,'FontSize',16);
legend(axes1,'show');

clear madeit; clear new;

%% OUTSIDE FIDUCIAL Z WITH R CUT ( > 172 MM)
close all;
madeit = zeros(length(abs_zs),3);
for i = 1:length(abs_zs)
    if ((abs_zs(i)<20.0) || (abs_zs(i)>172.0)) && (Rs(i) > 172.0)
        madeit(i,:) = [xs(i),ys(i),zs(i)];
        fprintf('%s%% done\n',num2str(100*i/length(abs_zs)))
    end
end
new = madeit;
madeit(~any(new,2),:) = [];

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter3(madeit(:,1),madeit(:,2),madeit(:,3),[],madeit(:,3));
xlabel('x (mm)');
zlabel('z (mm)');
title('title');
ylabel('y (mm)');
view(axes1,[-37.5 30]);
grid(axes1,'on');
set(axes1,'FontSize',16);
legend(axes1,'show');

clear madeit; clear new;

%% INSIDE FIDUCIAL Z
close all;
madeit = zeros(length(abs_zs),3);
for i = 1:length(abs_zs)
    if (abs_zs(i)>20.0) && (abs_zs(i)<172.0)
        madeit(i,:) = [xs(i),ys(i),zs(i)];
        fprintf('%s%% done\n',num2str(100*i/length(abs_zs)))
    end
end

new = madeit;
madeit(~any(new,2),:) = [];

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter3(madeit(:,1),madeit(:,2),madeit(:,3),[],madeit(:,3));
xlabel('x (mm)');
zlabel('z (mm)');
title('title');
ylabel('y (mm)');
view(axes1,[-37.5 30]);
grid(axes1,'on');
set(axes1,'FontSize',16);
legend(axes1,'show');

%clear madeit;
clear new;

%%

histogram(madeit(:,3),100)
xlabel('z (mm)');
ylabel('0CC events');
title('title');

clear madeit; clear new;

%% INSIDE FIDUCIAL Z WITH R CUT
close all;
madeit = zeros(length(abs_zs),3);
for i = 1:length(abs_zs)
    if (abs_zs(i)>20.0) && (abs_zs(i)<172.0) && (Rs(i) > 172.0)
        madeit(i,:) = [xs(i),ys(i),zs(i)];
        fprintf('%s%% done\n',num2str(100*i/length(abs_zs)))
    end
end

new = madeit;
madeit(~any(new,2),:) = [];

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter3(madeit(:,1),madeit(:,2),madeit(:,3),[],madeit(:,3));
xlabel('x (mm)');
zlabel('z (mm)');
title('title');
ylabel('y (mm)');
view(axes1,[-37.5 30]);
grid(axes1,'on');
set(axes1,'FontSize',16);
legend(axes1,'show');

clear madeit; clear new;





