%% Load Files
load('1CentroidArray2h.mat')

load('1TimeArray2h_u.mat')





%%
%% Load Files
load('1CentroidArray24h.mat')
load('1TimeArray24h_u.mat')




%% Test minimal centroid stuff
scale=10^6;
start=0*scale+1;
numdatapoints=min(16*scale, length(Centroidarray));

centroidtrunc=Centroidarray(start:numdatapoints, :, 1, 1);
timetrunc=timearray(start:numdatapoints, 1, 1, 1);
% plot(centroidtrunc(:, 1, 1, 1), centroidtrunc(:,2,1,1));

[~,b]=AngleArrays(centroidtrunc, timetrunc, false);
%% Check histogram
b
subplot(2,1,1)
histogram(real((b.turning)), 200)
subplot(2,1,2)
plot(b.direction)
%% Plot all data

load('1CentroidArray2h.mat')

load('1TimeArray2h_u.mat')




f=27
f=19
[a1,b1]=AngleArrays(Centroidarray(:, :, f, 1), timearray(:, 1, 1, 1), false);
figure(2); clf; histogram(log(b1.turning))

% load('1CentroidArray24h.mat')
% load('1TimeArray24h_u.mat')
% 
% [a1,b1]=AngleArrays(Centroidarray(:, :, f, 1), timearray(:, 1, 1, 1), false);
% hold on;
% histogram(log(b1.speed))

%% Check all
figure(2)
clf
subplot(8,1,1)
plot(b1.inx)
title('x')
subplot(8,1,2)
plot(b1.iny)
title('y')

subplot(8,1,3)
plot(b1.theta)
title('theta')

subplot(8,1,4)
plot(b1.r)
title('r')

subplot(8,1,5)
plot(b1.direction)
title('direction')

subplot(8,1,6)
plot(log(b1.speed))
title('log speed')

subplot(8,1,7)
plot(b1.turning)
title('turning')

subplot(8,1,8)
plot(b1.angle)
title('angle')




%%
s_i=0.20*10^5
% s_i=1
e_i=1*10^6

% s_i=4800
% e_i=5050;
e_i=min(length(Centroidarray(:,1,f,1)), e_i)
figure(3)
clf
% plot(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1))
q=quiver(Centroidarray(s_i:e_i-1,1,f,1), Centroidarray(s_i:e_i-1,2,f,1), diff(Centroidarray(s_i:e_i,1,f,1)), diff(Centroidarray(s_i:e_i,2,f,1)), 'off', LineWidth=.1)
q.Color=[0.8 0.8 0.8]
hold on 
% colormap("turbo
colormap(brewermap([],"RdBu")) % For diverging
% colormap("parula")
% maxval=50/max(cumsum(timearray(s_i:e_i, 1, 1, 1)));
% scatter(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1),1,b1.turning(s_i:e_i), "filled")
scatter(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1),20,(b1.direction(s_i:e_i)), "filled")
% 
hold on
colorbar
clim([-pi, pi])
xlim([0,95])
ylim([0,105])
% clim([5-(5-log(4))*2, 5])
% title('Speed')
% clabel('Speed')

%%
figure(3)
clf
% plot(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1))
q=quiver(Centroidarray(s_i:e_i-1,1,f,1), Centroidarray(s_i:e_i-1,2,f,1), diff(Centroidarray(s_i:e_i,1,f,1)), diff(Centroidarray(s_i:e_i,2,f,1)), 'off', LineWidth=.1)
q.Color=[0.8 0.8 0.8]
hold on 
% colormap("turbo
colormap(brewermap([],"RdBu")) % For diverging
% colormap("parula")
% maxval=50/max(cumsum(timearray(s_i:e_i, 1, 1, 1)));
scatter(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1),1,log(b1.speed(s_i:e_i)), "filled")
% scatter(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1),20,(b1.direction(s_i:e_i)), "filled")
% 
hold on
colorbar
% clim([-pi, pi])
xlim([0,95])
ylim([0,105])
clim([5-(5-log(4))*2, 5])
% title('Speed')
% clabel('Speed')

%%
figure(3)
clf
% plot(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1))
q=quiver(Centroidarray(s_i:e_i-1,1,f,1), Centroidarray(s_i:e_i-1,2,f,1), diff(Centroidarray(s_i:e_i,1,f,1)), diff(Centroidarray(s_i:e_i,2,f,1)), 'off', LineWidth=.1)
q.Color=[0.8 0.8 0.8]
hold on 
% colormap("turbo
colormap(brewermap([],"RdBu")) % For diverging
% colormap("parula")
% maxval=50/max(cumsum(timearray(s_i:e_i, 1, 1, 1)));
scatter(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1),1,(b1.angle(s_i:e_i)), "filled")
% scatter(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1),20,(b1.direction(s_i:e_i)), "filled")
% 
hold on
colorbar
clim([-1, 1])
xlim([0,95])
ylim([0,105])
% clim([5-(5-log(4))*2, 5])
% title('Speed')
% clabel('Speed')

%%
figure(3)
clf
% plot(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1))
q=quiver(Centroidarray(s_i:e_i-1,1,f,1), Centroidarray(s_i:e_i-1,2,f,1), diff(Centroidarray(s_i:e_i,1,f,1)), diff(Centroidarray(s_i:e_i,2,f,1)), 'off', LineWidth=.1);
q.Color=[0.8 0.8 0.8];
hold on 
% colormap("turbo
% colormap(brewermap([],"RdBu")) % For diverging
colormap("parula")
% maxval=50/max(cumsum(timearray(s_i:e_i, 1, 1, 1)));
scatter(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1),1,b1.r(s_i:e_i), "filled")
% scatter(Centroidarray(s_i:e_i,1,f,1), Centroidarray(s_i:e_i,2,f,1),20,(b1.direction(s_i:e_i)), "filled")
% 
hold on
colorbar
clim([0, .5])
xlim([0,95])
ylim([0,105])
% clim([5-(5-log(4))*2, 5])
% title('Speed')
% clabel('Speed')

%% Plot direction over time
figure(5)
disp('test')
% s_i=4800
% e_i=5050;
t=timearray(s_i:e_i, 1, 1, 1);
% disp()
disp(length(t))
t=t-t(1);
disp(t(2)-t(1))
plot(t,a1(s_i:e_i), LineWidth=1)
xlabel("Time(s)")

% whos
disp([s_e, s_i, s_e-s_i])
disp(s_e-s_i)
%% Plot direction over all time
figure()
t=timearray(:, 1, 1, 1);
t=t-t(1);
disp(t(2)-t(1))
plot(t,a1(:), LineWidth=1)
xlabel("Time(s)")


%% Look for plot with least distortion
figure(2)
clf
hold on
subplot(6,7,1)
for i=1:42
    subplot(6,7,i)
    plot(Centroidarray(:, 1, i, 1), Centroidarray(:, 2, i, 1))
end

%% Show directions
figure(3)
i=23
x=Centroidarray(1:end-1, 1, i, 1);
y=Centroidarray(1:end-1, 2, i, 1);
u=diff(Centroidarray(1:end, 1, i, 1));
v=diff(Centroidarray(1:end, 2, i, 1));
q=quiver(Centroidarray(1:end-1, 1, i, 1), Centroidarray(1:end-1, 2, i, 1), diff(Centroidarray(1:end, 1, i, 1)), diff(Centroidarray(1:end, 2, i, 1)),'off')
% quiver(x(1:25),y(1:25),u(1:25),v(1:25), "off")
q


%% Count transitions
binarized_behavior=round(a1(s_i:e_i)/2+.5);
figure(4)
clf
% plot(t, binarized_behavior)
% hold on
binarized_behavior_nonan=binarized_behavior(~isnan(binarized_behavior))
t_trunc=t(~isnan(binarized_behavior));
scatter(t_trunc, binarized_behavior_nonan)
figure(5)
clf
plot(t_trunc(1:end-1),abs(diff(binarized_behavior_nonan)))
changes_in_circling=abs(diff(binarized_behavior_nonan));
scatter(t_trunc(changes_in_circling>0),changes_in_circling(changes_in_circling>0))
% histogram(diff(t(binarized_behavior_nonan>0)), 10)