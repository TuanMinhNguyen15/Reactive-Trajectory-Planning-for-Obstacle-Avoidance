%% 1 Moving Obstacle - Linear Velocity
clear
close all
clc

syms x y

cushion = 1;
margin = 1;

x0 = 3;
y0 = 3;

xc = -2;
yc = -2;
a = 1;
b = 1;
xr = xc+0;
yr = yc+0.1;

xa = 0;
ya = 0;

vx_obs = 1;
vy_obs = 1;
v_obs = [vx_obs;vy_obs];

B = (((x-xc)/a)^2) +(((y-yc)/b)^2); 

n = [diff(B,x);diff(B,y)];
r = [x-xr ; y-yr];
e = [-diff(B,y);diff(B,x)];
E = [r,e];

D = [1-1/(B^(1/cushion)) , 0;
         0         , 1+1/(B^(1/cushion))];
     
M = E*D*inv(E);
M = subs(M,[x,y],[x/margin,y/margin]);

f = [-(x-xa);-(y-ya)];
% f = 0;   % TUAN
% fm_wrt_obs = M*(f-vo);

tstep = 0.01;
steps = 1000;

pos = [x0,y0];
pos_wrt_obs = pos;
pos_dot_sim = [];
for i = 1:steps
    pos_current = pos(end,:);
    pos_wrt_obs_current = pos_wrt_obs(end,:);
    pos_dot = double(subs(f,[x,y],pos_current));
    pos_dot_wrt_obs = pos_dot - v_obs;
    pos_dot_wrt_obs = double(subs(M,[x,y],pos_wrt_obs_current))*pos_dot_wrt_obs;
    pos_dot = pos_dot_wrt_obs + v_obs;
    
    pos_wrt_obs_next = pos_wrt_obs_current' + pos_dot_wrt_obs*tstep;
    pos_next = pos_current' + pos_dot*tstep;  % Assume v_obs is constant
    pos = [pos;pos_next'];
    pos_wrt_obs = [pos_wrt_obs; pos_wrt_obs_next'];
    pos_dot_sim = [pos_dot_sim;pos_dot'];
end
%% 1 Moving Obstacle - Linear Velocity - With Respect To Obstacle - Plot
close all
clc

theta = [0:0.01:2*pi];
r_theta  = [];
for i = 1:length(theta)
    r_theta = [r_theta,(a*b)/sqrt((b*cos(theta(i)))^2 + (a*sin(theta(i)))^2)];
end

x_ellipse = r_theta .* cos(theta);
y_ellipse = r_theta .* sin(theta);

plot(pos_wrt_obs(:,1),pos_wrt_obs(:,2))
hold on;
plot(x_ellipse+xc,y_ellipse+yc);
hold off;

%% 1 Moving Obstacle - Linear Velocity - Plot
figure

plot(pos(:,1),pos(:,2))
hold on;
plot(x_ellipse+xc,y_ellipse+yc);
hold off;

%% 1 Moving Obstacle - Linear Velocity - Movie
close all
clc

for k=1:length(pos)
    plot(pos(k,1),pos(k,2),'x')
    hold on
    quiver(pos(k,1),pos(k,2),pos_dot_sim(k,1),pos_dot_sim(k,2),0,'r');
    plot(x_ellipse+xc+vx_obs*tstep*(k-1),y_ellipse+yc+vy_obs*tstep*(k-1),'-');
    axis([-5,5,-5,5])
    hold off
    pause(tstep)
end
