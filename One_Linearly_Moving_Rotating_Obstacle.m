%% 1 Moving Obstacle - Linear & Angular Velocity
clear
close all
clc

syms x y w

Rotate = @(vec,theta) [cos(theta),-sin(theta);sin(theta),cos(theta)]*vec;

cushion = 1;
margin = 1;

x0 = -4;
y0 = -4;

xc = 3;
yc = 3;
a = 3;
b = 0.2;

xr = xc+0;
yr = yc-0.1;

xa = 4;
ya = 4;

B = ((x-xc)/a)^2+((y-yc)/b)^2;  

n = [diff(B,x);diff(B,y)];
r = [x-xr ; y-yr];
e = [-diff(B,y);diff(B,x)];
E = [r,e];

D = [1-1/(B^(1/cushion)) , 0;
         0         , 1+1/(B^(1/cushion))];
     
M = E*D*inv(E);
M = subs(M,[x,y],[x/margin,y/margin]);

f = [-(x-xa);-(y-ya)];

vx_obs = -1;
vy_obs = -1;
v_obs = [vx_obs;vy_obs];

w_obs = 0.5;
fw = [-w*(y-yc);w*(x-xc)];

tstep = 0.01;
steps = 500;

pos = [x0,y0];
pos_wrt_obs = pos;
pos_dot_sim = [];
for i = 1:steps
    pos_current = pos(end,:);
    pos_wrt_obs_current = pos_wrt_obs(end,:);
    pos_dot = double(subs(f,[x,y],pos_current));
    pos_dot_wrt_obs = pos_dot - v_obs;
    pos_dot_wrt_obs = Rotate(pos_dot_wrt_obs,-w_obs*tstep*(i-1));
    vw_obs = double(subs(fw,[x,y,w],[pos_wrt_obs_current,w_obs]));
    pos_dot_wrt_obs = pos_dot_wrt_obs - vw_obs;
    pos_dot_wrt_obs = double(subs(M,[x,y],pos_wrt_obs_current))*pos_dot_wrt_obs;
    pos_wrt_obs_next = pos_wrt_obs_current' + pos_dot_wrt_obs*tstep;
    pos_dot_wrt_obs = pos_dot_wrt_obs + vw_obs;
    pos_dot_wrt_obs = Rotate(pos_dot_wrt_obs,w_obs*tstep*(i-1));
    pos_dot = pos_dot_wrt_obs + v_obs;
    pos_next = pos_current' + pos_dot*tstep;
    
    pos = [pos;pos_next'];
    pos_wrt_obs = [pos_wrt_obs; pos_wrt_obs_next'];
    pos_dot_sim = [pos_dot_sim;pos_dot'];
end

%% 1 Moving Obstacle - Linear & Angular Velocity - With Respect To Obstacle - Plot
close all
clc

theta = [0:0.01:2*pi];
r_theta  = [];
for i = 1:length(theta)
    r_theta = [r_theta,(a*b)/sqrt((b*cos(theta(i)))^2 + (a*sin(theta(i)))^2)];
end
x_ellipse = r_theta .* cos(theta);
y_ellipse = r_theta .* sin(theta);
ellipse = [x_ellipse',y_ellipse'];

plot(pos_wrt_obs(:,1),pos_wrt_obs(:,2))
hold on;
plot(x_ellipse+xc,y_ellipse+yc);
hold off;

%% 1 Moving Obstacle - Linear & Angular Velocity - Plot
figure

plot(pos(:,1),pos(:,2))
hold on;
plot(x_ellipse+xc,y_ellipse+yc);
hold off;

%% 1 Moving Obstacle - Linear & Angular Velocity - Movie
close all
clc


for k=1:length(pos)-1
    plot(pos(k,1),pos(k,2),'x')
    hold on
    plot(xa,ya,'p')
    quiver(pos(k,1),pos(k,2),pos_dot_sim(k,1),pos_dot_sim(k,2),0,'r');
    ellipse_rotate = Rotate(ellipse',w_obs*tstep*(k-1))';
    plot(ellipse_rotate(:,1)+xc+vx_obs*tstep*(k-1),ellipse_rotate(:,2)+yc+vy_obs*tstep*(k-1));
    axis([-5,5,-5,5])
    hold off
    pause(tstep)
end
