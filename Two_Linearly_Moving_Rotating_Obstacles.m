%% 2 Moving Obstacles - Linear & Angular Velocity
clear
close all
clc

syms x y w B1_ B2_

Rotate = @(vec,theta) [cos(theta),-sin(theta);sin(theta),cos(theta)]*vec;

cushion = 1.5;  %1.5
margin = 1;

% Robot's initial position
x0 = -4;
y0 = -4;

% Goal point
xa = 4;
ya = 4;
f = [-(x-xa);-(y-ya)];

% Obstacle 1
xc1 = 3;  %0
yc1 = 3;  %0
a1 = 0.2;   %2
b1 = 3;   %1
xr1 = xc1+0;
yr1 = yc1+0.1;   %+0.1
vx_obs1 = -1;
vy_obs1 = -1;
v_obs1 = [vx_obs1;vy_obs1];
w_obs1 = 0;
fw1 = [-w*(y-yc1);w*(x-xc1)];
B1 = ((x-xc1)/a1)^2+((y-yc1)/b1)^2;  
n1 = [diff(B1,x);diff(B1,y)];
r1 = [x-xr1 ; y-yr1];
e1 = [-diff(B1,y);diff(B1,x)];
E1 = [r1,e1];
D1 = [1-1/(B1^(1/cushion)) , 0;
         0         , 1+1/(B1^(1/cushion))];     
M1 = E1*D1*inv(E1);
M1 = subs(M1,[x,y],[x/margin,y/margin]);

% Obstacle 2
xc2 = 0;  %-5
yc2 = 0;   % 5
a2 = 1;    % 1
b2 = 3;    % 3
xr2 = xc2+0;
yr2 = yc2-0.1;  %-0.1
vx_obs2 = -1;
vy_obs2 = -1;
v_obs2 = [vx_obs2;vy_obs2];
w_obs2 = 0.5;
fw2 = [-w*(y-yc2);w*(x-xc2)];
B2 = ((x-xc2)/a2)^2+((y-yc2)/b2)^2;  
n2 = [diff(B2,x);diff(B2,y)];
r2 = [x-xr2 ; y-yr2];
e2 = [-diff(B2,y);diff(B2,x)];
E2 = [r2,e2];
D2 = [1-1/(B2^(1/cushion)) , 0;
         0         , 1+1/(B2^(1/cushion))];     
M2 = E2*D2*inv(E2);
M2 = subs(M2,[x,y],[x/margin,y/margin]);

w1 = (B2_-1)/((B1_-1)+(B2_-1));
w2 = (B1_-1)/((B1_-1)+(B2_-1));

% Modulation Process
tstep = 0.01;
steps = 1000;

pos = [x0,y0];
pos_wrt_obs1 = pos;
pos_wrt_obs2 = pos;
pos_dot_avg = [];
pos_dot_v1 = [];
pos_dot_v2 = [];
pos_dot_wrt_obs1_sim = [];
pos_dot_wrt_obs2_sim = [];
B = [];
w_sim = [];

for i = 1:steps
    pos_current = pos(end,:);
    pos_dot = double(subs(f,[x,y],pos_current));
    
    pos_wrt_obs1_current = pos_wrt_obs1(end,:);
    pos_dot_wrt_obs1 = pos_dot - v_obs1;
    pos_dot_wrt_obs1 = Rotate(pos_dot_wrt_obs1,-w_obs1*tstep*(i-1));
    vw_obs1 = double(subs(fw1,[x,y,w],[pos_wrt_obs1_current,w_obs1]));
    pos_dot_wrt_obs1 = pos_dot_wrt_obs1 - vw_obs1;
    pos_dot_wrt_obs1 = double(subs(M1,[x,y],pos_wrt_obs1_current))*pos_dot_wrt_obs1;
    pos_dot_wrt_obs1 = pos_dot_wrt_obs1 + vw_obs1;
    pos_dot_wrt_obs1 = Rotate(pos_dot_wrt_obs1,w_obs1*tstep*(i-1));
    pos1_dot = pos_dot_wrt_obs1 + v_obs1;
    
    pos_wrt_obs2_current = pos_wrt_obs2(end,:);
    pos_dot_wrt_obs2 = pos_dot - v_obs2;
    pos_dot_wrt_obs2 = Rotate(pos_dot_wrt_obs2,-w_obs2*tstep*(i-1));
    vw_obs2 = double(subs(fw2,[x,y,w],[pos_wrt_obs2_current,w_obs2]));
    pos_dot_wrt_obs2 = pos_dot_wrt_obs2 - vw_obs2;
    pos_dot_wrt_obs2 = double(subs(M2,[x,y],pos_wrt_obs2_current))*pos_dot_wrt_obs2;
    pos_dot_wrt_obs2 = pos_dot_wrt_obs2 + vw_obs2;
    pos_dot_wrt_obs2 = Rotate(pos_dot_wrt_obs2,w_obs2*tstep*(i-1));
    pos2_dot = pos_dot_wrt_obs2 + v_obs2;
    
    B1_current = double(subs(B1,[x,y],pos_wrt_obs1_current));
    B2_current = double(subs(B2,[x,y],pos_wrt_obs2_current));
    w1_i = double(subs(w1,[B1_,B2_],[B1_current,B2_current]));
    w2_i = double(subs(w2,[B1_,B2_],[B1_current,B2_current]));
%     w1_i = 0;
%     w2_i = 1;
    pos1_dot_mag = norm(pos1_dot);
    pos2_dot_mag = norm(pos2_dot);
    pos1_dot_dir = pos1_dot/pos1_dot_mag;
    pos2_dot_dir = pos2_dot/pos2_dot_mag;
    pos_dot_mag = w1_i*pos1_dot_mag + w2_i*pos2_dot_mag;
    n_f = pos_dot/norm(pos_dot);
    e_f = [-n_f(2);n_f(1)];
    R_f = [n_f,e_f];
    v1 = R_f'*pos1_dot_dir;
    v2 = R_f'*pos2_dot_dir;
    k1 = acos(v1(1))*v1(2:end)/norm(v1(2:end));
    k2 = acos(v2(1))*v2(2:end)/norm(v2(2:end));
    k = w1_i*k1 + w2_i*k2;
    v = [cos(norm(k)) ; sin(norm(k))*k/norm(k)];
    pos_dot_dir = R_f*v;
    pos_dot = pos_dot_mag*pos_dot_dir;
    
    pos_dot_avg = [pos_dot_avg;pos_dot'];
    pos_dot_v1 = [pos_dot_v1;pos1_dot'];
    pos_dot_v2 = [pos_dot_v2;pos2_dot'];
    B = [B;B1_current,B2_current];
    w_sim = [w_sim;w1_i,w2_i];
    
    pos_dot_wrt_obs1 = pos_dot - v_obs1;
    pos_dot_wrt_obs1 = Rotate(pos_dot_wrt_obs1,-w_obs1*tstep*(i-1));
    vw_obs1 = double(subs(fw1,[x,y,w],[pos_wrt_obs1_current,w_obs1]));
    pos_dot_wrt_obs1 = pos_dot_wrt_obs1 - vw_obs1;
    
    pos_dot_wrt_obs2 = pos_dot - v_obs2;
    pos_dot_wrt_obs2 = Rotate(pos_dot_wrt_obs2,-w_obs2*tstep*(i-1));
    vw_obs2 = double(subs(fw2,[x,y,w],[pos_wrt_obs2_current,w_obs2]));
    pos_dot_wrt_obs2 = pos_dot_wrt_obs2 - vw_obs2;
    
    pos_next = pos_current' + pos_dot*tstep;
    pos_wrt_obs1_next = pos_wrt_obs1_current' + pos_dot_wrt_obs1*tstep;
    pos_wrt_obs2_next = pos_wrt_obs2_current' + pos_dot_wrt_obs2*tstep;
    
    pos = [pos;pos_next'];
    pos_wrt_obs1 = [pos_wrt_obs1; pos_wrt_obs1_next'];
    pos_wrt_obs2 = [pos_wrt_obs2; pos_wrt_obs2_next'];
    pos_dot_wrt_obs1_sim = [pos_dot_wrt_obs1_sim;pos_dot_wrt_obs1'];
    pos_dot_wrt_obs2_sim = [pos_dot_wrt_obs2_sim;pos_dot_wrt_obs2'];
end

%% 2 Moving Obstacle - Linear & Angular Velocity - Plot
close all

plot(pos(:,1),pos(:,2))

%% 2 Moving Obstacles - Linear & Angular Velocity - Movie
close all
clc

% Rotate = @(vec,theta) [cos(theta) -sin(theta);sin(theta) cos(theta)]*vec;

theta = [0:0.01:2*pi];

r1_theta  = [];
for i = 1:length(theta)
    r1_theta = [r1_theta,(a1*b1)/sqrt((b1*cos(theta(i)))^2 + (a1*sin(theta(i)))^2)];
end
x_ellipse1 = r1_theta .* cos(theta);
y_ellipse1 = r1_theta .* sin(theta);
ellipse1 = [x_ellipse1',y_ellipse1'];

r2_theta  = [];
for i = 1:length(theta)
    r2_theta = [r2_theta,(a2*b2)/sqrt((b2*cos(theta(i)))^2 + (a2*sin(theta(i)))^2)];
end
x_ellipse2 = r2_theta .* cos(theta);
y_ellipse2 = r2_theta .* sin(theta);
ellipse2 = [x_ellipse2',y_ellipse2'];

for k=1:length(pos)-1
%     subplot(3,1,1)
    plot(pos(k,1),pos(k,2),'x')
    hold on
    quiver(pos(k,1),pos(k,2),pos_dot_avg(k,1),pos_dot_avg(k,2),'k');
    quiver(pos(k,1),pos(k,2),pos_dot_v1(k,1),pos_dot_v1(k,2),'r');
    quiver(pos(k,1),pos(k,2),pos_dot_v2(k,1),pos_dot_v2(k,2),'b');
    ellipse1_rotate = Rotate(ellipse1',w_obs1*tstep*(k-1))';
    plot(ellipse1_rotate(:,1)+xc1+vx_obs1*tstep*(k-1),ellipse1_rotate(:,2)+yc1+vy_obs1*tstep*(k-1));
    ellipse2_rotate = Rotate(ellipse2',w_obs2*tstep*(k-1))';
    plot(ellipse2_rotate(:,1)+xc2+vx_obs2*tstep*(k-1),ellipse2_rotate(:,2)+yc2+vy_obs2*tstep*(k-1));
    axis([-5,5,-5,5])
    title(sprintf("B1 = %.2f   and   B2 = %.2f \n w1 = %.2f   and   w2 = %.2f",B(k,1),B(k,2),w_sim(k,1),w_sim(k,2)));
    hold off
    
%     subplot(3,1,2)
%     plot(pos_wrt_obs1(k,1),pos_wrt_obs1(k,2),'x')
%     hold on;
%     plot(x_ellipse1+xc1,y_ellipse1+yc1);
%     quiver(pos_wrt_obs1(k,1),pos_wrt_obs1(k,2),pos_dot_wrt_obs1_sim(k,1),pos_dot_wrt_obs1_sim(k,2),'r');
%     hold off;
%     
%     subplot(3,1,3)
%     plot(pos_wrt_obs2(k,1),pos_wrt_obs2(k,2),'x')
%     hold on;
%     plot(x_ellipse2+xc2,y_ellipse2+yc2);
%     quiver(pos_wrt_obs2(k,1),pos_wrt_obs2(k,2),pos_dot_wrt_obs2_sim(k,1),pos_dot_wrt_obs2_sim(k,2),'b');
%     hold off;
    
    pause(tstep)
end