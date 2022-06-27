%% 1 Static Obstacle 
clear
close all
clc

syms x y

x0 = 5;
y0 = -5;

xc = 3;
yc = -3;
a = 2;
b = 1;
xr = xc+0;
yr = yc+0.1;

xa = 0;
ya = 0;

B = (((x-xc)/a)^2) +(((y-yc)/b)^2); 

n = [diff(B,x);diff(B,y)];
r = [x-xr ; y-yr];
e = [-diff(B,y);diff(B,x)];
E = [r,e];

D = [1-1/B , 0;
         0         , 1+1/B];
     
M = E*D*inv(E);

f = [-(x-xa);-(y-ya)];
fm = M*f;


tstep = 0.01;
steps = 500;

pos_sim = [x0 , y0];
pos_dot_sim = [];
B_sim = [];
for i = 1:steps
    pos_current = pos_sim(end,:);
    B_sim = [B_sim;double(subs(B,[x,y],pos_current))];
    pos_dot = double(subs(fm,[x,y],pos_current));
    pos_new = pos_current' + pos_dot*tstep;
    pos_sim = [pos_sim;pos_new'];
    pos_dot_sim = [pos_dot_sim;pos_dot'];
end

%% 1 Static Obstacle  - Plot
close all
clc

theta = [0:0.01:2*pi];
r_theta  = [];
for i = 1:length(theta)
    r_theta = [r_theta,(a*b)/sqrt((b*cos(theta(i)))^2 + (a*sin(theta(i)))^2)];
end

x_ellipse = r_theta .* cos(theta);
y_ellipse = r_theta .* sin(theta);
plot(pos_sim(:,1),pos_sim(:,2))
hold on;
plot(x_ellipse+xc,y_ellipse+yc);
hold off;

%% 1 Static Obstacle  - Movie
close all
clc

for k=1:length(pos_sim)-1
    plot(pos_sim(k,1),pos_sim(k,2),'x')
    hold on
    plot(x_ellipse+xc,y_ellipse+yc);
    hold on
    quiver(pos_sim(k,1),pos_sim(k,2),pos_dot_sim(k,1),pos_dot_sim(k,2),0)
    axis([-5,5,-5,5])
    title(sprintf("B = %.2f \n x = %.2f   and   y = %.2f",B_sim(k,1),pos_sim(k,1),pos_sim(k,2)));
    hold off
    pause(0.05)
end
