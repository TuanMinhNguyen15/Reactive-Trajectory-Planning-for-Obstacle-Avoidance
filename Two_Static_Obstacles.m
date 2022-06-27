%% 2 Static Obstacles
clear 
close all
clc

syms x y

x0 = -2;
y0 = 0.1;

xc1 = 0.5;
yc1 = 0;
a1 = 1;
b1 = 1;

xc2 = -0.5;
yc2 = 0;
a2 = 1;
b2 = 1;

xa = 4;
ya = 0;

B1 = ((x-xc1)/a1)^2+((y-yc1)/b1)^2;  
B2 = ((x-xc2)/a2)^2+((y-yc2)/b2)^2;  

r1 = [x-0;y-0];
e1 = [-diff(B1,y);diff(B1,x)];
r2 = [x-0;y-0];
e2 = [-diff(B2,y);diff(B2,x)];

E1 = [r1,e1];
D1 = [1-1/B1 , 0;
         0         , 1+1/B1];
M1 = E1*D1*inv(E1);

E2 = [r2,e2];
D2 = [1-1/B2 , 0;
         0         , 1+1/B2];
M2 = E2*D2*inv(E2);

f = [-(x-xa);-(y-ya)];
fm1 = M1*f;
fm2 = M2*f;

w1 = (B2-1)/((B1-1)+(B2-1));
w2 = (B1-1)/((B1-1)+(B2-1));


tstep = 0.01;
steps = 500;

pos_sim = [x0 , y0];
pos_dot_sim = [];
for i = 1:steps
    pos_current = pos_sim(end,:);
    
    w1_i = double(subs(w1,[x,y],pos_current));
    w2_i = double(subs(w2,[x,y],pos_current));
    pos1_dot = double(subs(fm1,[x,y],pos_current));
    pos2_dot = double(subs(fm2,[x,y],pos_current));
    pos1_dot_mag = norm(pos1_dot);
    pos2_dot_mag = norm(pos2_dot);
    pos1_dot_dir = pos1_dot/pos1_dot_mag;
    pos2_dot_dir = pos2_dot/pos2_dot_mag;
    
    pos_dot_mag = w1_i*pos1_dot_mag + w2_i*pos2_dot_mag;
    
    pos_dot_f = double(subs(f,[x,y],pos_current));
    n_f = pos_dot_f/norm(pos_dot_f);
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
    
    pos_new = pos_current' + pos_dot*tstep;
    pos_sim = [pos_sim;pos_new'];
    pos_dot_sim = [pos_dot_sim;pos_dot'];
end

%% 2 Static Obstacles - Plot
close all
clc

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

plot(pos_sim(:,1),pos_sim(:,2))
hold on;
plot(x_ellipse1+xc1,y_ellipse1+yc1);
plot(x_ellipse2+xc2,y_ellipse2+yc2);
hold off;

%% 2 Static Obstacles  - Movie
close all
clc

theta = [0:0.01:2*pi];
x_circle = cos(theta);
y_circle = sin(theta);



for k=1:length(pos_sim)-1
    plot(pos_sim(k,1),pos_sim(k,2),'x')
    hold on
    plot(xa,ya,'p')
    plot(x_ellipse1+xc1,y_ellipse1+yc1);
    plot(x_ellipse2+xc2,y_ellipse2+yc2);
    hold on
    quiver(pos_sim(k,1),pos_sim(k,2),pos_dot_sim(k,1),pos_dot_sim(k,2),0,'r')
    axis([-5,5,-5,5])
    hold off
    pause(tstep)
end
