%parameters
alpha = 0.8;
g = 9.81;
omega = 5500;
% v_b = 5;
% a = v_b/omega;
% sigma = a*omega^2/(pi*9.81);

%红点的个数,为一整数
point_number = 50;

%速度的计算范围
v_range = 50;

sigma = 1154;
a = sigma*9.81*pi/(omega^2);
v_b = a*omega;

% filename = 'v_b_5_alpha_0.8.xlsx';
% A = xlsread(filename);
% 
% t_0 = A(2,:);
t = mod(t,1);
% w = A(1,:);

v = pi*g.*(w+ sigma.*cos(2.*pi.*t))/omega;
p_1 = [0];

v_0 = 0: 0.001: v_range;
%p_v = @(v,v_b,alpha) (2.*((2.*(1-alpha)./((1+alpha).^2.*(v_b.^2))).^((3+alpha)/(2.*(1+alpha)))./(gamma((3+alpha)./(2.*(1+alpha))))).*(v.^(2./(1+alpha))).*exp(-2.*(1-alpha).*(v.^2)./(((1+alpha).^2).*(v_b.^2))));
p_0 = p_v(v_0,v_b,0.56);
count = 0;

%the distribution from the simulation
for i = 1:1000
    if v(i) < 0 
        v(i) = 0;
        count = count + 1;
    end

    temp = fix(v(i)/(v_range/point_number)) + 1;
    num = size(p_1);
    if num(2) < temp
        p_1(temp) = 0;
    end
    p_1(temp) = p_1(temp)+1;
end

v_1 = 0:50/50:50;


p_1 = p_1./(1000*v_range/point_number);

p_1(point_number+1) = 0;


plot(v_0,p_0);
hold on;
plot(v_1,p_1,'ro','MarkerSize',3);
title("I=2.8A, alpha=0.56, omega=2950/s")
legend('Theoretical Result','Matlab Simulation')
xlabel('v/(m.s^-1)')
ylabel('P(v)')

%Boltzmann Distribution function
function prob = p_v(v, v_b, alpha)
    prob = 2.*((2.*(1-alpha)./((1+alpha).^2.*(v_b.^2))).^((3+alpha)/(2.*(1+alpha)))./(gamma((3+alpha)./(2.*(1+alpha))))).*(v.^(2./(1+alpha))).*exp(-2.*(1-alpha).*(v.^2)./(((1+alpha).^2).*(v_b.^2)));
end