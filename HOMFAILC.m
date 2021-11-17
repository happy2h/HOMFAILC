% 期望轨迹
for t = 1:1:100
	if t <= 30
		yd(t+1) = 0.5*(-1)^round(t/10);
	elseif t > 30 && t <= 70
		yd(t+1) = 0.5*sin((t*pi)/10) + 0.3*cos((t*pi)/10);
	else
		yd(t+1) = 0.5*(-1)^round(t/10);
    end
end

% 参数设置
epsilon = 0.01; % 拟伪偏导重置阈值
lambda = 0.6; % 控制律步长
rho = 1;  % 控制律增益
mu = 0.1;    % 参数估计步长
eta = 0.6;  % 参数估计增益
% 最佳高阶参数
alpha_1 = 1.2;
alpha_2 = -0.1;
alpha_3 = -0.1;

% 迭代过程
i_n = 100; % 迭代次数
y(1:i_n,1:100) = 0; % 设置输出
for k = 1:1:i_n
	for t = 1:1:100 % 每次迭代学习所有时刻点
		if k == 1
			phi(k,t) = 10; % 初始拟伪偏导
		elseif k == 2
			phi(k,t) = phi(k-1,t) + (eta*(u(k-1,t) - 0)/(mu + norm(u(k-1,t) - 0)^2))*(y(k-1,t+1) - 0 - phi(k-1,t)*(u(k-1,t) - 0));
		elseif k == 3
			phi(k,t) = phi(k-1,t) + (eta*(u(k-1,t) - 0)/(mu + norm(u(k-1,t) - 0)^2))*(y(k-1,t+1) - 0 - phi(k-1,t)*(u(k-1,t) - 0));
		else 
			phi(k,t) = (u(k-1,t) - u(k-2,t))*(y(k-1,t+1) - y(k-2,t+1))/(mu + norm(u(k-1,t) - u(k-2,t))^2) + (mu*eta/(mu + norm(u(k-1,t) - u(k-2,t))^2))*(alpha_1*phi(k-1,t) + alpha_2*phi(k-2,t) + alpha_3*phi(k-3,t));
        end
		if k == 1
			u(k,t) = 0; % 初始控制信号
		else
			u(k,t) = u(k-1,t) + (rho*phi(k,t)/(lambda + norm(phi(k,t))^2))*e(k-1,t+1);
        end
	
	%% 判断拟伪偏导是否需要重置
        if k >= 2 && (phi(k,t) <= epsilon || (abs(u(k,t) - u(k-1,t) <= epsilon)))
			phi(k,t) = phi(1,t);
		end

		% 仿真的非线性系统函数，对于实际控制来说，系统函数是隐式未知的
		if t <= 50
			y(k,t+1) = y(k,t)/(1 + norm(y(k,t))^2) + (u(k,t))^3;
		else
			y(k,t+1) = (y(k,t)*y(k,t-1)*y(k,t-2)*u(k,t-1)*(y(k,t-2) - 1) + (1+round(t/50))*u(k,t))/(1 + norm(y(k,t-1))^2 + norm(y(k,t-2))^2);
		end
		e(k,t+1) = yd(t+1) - y(k,t+1);
	end
end
for i = 1:1:i_n
	e_max(i) = max(abs(e(i,:)));
end

figure(1) 
plot(yd,'r'); hold on;
plot(y(i_n,:),'b');
figure(2);
plot(e_max);title('error of time k');
