%????
for t = 1:1:100
	if t <= 30
		yd(t+1) = 0.5*(-1)^round(t/10);
	elseif t > 30 && t <= 70
		yd(t+1) = 0.5*sin((t*pi)/10) + 0.3*cos((t*pi)/10);
	else
		yd(t+1) = 0.5*(-1)^round(t/10);
    end
end

%?????
epsilon = 0.05;
lambda = 0.6; %0.6
rho = 1;  %1
mu = 0.1;    %1
eta = 0.6;  %0.6
alpha_1 = 0.4;
alpha_2 = 0.4;
alpha_3 = 0.2;

%??ƹ???
i_n = 100; %???????
y(1:i_n,1:100) = 0;
for k = 1:1:i_n
	for t = 1:1:100 %ʱ????
		if k == 1
			phi(k,t) = 10;
		elseif k == 2
			phi(k,t) = phi(k-1,t) + (eta*(u(k-1,t) - 0)/(mu + norm(u(k-1,t) - 0)^2))*(y(k-1,t+1) - 0 - phi(k-1,t)*(u(k-1,t) - 0));
		elseif k == 3
			phi(k,t) = phi(k-1,t) + (eta*(u(k-1,t) - 0)/(mu + norm(u(k-1,t) - 0)^2))*(y(k-1,t+1) - 0 - phi(k-1,t)*(u(k-1,t) - 0));
		else 
			phi(k,t) = (u(k-1,t) - u(k-2,t))*(y(k-1,t+1) - y(k-2,t+1))/(mu + norm(u(k-1,t) - u(k-2,t))^2) + (mu*eta/(mu + norm(u(k-1,t) - u(k-2,t))^2))*(alpha_1*phi(k-1,t) + alpha_2*phi(k-2,t) + alpha_3*phi(k-3,t));
        end
		if k == 1
			u(k,t) = 0;
		else
			u(k,t) = u(k-1,t) + (rho*phi(k,t)/(lambda + norm(phi(k,t))^2))*e(k-1,t+1);
        end
        if k >= 2 && (phi(k,t) <= epsilon || (abs(u(k,t) - u(k-1,t) <= epsilon)))
			phi(k,t) = phi(1,t);
		end

		%ϵͳ????
		if t <= 50
			y(k,t+1) = y(k,t)/(1 + norm(y(k,t))^2) + (u(k,t))^3;
		else
			y(k,t+1) = (y(k,t)*y(k,t-1)*y(k,t-2)*u(k,t-1)*(y(k,t-2) - 1) + 0.1*round(t/50)*u(k,t))/(1 + norm(y(k,t-1))^2 + norm(y(k,t-2))^2);
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
