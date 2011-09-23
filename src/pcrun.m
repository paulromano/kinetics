%% Kinetics Parameters
beta_sum = 0.00765; L = 1e-2;
lambda = [3.0100 1.1400 0.3010 0.1110 0.0305 0.0124];
beta = beta_sum*[0.041 0.115 0.396 0.196 0.219 0.033];

%% Initial Conditions
m = length(lambda);
init_cond = zeros(m+1,1);
init_cond(1) = 6e-4;
for i = 1:m
    init_cond(i+1) = init_cond(1) * beta(i) / (L* lambda(i));
end

%% Run Code
step = 0.1; target = 400;
insert_time = 150; velocity = 12; rho_ex = 0.40;
z = solvePKE(lambda, beta, beta_sum, L, target, step, 1, ...
    init_cond, velocity, insert_time, rho_ex);
init_cond = [0;init_cond];
z = [init_cond,z];

%% Plot Results
plot(z(1,:),z(2,:),'k');
box on
xlabel('Time (s)')
ylabel('Power (W)')
xlim([0 target])
set(gca,'YScale','log')

% Determine insertion power and maximum power
max_power = 0;
for i = 1 : size(z,2)
    if z(2,i) > max_power
        max_power = z(2,i);
        index = i;
    end
end
insert_power = z(2,1+round(insert_time/step));
display(insert_power)
display(max_power)
line([0,round(index*step)],[max_power,max_power],'Color','g','LineStyle',':')
line([round(index*step),round(index*step)],[1e-4,max_power],...
    'Color','g','LineStyle',':')
line([insert_time,insert_time],[1e-4,insert_power],...
    'Color','r','LineStyle','-.')
line([0,insert_time],[insert_power,insert_power],...
    'Color','r','LineStyle','-.')