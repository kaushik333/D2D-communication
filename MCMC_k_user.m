clc;
clear;
close;
Pd_intended = input('What is the power intended at the base station by each cellular user?');
P_d2d_intended = input('What is the total intended interference at the BS due to all D2D users?');
Pb_intended = input('What is the intended power at the cell edge? (in watts please)');

%num = [5 6 7 8 9 10 15 20 25 30 35 40 45 50 60 70 80 100 150 200];
%num = [1:9 10:5:25 30:10:70 100];
num = [80, 90];
side_len = 400; % units: meters

N = 1; % noise power

eta = 4;

d0 = 10; % meters
%sum_rates = zeros(size(num,2), 100);
for s = 1:size(num,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting the data points
%%% Red is transmitter and blue is Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(121);
num_users = num(s);
points_required = 2*(num_users) + 1;
%side_len = 400; % units: meters
for iter = 1:100
t = side_len*rand(num_users,2);
r = side_len*rand(num_users,2);
t = [t; side_len/2 side_len/2]; %% Base station
r = [r; side_len/2 side_len/2]; %% Base station

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate distance from every Tx or every other Rx and the BS. Note that
%%% the last element of the array will be 0, since its the distance between
%%% the base station and itself. Then calculate the channel gains too
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:num_users+1
    for j=1:num_users+1
        dis(i,j) = norm(t(i,:)-r(j,:));
        h(i,j) = min(1,((dis(i,j)/d0)^(-eta)));
    end
end

%rand_vec = round(rand(1,num_users));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MCMC Optimisation Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% while(size(find(rand_vec==1),2)<=2)
%     rand_vec = round(rand(1,num_users));
% end
% disp(rand_vec);

tic;

%rand_vec = ones(1,num_users);
rand_vec = round(rand(1,num_users));
%rand_vec = zeros(1,num_users);
disp(rand_vec);
d2d_pair_first = find(rand_vec==1) % All d2d users
cell_pair_first = find(rand_vec==0) % All cellular users
[R_cell_up_first, R_cell_down_first, sum_rate] = sumrate(d2d_pair_first, cell_pair_first, dis, h, num_users, Pd_intended, Pb_intended, P_d2d_intended, d0, eta, side_len)
a_old = 0;
%a_new = randi(num_users);
a_new = 0; 

if (2^num_users) < 10000
    iterations = 2^(num_users-1);
else
    iterations = 200*num_users;
end

for i = 1:iterations

rand_vec_new = rand_vec;

%To pick a random element and toggle its state

%1 means d2d user and 0 means cellular user


while(a_old==a_new)
    a_new = randi(num_users);
end
%a_new
if(rand_vec_new(a_new)==0)
    rand_vec_new(a_new)=1;
elseif(rand_vec_new(a_new)==1)
    rand_vec_new(a_new)=0;
end


a_old = a_new;
d2d_pair = find(rand_vec_new==1); % All d2d users
cell_pair = find(rand_vec_new==0); % All cellular users
%disp('iteration is');disp(i);
% cell_pair_first
% R_cell_up_first
[delta_R, R_cell_up, R_cell_down] = sumrate1(d2d_pair, d2d_pair_first, cell_pair, cell_pair_first, R_cell_up_first, R_cell_down_first, h, dis, num_users, Pd_intended, Pb_intended, P_d2d_intended, d0, eta, side_len);
disp(['delta R is ', num2str(delta_R)]);
disp(num(s));
disp(iter);
if(delta_R > 0)
    rand_vec = rand_vec_new;
    d2d_pair_first = d2d_pair;
    cell_pair_first = cell_pair;
    R_cell_up_first = R_cell_up;
    R_cell_down_first = R_cell_down;
    %sum_rate = sum_rate_new;
end



end

time = toc; 
disp('The time taken is : '); disp(time);

disp(rand_vec);

d2d_pair = find(rand_vec==1); % All d2d users
cell_pair = find(rand_vec==0); % All cellular users

[R_up, R_down, sum_rates(s,iter)] = sumrate(d2d_pair, cell_pair, dis, h, num_users, Pd_intended, Pb_intended, P_d2d_intended, d0, eta, side_len);

end
end
sum_rates_new = sum_rates';
sum = mean(sum_rates_new);
plot(num, sum);