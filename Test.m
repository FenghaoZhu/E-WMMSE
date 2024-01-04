% 此脚本用来运行WMMSE和E-WMMSE这两个算法的性能差距，指标包括运行时间，迭代次数以及最终速率
% 暂时此脚本只支持单基站的仿真情景

% 定义运行参数
clc;clear;
rng(1);
load('processed_channel.mat','I','N','T','ite','HH','Theta',...
    'GG', 'snr','omega'); % 加载信道，用户信息等参数

K=1; % 基站个数，目前由于R-WMMSE程序只支持一个基站，故只能固定为1
% T=64; % 基站天线数量
R=1; % 每个用户天线数量
epsilon=0.001; % 收敛设定的限制
sigma2=1; % 噪声功率
% snr=10; % 用户的信噪比
% I=16; % 用户数量
alpha1=ones(I,K); % 用户权重
alpha1(:,1) = omega;
d=1; % 每个用户流数
max_iter=100; % 最大的迭代次数
num_sample = 100; % 信道样本数量
beta = 0.9; % E-WMMSE里面的超参数

% % 生成信道，以便于在相同信道的条件下测试R-WMMSE和WMMSE的性能差别
% H_WMMSE = cell(I,K,K); % 生成原版WMMSE所需的信道
% for i=1:I
%     for k = 1:K
%         for j=1:K
%            H_WMMSE{i,k,j}=sqrt(1/2)*(randn(R,T)+1i*randn(R,T)); % 圆复高斯信道
%         end
%     end
% end
% 
% H_E_WMMSE = zeros(R,T,I); % 生成R-WMMSE所需的信道
% for i=1:I
%     H_E_WMMSE(:,:,i) = H_WMMSE{i};
% end

H = zeros(N,R,I); % RIS到每个用户的信道 
theta = exp(1j*zeros(N,1)); % 随机初始化RIS相位，每个元素模长为1
G = squeeze(GG(:,:,1));
for i=1:I
    H(: , 1, i)=HH(i,:,1); % 从HH和GG中取第一个样本
end

% 计算级联信道
cascaded_chanel = zeros(R, T, I);
for i=1:I
cascaded_chanel(:,:,i) = H(:,:,i)'*diag(theta)*G; 
end

H_E_WMMSE = cascaded_chanel;

% 运行对比测试,生成前两张图
[iter_E_WMMSE, time_E_WMMSE, rate_E_WMMSE] = Test_E_WMMSE(H_E_WMMSE, K,T,R,epsilon,sigma2,snr,I,alpha1,d,max_iter,beta);
[iter_WMMSE, time_WMMSE, rate_WMMSE] = Test_WMMSE(H_E_WMMSE, K,T,R,epsilon,sigma2,snr,I,alpha1,d,max_iter);

figure(1);
plot(0:iter_E_WMMSE-1,rate_E_WMMSE, '-sb')
hold on
plot(0:iter_WMMSE-1,rate_WMMSE, '-*r')
grid on
xlabel('Iterations')
ylabel('Sum rate (bits per channel use)')
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
title(['Compare WMMSE with E-WMMSE, K=',num2str(K), ',', 'T=', num2str(T), ',', 'R=', num2str(R), ',','\epsilon=', num2str(epsilon)])
legend('E-WMMSE','WMMSE')
savefig(['./figs/T=',num2str(T),', I=',num2str(I),', d=',num2str(d),',',num2str(snr),'dB, sumrate.fig'])

figure(2);
plot(0:iter_E_WMMSE-1,time_E_WMMSE, '-bs')
hold on
plot(0:iter_WMMSE-1,time_WMMSE,  'd-r', 'MarkerFaceColor', 'r')
grid on
xlabel('Iterations')
ylabel('Time (s)')
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
title(['Compare WMMSE with E-WMMSE, K=',num2str(K), ',','T=', num2str(T), ',','R=', num2str(R), ',','\epsilon=', num2str(epsilon)])
legend('E-WMMSE','WMMSE')
savefig(['./figs/T=',num2str(T),', I=',num2str(I),', d=',num2str(d),',',num2str(snr),'dB, runtime.fig'])

% 运行对比测试，生成运行时间关于用户数量变化的图
% begin_user=5;
% end_user=45;
% len=end_user-begin_user+1;
% time1 = zeros(len,1); % WMMSE的不同用户数量的运行时间
% time2 = zeros(len,1); % R-WMMSE的不同用户数量的运行时间
% 
% bar = waitbar(0,'开始测试');    % waitbar显示进度条
% for num_user=begin_user:end_user
%     for f = 1:num_sample
%     [iter_R_WMMSE, time_R_WMMSE, rate_R_WMMSE] = Test_R_WMMSE(K,T,R,epsilon,sigma2,snr,num_user,ones(num_user,1),d,max_iter);
%     [iter_WMMSE, time_WMMSE, rate_WMMSE] = Test_WMMSE(K,T,R,epsilon,sigma2,snr,num_user,ones(num_user,1),d,max_iter);
%     time1(num_user-begin_user+1)=time1(num_user-begin_user+1)+time_WMMSE(iter_WMMSE);
%     time2(num_user-begin_user+1)=time2(num_user-begin_user+1)+time_R_WMMSE(iter_R_WMMSE);
%     str=['计算中...',num2str(100*f/100),'%'];% 百分比形式显示处理进程
%     waitbar((num_user-begin_user+1)/len,bar,str) % 更新进度条bar
%     end
%     time1(num_user-begin_user+1) = time1(num_user-begin_user+1) / (num_sample);
%     time2(num_user-begin_user+1) = time2(num_user-begin_user+1) / (num_sample);
% end
% close(bar); % 循环结束关闭进度条
% 
% figure(1);
% plot(begin_user:end_user,time2, '-sb') % R-WMMSE的不同用户数量的运行时间
% hold on
% plot(begin_user:end_user,time1, '-*r') % WMMSE的不同用户数量的运行时间
% grid on
% xlabel('Number of users')
% ylabel('Average CPU Time (s)')
% set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
% title(['Compare WMMSE with R-WMMSE, K=',num2str(K), ',', 'T=', num2str(T), ',', 'R=', num2str(R), ',','\epsilon=', num2str(epsilon)])
% legend('R-WMMSE','WMMSE')
% savefig('./figs/T=128, I=12, d=2, 10dB, runtime.fig')

% 运行对比测试，生成运行时间关于基站天线数量变化的图
% antenna_num_pool = [32, 64, 128, 256, 512];
% len = length(antenna_num_pool);
% time1 = zeros(len,1); % WMMSE的不同用户数量的运行时间
% time2 = zeros(len,1); % R-WMMSE的不同用户数量的运行时间
% bar = waitbar(0,'开始测试');    % waitbar显示进度条
% for num_antenna_index=1:len
%     for f = 1:num_sample
%     [iter_R_WMMSE, time_R_WMMSE, rate_R_WMMSE] = Test_R_WMMSE(K,antenna_num_pool(num_antenna_index),R,epsilon,sigma2,snr,I,ones(I,1),d,max_iter);
%     [iter_WMMSE, time_WMMSE, rate_WMMSE] = Test_WMMSE(K,antenna_num_pool(num_antenna_index),R,epsilon,sigma2,snr,I,ones(I,1),d,max_iter);
%     time1(num_antenna_index)=time1(num_antenna_index)+time_WMMSE(iter_WMMSE);
%     time2(num_antenna_index)=time2(num_antenna_index)+time_R_WMMSE(iter_R_WMMSE);
%     str=['计算中...',num2str(100*f/100),'%'];% 百分比形式显示处理进程
%     waitbar((num_antenna_index)/len,bar,str) % 更新进度条bar
%     end
%     time1(num_antenna_index) = time1(num_antenna_index) / (num_sample);
%     time2(num_antenna_index) = time2(num_antenna_index) / (num_sample);
% end
% close(bar); % 循环结束关闭进度条
% 
% figure(1);
% semilogy(antenna_num_pool,time2, '-sb') % R-WMMSE的不同用户数量的运行时间
% hold on
% semilogy(antenna_num_pool,time1, '-*r') % WMMSE的不同用户数量的运行时间
% grid on
% xlabel('Number of BS Antennas')
% ylabel('Average CPU Time (s)')
% set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
% title(['Compare WMMSE with R-WMMSE, K=',num2str(K), ',', 'T=', num2str(T), ',', 'R=', num2str(R), ',','\epsilon=', num2str(epsilon)])
% legend('R-WMMSE','WMMSE')
% savefig('./figs/I=16, d=4, 10dB, runtime1.fig')

% 运行遍历测试，测试所有的可能组合并且保存数据以及图
% T_pool = [32, 64, 128, 256, 512, 1024]; % 天线数量范围
% I_num_pool = 5:5:45; % 用户数量范围
% d_num_pool = 1:8;  % 用户数据流数数量范围
% snr_num_pool = -10:5:30; % 信噪比范围
% [iter_R_WMMSE, time_R_WMMSE, rate_R_WMMSE] = Test_R_WMMSE(K,T,R,epsilon,sigma2,snr,I,alpha1,d,max_iter);
% [iter_WMMSE, time_WMMSE, rate_WMMSE] = Test_WMMSE(K,T,R,epsilon,sigma2,snr,I,alpha1,d,max_iter);
% 
% figure(1);
% plot(0:iter_R_WMMSE-1,rate_R_WMMSE, '-sb')
% hold on
% plot(0:iter_WMMSE-1,rate_WMMSE, '-*r')
% grid on
% xlabel('Iterations')
% ylabel('Sum rate (bits per channel use)')
% set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
% title(['Compare WMMSE with R-WMMSE, K=',num2str(K), ',', 'T=', num2str(T), ',', 'R=', num2str(R), ',','\epsilon=', num2str(epsilon)])
% legend('R-WMMSE','WMMSE')
% savefig(['./figs/T=',num2str(T),', I=',num2str(I),', d=',num2str(d),',',num2str(snr),'dB, sumrate.fig'])
% 
% figure(2);
% plot(0:iter_R_WMMSE-1,time_R_WMMSE, '-bs')
% hold on
% plot(0:iter_WMMSE-1,time_WMMSE,  'd-r', 'MarkerFaceColor', 'r')
% grid on
% xlabel('Iterations')
% ylabel('Time (s)')
% set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
% title(['Compare WMMSE with R-WMMSE, K=',num2str(K), ',','T=', num2str(T), ',','R=', num2str(R), ',','\epsilon=', num2str(epsilon)])
% legend('R-WMMSE','WMMSE')
% savefig(['./figs/T=',num2str(T),', I=',num2str(I),', d=',num2str(d),',',num2str(snr),'dB, runtime.fig'])





