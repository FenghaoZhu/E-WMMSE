function [iter1, time, rate] = Test_E_WMMSE(H_E_WMMSE, K, T, R, epsilon, sigma2, snr, I, alpha1, d, max_iter, beta)
% TEST_E_WMMSE 是用来测试WMMSE性能的函数
% 输入函数运行参数，返回迭代次数，运行时间与速率信息
P = db2pow(snr)*sigma2; % 发射功率
rate = []; % 初始化一个空向量记录rate
time = []; % 初始化一个空运行时间记录rate
tic; %开始计时
begin_time = toc; % 标记开始时间，结束时间减去开始时间就是使用的时间
time = [time 0];

% 初始化信道向量
H = H_E_WMMSE;

rate = []; % 初始化一个空向量记录rate

% 初始化W和U矩阵
U =randn(R,d,I) + 1j*randn(R,d,I);
W = zeros(d,d,I);
for i=1:I
    W(:,:,i)=eye(d,d);
end

% 初始化波束赋形矩阵
V = zeros(T,d,I); % 每个基站的波束赋形矩阵
for i=1:I
        v = sqrt(1/2)*(randn(T,d)+1i*randn(T,d));
        V(:,:, i) = sqrt(P/(I*trace(v*v')))*v;
end 

% 求初始化发射波束V后求系统和速率
rate_old = sum_rate(H,V,sigma2,R,I,alpha1);
rate = [rate rate_old];

V_list = []; % 初始化一个空向量记录V_list
iter1 = 1; % 迭代次数
while(1)
    % E-WMMSE的EXTRAPOLATION TECHNIQUE，特殊的初始化
    if iter1 >=3
        Q = normalize_V(V_list(:,(iter1-2)*d + 1:(iter1-1)*d,:) + beta*(V_list(:,(iter1-2)*d + 1:(iter1-1)*d,:)-V_list(:,(iter1-3)*d + 1:(iter1-2)*d,:)),I, P);
    elseif iter1 == 2
        Q = V_list(:,(iter1-2)*d+1:(iter1-1)*d,:);
    elseif iter1 == 1
        Q = V;
    end
    U = find_U(H,Q,sigma2, P, R,I,d); 
    W = find_W(U,H,Q,R,I,d); 
    V = find_V(alpha1,H,sigma2,U,W,T, R, I,d ,P); 
    rate_new = sum_rate(H,V,sigma2,R,I,alpha1); % 计算和速率
    V_list = [V_list V]; % 初始化一个空向量记录V_list
    rate = [rate rate_new];
    elapsed_time = toc;
    elapsed_time = elapsed_time - begin_time;
    time = [time elapsed_time];
    iter1 = iter1 + 1;
    if abs(rate_new-rate_old) / rate_old < epsilon || iter1 > max_iter
        break;
    end
    rate_old = rate_new;
end
end

