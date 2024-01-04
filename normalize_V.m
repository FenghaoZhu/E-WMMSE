function [V_normalized] = normalize_V(V, I, P)
%此函数用来做波束赋形向量的功率约束
sum_power = 0; % 现有的波束赋形向量功率和
for i=1:I
        v = V(:,:, i);
        sum_power = sum_power + trace(v*v');
end 
V = sqrt(P/sum_power) * V;
V_normalized = V;
end

