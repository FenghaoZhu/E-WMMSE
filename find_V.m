function V = find_V(alpha1, H,sigma2, U, W, T , R ,I ,d ,P )

    % eta = alpha1(1,1) * trace(U(:,:,1)*W(:,:,1)*(U(:,:,1)'));
    % W_hat = alpha1(1,1) * W(:,:,1);
    % U_hat = U(:,:,1);
    % % 计算分块矩阵，为了进一步降低复杂度
    % for i = 2:I
    %     W_hat = blkdiag(W_hat, alpha1(i,1) * W(:,:,i));
    %     U_hat = blkdiag(U_hat, U(:,:,i));
    %     eta = eta + alpha1(i,1) * trace(U(:,:,i)*W(:,:,i)*(U(:,:,i)'));
    % end
    % 
    % eta = eta * sigma2 / P;
    % % 计算低维度的坍塌矩阵X_hat
    % X_hat = U_hat / (eta * inv(W_hat) + U_hat' * H_hat * U_hat);
    % 
    % % 将低维度的坍塌矩阵X_hat由二维矩阵形式转换为三维矩阵形式以兼容代码数据格式
    % X = zeros(R*I,d,I); % 每个基站的坍塌波束赋形矩阵
    % for i = 1:I
    %     X(:,:,i) = X_hat(:, d*(i-1) + 1:d*i);
    % end
    
 
    J=zeros(T, T);
    V=zeros(T,d, I);
   
            for l=1:I
                J = J + alpha1(l, 1) * H(:,:,l)'*U(:,:,l)*W(:,:,l)*(U(:,:,l)')*(H(:,:,l));   
            end

            max_iter = 100; % 二分法查找最优对偶变量\mu
            mu = zeros(1,1);
            mu_min = 0;
            mu_max = 10;
            iter = 0;
            while(1)
                mu1 = (mu_max+mu_min) / 2;
                P_tem = 0;

                for i=1:I % 计算功率和
                    V_tem = ((J+mu1*eye(T))) \ (alpha1(i,1)*(H(:,:,i)'*U(:,:,i)*W(:,:,i))); % 公式15
                    P_tem = P_tem + real(trace(V_tem*V_tem'));
                end

                if P_tem > P
                    mu_min = mu1;
                else
                    mu_max = mu1;
                end
                iter = iter + 1;

                if abs(mu_max - mu_min) < 1e-5 || iter > max_iter
                    break
                end

            end

            mu = mu1;

            for l=1:I
                V(:,:,l) = (J+mu*eye(T, T)) \ (alpha1(l, 1) * (H(:,:,l)'*U(:,:,l)*W(:,:,l))); 
            end

end































%     V = cell(I,K);
%     A = cell(K);
%     A(:) = {zeros(T,T)};
% 
% 
%     for k=1:K   % 公式15括号内求和部分     
%         for j=1:K
%             for l=1:I
%                 A{k} = A{k} + alpha1(l,j)*H{l,j,k}'*U{l,j}*W{l,j}*(U{l,j}')*H{l,j,k};
%             end
%         end   
%     end 
% 
%     max_iter = 100; % 二分法查找最优对偶变量\mu
%     mu = zeros(K,1);
%     for k=1:K % 对每个基站迭代寻找最优\mu
%         mu_min = 0;
%         mu_max = 10;
%         iter = 0;
%         while(1)
%             mu1 = (mu_max+mu_min) / 2;
%             P_tem = 0;
%             for i=1:I % 计算功率和
%                 V_tem = inv((A{k}+mu1*eye(T)))*alpha1(i,k)*((H{i,k,k}')*U{i,k}*W{i,k}); % 公式15
%                 P_tem = P_tem + real(trace(V_tem*V_tem'));
%             end
%             if P_tem > P
%                 mu_min = mu1;
%             else
%                 mu_max = mu1;
%             end
%             iter = iter + 1;
% 
%             if abs(mu_max - mu_min) < 1e-5 || iter > max_iter
%                 break
%             end
%         end
%         mu(k) = mu1;
%     end
% 
%     for i=1:I
%         for k=1:K
%             V{i,k} =  inv((A{k}+mu(k)*eye(T)))*alpha1(i,k)*((H{i,k,k}')*U{i,k}*W{i,k}); % 公式15
%         end 
%     end 
% end