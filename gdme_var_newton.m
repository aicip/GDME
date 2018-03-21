
function s = gdme_var_newton(x,A,SNR);

% Gradient descent maximum entropy using variational Newton's method

delta = 20;
x = [x; delta];
A = [A; ones(1,size(A,2))*delta];
[BandNum,ClsNum] = size(A);

%Initialize
s_old = zeros(ClsNum,1);
s = 1/ClsNum*ones(ClsNum,1);
lambda = zeros(BandNum,1);

if SNR == 5
    th1 = 1e-9; % error tolerance
    th2 = 50;  % maximum iteration number 
    eta = 0.05; % learning rate
elseif SNR == 10
    th1 = 1e-9;
    th2 = 80;
    eta = 0.08;
elseif SNR == 15
    th1 = 1e-9;
    th2 = 100;
    eta = 0.1;
else
    th1 = 1e-9; 
    th2 = 50; 
    eta = 0.15; 
end

ck = 1;
k = 1;
Nk = A';

while norm(s-s_old)>th1 
    
    ck = ck*10;
    
    s_old = s;
    ss(:,k) = s;
    
    % stop is one of s components is zero
    idx = find(s==0);
    if ~isempty(idx)
        break;
    end
    
    hs = A*s-x;
    Hk = diag(1./s);
    grad = 1+log(s);
    gradLc = grad+A'*(lambda+ck*hs);
        
%     % Update of s
%     % It is better when no pure pixels present in the scene
%     h = -inv(Hk+ck*Nk*Nk')*gradLc;
%     beta = min(1, 0.9*min(s./abs(h)));
%     s = s + beta*h;
    
    % it is better when there exist pure pixels
    s = max(s - inv(Hk+ck*Nk*Nk')*gradLc, eps);
    
    % Update of lambda
    lambda = lambda+ck*(hs+Nk'*(s-s_old));  

    
    k = k+1;
    if k>th2
        break;
    end
   
end


