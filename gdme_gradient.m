
function s = gdme_gradient(x,A,SNR);

% Gradient descent maximum entropy using steepest descent

[BandNum,ClsNum] = size(A);

%Initialize
s_old = zeros(ClsNum,1);
s = 1/ClsNum*ones(ClsNum,1);
mu = zeros(BandNum,1);
dmu = zeros(BandNum,1);
err = 100;

k = 1;
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
    th2 = 150; 
    eta = 0.15; 
end

while norm(s-s_old)>th1 
    s_old = s;
    ss(:,k) = s;
    
    % stop is one of s components is zero
    idx = find(s==0);
    if ~isempty(idx)
        break;
    end
    
    % Update of Lagrange multiplier
    dmu = eta*(x-A*s);
    mu = mu-dmu;       
    
    % Calculate s 
    mu0 = log(sum(exp(-A'*mu)));
    s = exp(-A'*mu-mu0);
      
    k = k+1;
    if k>th2
        break;
    end
   
end


