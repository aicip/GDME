
function [Aest, sest, time] = gdme(mixed, SNR, algo, c, A)

% Gradient descent maxent (GDME) algorithm
% 
% mixed - mixture data L-by-N, L: number of bands, N: number of pixes
% SNR   - signal to noise ratio
% algo  - type of learning algorithm
%       - 1: lagrange steepest descent, 2: variation of Newton's method
% c     - Estimated number of endmembers
% A     - actural endmembers, for testing abundance estimation only
%
%
% sest  - estimated abundances
% Aest  - estimated endmembers
% time  - running time
%



t0 = cputime;
[L, N] = size(mixed);

% Denoise using svd
[UU, SS, VV] = svds(mixed*mixed'/N,c);
Lowmixed = UU'*mixed;
mixed = UU*Lowmixed;

% Calculate the first endmember signature
[tmp, idx] = max(sum(mixed.*mixed,1));
a0 = mixed(:,idx);
EndIdx = idx;

% Find another pixel which is the furthest away from a0
merr = [];
for i=1:N
    x = mixed(:,i);
    tmp = sum((x-a0).*(x-a0));
    merr = [merr tmp];
end
[tmp idx] = sort(merr);

% Main loop
Aest = a0;
mlse(1) = 1000;
R = 150;
SNR_th = 30 + 10*log10(c);
k = 1;
if SNR < 25
    theta = 6/180*pi;
else
    theta = 2/180*pi;
end

while k<c
    k = k+1
    
    % If the snr high, use single pure pixel as endmember estimation;
    if SNR > SNR_th
        Aest = [Aest mixed(:,idx(end))];
        EndIdx = [EndIdx idx(end)];
        
        % otherwise, use the average of set of pure pixels 
    else 
        a1 = mixed(:,idx(end));
        EndIdx = [EndIdx idx(end)];
        ta = a1;
        for i=length(idx)-1:-1:length(idx)-R
            x = mixed(:,idx(i));
            angle = acos(x'*a1/(norm(x)*norm(a1)+1e-9));
            if angle<theta
                ta = [ta x];
                EndIdx = [EndIdx idx(i)];
            end
        end
        Aest = [Aest mean(ta,2)];
    end
    
    % Use gdme to unmix all image pixels
    sest = zeros(k, N);
    switch algo
        case 1 % sdme
            for i=1:N
                sest(:,i) = gdme_gradient(mixed(:,i), Aest, SNR);
            end
            
        case 2 % variational newton
            for i=1:N
                sest(:,i) = gdme_var_newton(mixed(:,i), Aest, SNR);
            end
    end
    
    % Find the largest lse   
    lserr = sum((mixed-Aest*sest).*(mixed-Aest*sest),1); % LSE measure
    [tlse, idx] = sort(lserr);
end

time = cputime-t0;
