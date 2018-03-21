function [mixed, abf] = getSynData(A, type, win, s, b, pure)

%
% Generate synthetic data.
% The spectra of simulated data is obtained from the USGS library "signatures"
%
% Input
%   - A: matrix of reflectances
%   - type: the type of synthetic data
%   - win: block size
%   - s: signal space
%   - b: background space
%   - pure: 0 - no pure pixels, 1 - exist pure pixel
%
% Output
%   - mixed: generated synthetic mixed data
%   - abf: actual abundance fractions
%
% The pure pixels can be removed by adding the following two lines
%        ----Index = ceil(find(abf>0.8)/c);
%        ----abf(:,Index) = 1/c*ones(c,1)*ones(1,length(Index));
%


[band, c] = size(A);
dim = 64;

switch type
    
    case 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % syndata 1: the image is divided into 9 blocks and 4 end-members are used.
        % There exist pure pixels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('The used synthetic data is type 1: block distribution!');
        abf = zeros(floor(dim/3),floor(dim/3));
        maxd = ceil(sqrt((1-dim/3)^2+(1-dim/3)^2));
        for i=1:floor(dim/3)*2-1
            for j=1:floor(dim/3)*2-1
                dist = sqrt((i-dim/3)^2+(j-dim/3)^2);
                dist = (maxd-dist)/maxd;
                abf(i,j) = dist;
            end
        end
        
        [M,N] = size(abf);
        abf1 = [abf zeros(M,dim-N)];
        abf1 = [abf1; zeros(dim-M,dim)];
        abf2 = [zeros(M,dim-N) abf];
        abf2 = [abf2; zeros(dim-M,dim)];
        abf3 = [zeros(dim-M,N); abf];
        abf3 = [abf3 zeros(dim,dim-N)];
        abf4 = [zeros(M,dim-N) abf];
        abf4 = [zeros(dim-M,dim);abf4];
        total = abf1+abf2+abf3+abf4;
        abf1 = reshape(abf1./total,1,dim*dim);
        abf2 = reshape(abf2./total,1,dim*dim);
        abf3 = reshape(abf3./total,1,dim*dim);
        abf4 = reshape(abf4./total,1,dim*dim);
        abf = [abf1;abf2;abf3;abf4];
        clear abf1 abf2 abf3 abf4 total;
        
        %         Index = ceil(find(abf>0.8)/c);
        %         abf(:,Index) = 1/c*ones(c,1)*ones(1,length(Index));
        mixed = reshape((A*abf)',dim,dim,band);
        
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % syndata 2: 4 end-members
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('The used synthetic data is type 2: circular shape!');
        abf = zeros(dim,dim);
        
        % 		maxd = ceil(sqrt(dim^2+dim^2));
        % 		for i=1:dim
        % 		    for j=1:dim
        % 		        dist = sqrt(i^2+j^2);
        % 		        dist = (maxd-dist)/maxd;
        % 		        abf(i,j) = dist;
        % 		    end
        % 		end
        
        %%%%% there exist pure pixels
        for i=1:dim
            for j=1:dim
                dist = sqrt(i^2+j^2);
                if dist < dim
                    dist = (dim-dist)/dim;
                    abf(i,j) = dist;
                end
            end
        end
        
        abf1 = abf;
        abf2 = abf(end:-1:1,:);
        abf3 = abf(:,end:-1:1);
        abf4 = abf(end:-1:1,end:-1:1);
        total = abf1+abf2+abf3+abf4;
        abf1 = reshape(abf1./total,1,dim*dim);
        abf2 = reshape(abf2./total,1,dim*dim);
        abf3 = reshape(abf3./total,1,dim*dim);
        abf4 = reshape(abf4./total,1,dim*dim);
        abf = [abf1;abf2;abf3;abf4];
        clear abf1 abf2 abf3 abf4 total;
        
        Index = ceil(find(abf>0.8)/c);
        abf(:,Index) = 1/c*ones(c,1)*ones(1,length(Index));
        mixed = reshape((A*abf)',dim,dim,band);
        
        
    case 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % syndata 3: randomly generate abundance fractions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('The used synthetic data is type 3: randomly generate abundance!');
        
        % uniform distribution, the mean is 1/c
        abf = rand(c,dim*dim);
        abf = abf./(ones(c,1)*sum(abf,1));
        mixed = reshape((A*abf)',dim,dim,band);
        
        %         % nonuniform distribution
        %         abf = zeros(c,dim*dim);
        % 		tmps = ones(1,dim*dim);
        % 		ss = 0;
        % 		for i=1:c-1
        %             s = rand(1,dim*dim).*tmps;
        %             abf(i,:) = s;
        %             ss = ss+s;
        %             tmps = 1-ss;
        % 		end
        % 		abf(c,:) = 1-ss;
        %         mixed = reshape((A*abf)',dim,dim,band);
        
    case 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % syndata 4: there exist anomaly pixels arrayed in a 7-by-7 array, 
        % 5 end-members are used.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('The used synthetic data is type 4: there exist anomoly pixels!');
        abf = rand(2,dim*dim);
        abf = abf./(ones(2,1)*sum(abf,1));
        abf = [abf;zeros(1,dim*dim)];
        
        abf = reshape(abf',dim,dim,3);
        row = floor(dim/8);
        col = floor(dim/8);
        per = 0.8;
        for i=row-4:row:dim
            abf(i,col-4:col:dim,1:2) = abf(i,col-4:col:dim,1:2)*(1-per);
            abf(i,col-4:col:dim,3) = per;
            per = per-0.1;
        end
        
        abf = reshape(abf,dim*dim,3)';
        mixed = reshape((A*abf)',dim,dim,band);
        
        
    case 5
        disp('The used synthetic data is type 5: small blocks with pure classes!');
        label = ones((dim/8)^2,1);
        num = floor(length(label)/c);
        
        for i=1:c-1
            label((i-1)*num+1:i*num) = (i+1); 
        end
        
        ridx = randperm(length(label));
        label = label(ridx)';
        label = reshape(label,dim/8,dim/8);
        abf = zeros(dim,dim,c);
        img = zeros(dim,dim);
        for i=1:dim
            for j=1:dim
                for cls = 1:c
                    if label(floor((i-1)/8)+1,floor((j-1)/8)+1) == cls
                        tmp = zeros(c,1);
                        tmp(cls) = 1;
                        abf(i,j,:) = tmp;
                        img(i,j) = c;
                    end
                end
            end
        end
        
        %low pass filter
        H = ones(win,win)/(win*win);
        img_fil = filter2(H,img);
        for i=1:c
            abf(:,:,i) = filter2(H,abf(:,:,i));
        end
        abf = abf(ceil(win/2):end-floor(win/2),ceil(win/2):end-floor(win/2),:);
        
        
        % generate mixtures
        [M,N,c] = size(abf);
        abf = reshape(abf,M*N,c)';
        
        % remove pure pixels
        if pure == 0
            Index = ceil(find(abf>0.8)/c);
            abf(:,Index) = 1/c*ones(c,1)*ones(1,length(Index));
        end
        
        mixed = reshape((A*abf)',M,N,band);
        
    case 6
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % syndata 6: subspace signals
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp('The used synthetic data is type 6: subspace signal!');
        abf = rand(b,dim*dim);
        abf = abf./(ones(b,1)*sum(abf,1));
        abf = [abf;zeros(s,dim*dim)];
        
        abf = reshape(abf',dim,dim,s+b);
        row = floor(dim/8);
        col = floor(dim/8);
        per = 0.8;
        for i=row-4:row:dim
            abf(i,col-4:col:dim,1:b) = abf(i,col-4:col:dim,1:b)*(1-per);
            tmp = rand(s,1);
            tmp = tmp./sum(tmp);
            abf(i,col-4:col:dim,b+1:end) = ones(8,1)*per*tmp';
            per = per-0.1;
        end
        
        abf = reshape(abf,dim*dim,s+b)';
        mixed = reshape((A*abf)',dim,dim,band);  
        
    case 7 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % simulate random probe distribution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        abf = zeros(dim, dim, c);
        [M,N,c] = size(abf);
        fidx_x = [];
        fidx_y = [];
        for i=1:c
            for j=1:200
                idx_x = ceil(rand*(dim-3))+2;
                idx_y = ceil(rand*(dim-3))+2;
                abf(idx_x-1:idx_x+1, idx_y-1:idx_y+1, i) = rand(3,3);
            end
        end       
        
        abf = reshape(abf, dim*dim, c);
        abf = abf';
        idx = find(sum(abf,1)>0); 
        abf(:,idx) = abf(:,idx)./(ones(c,1)*sum(abf(:,idx),1));
        mixed = reshape((A*abf)',dim,dim,band);  
        
    otherwise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % random distribution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        abf = rand(c,dim*dim);
        abf = abf./(ones(c,1)*sum(abf,1));
        Index = ceil(find(abf>0.8)/c);
        abf(:,Index) = 1/c*ones(c,1)*ones(1,length(Index));
        mixed = reshape((A*abf)',dim,dim,band);
        
end
