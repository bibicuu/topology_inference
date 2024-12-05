function [S_out] = mbinarize(S_in,option)
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
% S_in => Input matrix (GSO) to binarize
% option => Option to binarize. 


%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
% S_out => Binarized matrix
% figure(40+option)
% subplot(211)
% imagesc(S_in); colorbar;
% title('S_{in}')


    O = size(S_in,1);
    if sum(sum(isnan(S_in))) ~= 0 % If any NaN => Input is all 1
       S_in = ones(size(S_in)); 
    end
    if option == 1   % Below/above the mean of the min&max value
        th =(max(S_in(:))-min(S_in(:)))/2; % Mean of the min&max value
        S_in(S_in>=th)=1;
        S_in(S_in < th)=0; 
        S_out = S_in;
		

    elseif option == 2 % k-means clustering in 2 groups
        [idx,~] = kmeans(S_in(:),2);
        if idx(1)==2
            idx(idx==2) = 0;
        else
            idx(idx==1) = 0;
            idx(idx==2) = 1;
        end
        S_out = reshape(idx, [O O]); 
		
		
    elseif option == 3   % Below/above a fixed threshold
        th = 0.1;
        S_in(S_in>=th)=1;
        S_in(S_in < th)=0;
        S_out = S_in;

    elseif option == 4 % Histogram of S_in finds the range with more values 
        figure(50)
        aux = histogram(S_in);
        valores = aux.Values;
        intervalo = aux.BinEdges;
        [~,v2]= max(valores); 
        th_inf = intervalo(v2);
        th_sup = intervalo(v2+1);

        lgc = logical((S_in>=th_inf).*(S_in<=th_sup));
        S_in(lgc) = 0;
        S_in(~lgc) = 1;
        S_out = S_in;
        close(50)
		
    elseif option == 5 % Mean of the values 
        th = mean(S_in(:)); 
        S_in(S_in>=th)=1;
        S_in(S_in < th)=0; 
        S_out = S_in;
		
    elseif option == 6
        th = 0.8*mean(S_in(:)); 
        S_in(S_in>=th)=1;
        S_in(S_in < th)=0; 
        S_out = S_in;
		
    elseif option == 7
        th = 0.5 * mean(S_in(:)); 
        S_in(S_in>=th)=1;
        S_in(S_in < th)=0; 
        S_out = S_in;
    elseif option == 8 
        figure(50)
        aux = histogram(S_in);
        valores = aux.Values;
        intervalo = aux.BinEdges;
        [~,v2]= max(valores); 
        th_inf = intervalo(v2);
        th_sup = intervalo(v2+1);

        lgc = logical((S_in>=th_inf));
        S_in(lgc) = 0;
        S_in(~lgc) = 1;
        S_out = S_in;
        close(50)
    elseif option == 9
        th = 0.5 * max(S_in(:)); 
        S_in(S_in>=th)=1;
        S_in(S_in < th)=0; 
        S_out = S_in;

    end


end
