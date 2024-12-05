function [metrics] = graph_learning_perf_eval(A,S)
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
% A => Groundtrue adjacency
% S => Estimated adjacency

%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
% [precision,recall,f,NMI,num_of_edges, err_amb, err_no_amb, new_err_amb, new_err_no_amb]
% metrics => Includes precision, recall, f, NMI, num_of_edges, err_amb, err_no_amb, err_norm1, err_noamb_norm1

n_options = 9;

precision = zeros(1,n_options);
recall = zeros(1,n_options);
f = zeros(1,n_options);
NMI = zeros(1,n_options);
num_of_edges_true = zeros(1,n_options);
num_of_edges_est = zeros(1,n_options);
shamming_dist = zeros(1,n_options);
TPR = zeros(1, n_options);
FDR = zeros(1, n_options);

for option = 1:n_options
    % To binarize Laplacian
    % L_0 = diag(sum(A,2)) - A; % Real Laplacian (should be binarized) 
    % L = mbinarize(abs(diag(sum(S,2)) - S),option, A); % Estimated Laplacian binarized
    
	% To binarize GSO (Adjacency)
	L_0 = A>0;
	L = mbinarize(S, option);

    shamming_dist(option) = sum(sum((L_0 ~= 0) ~= (L ~= 0)));
    
    
    %% If there's no graph
    % if sum(isnan(L))>0
    %     metrics = NaN;
    %     % precision = NaN;
    %     % recall = NaN;
    %     % f = NaN;
    %     % NMI = NaN;
    %     % num_of_edges_true = NaN;
    %     % num_of_edges_est = NaN;
    %     return
	% end
	%% Calculate TPR y FDR
	TP = sum((L == 1) & (L_0 == 1), 'all');  % Verdaderos positivos
	FN = sum((L == 0) & (L_0 == 1), 'all');  % Falsos negativos
	FP = sum((L == 1) & (L_0 == 0), 'all');  % Falsos positivos
	
	
	% Calcular TPR (True Positive Ratio) o Recall
	TPR(option) = TP / (TP + FN); % NaN solo si no hay grafo real
	
	% Calcular FDR (False Discovery Ratio)
	FDR(option) = FP / (TP + FP); % NaN si no he estimado nada

    %% Calculate num_of_edges true
    edges_groundtruth = squareform(L_0-diag(diag(L_0)))~=0; % Binary L_0 to see if there's an edge => Adjacency binary (ground truth)
    num_of_edges_true(option) = sum(edges_groundtruth); % total nr of edges in the learned matrix

    %% Calculate num_of_edges estimated
    edges_learned = squareform(L-diag(diag(L)))~=0; % Binary L_0 to see if there's an edge => Adjacency binary (learned)
    num_of_edges_est(option) = sum(edges_learned); % total nr of edges in the learned matrix

    %% Calculate fscores
    [precision(option),recall(option)] = perfcurve(double(edges_groundtruth),double(edges_learned),1,'Tvals',1,'xCrit','prec','yCrit','reca');

    if precision(option) == 0 && recall(option) == 0 % Avoid NaN in the formula 
        f(option) = 0;
    else
        f(option) = 2*precision(option)*recall(option)/(precision(option)+recall(option)); % Formula to calculate f score
    end
    NMI(option) = perfeval_clus_nmi(double(edges_groundtruth),double(edges_learned)); % Normalized Mutual Information

	if num_of_edges_est(option) == 0
		NMI(option) = NaN;
		precision(option) = NaN;
		recall(option) = NaN;
		fscore(option) = NaN;
		TPR(option) = NaN;
		FDR(option) = NaN;
		shamming_dist(option) = NaN;
	end
end


%% Calculate errors  err_amb, new_err_amb, err_no_amb & new_err_no_amb
A_true = A;
A_est = S;

err_amb = norm(A_true-A_est, 'fro')^2 / norm(A_true,'fro')^2;
new_err_amb = norm(A_true - A_est, 1)^2 / norm(A_true,1)^2;
if (norm(A_est, 'fro') ~= 0)
    err_no_amb = norm((A_true ./ norm(A_true, 'fro')) - (A_est ./ norm(A_est, 'fro')),'fro')^2;
    aux = (A_true ./ norm(A_true(:), 1)) - (A_est ./ norm(A_est(:), 1));
    new_err_no_amb = norm(aux(:), 1)^2;
else
    err_no_amb = norm((A_true ./ norm(A_true, 'fro')),'fro')^2;
    aux = (A_true ./ norm(A_true(:), 1)) ;
    new_err_no_amb = norm(aux(:), 1)^2;
    precision = 'CERO'; % Actuar como flag para ver que penalizacion demasiado grande
end

% metrics => Includes precision, recall, f, NMI, num_of_edges, err_amb, err_no_amb, err_norm1, err_noamb_norm1

metrics.fscore = f; 
metrics.precision = precision;
metrics.recall = recall;
metrics.NMI = NMI;
metrics.num_of_edges_true = num_of_edges_true;
metrics.num_of_edges_est = num_of_edges_est;
metrics.shamming_dist = shamming_dist;
metrics.TPR = TPR;
metrics.FDR = FDR;


metrics.err_amb = err_amb;
metrics.new_err_amb = new_err_amb;
metrics.err_no_amb = err_no_amb;
metrics.new_no_err_amb = new_err_no_amb;
