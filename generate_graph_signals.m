function [X, C, snr, snr_grafo, noise] = generate_graph_signals(S, params, graph_params)
% Only SEM signals
    
    N = size(S,1); % Number of nodes
    M = params.M; % Number of samples to generate
    cov_matrix = params.cov_matrix; % std of the noise
    h_mu = zeros(N,1); % Means of the multivariant distribution
	% Generate signals
    noise = mvnrnd(h_mu, cov_matrix, M)';

	% 
	% if (graph_params.weighted == 0) % Fix matrix to avoid singularity
	% 	[V, eignv] = eig(S);
	% 	post_eig = eignv(eignv > 0);
	% 	max_eig = max(post_eig);
	% 	min_eig = min(post_eig);
	% 
	% 	if (max_eig > 1) && (graph_params.scale_graph > 1)
		% 	S_scaled = S ./ (graph_params.scale_graph * max_eig);
		% 	S = S_scaled;
	% 	elseif (min_eig < 1) && (graph_params.scale_graph < 1)
		% 	if (min_eig > 1e-3)
			% 	S_scaled = S ./ (graph_params.scale_graph * min_eig);
			% 	S = S_scaled;
		% 	else
			% 	eignv(eignv > 0 & eignv < 1e-3) = 0;
			% 	post_eig = eignv(eignv > 0);
			% 	min_eig = min(post_eig);
			% 	S = V * eignv * inv(V);
			% 	S_scaled = S ./ (graph_params.scale_graph * min_eig);
			% 	S = S_scaled;
		% 	end
	% 
	% 	else
		% 	disp('Graph scale not set properly')
	% 	end
	% 
	% end






	% S = S ./ max(S(:));
    X = (eye(N)-S) \ noise; % Signals in the graph according to SEM
	if (graph_params.weighted == 1)
		norms = sqrt(sum(X.^2));  % Esto calcula la norma de cada columna
		X_normalized = X ./ norms;
		X = X_normalized;
	end


	C = X*X'/M; % Covariance matrix of X
	% snr
	p_x = norm(X, 'fro')^2/M; % Power of the signal generated
	p_x_minusnoise = norm(X - N, 'fro')^2/M;
	p_n = norm(noise, 'fro')^2/M; % Noise power
	snr = p_x/p_n;
	snr_grafo = p_x_minusnoise/p_n;
        
end
