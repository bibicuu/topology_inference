function [S, exec_time, term1, term2] = homocedastic_Cov_sqrt_Estimate(lambda_scale, C, w, M)  

	tic
	N = size(C,1);
	lambda = lambda_scale .* w .* sqrt(log(N)/M);
	
	
	%% Optimization problem
	cvx_begin quiet
	cvx_solver mosek;
	
	variable S(N,N) symmetric nonnegative % Non directed graph

	
	frob_term = norm((eye(N) - S) * sqrtm(C), 'fro');
	
	term1 = 1/(2*w) * square_pos(frob_term);
	term2 = norm(S(:), 1);
	obj = term1 + lambda .* term2;
	fprintf('Lambda true %d\n', lambda);
	minimize(obj)  
	% Constraints
	subject to
    	diag(S) == 0;
	cvx_end
	
	exec_time = toc;


end
