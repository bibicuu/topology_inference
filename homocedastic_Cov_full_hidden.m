function [S_O, exec_time, K, term_fro, term_norm1, term_normnuc, term_normlasso] = homocedastic_Cov_full_hidden(lambda_scale, C_O, w, alpha, beta, M)   

% Formulación 1 de lo que hablé con Antonio el 02/09 = > Agrupar S_OH*C_HO,
% sin alternating, low rank y row sparse

    %% Start of the function
    
    tic
    O = size(C_O,1);
    lambda = lambda_scale .* w .* sqrt(log(O)/M); 
	beta_tmp = beta .* w .* sqrt(log(O)/M);
	beta = beta_tmp;


    %% Optimization problem

	%% Resolver S_O y K
    cvx_begin quiet
    cvx_solver mosek;
    
    % Variables
    variable S_O(O,O) symmetric nonnegative
    variable K(O,O) 

    frobenius_term = norm(((eye(O) - S_O) * sqrtm(C_O) - K * inv(sqrtm(C_O)) ), 'fro');
    K_row_norms = sum(norms(K, 2, 2));
	

	term_fro = 1/(2*w) * square_pos(frobenius_term);
	term_norm1 = norm(S_O(:), 1);
	term_normnuc = norm_nuc(K);
	term_normlasso = sum(K_row_norms);
	
    obj = term_fro + lambda * term_norm1 + ...
         beta * alpha * term_normnuc +  beta * (1-alpha) * term_normlasso;
            
   
    minimize(obj)  

    
    % Constraints
    subject to
        diag(S_O) == 0;
		% K(1:13) == 0;
		% K(15:end) == 0;
    cvx_end

    exec_time = toc;
%    fprintf('Beta %d, Term frob %d \n\n\n', beta, term_fro)
% figure; imagesc(K); colorbar; title(sprintf('Alpha %d, beta %d', alpha, beta))
end
