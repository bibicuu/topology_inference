function [A] = generate_graph(params)  
% WEIGHTED NOW
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
% g_type => Graph type
% params.N => Number of nodes
% params.p => Case g_type='ER'
% params.symmetric => Case g_type = 'ER'

% params.m => Case g_type='BA'
% params.connected => Case g_type='RBF'
% params.K => Case g_type='RING' or g_type='SW'
% params.beta => Case g_type='SW'


%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
% A => Adjacency matrix
	A = [];
    N = params.N;
	g_type = params.g_type;
    
    
	switch g_type
        case 'ER' % Erdos-Renyi graph
            p = params.p; % probability to connect an edge (indep of the rest)
            symmetric = params.symmetric;
            A = generate_connected_ER(N, p, symmetric);
			weights = zeros(N); 
			if (params.weighted == 0)
				[V, eignv] = eig(A);
				post_eig = eignv(eignv > 0);
				max_eig = max(post_eig);
				min_eig = min(post_eig);

				if (max_eig >= 1) && (params.scale_graph > 1)
					A_scaled = A ./ (params.scale_graph * max_eig);
					A = A_scaled;
				elseif (min_eig < 1) && (params.scale_graph < 1)
					if (min_eig > 1e-3)
						A_scaled = A ./ (params.scale_graph * min_eig);
						A = A_scaled;
					else
						eignv(eignv > 0 & eignv < 1e-3) = 0;
						post_eig = eignv(eignv > 0);
						min_eig = min(post_eig);
						A = V * eignv * inv(V);
						A_scaled = A ./ (params.scale_graph * min_eig);
						A = A_scaled;
					end

				else
					disp('Graph scale not set properly')
				end
		
			else
				if symmetric == 1
    				for i = 1:N
        				for j = i+1:N
            				if A(i,j) == 1
                				weight_value = 0.5 + (rand * 0.5); % Valores U([0.5, 1])
                				weights(i,j) = weight_value;
                				weights(j,i) = weight_value;
            				end
        				end
    				end
				else
    				% Si no es simétrico, asignar de manera normal
    				weights(A == 1) = 0.5 + (rand(nnz(A == 1), 1) * 0.5); % Valores U([0.5, 1])
				end
				A = params.scale_graph .* weights;
			end

			
			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
        case 'TREE' % Tree graph
            A = preferential_attachment_graph(N,1);
        case 'BA' % Barabasi Albert graph => Prefer a connection if a node has already many connections
            m = params.m;
            %A = scale_free(N,m,m0);
            A = preferential_attachment_graph(N,m);
        case 'RBF' % Radial based graph
			max_tries = 10;
            connected = 1;%params.connected;
            T = 0.75;
            s = 0.5;
            A = gaussian_graph(N,T,s);
			% connected = 0;
            if connected %generate connected graph
                k = 1;
                lambdas = eig(diag(sum(A))-A);
                while((abs(lambdas(2)) < 1e-6) && (k < max_tries))
                    A = gaussian_graph(N,T,s);
                    k = k+1;
                    lambdas = eig(diag(sum(A))-A);
                end
                if k >= max_tries
                    disp('ERR: Generated unconnected RBF graph')
                end
            end
        case 'RING' % Ring graph
            K = params.K;
            A = small_world(N,K,0);
        case 'SW' % Small world graph => Arrive in a few steps anywhere
            K = params.K;
            beta = 0.1;
            A = small_world(N,K,beta);

        case 'SBM' % Stochastic 
            disp('NOT AVAILABLE YET')
        otherwise
            disp('ERR: Unkown graph type')
            return
	end
end    
