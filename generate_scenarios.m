function generate_scenarios(scenarios_params, nr_scenarios, signal_params, nr_signals, graph_params, nr_graphs)


matrix_random  = [
     7     4    11     7;
     6    12     5     1;
     3    14    12    13;
    20    13     8    11;
     2    19    10     9;
    15    16     9    17;
     7     5    19     4;
    14    13    19    16;
    13    15    17    10;
     3    18     3     9;
     9    10    17    14;
     6    18     6    20;
     3    14    11    20;
    11     3     3    20;
     7    19    10    11;
     2     1     7    14;
     2     4    13    17;
     9    13     8     5;
    15     5    13    14;
    17    14    18    11;
     9    12    18     5;
    14    16    11     9;
     9     6    13    13;
     5     5     7    14;
    15    19     2    14;
     9     3     1    19;
    19     5    13    13;
    13    12    10     2;
     8    14     1     9;
     8     6    12    10;
];

	type_of_hidden = scenarios_params.type_of_hidden;
	H = scenarios_params.H;

	% Porque solo puedo generar varios escenarios si estoy en random
	if strcmp(type_of_hidden, 'min') || strcmp(type_of_hidden, 'max')
		nr_scenarios = 1;
	end


	for ID_graph = 1:nr_graphs
		% Load graph
		N = graph_params.N;
		g_type = graph_params.g_type;
		prob = graph_params.p;
		symmetric = graph_params.symmetric;
		scale_graph = graph_params.scale_graph;
		weighted = graph_params.weighted;
		A_true = load(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d.mat', ...
				ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_graph, N, scale_graph,g_type, prob, symmetric, weighted)).A;

			
		for ID_signal = 1:nr_signals
			% Load signals
			M = signal_params.M;
			w = signal_params.w;
			data_signal = load(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/signal%d_N%d_M%d_w%d.mat', ...
				ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_signal, N, M, w));
			X_true = data_signal.X;
			C_true = data_signal.C;
			Noise = data_signal.Noise;
			% Load data
			A_scenario_cell = cell(1,nr_scenarios);
			X_scenario_cell = cell(1,nr_scenarios); 
			C_scenario_cell = cell(1,nr_scenarios);
			s_n_cell = cell(1,nr_scenarios);
			s_h_cell = cell(1,nr_scenarios);
			noise_power_cell = cell(1,nr_scenarios);
			degrees_cell = cell(1,nr_scenarios);

			parfor ID_scenario = 1:nr_scenarios
				if ~exist(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/scenario%d_H%d_%s.mat', ...
						ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ID_scenario, H, type_of_hidden), 'file')	
					% Make scenario
					O = N - H;

					if strcmp(type_of_hidden, 'rand') && (H>0)
						s_h = matrix_random(ID_graph, 1:H);
						s_n = setdiff(1:20, s_h);
					else
					[s_n, s_h] = select_hidden_nodes(type_of_hidden, O, A_true);
					end
					s_n_cell{ID_scenario} = s_n;
					s_h_cell{ID_scenario} = s_h;
					X_scenario_cell{ID_scenario} = X_true(s_n, :);
					C_scenario_cell{ID_scenario} = C_true(s_n, s_n);
					A_scenario_cell{ID_scenario} = full(A_true(s_n,s_n));
					noise_power_cell{ID_scenario} = sum(Noise.^2, 2) / size(Noise, 2);
					degrees_cell{ID_scenario} = sum(A_true);

				else
					A_scenario_cell{ID_scenario} = [];
					X_scenario_cell{ID_scenario} = []; 
					C_scenario_cell{ID_scenario} = [];
					s_n_cell{ID_scenario} = [];
					s_h_cell{ID_scenario} = [];
					noise_power_cell{ID_scenario} = [];
					degrees_cell{ID_scenario} = [];
				end
			end
		
		
		
			for ID_scenario = 1:nr_scenarios
				if ~isempty(A_scenario_cell{ID_scenario})
					mkdir(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/', ...
						ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden))
					filename = sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/scenario%d_H%d_%s.mat', ...
						ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ID_scenario, H, type_of_hidden);
					A = A_scenario_cell{ID_scenario};
					X = X_scenario_cell{ID_scenario}; 	
					% Normalizar cada columna dividiendo por su respectiva norma
					

					C = C_scenario_cell{ID_scenario};
					s_n = s_n_cell{ID_scenario};
					s_h = s_h_cell{ID_scenario};
					degrees = degrees_cell{ID_scenario};
					noise_power = noise_power_cell{ID_scenario};
					
					save(filename, 'A', 'X', 'C', 's_n', 's_h', 'graph_params', 'signal_params', 'scenarios_params', 'degrees', 'noise_power');
				end
			end
		end
	end

	
end



