function generate_signals(signal_params, nr_signals, graph_params, nr_graphs)
% signal<ID_signal>_N<N>_M<M>_w<w>.mat
% sprintf('signal%d_N%d_M%d_w%d.mat', ID_signal, N, M, w)

	for ID_graph = 1:nr_graphs
		% Prepare data to save
		X_cell = cell(1, nr_signals);
		C_cell = cell(1, nr_signals);
		snr_cell = cell(1, nr_signals);
		snr_grafo_cell = cell(1, nr_signals);
		noise_cell = cell(1, nr_signals);

		% Create data
		N = graph_params.N;
		g_type = graph_params.g_type;
		prob = graph_params.p;
		symmetric = graph_params.symmetric;
		scale_graph = graph_params.scale_graph;
		weighted = graph_params.weighted;
		A = load(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d.mat', ...
				ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_graph, N, scale_graph, g_type, prob, symmetric, weighted)).A;
		% Data for signals
		M = signal_params.M;
		w = signal_params.w;
		parfor ID_signal = 1:nr_signals
			if ~exist(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/signal%d_N%d_M%d_w%d.mat', ...
				ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_signal, N, M, w), 'file')
				[X_cell{ID_signal}, C_cell{ID_signal}, snr_cell{ID_signal}, snr_grafo_cell{ID_signal}, noise_cell{ID_signal}] = ... 
				generate_graph_signals(A, signal_params, graph_params);
			else
				X_cell{ID_signal} = [];
				C_cell{ID_signal} = [];
				snr_cell{ID_signal} = [];
				snr_grafo_cell{ID_signal} = [];
				noise_cell{ID_signal} = [];
			end
		end

		% Save data
		for ID_signal = 1:nr_signals
			if ~isempty(X_cell{ID_signal})
				mkdir(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/', ...
				ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w));
				filename = sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/signal%d_N%d_M%d_w%d.mat', ...
				ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_signal, N, M, w);
				X = X_cell{ID_signal};
				C = C_cell{ID_signal};
				Noise = noise_cell{ID_signal};
				snr = snr_cell{ID_signal};
				snr_grafo = snr_grafo_cell{ID_signal};
				save(filename, 'X', 'Noise', 'A', 'snr', 'snr_grafo', 'C', 'graph_params', 'signal_params');
			end
		end

	end
end
