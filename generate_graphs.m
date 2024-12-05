function generate_graphs(graph_params, nr_graphs)
% Este codigo generará todos los grafos (cada hilo hará 1). Luego se guardarán 

% /simulations/graph<ID_graph>_N<N>_sc<scale>_<type_graph>_p<.2f>_<sim>

	N = graph_params.N;
	g_type = graph_params.g_type;
	% if strcmp(graph_params, )
	prob = graph_params.p;
	symmetric = graph_params.symmetric;
	scale_graph = graph_params.scale_graph;
	weighted = graph_params.weighted;

	graph_params.symmetric = strcmp(symmetric, 'sim');

	A_cell = cell(1, nr_graphs);

	% Create graphs
	parfor ID_graph = 1:nr_graphs
		if ~exist(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d.mat', ...
				ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_graph, N, scale_graph, g_type, prob, symmetric, weighted), 'file')
			[A_cell{ID_graph}] = generate_graph(graph_params);
		else
			A_cell{ID_graph} = [];
		end
	end

	% Save graphs
	for ID_graph = 1:nr_graphs
		if ~isempty(A_cell{ID_graph})
			mkdir(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d', ID_graph, N, scale_graph, g_type, prob, symmetric, weighted))
			filename = sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d.mat', ...
				ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_graph, N, scale_graph, g_type, prob, symmetric, weighted);
    		A = A_cell{ID_graph};
    		save(filename, 'A', 'graph_params');
		end
	end
	
end
