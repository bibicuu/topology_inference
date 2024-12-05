
function normal_executions(M_list, H_list, lambda_scale_list, beta_list, alpha_list, omega, wei, sim, type_hid, tgraph)
% clear all; close all; clc
%% 
% M_list = [40, 1e3];
% H_list = [1,2,3];
% lambda_scale_list
% beta_list
% alpha_list

% Generate graphs
nr_graphs = 20; % Generated randomly ER // diff realizations
graph_params.N = 20;
if tgraph == 1
	graph_params.g_type = 'ER';
	graph_params.p = 0.2;
elseif tgraph == 2
	graph_params.g_type = 'ER';
	graph_params.p = 0.1;
elseif tgraph == 3
	graph_params.g_type = 'BA';
	graph_params.m = 2;
	graph_params.p = 2; % To put m in name
elseif tgraph == 4
	graph_params.g_type = 'BA';
	graph_params.m = 4;
	graph_params.p = 4; % To put m in name
elseif tgraph == 5
	graph_params.g_type = 'RBF';
	graph_params.p = 0; % To put m in name
elseif tgraph == 6
	graph_params.g_type = 'SW';
	graph_params.K = 2;
	graph_params.p = 2;
elseif tgraph == 7
	graph_params.g_type = 'SW';
	graph_params.K = 4;
	graph_params.p = 4;
else
	err('Ya la has liado')
end

if (sim == 1)
graph_params.symmetric = 'sim';

else
graph_params.symmetric = 'nosim';

end

% Es buena idea que sea scale_graph > 3 si no es weighted
graph_params.scale_graph = 1; % If < 1 then scale with min_eig. Else with max_eig
if (wei == 0)
	graph_params.scale_graph = 2;
end
graph_params.weighted = wei; % IF It is 0, scale graph cannot be 1!!!

generate_graphs(graph_params, nr_graphs);
disp('Finished generating graphs')

% Generate data for each graph
nr_signals = 3; % Generated random gaussian // diff realizations (for each graph)
signal_params.w = omega; 
signal_params.M = 1e3; 
signal_params.cov_matrix = signal_params.w * eye(graph_params.N);

for M = M_list
	signal_params.M = M;
generate_signals(signal_params, nr_signals, graph_params, nr_graphs);
disp('Finished generating signals')
end

% Generate hidden scenario for each graph/data
nr_scenarios = 1; % Gen

scenarios_params.type_of_hidden = type_hid; % rand min o max
scenarios_params.H = 2;
for M = M_list
	signal_params.M = M;
for H = H_list
	scenarios_params.H = H;
	generate_scenarios(scenarios_params, nr_scenarios, signal_params, nr_signals, graph_params, nr_graphs);
	disp('Finished generating scenarios')
end
end

%% Solve scenarios
graphs_to_solve = 1:1:nr_graphs; % Graphs to solve. Default all created
signals_to_solve = 1:1:nr_signals; % Signals to solve. Default all created
scenarios_to_solve = 1:1:nr_scenarios; % Scenarios to solve. Default all created

% graphs_to_solve = 1:1:1; % Graphs to solve. Default all created
% signals_to_solve = 1:1:1; % Signals to solve. Default all created
% scenarios_to_solve = 1:1:1; % Scenarios to solve. Default all created


% Data to solve (default created)
N = graph_params.N;
g_type = graph_params.g_type;
prob = graph_params.p;
symmetric = graph_params.symmetric;
weighted = graph_params.weighted;
scale_graph = graph_params.scale_graph;
M = signal_params.M;
w = signal_params.w;
H = scenarios_params.H;
type_of_hidden = scenarios_params.type_of_hidden;

%% ---------------------------------------------------------------------------------------------------------------
%% Solve no hidden for each scenario
% exps_lambda = linspace(-6,6,5);
% lambda_scale_list = [1e-3];

if (strcmp(type_of_hidden, 'min')) || (strcmp(type_of_hidden, 'max'))
	scenarios_to_solve = 1;
end 

for ID_graph = graphs_to_solve
	for ID_signal = signals_to_solve
		for M = M_list % DELETE THIS
		for ID_scenario = scenarios_to_solve
			for H = H_list % DELETE THIS
				% Load graph signal y scenario

				scenario_data = load(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/scenario%d_H%d_%s.mat', ...
						ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ID_scenario, H, type_of_hidden));



				if ~exist(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/', ...
							ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden), 'dir')
					mkdir(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/', ...
							ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden));
				end
				O = N - H;
				C = scenario_data.C;
				A = scenario_data.A; % This is the one that has hidden inside. It's not ground truth
				X = scenario_data.X;
				degrees = scenario_data.degrees;
				noise_power = scenario_data.noise_power;
				s_h = scenario_data.s_h;

				if ~isempty(s_h)
					noise_power_hidden = mat2str(noise_power(s_h),2);
					degrees_hidden = mat2str(degrees(s_h),2);
				else
					noise_power_hidden = [];
					degrees_hidden = [];
				end
				mean_degrees = mean(degrees);


				% frob_A = norm(A, 'fro')^2;
				% frob_X = norm(X, 'fro')^2;
				frob_dif = norm((X-A*X), 'fro')^2/(2*w*M);

				% Create cells to save output for parfor
				metrics_cell = cell(1, length(lambda_scale_list));
				S_cell = cell(1, length(lambda_scale_list));
				term1_cell = cell(1, length(lambda_scale_list));
				term2_cell = cell(1, length(lambda_scale_list));
				exec_time_cell = cell(1, length(lambda_scale_list));
				parfor lambda_scale_idx = 1:length(lambda_scale_list)
					lambda_scale = lambda_scale_list(lambda_scale_idx);
					% if ~exist(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/Covsqrt_Estimate_N%d_M%d_H%d_l%d_w%d.mat', ...
						% 	ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ...
						% 	N, M, H, lambda_scale, w), 'file')
						if (sim == 1)
						[S, exec_time, term1, term2] = homocedastic_Cov_sqrt_Estimate(lambda_scale, C, w, M);  
						else
						[S, exec_time, term1, term2] = homocedastic_Cov_sqrt_Estimate_nosim(lambda_scale, C, w, M);  
						end

						[metrics] = graph_learning_perf_eval(A,S);

						fscores_disp = metrics.fscore([2,4,5,6]);
						fscore_str = mat2str(fscores_disp,3);
						errores_str = mat2str([metrics.err_amb, metrics.err_no_amb, metrics.new_err_amb, metrics.new_no_err_amb],3);

fprintf('H = %d; (1/2wM)||X-AX||_F^2 = %.1f <<< = >>>> fscore %s ; errores %s\n', ...
							H,  frob_dif, fscore_str, errores_str);
% fprintf('No hidden lambda 10^{%.1f}: fscore %s ; errores %s\n\n', log10(lambda_scale),fscore_str, errores_str);
fprintf('-----------------------\n');
						exec_time_cell{lambda_scale_idx} = exec_time;
						metrics_cell{lambda_scale_idx} = metrics;
						S_cell{lambda_scale_idx} = S;
						term1_cell{lambda_scale_idx} = term1;
						term2_cell{lambda_scale_idx} = term2;
					% end
				end

				% Save the scenarios
				for lambda_scale_idx = 1:length(lambda_scale_list)
					lambda_scale = lambda_scale_list(lambda_scale_idx);
					if ~isempty(S_cell{lambda_scale_idx})

						metrics = metrics_cell{lambda_scale_idx};
						S = S_cell{lambda_scale_idx};
						term1 = term1_cell{lambda_scale_idx};
						term2 = term2_cell{lambda_scale_idx};
						exec_time = exec_time_cell{lambda_scale_idx};
fscores_disp = metrics.fscore(5);
fscore_str = mat2str(fscores_disp,3);
errores_str = mat2str([metrics.err_amb, metrics.err_no_amb],3);
fprintf('No hidden lambda 10^{%.1f}: fscore %s ; errores %s\n\n', log10(lambda_scale),fscore_str, errores_str);

						filename = sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/Covsqrt_Estimate_N%d_M%d_H%d_l%d_w%d.mat', ...
							ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ...
							N, M, H, lambda_scale, w);
						save(filename, 'metrics', 'exec_time', 'A', 'S', 'N', 'M', 'H', 'lambda_scale', 'w', 'term1', 'term2');
					end
				end
			end
			disp('------Next scen---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
		end
		end
	end

fprintf('--------------------- FINISHED GRAPH %d (--------------------- \n', ID_graph)
end
disp('Finished executions no hidden')

%% ---------------------------------------------------------------------------------------------------------------
% %% Solve hidden for each scenario
% lambda_scale_list = lambda_scale_list ./ 1e3;
% exps_beta = 6.5;
% beta_list = 10.^exps_beta;
% beta_list = 0;
% alpha_list = [0];

if (strcmp(type_of_hidden, 'min')) || (strcmp(type_of_hidden, 'max'))
	scenarios_to_solve = 1;
end 

for ID_graph = graphs_to_solve
	for ID_signal = signals_to_solve
		for M = M_list
		for ID_scenario = scenarios_to_solve
			for H = H_list % DELETE THIS
			% Load graph signal y scenario

				scenario_data = load(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/scenario%d_H%d_%s.mat', ...
						ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ID_scenario, H, type_of_hidden));


				if ~exist(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/', ...
							ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden), 'dir')
					mkdir(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/', ...
							ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden));
				end
				O = N - H;
				C = scenario_data.C;
				A = scenario_data.A; % This is the one that has hidden inside. It's not ground truth
				X = scenario_data.X;
				degrees = scenario_data.degrees;
				noise_power = scenario_data.noise_power;
				s_h = scenario_data.s_h;

				if ~isempty(s_h)
					noise_power_hidden = mat2str(noise_power(s_h),2);
					degrees_hidden = mat2str(degrees(s_h),2);
				else
					noise_power_hidden = [];
					degrees_hidden = [];
				end

				mean_degrees = mean(degrees);
				% frob_A = norm(A, 'fro')^2;
				% frob_X = norm(X, 'fro')^2;
				frob_dif = norm((X-A*X), 'fro')^2/(2*w*M);

				% Create cells to save output for parfor


				for lambda_scale = (lambda_scale_list)
					for alpha = alpha_list
						metrics_cell = cell(1, length(beta_list));
						S_cell = cell(1, length(beta_list));
						term_fro_cell = cell(1, length(beta_list));
						term_norm1_cell = cell(1, length(beta_list));
						term_normnuc_cell = cell(1, length(beta_list));
						term_normlasso_cell = cell(1, length(beta_list));
						exec_time_cell = cell(1, length(beta_list));
						K_cell = cell(1, length(beta_list));
						for beta_idx = 1:length(beta_list)
							beta = beta_list(beta_idx);
							if ~exist(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/Cov_full_N%d_M%d_H%d_l%d_w%d_alpha%d_beta%d.mat', ...
									ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ...
									N, M, H, lambda_scale, w, alpha, beta), 'file')
								if (sim == 1)
								[S, exec_time, K, term_fro, term_norm1, term_normnuc, term_normlasso] = homocedastic_Cov_full_hidden(lambda_scale, C, w, alpha, beta, M) ; 						
								else
									[S, exec_time, K, term_fro, term_norm1, term_normnuc, term_normlasso] = homocedastic_Cov_full_hidden_nosim(lambda_scale, C, w, alpha, beta, M) ; 						
								end

								[metrics] = graph_learning_perf_eval(A,S);


fscores_disp = metrics.fscore([2,4,5,6]);
fscore_str = mat2str(fscores_disp,3);
errores_str = mat2str([metrics.err_amb, metrics.err_no_amb, metrics.new_err_amb, metrics.new_no_err_amb], 3);
% fprintf('H = %d; mean_deg = %.2f ; deg_hidden = %s; noise_hidden = %s; (1/2wM)||X-AX||_F^2 = %.1f <<< = >>>> l: %d ; b %d; alpha %.1f<<< = >>>> fscore %s ; errores %s\n', ...
% 	H, mean_degrees, degrees_hidden, noise_power_hidden, frob_dif, lambda_scale, beta, alpha, fscore_str, errores_str);
fprintf('Hidden beta 10^{%.1f}  <<< = >>>> fscore %s ; errores %s\n', log10(beta), fscore_str, errores_str);
								exec_time_cell{beta_idx} = exec_time;
								metrics_cell{beta_idx} = metrics;
								S_cell{beta_idx} = S;
								term_fro_cell{beta_idx} = term_fro;
								term_norm1_cell{beta_idx} = term_norm1;
								term_normnuc_cell{beta_idx} = term_normnuc;
								term_normlasso_cell{beta_idx} = term_normlasso;
								K_cell{beta_idx} = K; 
							end
						end

						for beta_idx = 1:length(beta_list)
							beta = beta_list(beta_idx);
							if ~isempty(S_cell{beta_idx})
								exec_time = exec_time_cell{beta_idx};
								metrics = metrics_cell{beta_idx} ;
								S = S_cell{beta_idx};
								term_fro = term_fro_cell{beta_idx};
								term_norm1 = term_norm1_cell{beta_idx};
								term_normnuc = term_normnuc_cell{beta_idx};
								term_normlasso = term_normlasso_cell{beta_idx};
								K = K_cell{beta_idx};
								filename = sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/Cov_full_N%d_M%d_H%d_l%d_w%d_alpha%d_beta%d.mat', ...
									ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ...
									N, M, H, lambda_scale, w, alpha, beta);
								save(filename, 'metrics', 'exec_time', 'A', 'S', 'N', 'M', 'H', 'lambda_scale', 'w', 'term_fro', 'term_norm1', 'term_normnuc', 'term_normlasso', 'K');
							end
						end
					end
				end

			end
			disp('---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')

		end
		end
		fprintf('Executions hidden: Finished signal %d\n', ID_signal)
	end
	fprintf('Executions hidden: Finished graph %d\n', ID_graph)
	fprintf('--------------------- FINISHED GRAPH %d (--------------------- \n', ID_graph)

end
disp('Finished executions hidden')

disp('FINISHED EVERYTHING')
end

