% for type_binarize = 1:8
close all; clc;

M_plot = [1e6];
H_plot = 0:1:4;
nr_graphs = 30;
nr_signals = 3; % Generated random gaussian // diff realizations (for each graph)
nr_scenarios = 1;

% Graph
graph_params.N = 20;
graph_params.g_type = 'ER';
graph_params.p = 0.2;
graph_params.symmetric = 'nosim';
graph_params.scale_graph = 1; % If < 1 then scale with min_eig. Else with max_eig
graph_params.weighted = 1; % IF It is 0, scale graph cannot be 1!!!
% Signal
signal_params.w = 1e-2; 
signal_params.cov_matrix = signal_params.w * eye(graph_params.N);
% Scenario
scenarios_params.type_of_hidden = 'min'; % rand min o max

N = graph_params.N;
g_type = graph_params.g_type;
prob = graph_params.p;
symmetric = graph_params.symmetric;
weighted = graph_params.weighted;
scale_graph = graph_params.scale_graph;
w = signal_params.w;
type_of_hidden = scenarios_params.type_of_hidden;


alpha_plot = [0,1];

% alpha = 0;
beta = 10^(6.5);
lambda_scale = 1e-4;

% metrics = {'fscore', 'shamming_dist', 'TPR', 'FDR', 'err_amb', 'new_err_amb', 'err_no_amb', 'new_no_err_amb'};
metrics = {'fscore', 'shamming_dist', 'err_amb', 'err_no_amb'};
% for type_binarize = 1:8
type_binarize = 5;
metrics_name = cell(1, length(metrics));

for i = 1:length(metrics)
    metric_name = metrics{i};
    tensor_name = [metric_name, '_nohid_tensor'];
    eval([tensor_name, ' = zeros(nr_graphs, nr_signals, nr_scenarios, length(H_plot), length(M_plot));']);
    metrics_name{i} = tensor_name; % Guarda el nombre en la lista

	tensor_name = [metric_name, '_hid_tensor'];
    eval([tensor_name, ' = zeros(nr_graphs, nr_signals, nr_scenarios, length(H_plot), length(M_plot), length(alpha_plot));']);
    metrics_name{i+length(metrics)} = tensor_name;
end



for i = 1:length(metrics_name)
    tensor_name = metrics_name{i};
	for M_idx = 1:length(M_plot)
		M = M_plot(M_idx);
		for H_idx = 1:length(H_plot)
			H = H_plot(H_idx);
			for ID_graph = 1:nr_graphs
				for ID_signal = 1:nr_signals
					for ID_scenario = 1:nr_scenarios
						if (i <= length(metrics))
							data_nohid = load(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/Covsqrt_Estimate_N%d_M%d_H%d_l%d_w%d.mat', ...
													ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ...
													N, M, H, lambda_scale, w)).metrics;
							metric_idx = mod(i,length(metrics));
							if metric_idx == 0
								metric_idx = length(metrics);
							end
							metric = metrics(metric_idx);

							if strcmp(metric, 'fscore')
    							eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, H_idx, M_idx) = data_nohid.fscore(', num2str(type_binarize), ');']);
							elseif strcmp(metric, 'shamming_dist')
    							eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, H_idx, M_idx) = data_nohid.shamming_dist(', num2str(type_binarize), ');']);
							else
    							eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, H_idx, M_idx) = data_nohid.', metric{1}, ';']);
							end

							

							
						else
							for alpha_idx = 1:length(alpha_plot)
								alpha = alpha_plot(alpha_idx);
								data_hid = load(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/Cov_full_N%d_M%d_H%d_l%d_w%d_alpha%d_beta%d.mat', ...
																ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ...
																N, M, H, lambda_scale, w, alpha, beta)).metrics;
		
								metric_idx = mod(i,length(metrics));
								if metric_idx == 0
									metric_idx = length(metrics);
								end
								metric = metrics(metric_idx);
								
								if strcmp(metric, 'fscore')
    								eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, H_idx, M_idx, alpha_idx) = data_hid.fscore(', num2str(type_binarize), ');']);
								elseif strcmp(metric, 'shamming_dist')
    								eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, H_idx, M_idx, alpha_idx) = data_hid.shamming_dist(', num2str(type_binarize), ');']);
								else
    								eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, H_idx, M_idx, alpha_idx) = data_hid.', metric{1}, ';']);
								end
							end

						end
					end
				end
			end
		end
	end
end

for i = 1:length(metrics_name)
    tensor_name = metrics_name{i}; % Nombre del tensor actual
    base_name = erase(tensor_name, '_tensor'); % Elimina el sufijo '_tensor'
    
    % Calcula la media y la mediana a lo largo de la cuarta dimensión (H_idx)
    mean_value = squeeze(squeeze(mean(mean(mean(eval(tensor_name), 1), 2), 3)));
    median_value = squeeze(squeeze(median(median(median(eval(tensor_name), 1), 2), 3)));
    
    % Asigna los resultados a nuevas variables con los nombres deseados
    eval([base_name, '_mean = mean_value;']);
    eval([base_name, '_median = median_value;']);
end



%% Define la métrica que deseas visualizar
% fscore, err_amb, err_no_amb, shamming_dist, new_err_amb, err_new_no_amb
metric_figure = 'fscore'; % Cambia 'err_no_amb' por cualquier métrica, como 'fscore'
metric_name = 'F-score';
logaritmic = 0;


% Genera los nombres dinámicos para la media y mediana de HID y NO HID
hid_mean_name = [metric_figure, '_hid_mean'];
nohid_mean_name = [metric_figure, '_nohid_mean'];
hid_median_name = [metric_figure, '_hid_median'];
nohid_median_name = [metric_figure, '_nohid_median'];

% Gráfico de la media
close all;
figure;
plot(H_plot, eval(nohid_mean_name), '-s', 'DisplayName', ['No hid'], 'LineWidth', 1.5);
hold on;
for alpha_idx = 1:length(alpha_plot)
	plot(H_plot, eval([hid_mean_name, '(:,alpha_idx)']), '-o', 'DisplayName', ['Hid ', sprintf(' \\(\\alpha = %.1f\\)', alpha_plot(alpha_idx))], 'LineWidth', 1.5);
end
xlabel('H','Interpreter', 'latex');
ylabel(sprintf('%s', metric_name), 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex', 'Location', 'northeast');
% title(['Comparación de ', metric_figure, ' para HID y NO HID - Media']);
grid on; box on;
if logaritmic == 1
set(gca, 'YScale', 'log');
end
hold off;



% Gráfico de la mediana
figure;
plot(H_plot, eval(nohid_median_name), '-s', 'DisplayName', ['No hid'], 'LineWidth', 1.5);
hold on;

for alpha_idx = 1:length(alpha_plot)
	plot(H_plot, eval([hid_median_name, '(:,alpha_idx)']), '-o', 'DisplayName',  ['Hid', sprintf(' \\(\\alpha = %.1f\\)', alpha_plot(alpha_idx))], 'LineWidth', 1.5);
end
xlabel('H', 'Interpreter', 'latex');
ylabel(sprintf('%s', metric_name), 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex', 'Location', 'northeast');
% title(['Comparación de ', metric_figure, ' para HID y NO HID - Mediana']);
grid on; box on;
if logaritmic == 1
set(gca, 'YScale', 'log');
end
hold off;
% end





%% mkdir(sprintf('figures_exp1'))
if (type_binarize == 2) && strcmp(metric_figure, 'fscore')
	metric_name = 'fscore2';
end
if logaritmic == 1
figure(1)
exportgraphics(gcf, sprintf('figures_exp1/%s_plot_media.png', metric_name), 'Resolution', 600);
figure(2)
exportgraphics(gcf, sprintf('figures_exp1/%s_plot_mediana.png', metric_name), 'Resolution', 600);
else
figure(1)
exportgraphics(gcf, sprintf('figures_exp1/%s_semilogy_media.png', metric_name), 'Resolution', 600);
figure(2)
exportgraphics(gcf, sprintf('figures_exp1/%s_semilogy_mediana.png', metric_name), 'Resolution', 600);
end
