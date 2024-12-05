close all; clear all; clc;

% Hacer tambien con H = 1e6
M = [1e6];
H = [2];

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


metrics = {'fscore', 'shamming_dist', 'err_amb', 'err_no_amb'};
type_binarize = 2;
metrics_name = cell(1, length(metrics));

% %% Analyze lambda
% 
% 
% beta = 10.^(6.5);
% alpha = 0;
% lambda_exps = -10:0.1:2;
% 
% 
% for i = 1:length(metrics)
%     metric_name = metrics{i};
%     tensor_name = [metric_name, '_nohid_tensor'];
%     eval([tensor_name, ' = zeros(nr_graphs, nr_signals, nr_scenarios, length(lambda_exps));']);
%     metrics_name{i} = tensor_name; % Guarda el nombre en la lista
% 
% 	tensor_name = [metric_name, '_hid_tensor'];
%     eval([tensor_name, ' = zeros(nr_graphs, nr_signals, nr_scenarios, length(lambda_exps));']);
%     metrics_name{i+length(metrics)} = tensor_name;
% end
% 
% 
% 
% for i = 1:length(metrics_name)
%     tensor_name = metrics_name{i};
% 	for ID_graph = 1:nr_graphs
% 		for ID_signal = 1:nr_signals
% 			for ID_scenario = 1:nr_scenarios
% 				for lambda_idx = 1:length(lambda_exps)
% 					lambda_scale = 10.^(lambda_exps(lambda_idx));
% 
% 
% 					if (i <= length(metrics))
% 						data_nohid = load(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/Covsqrt_Estimate_N%d_M%d_H%d_l%d_w%d.mat', ...
% 												ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ...
% 												N, M, H, lambda_scale, w)).metrics;
% 						metric_idx = mod(i,length(metrics));
% 						if metric_idx == 0
% 							metric_idx = length(metrics);
% 						end
% 						metric = metrics(metric_idx);
% 
% 						if strcmp(metric, 'fscore')
% 							eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, lambda_idx) = data_nohid.fscore(', num2str(type_binarize), ');']);
% 						elseif strcmp(metric, 'shamming_dist')
% 							eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, lambda_idx) = data_nohid.shamming_dist(', num2str(type_binarize), ');']);
% 						else
% 							eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, lambda_idx) = data_nohid.', metric{1}, ';']);
% 						end
% 
% 					else
% 
% 						data_hid = load(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/Cov_full_N%d_M%d_H%d_l%d_w%d_alpha%d_beta%d.mat', ...
% 														ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ...
% 														N, M, H, lambda_scale, w, alpha, beta)).metrics;
% 
% 						metric_idx = mod(i,length(metrics));
% 						if metric_idx == 0
% 							metric_idx = length(metrics);
% 						end
% 						metric = metrics(metric_idx);
% 
% 						if strcmp(metric, 'fscore')
% 							eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, lambda_idx) = data_hid.fscore(', num2str(type_binarize), ');']);
% 						elseif strcmp(metric, 'shamming_dist')
% 							eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, lambda_idx) = data_hid.shamming_dist(', num2str(type_binarize), ');']);
% 						else
% 							eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, lambda_idx) = data_hid.', metric{1}, ';']);
% 						end
% 
% 					end
% 				end
% 			end
% 		end
% 	end
% end
% 
% 
% 
% for i = 1:length(metrics_name)
%     tensor_name = metrics_name{i}; % Nombre del tensor actual
%     base_name = erase(tensor_name, '_tensor'); % Elimina el sufijo '_tensor'
% 
%     % Calcula la media y la mediana a lo largo de la cuarta dimensión (H_idx)
%     mean_value = squeeze(squeeze(mean(mean(mean(eval(tensor_name), 1), 2), 3)));
%     median_value = squeeze(squeeze(median(median(median(eval(tensor_name), 1), 2), 3)));
% 
%     % Asigna los resultados a nuevas variables con los nombres deseados
%     eval([base_name, '_mean = mean_value;']);
%     eval([base_name, '_median = median_value;']);
% end
% 
% 
% 
% 
% % Define la métrica que deseas visualizar
% % fscore, err_amb, err_no_amb, shamming_dist, new_err_amb, err_new_no_amb
% %%
% metric_figure = 'err_amb'; % Cambia 'err_no_amb' por cualquier métrica, como 'fscore'
% metric_name = 'nse';
% logaritmic = 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Genera los nombres dinámicos para la media y mediana de HID y NO HID
% hid_mean_name = [metric_figure, '_hid_mean'];
% nohid_mean_name = [metric_figure, '_nohid_mean'];
% hid_median_name = [metric_figure, '_hid_median'];
% nohid_median_name = [metric_figure, '_nohid_median'];
% 
% % Gráfico de la media
% close all;clc;
% figure;
% plot(10.^lambda_exps, eval(nohid_mean_name), '-s', 'DisplayName', ['No hid'], 'LineWidth', 1.5);
% hold on;
% plot(10.^lambda_exps, eval([hid_mean_name]), '-o', 'DisplayName', ['Hid '], 'LineWidth', 1.5);
% 
% xlabel('$log(\lambda_{sc})$', 'Interpreter', 'latex');
% ylabel(sprintf('%s', metric_name), 'Interpreter', 'latex');
% legend('show', 'Interpreter', 'latex', 'Location', 'northwest');
% % title(['Comparación de ', metric_figure, ' para HID y NO HID - Media']);
% grid on; box on;
% set(gca, 'XScale', 'log');
% if logaritmic == 1
% set(gca, 'YScale', 'log');
% end
% hold off;
% 
% 
% 
% % Gráfico de la mediana
% figure;
% plot(10.^lambda_exps, eval(nohid_median_name), '-s', 'DisplayName', ['No hid'], 'LineWidth', 1.5);
% hold on;
% plot(10.^lambda_exps, eval([hid_median_name]), '-o', 'DisplayName',  ['Hid'], 'LineWidth', 1.5);
% 
% xlabel('$log(\lambda_{sc})$', 'Interpreter', 'latex');
% ylabel(sprintf('%s', metric_name), 'Interpreter', 'latex');
% legend('show', 'Interpreter', 'latex', 'Location', 'northwest');
% % title(['Comparación de ', metric_figure, ' para HID y NO HID - Mediana']);
% grid on; box on;
% set(gca, 'XScale', 'log');
% if logaritmic == 1
% set(gca, 'YScale', 'log');
% end
% hold off;
% 
% 
% 
% %%
% % if (type_binarize == 2) && strcmp(metric_figure, 'fscore')
% % 	metric_name = 'fscore2';
% % end
% % if logaritmic == 0
% % figure(1)
% % exportgraphics(gcf, sprintf('figures_hiperpams/%s_plot_media.png', metric_name), 'Resolution', 600);
% % figure(2)
% % exportgraphics(gcf, sprintf('figures_hiperpams/%s_plot_mediana.png', metric_name), 'Resolution', 600);
% % else
% % figure(1)
% % exportgraphics(gcf, sprintf('figures_hiperpams/%s_semilogy_media.png', metric_name), 'Resolution', 600);
% % figure(2)
% % exportgraphics(gcf, sprintf('figures_hiperpams/%s_semilogy_mediana.png', metric_name), 'Resolution', 600);
% % end




%% Analyze beta

beta_exps = -2:0.1:8;
alpha = 0;
lambda_scale = 10.^(-6);


for i = 1:length(metrics)
    metric_name = metrics{i};
    tensor_name = [metric_name, '_hid_tensor'];
    eval([tensor_name, ' = zeros(nr_graphs, nr_signals, nr_scenarios, length(beta_exps));']);
    metrics_name{i} = tensor_name; % Guarda el nombre en la lista
end


for i = 1:length(metrics_name)
    tensor_name = metrics_name{i};
	for ID_graph = 1:nr_graphs
		for ID_signal = 1:nr_signals
			for ID_scenario = 1:nr_scenarios
				for beta_idx = 1:length(beta_exps)
					beta = 10.^(beta_exps(beta_idx));


					data_hid = load(sprintf('./simulations/graph%d_N%d_sc%d_%s_p%.2f_%s_weigh%d/signal%d_N%d_M%d_w%d/scenario%d_H%d_%s/sols/Cov_full_N%d_M%d_H%d_l%d_w%d_alpha%d_beta%d.mat', ...
													ID_graph, N, scale_graph, g_type, prob, symmetric, weighted, ID_signal, N, M, w, ID_scenario, H, type_of_hidden, ...
													N, M, H, lambda_scale, w, alpha, beta)).metrics;

					metric_idx = mod(i,length(metrics));
					if metric_idx == 0
						metric_idx = length(metrics);
					end
					metric = metrics(metric_idx);

					if strcmp(metric, 'fscore')
						eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, beta_idx) = data_hid.fscore(', num2str(type_binarize), ');']);
					elseif strcmp(metric, 'shamming_dist')
						eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, beta_idx) = data_hid.shamming_dist(', num2str(type_binarize), ');']);
					else
						eval([tensor_name, '(ID_graph, ID_signal, ID_scenario, beta_idx) = data_hid.', metric{1}, ';']);
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





% Define la métrica que deseas visualizar
% fscore, err_amb, err_no_amb, shamming_dist, fscore2
%%
metric_figure = 'err_amb'; % Cambia 'err_no_amb' por cualquier métrica, como 'fscore'
metric_name = 'nse';
logaritmic = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Genera los nombres dinámicos para la media y mediana de HID y NO HID
hid_mean_name = [metric_figure, '_hid_mean'];
% nohid_mean_name = [metric_figure, '_nohid_mean'];
hid_median_name = [metric_figure, '_hid_median'];
% nohid_median_name = [metric_figure, '_nohid_median'];

% Gráfico de la media
close all; clc;
figure;
% plot(10.^beta_exps, eval(nohid_mean_name), '-s', 'DisplayName', ['No hid'], 'LineWidth', 1.5);
% hold on;
plot(10.^beta_exps, eval([hid_mean_name]), '-o', 'DisplayName', ['Hid '], 'LineWidth', 1.5);

xlabel('$log(\beta_{sc})$', 'Interpreter', 'latex');
ylabel(sprintf('%s', metric_name), 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex', 'Location', 'northeast');
% title(['Comparación de ', metric_figure, ' para HID y NO HID - Media']);
grid on; box on;
set(gca, 'XScale', 'log');
if logaritmic == 1
set(gca, 'YScale', 'log');
end
hold off;



% Gráfico de la mediana
figure;
% plot(10.^beta_exps, eval(nohid_median_name), '-s', 'DisplayName', ['No hid'], 'LineWidth', 1.5);
% hold on;
plot(10.^beta_exps, eval([hid_median_name]), '-o', 'DisplayName',  ['Hid'], 'LineWidth', 1.5);

xlabel('$log(\beta_{sc})$', 'Interpreter', 'latex');
ylabel(sprintf('%s', metric_name), 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex', 'Location', 'northeast');
% title(['Comparación de ', metric_figure, ' para HID y NO HID - Mediana']);
grid on; box on;
set(gca, 'XScale', 'log');
if logaritmic == 1
set(gca, 'YScale', 'log');
end
hold off;

%%
% if (type_binarize == 2) && strcmp(metric_figure, 'fscore')
% 	metric_name = 'fscore2';
% end
% if logaritmic == 0
% figure(1)
% exportgraphics(gcf, sprintf('figures_hiperpams/%s_plot_media.png', metric_name), 'Resolution', 600);
% figure(2)
% exportgraphics(gcf, sprintf('figures_hiperpams/%s_plot_mediana.png', metric_name), 'Resolution', 600);
% % 
% else
% figure(1)
% exportgraphics(gcf, sprintf('figures_hiperpams/%s_semilogy_media.png', metric_name), 'Resolution', 600);
% figure(2)
% exportgraphics(gcf, sprintf('figures_hiperpams/%s_semilogy_mediana.png', metric_name), 'Resolution', 600);
% end







