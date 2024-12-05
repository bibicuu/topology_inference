function [s_n, s_h] = select_hidden_nodes(sel_type, O, A)
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
% sel_type => What nodes to select (min, max, min_smooth, max_smooth,
% max_diff_sm, rand)
% O => Observed nodes
% L => Laplacian
% X => Signals in the graph


%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
% s_h => Position of the selected hidden nodes
% s_n =>  Position of the selected observed nodes 

    N = size(A,1);
    s_n = [];
    s_h = [];
    switch sel_type
        % select nodes with less links
        case 'min' 
            connections = sum(A);
            [~,node_pos] = sort(connections,'descend');
        % select nodes with more links
        case 'max'
            connections = sum(A);
            [~,node_pos] = sort(connections,'ascend');
        % % select nodes which when removed the smoothness is smaller
        % case 'min_smooth'
        %     smooth = node_smoothness(N, L, X);
        %     [~,node_pos] = sort(smooth,'descend');
        % % select nodes which when removed the smoothness is bigger
        % case 'max_smooth'
        %     smooth = node_smoothness(N, L, X);
        %     [~,node_pos] = sort(smooth,'ascend');
        % % select nodes which when removed the variation of the smoothness
        % % is bigger
        % case 'max_diff_sm'
        %     smooth = node_smoothness(N, L, X);
        %     diff_sm = abs(smooth-trace(X'*L*X));
        %     [~,node_pos] = sort(diff_sm,'ascend');
        case 'rand'
            [~,node_pos] = sort(rand(N,1));
        % select nodes at random
        otherwise
            disp('ERR: unknown method for selecting the hidden node')
            return
    end
    s_n = sort(node_pos(1:O));
    s_h = sort(node_pos(O+1:N));
end

