function fsc = fscore(Strue,Sest)
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
% Strue => True adjacency matrix
% Sest => Estimated adjacency matrix


%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
% fsc => fscore

    true_pos = sum(Strue(Sest==1)==1);
    total_pos = sum(sum(Sest==1));
    precision= true_pos/total_pos;
    relevant_elems = sum(sum(Strue==1));
    recall = true_pos/relevant_elems;
    fsc= 2*precision*recall/(precision + recall);
end
