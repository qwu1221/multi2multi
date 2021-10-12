function new_clu = Refine_spatial(t_in,distW0, node_weights)
% graph to indicate the connectivity for the detected voxels

distW = (distW0<2)-eye(size(distW0,2));
distWt = distW(t_in,t_in);
t_dist=ones(size(t_in,2),size(t_in,2));
t_dist = distWt.*t_dist;
%figure;imagesc(t_dist);colorbar;

% find maximal connected component
t_graph = graph(t_dist);
bins = conncomp(t_graph);
binnodes = accumarray(bins', 1:numel(bins), [], @(v) {sort(v')});
val = unique(bins);
cnt = histc(bins,val);
%large_clu =find(cnt>10);
[~,cid] =sort(cnt,'descend');
num_com = size(val,2);
for i=1:num_com
    clu{i} = t_in(find(bins==cid(i)));
end


% construct the node-weighted graph for searching the shortest path
% utilize the sum of SNPs in IGDB
W_weight = (node_weights'*ones(1,size(distW0,2)))';
W_weight = (1./W_weight).*(distW0<2);


com_seq = 1:num_com;
for k=1:(num_com-1)
    %[k size(com_seq,2)]
    ifpath = zeros(size(com_seq,2)-1,size(com_seq,2));
    for i=1:(size(com_seq,2)-1)
     for j=(i+1):size(com_seq,2)
        clui = clu{com_seq(i)};
        cluj = clu{com_seq(j)};
        %[i j]
        %[size(clui,2) size(cluj,2)]
        distij =  min(min(distW0(clui,cluj)));
        [r,c] = find(distW0(clui,cluj)==distij);
        r_raw = clui(r);
        c_raw = cluj(c);
        [e L] = dijkstra(W_weight,r_raw(1),c_raw(1));
        if size(L,2)-2<(size(clui,2)+size(cluj,2))/2
            ifpath(i,j) = e;
        end
        %size(L)
        %path = [path L];
        end
    end
    if any(ifpath>0)==0 
        break
    else
        eMin = min(ifpath(ifpath>0));
        [idx1 idx2] = find(ifpath==eMin); 
        clui = clu{com_seq(idx1)};
        cluj = clu{com_seq(idx2)};
        [r,c] = find(distW0(clui,cluj)==distij);
        r_raw = clui(r);
        c_raw = cluj(c);
        [e L] = dijkstra(W_weight,r_raw(1),c_raw(1));
 
        clu{max(com_seq)+1} = unique([clui cluj L]);
        com_seq = [com_seq max(com_seq)+1];
        com_seq = setdiff(com_seq, [com_seq(idx1) com_seq(idx2)]);
    end
end
new_clu = clu{com_seq};
%t_in_new = unique([new_clu{:}]);
end



