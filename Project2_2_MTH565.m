%Vikram Vijayakumar (02068559)
%MTH 565 Project 2_2
% https://github.com/ivanbrugere/matlab-networks-toolbox.
%I had reviewed the above github link to obtain examples for clustering
%coefficients

n = 100;  %Number of vertices
e = 495;  %Number of edges

for i = 1:3
    A = zeros(n);  %Start with an empty adjacency matrix
    edges_added = 0;  %Initialize the number of edges added
    
    while edges_added < e
        %Randomly pick two vertices
        u = randi(n);
        v = randi(n);
        
        if u ~= v && A(u, v) == 0  %Ensure they are distinct and no edge exists yet
            A(u, v) = 1;
            A(v, u) = 1;  %Create a symmetric graph
            edges_added = edges_added + 1;
        end
    end
    
    %Plot the graph
    G = graph(A);
    figure;
    plot(G, 'Layout', 'force');
    title(['Graph Sample ', num2str(i), ' from Gt(100,0.1)and Edges = 495']);

    %Calculate density of the graph
    num_edges = numedges(G);  %nnumber of edges
    density = 2 * num_edges / (n * (n - 1));
    disp(['Density of Graph ', num2str(i), ': ', num2str(density)]);

    %Calculate mean clustering coefficient
    A = adjacency(G);  % Adjacency matrix of your graph
    local_clustering = local_clustering_coeff(A);  % Compute local clustering coefficient
    mean_clustering = mean(local_clustering);  % Mean clustering coefficient
    disp(['Mean Clustering Coefficient: ', num2str(mean_clustering)]);

    %Calculate global clustering coefficient
    global_clustering = global_clustering_coeff(A);  % Compute global clustering coefficient
    disp(['Global Clustering Coefficient: ', num2str(global_clustering)]);

    %Calculate the diameter of the graph
    D = distances(G); %path between two vertices
    diameter = max(D(:)); %longest path between two vertices
    disp(['Diameter = ', num2str(diameter)]);

end

%Mean and global clustering fucntions
function global_C = global_clustering_coeff(A)
    triangles = trace(A^3) / 6;  % Number of triangles in the graph
    degree_sum = sum(sum(A));  % Sum of degrees of all vertices
    triples = 0.5 * sum(sum(A * A));  % Total number of connected triplets
    global_C = (6 * triangles) / triples;  % Global clustering coefficient (transitivity)
end

function C = local_clustering_coeff(A)
    n = size(A, 1);  %Number of vertices
    C = zeros(n, 1);  %Initialize local clustering coefficients
    for i = 1:n
        neighbors = find(A(i, :));  %Find the neighbors of a vertex i
        k_i = length(neighbors);  %Degree of a vertex i
        if k_i >= 2
            subgraph = A(neighbors, neighbors);  %Subgraph of neighbors of v
            num_triangles = sum(subgraph(:)) / 2;  % Number of triangles involving vertex i
            C(i) = (2 * num_triangles) / (k_i * (k_i - 1));  % Local clustering coefficient
        end
    end
end