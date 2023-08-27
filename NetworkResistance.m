%% Inverse Resistance Problem - Connor Blake 2023
clear
clc
close all
% To play with this model, tweak the Hyperparameters and RHO sections
%% Hyperparameters
l = 10; % length of wire
n = 500; % number of intervals, 1000 takes ~30 seconds on my laptop, 500 ~2 seconds

doTadpole31 = false; % ignore - for verifying R2G and G2R on a small, branched case
doCircle = true;

printGraph = false; % prints out the cycle graph to verify it is, in fact, a circle
printMatrices = false; % don't do this for large matrices
printOutputs = true; % error norms, times
printPlots = true; % graph of original and reconstructed functions
%% Define Graph
% tadpole graph
% if doTadpole31
%     n = 4;
%     g = [1 1 1 1]; % conductances ie 1/R_edge
%     a1 = [-1;0;1;0];
%     a2 = [-1;1;0;0];
%     a3 =[0;-1;1;0];
%     a4=[0;0;-1;1];
%     A = [a1 a2 a3 a4];
% end
% circle graph
if doCircle
    g = [1 1 1 1 1 1 1 1 1 1];
    g = 1./rho(n,l);
    A = zeros(n,n);
    for i = 1:n
        if i == n
            A(n,n) = 1;
            A(1,n) = -1;
        else
            A(i,i) = -1;
            A(i+1,i) = 1;
        end
    end
end

%% Graph Laplacian and Adjacency generation
Adj = A*A';
Adj = Adj - diag(diag(Adj));
Adj = abs(Adj);
G = A*diag(g)*A';
x = linspace(0,l,n);
rho_exact = 1./g; % exact discrete
%% Linear Algebra computation
% forward R computation
if printOutputs
    tic;
end
R = G2R(G);
if printOutputs
    disp("G-->R")
    toc
end
% inverse G computation
if printOutputs
    tic;
end
Gi = R2G(R,Adj);
if printOutputs
    disp("R-->G")
    toc
end
gi = diag(pinv(A)*Gi*pinv(A'))';
rho_ila = 1./gi; % inverse via inverse method


%% Functional Equation computation
% https://mathoverflow.net/questions/453133/how-to-find-the-inverse-of-a-product-of-two-integral-equations
R_i = R(1,:);
[R_total,Rindex] = max(R_i);
R_total = R_total*4;
rho_ifa = abs(discrete_deriv(x,.5*(R_total+sqrt(R_total*R_total-4*R_total*R_i))));
rho_ifa(Rindex) = .5*(rho_ifa(Rindex+1)+rho_ifa(Rindex-1));
rho_ifa = rho_ifa*l/n;


%% OUTPUTS
if printGraph
    Graph = graph(Adj);
    Graph.plot();
end
if printMatrices
    disp("A")
    disp(A)
    disp("G")
    disp(G)
    disp("FORWARD R")
    disp(R);
    disp("INVERTED G")
    disp(Gi)
end
if printOutputs
    disp("n")
    disp(n)
    disp("l")
    disp(l)
    % LA
    disp("R_i LA Max Gap")
    disp(max(abs(rho_exact-rho_ila))) % what i really want - computed inverted resistance via linalg
    disp("Linalg r-inversion error norm")
    disp(norm(rho_exact-rho_ila))
    % FA
    disp("R_i FA Max Gap")
    disp(max(abs(rho_exact-rho_ifa)))
    disp("FA r-inversion error norm")
    disp(norm(rho_exact-rho_ifa))
end
if printPlots
    plot(x, rho_ila);
    hold on;
    plot(x, rho_exact);
    hold on;
    plot(x, rho_ifa);
    hold off;
    ylim([0,1.2*max([rho_ifa,rho_exact,rho_ila])])
    xlabel('X-position');
    ylabel('Rho');
    legend('Rho_{i-linalg}', 'Rho','Rho_{i-fa}');
end

%% Utility Functions
function Rout = G2R(g)
    % https://web.stanford.edu/~boyd/papers/pdf/eff_res.pdf
    n = size(g,1);
    Rout = zeros(n,n);
    gi = pinv(g);
    for i = 1:n
        for j = 1:i
            e_i = zeros(n, 1);
            e_j = zeros(n, 1);
            e_i(i) = 1;
            e_j(j) = 1;
            v1 = gi*(e_i-e_j);
            Rout(i,j) = abs(v1(i)-v1(j));
        end
    end
    Rout = Rout + Rout';
end
function Gout = R2G(r,adj)
    % https://mathoverflow.net/questions/303120/is-it-possible-to-compute-a-valid-laplacian-matrix-from-an-effective-resistance
    n = size(r,1);
    mu = 2*ones(n,1);
    for i = 1:n
        indices = find(adj(:,i) == 1);
        for j = 1:length(indices)
            currentIndex = indices(j);
            mu(i) = mu(i) - r(i,currentIndex);
        end
    end
    Gout = -2*pinv(r)+2*(mu*mu')/(mu'*r*mu); % inverts R to find G using adjacency matrix
end

function deriv = discrete_deriv(x,v)
    %https://math.stackexchange.com/questions/302160/correct-way-to-calculate-numeric-derivative-in-discrete-time
    deriv = zeros(size(v));
    for i = 1:length(v)
        if i == 1
            deriv(1) = v(2)-v(1);
        elseif i == length(v)
            deriv(i) = v(i)-v(i-1);
        else
            deriv(i) = .5*(v(i+1)-v(i-1));
        end
    end
    deriv = deriv*length(x)/(x(end)-x(1));
end
%% RHO
function rho = rho(num,l)
    x = linspace(0,l,num);
    rho = (l-x).*x+.2;
    rho = cos(4*pi/l*x) + 2;
    %rho = heaviside(x-5)+1+.1*x; % I know this causes weird problem, but
    %that's on my janky discrete derivative operator not the math
    plot(x,rho)
end