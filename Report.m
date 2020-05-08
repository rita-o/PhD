%% Add paths and load data
close all
clc 
clear all

%cd('C:\Users\biol-etupol204-2\Desktop\2020_HandsOn_Connectomics')
cd('C:\Users\Rita\Documents\Rita\Work\_PhD CHUV\Diffusion_Course\2020_HandsOn_Connectomics')
addpath(genpath('Utilities'))

load('Data\Connectome_HandsOn_Dataset.mat')

%%------------------------------------------------------------------------
%% Controls - get data
close all

% Get matrix for one subject and binarize
A = SC_ctrl(:,:,1);
th = 0;
B = double(A>th); 
nn = size(B,1);

% Plot
figure, imagesc(log(A)), axis equal tight, colorbar, title('Connectivity matrix')
figure, imagesc(B), axis equal tight, title('Connectivity matrix binarize')

% Get consensus matrix
ns = size(SC_ctrl,3);
for i=1:ns
    SC_ctrl_bin(:,:,i)=double(SC_ctrl(:,:,i)>th);
end
C = sum(SC_ctrl_bin,3)/ns;

% Threshold consensus matrix
th_perc = 0.5;
SC_controls = double(C>=th_perc);

% Plot
figure, imagesc(SC_controls),axis equal tight, title('Connectivity matrix binarize after consensus')

%% Patients  - get data
close all

% Get matrix for one subject and binarize
A = SC_schz(:,:,1);
th = 0;
B = double(A>th);

% Plot
figure, imagesc(log(A)), axis equal tight, colorbar, title('Connectivity matrix')
figure, imagesc(B), axis equal tight, title('Connectivity matrix binarize')

% Get consensus matrix
ns = size(SC_schz,3);
for i=1:ns
    SC_schz_bin(:,:,i)=double(SC_schz(:,:,i)>th);
end
C = sum(SC_schz_bin,3)/ns;

% Threshold consensus matrix
th_perc = 0.5;
SC_patients = double(C>=th_perc);

% Plot
figure, imagesc(SC_patients),axis equal tight, title('Connectivity matrix binarize after consensus')

%%-------------------------------------------------------------------------
%% Controls - Compute metrics
% Controls network
D = distance_bin(SC_controls);
Cpl = charpath(D);
Cl = clustering_coef_bu(SC_controls); 
Cl_avg = mean(Cl);

% Random network
R = randmio_und(SC_controls,100); % 100 number of iterations

D_rand = distance_bin(R);
Cpl_rand = charpath(D_rand);
Cl_rand =clustering_coef_bu(R);
Cl_rand_avg= mean(Cl_rand);

% Ratios 
Cpl_ratio=Cpl/Cpl_rand;
Cl_ratio = Cl_avg/Cl_rand_avg;

% Small word index
SW_controls = Cl_ratio/Cpl_ratio;

%% Patients - Compute metrics
% Patient network
D = distance_bin(SC_patients);
Cpl = charpath(D);
Cl = clustering_coef_bu(SC_patients); 
Cl_avg = mean(Cl);

% Random network
R = randmio_und(SC_patients,100); % 100 number of iterations

D_rand = distance_bin(R);
Cpl_rand = charpath(D_rand);
Cl_rand =clustering_coef_bu(R);
Cl_rand_avg= mean(Cl_rand);

% Ratios 
Cpl_ratio=Cpl/Cpl_rand;
Cl_ratio = Cl_avg/Cl_rand_avg;

% Small word index
SW_patients = Cl_ratio/Cpl_ratio;

% Yes the small world is preserved in patients
% In a randomize network you need to preserve the degree distribution
% but not the density??

%% Network density and connectivity
close all

% Variables initialization
d_ctrl = zeros(ns,1);
d_schz = zeros(ns,1);
w_ctrl = zeros(ns,1);
w_schz = zeros(ns,1);

% Loop over subjects
for i = 1:ns
    % CTRL
    this_network = SC_ctrl(:,:,i);
    d_ctrl(i) = nnz(this_network)/(nn*(nn-1));
    w_ctrl(i) = sum(sum(triu(this_network)));
    % SCHZ
    this_network = SC_schz(:,:,i);
    d_schz(i) = nnz(this_network)/(nn*(nn-1));
    w_schz(i) = sum(sum(triu(this_network)));
end

% Statistical comparison
p_d = ranksum(d_ctrl,d_schz);
p_w = ranksum(w_ctrl,w_schz);
figure, boxplot([d_ctrl,d_schz],{'CTRL','SCHZ'}), title(['Density, p=',num2str(p_d)]);
figure, boxplot([w_ctrl,w_schz],{'CTRL','SCHZ'}), title(['Weight sum, p=',num2str(p_w)]);

% Patients and controls have different densities

% ranksum compares 2 groups (non parametric)
% discuss multiple comprisons

%% Efficiency
close all

% Variables initialization
Eff_ctrl = zeros(ns,1); Eff_schz = zeros(ns,1);
Cl_ctrl = zeros(ns,1); Cl_schz = zeros(ns,1);

% Loop over subjects
for i = 1:ns
    % CTRL
    this_network = SC_ctrl(:,:,i);
    this_network = this_network./(sum(this_network(:))/2);
    Eff_ctrl(i) = efficiency_wei(this_network);
    % SCHZ
    this_network = SC_schz(:,:,i);
    this_network = this_network./(sum(this_network(:))/2);
    Eff_schz(i) = efficiency_wei(this_network);
end

% Statistical comparison
p_Eff = ranksum(Eff_ctrl,Eff_schz);
figure, boxplot([Eff_ctrl,Eff_schz],{'CTRL','SCHZ'}),
title('Weighted network efficiency');

%  For ease of interpretation of the local efficiency it may be
%   advantageous to rescale all weights to lie between 0 and 1

% The efficient is increased

%% Nodal Closeness centralities
close all

% Variables initialization
ClCent_ctrl = zeros(ns,nn);
ClCent_schz = zeros(ns,nn);

% Loop over subjects
for i = 1:ns
    % CTRL
    this_network = SC_ctrl(:,:,i);
    this_network = this_network./(sum(this_network(:))/2);
    
    connection_lengths = 1 ./ this_network;
    connection_lengths(isinf(connection_lengths))= 0;
    [this_d,~] = distance_wei(connection_lengths);
    ClCent_ctrl(i,:) =  (nn-1)./sum(this_d);
    
    % SCHZ
    this_network = SC_schz(:,:,i);
    this_network = this_network./(sum(this_network(:))/2);
    connection_lengths = 1 ./ this_network;
    connection_lengths(isinf(connection_lengths))= 0;
    [this_d,~] = distance_wei(connection_lengths);
    ClCent_schz(i,:) = (nn-1)./sum(this_d);
    
    
end

% Plot
figure,subplot(1,2,1),hist(ClCent_ctrl(:)),title('CTRL');
subplot(1,2,2), hist(ClCent_schz(:)), title('SCHZ');

%% Compare centralities
close all

% Variables initialization
p_ClCent = zeros(nn,1);

% Loop for nodes
for i = 1:nn
    p_ClCent(i)=ranksum(ClCent_ctrl(:,i),ClCent_schz(:,i));
    
end

% Plot
figure, plot(p_ClCent,'o'), hold on;
plot([1 nn],[0.05 0.05],'--');
xlabel('node'), ylabel('p-values');
title('Nodal p-values CTRLvsSCHZ');

% FDR correction
ii_survive = FDR_main(p_ClCent, 0.05, 'bh95');
labels(ii_survive)

% Bonferoni correction
ii_survive2 = labels(p_ClCent<(0.05/nn));


%% Toolbox visuzalize
close all

% Create weighted matrix
SCw = zeros(nn);
for i=1:nn
    for j=1:nn
        this_connection =  SC_schz(i,j,:);
        SCw(i,j)=mean(this_connection(this_connection>0)); % filtering values larger than 0
    end
end
SCw(isnan(SCw))=0;
SCw = SCw.* SC_patients;

% Create edge file
dlmwrite('my_SCw.edge',SCw,'delimiter','\t')

% Create node file
labels(p_ClCent<0.05)

hubs = [9,12,15,28,54,56,68,69,82];
n_color=zeros(nn,1);
n_color(hubs)=1;
n_size = sum(SCw);
load('C:\Users\Rita\Documents\Rita\Work\_PhD CHUV\Diffusion_Course\2020_HandsOn_Connectomics\Data\BNV_Data\centroids_bert.mat');
%load('C:\Users\biol-etupol204-2\Desktop\2020_HandsOn_Connectomics\Data\BNV_Data\centroids_bert.mat')

filepath='my_nodes.node';
generate_node_file(filepath,coord,n_color,n_size);

BrainNet

%%

k = sum(SC_patients);
figure
subplot(2,1,1)
hist(k,10)
title('Degree distribution')

idx = [9,12,15,28,54,56,68,69,82];
degrees_impared = k(idx);
subplot(2,1,2)
hist(degrees_impared,10)
title('Degree distribution - impaired nodes')
