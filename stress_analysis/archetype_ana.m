%this script preforms PCA and then archetype analysis based on the genes
%that define each cell state, as described in the Baron et al. 2020.
%requirements are: the magic package
%(https://github.com/KrishnaswamyLab/MAGIC) and an example dataset
%PDAC_cancer_cells.mat. the format of the example dataset is raw counts (UMI) and
%genes as rows, cancer cells as columns.

pca_ana       = 1;%this part do the PCA and color expression based on gene sets
bootstrap_ana = 0;%this part gives a p-value for enrichment in the vertex

try
    A;
catch
    load PDAC_cancer_cells.mat %switch to your own dataset if intrested
    load cmap_color_blind2.mat %for colorblind
    set(0,'DefaultFigureColormap',flip(cmap_color_blind2));
end

if (pca_ana)
    %normalize, smooth, transform and z-score
    A_tpm = median(sum(A))*bsxfun(@rdivide,A,sum(A));
    A_smooth_tpm = run_magic(A',4)';
    A_ans = sqrt(A_smooth_tpm)+sqrt(A_smooth_tpm+1);
    A_z = zscore(A_smooth_tpm,0,2);
    A_z_no_smooth = zscore(A_tpm,0,2);
       
    %find informative genes for PCA
    max_read_per_cell = mean(A_tpm')';
    fano_factor       = var(A_tpm')'./mean(A_tpm,2);
    thresh = 0;
    informative_genes = intersect(find(fano_factor>nanmean(fano_factor)+thresh*nanstd(fano_factor)),find(max_read_per_cell>=nanmean(max_read_per_cell)+thresh*nanstd(max_read_per_cell)));
    
    %remove mito and ribo genes for the PCA but keep in the expression
    load ribo_mito_genes.mat
    [~,~,ribo_idx] = intersect(ribo_genes,gene_names);
    [~,~,mito_idx] = intersect(mit_genes,gene_names);
    informative_genes_new = setdiff(informative_genes,[ribo_idx;mito_idx]); length(informative_genes_new)
    
    %pca 
    [coeff, score, latent, tsquared, explained] = pca(A_ans(informative_genes_new,:)');   
    figure; scatter(score(:,1),score(:,2),30,'k','filled');
    
    %show some markers genes for the cancer cell states. can switch to any
    %gene of intrest
    figure; subplot(1,4,1);
    gene = 'JUN'; gene_ID = strmatch(gene,gene_names,'exact');
    scatter(score(:,1),score(:,2),30,A_z(gene_ID,:),'filled'); hold on; title(gene)
    axis off; caxis([-2 2]);
    subplot(1,4,2);
    gene = 'FOS'; gene_ID = strmatch(gene,gene_names,'exact');
    scatter(score(:,1),score(:,2),30,A_z(gene_ID,:),'filled'); hold on; title(gene)
    axis off; caxis([-2 2]);
    subplot(1,4,3);
    gene = 'UBC'; gene_ID = strmatch(gene,gene_names,'exact');
    scatter(score(:,1),score(:,2),30,A_z(gene_ID,:),'filled'); hold on; title(gene)
    axis off; caxis([-2 2]);
    
    %stress sign.
    try
        humanstressgenes;
    catch
        load human_stress_genes.mat
    end
    [~,~,stress_sig_genes] = intersect(humanstressgenes,gene_names);
    color_vec = zscore(sum(A_z(stress_sig_genes,:)),0,2);
    subplot(1,4,4);
    scatter(score(:,1),score(:,2),30,color_vec,'filled');
    axis off; title('stress sig'); caxis([-2 2]);
    colorbar('Position',...
            [0.95 0.106369414295492 0.017 0.4]);
end

if (bootstrap_ana)
    %stress sign.
    try
        humanstressgenes;
    catch
        load human_stress_genes.mat
    end
    [~,~,stress_sig_genes] = intersect(humanstressgenes,gene_names);
    stress_sig = sum(A_tpm(stress_sig_genes,:));
    
    genes_to_take = informative_genes_new;
    %creating random sets
    num = 10^3;
    random_sets = (randi(length(genes_to_take),[length(stress_sig_genes)],num));
    random_dis = nan(size(A,2),num);
    for r = 1:num
        random_dis(:,r) = sum(A_tpm(genes_to_take(random_sets(:,r)),:));
    end
    
    figure; histogram(mean(random_dis,2)); hold on;
    line([mean(stress_sig) mean(stress_sig)],[0 120]);
    %final p-value
    pval = (100-invprctile(mean(random_dis,2),mean(stress_sig)))./100;
    
end
