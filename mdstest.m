main()

function [mat_t1, mat_t2] = load_IMAGEN(mid_or_sst)
    addpath('G:/.shortcut-targets-by-id/1Nj5b1RhD0TcXoswrxiV5gkSPPq4fnGfu/IMAGEN_master_data/Matrices_Qinghao_new/matrices')
    
    disp("loading bsl")
    mat = load_IMAGEN_file(mid_or_sst, "bsl");
    mat_t1 = tanh(mat);

    disp("loading fu2")
    mat = load_IMAGEN_file(mid_or_sst, "fu2");
    mat_t2 = tanh(mat);
    
    for i = 1:size(mat_t1, 3)
        for j = 1:size(mat_t1, 1)
            mat_t1(j,j,i) = 1;
            mat_t2(j,j,i) = 1;
        end
    end
end

function mat = load_IMAGEN_file(mid_or_sst, bsl_or_fu2)
    load("mats_" + mid_or_sst + "_" + bsl_or_fu2 + ".mat", ...
        "mats_" + mid_or_sst)
    if mid_or_sst == "mid"
        mat = mats_mid;
    else
        mat = mats_sst;
    end
end

function ret = get_master_data()
    addpath('G:\.shortcut-targets-by-id\1Y42MQjJzdev5CtNSh2pJh51BAqOrZiVX\IMAGEN\behav_variables')
    ret = readtable("IMAGEN_master_data_sheet_2022-01-03.xlsx");
end

function ret=w_dist(xi, xj)
    ret = Metrics.Wasserstein.distance(xi,xj);
end

function write_dist_matrix()
    gcp; % create a default parallel pool if none exists
    [mat_t1, mat_t2] = load_IMAGEN("mid");
    m_data = get_master_data();
    
    good_i = [];
    for i = 1:size(mat_t1, 3)
        if ~(any(any(isnan(mat_t1(:,:,i))))) && ...
           ~(any(any(isnan(mat_t2(:,:,i))))) && ...
            m_data{i,"audit_c_valid_bsl"} == 1 && ...
            m_data{i,"audit_c_valid_fu2"} == 1
            good_i = [good_i i];
        end
    end
    
    %good_i = good_i(1:100);
    good_mat_1 = mat_t1(:,:,good_i);
    good_mat_2 = mat_t2(:,:,good_i);
    good_mat_all = cat(3,good_mat_1, good_mat_2);
    
    distances_ut = zeros(length(good_mat_all), length(good_mat_all));
    
    hundred_sixteen_time = 110;
    time_est = hundred_sixteen_time * (length(good_mat_all)/100)^2 * 16/gcp().NumWorkers;
    disp(["Wass. time ~ " num2str(round(time_est)) " seconds"]);
    for i = 1:length(good_mat_all)
        mat_i = good_mat_all(:,:,i);
        parfor j = (i+1):length(good_mat_all)
            mat_j = good_mat_all(:,:,j);
            if i < j
                distances_ut(i,j) = Metrics.Wasserstein.distance(mat_i, mat_j);
            end
        end
    end

    distances = distances_ut + tril(distances_ut',-1);
    writematrix(distances, "E:\\dev\\Repos\\matrix-metrics\\dm_temp.txt");

    good_audit_1 = m_data{good_i, "audit_c_total_score_bsl"};
    good_audit_2 = m_data{good_i, "audit_c_total_score_fu2"};
    good_audit = [good_audit_1 good_audit_2];
    
    writematrix(good_audit, "E:\\dev\\Repos\\matrix-metrics\\b_temp.txt")
end


function main()
    write_dist_matrix()
    distances = readmatrix("E:\\dev\\Repos\\matrix-metrics\\dm_temp.txt");
    good_audit = readmatrix("E:\\dev\\Repos\\matrix-metrics\\b_temp.txt");
    
    options = statset('MaxIter', 10000);

    n_is = 70;
    stresses = zeros(n_is);
    rsquareds = zeros(n_is);
    rsquareds_adj = zeros(n_is);
    is = 1:n_is;
    for i = is
        disp(i);
        [X,stress] = mdscale(distances, i, 'Criterion','metricstress', 'Options', options);
        stresses(i) = stress;
        
        b = fitlm(X, good_audit);
        rsquareds(i) = b.Rsquared.Ordinary;
        rsquareds_adj(i) = b.Rsquared.Adjusted;
        disp(b);
    end

    plot(is, stresses);
    hold on;
    plot(is, rsquareds);
    plot(is, rsquareds_adj);
    legend("stress", "r squared", "r squared (adj)");
    hold off;

end