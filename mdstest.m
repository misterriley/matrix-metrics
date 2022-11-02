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

    good_audit_1 = m_data{good_i, "audit_c_total_score_bsl"};
    good_audit_2 = m_data{good_i, "audit_c_total_score_fu2"};
    good_audit = cat(1, good_audit_1, good_audit_2);
    writematrix(good_audit, "b_temp.txt")
    
    good_mat_1 = mat_t1(:,:,good_i);
    good_mat_2 = mat_t2(:,:,good_i);
    good_mat_all = cat(3,good_mat_1, good_mat_2);
    
    distances_ut = zeros(length(good_mat_all), length(good_mat_all));
    
    gcp; % create a default parallel pool if none exists
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
    writematrix(distances, "dm_temp.txt");
end

function test_vector_embedding()
    %write_dist_matrix()
    distances = readmatrix("dm_temp.txt");
    good_audit = readmatrix("b_temp.txt");
    n_subs = length(good_audit)/2;
    
    options = statset('MaxIter', 10000);

    n_is = 10;
    stresses = zeros(n_is);
    rsquareds = zeros(n_is);
    rsquareds_adj = zeros(n_is);
    is = 1:n_is;
    t1_audits = good_audit(1:n_subs);
    t2_audits = good_audit(n_subs+1:2*n_subs);
    audit_diffs = t2_audits - t1_audits;
    for i = is
        dim = fibonacci(i + 1);
        disp(dim);
        [X,stress] = mdscale(distances, dim, 'Criterion','metricstress', 'Options', options);
        stresses(i) = stress;

        t1_embeds = X(1:n_subs,:);
        t2_embeds = X(n_subs+1:2*n_subs,:);

        vecs = t2_embeds - t1_embeds;
        
        b = fitlm(vecs, audit_diffs);
        rsquareds(i) = b.Rsquared.Ordinary;
        rsquareds_adj(i) = b.Rsquared.Adjusted;
        disp(b);
    end

    display_graph(dims, stresses, rsquareds, rsquareds_adj);
end

function test_point_embedding(distances, good_audit)
    %write_dist_matrix()
    %distances = readmatrix("dm_temp.txt");
    %good_audit = readmatrix("b_temp.txt");
    
    options = statset('MaxIter', 10000);

    n_is = 3;
    stresses = zeros(n_is);
    rs = zeros(n_is);
    rs_adj = zeros(n_is);
    is = 1:n_is;
    dims = zeros(n_is);
    for i = is
        disp(i);
        if i == min(is)
            dims(i) = i;
        else
            dims(i) = max(i+1, int16(dims(i-1)*1.25));
        end
        [X,stress] = mdscale(distances, dims(i), 'Criterion','metricstress', 'Options', options);
        stresses(i) = stress;
        
        b = fitlm(X, good_audit);
        rs(i) = sqrt(b.Rsquared.Ordinary);
        rs_adj(i) = sqrt(b.Rsquared.Adjusted);
        disp(b);
    end

    display_graph(dims, stresses, rs, rs_adj);
end

function display_graph(dims, stresses, rsquareds, rsquareds_adj)
    plot(dims, stresses);
    hold on;
    plot(dims, rsquareds);
    plot(dims, rsquareds_adj);
    legend("stress", "r squared", "r squared (adj)");
    hold off;
end

function build_t1_and_t2_matrices()
    distances = readmatrix("dm_temp_t1_t2.txt");
    good_audit = readmatrix("b_temp_t1_t2.txt");
    n_subs = length(good_audit)/2;

    sq_1 = 1:n_subs;
    sq_2 = n_subs + 1:2*n_subs;

    d_t1 = distances(sq_1,sq_1);
    d_t2 = distances(sq_2,sq_2);

    a_t1 = good_audit(sq_1);
    a_t2 = good_audit(sq_2);

    writematrix(d_t1, "dm_temp_t1.txt");
    writematrix(d_t2, "dm_temp_t2.txt");

    writematrix(a_t1, "b_temp_t1.txt");
    writematrix(a_t2, "b_temp_t2.txt");
end

function main()
    distances = readmatrix("dm_temp_t2.txt");
    good_audit = readmatrix("b_temp_t2.txt");
    %test_vector_embedding()
    test_point_embedding(distances, good_audit)
    %build_t1_and_t2_matrices();
end