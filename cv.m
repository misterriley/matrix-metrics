%main_cv()

function mat = load_IMAGEN_file(mid_or_sst, bsl_or_fu2)
    load("mats_" + mid_or_sst + "_" + bsl_or_fu2 + ".mat", ...
        "mats_" + mid_or_sst)
    if mid_or_sst == "mid"
        mat = mats_mid;
    else
        mat = mats_sst;
    end
end

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

function ret = get_master_data()
    addpath('G:\.shortcut-targets-by-id\1Y42MQjJzdev5CtNSh2pJh51BAqOrZiVX\IMAGEN\behav_variables')
    ret = readtable("IMAGEN_master_data_sheet_2022-01-03.xlsx");
end

function ret = get_good_audit_scores(timepoint)
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
    
    good_t1 = m_data{good_i, "audit_c_total_score_bsl"};
    good_t2 = m_data{good_i, "audit_c_total_score_fu2"};
    if timepoint == "t1"
        ret = good_t1;
    elseif timepoint == "t2"
        ret = good_t2;
    else
        ret = cat(1, good_t1, good_t2);
    end
end

function filepath = get_generic_filepath(folder, metricname, dataset, filename)
    sep = "/";
    filepath = strcat(".", sep, folder, ...
        sep, metricname, ...
        sep, dataset, ...
        sep, filename);
end

function filepath = get_distance_filepath(metricname, dataset, timepoint)
    filename = strcat(timepoint, ".txt");
    filepath = get_generic_filepath("distances", metricname, dataset, filename);
end

function ret = read_distance_matrix(metricname, dataset, timepoint)
    filepath = get_distance_filepath(metricname, dataset, timepoint);
    ret = readmatrix(filepath);
end

function ret = get_preds(i_train, i_test, distmatrix, good_audit, dim)
    distmatrix_train = distmatrix(i_train, i_train);
    [L, ~] = cmdscale(distmatrix_train, dim); % classical mds for embedding the landmarks
    L_sharp = pinv(L); % pseudoinverse of embedding matrix

    delta_as = distmatrix(i_train, i_test).^2;
    delta = distmatrix_train.^2;
    delta_mu = mean(delta, 2); % average of the columns
    x_as = 1/2 * L_sharp * (repmat(delta_mu, 1, sum(i_test)) - delta_as);
   
    lm = fitlm(L, good_audit(i_train));
    ret = lm.predict(x_as');
end

function graph_cv(metricname, dataset, timepoint, dim_seq, repetitions, folds)
    rhos = zeros(size(dim_seq));
    good_audit = get_good_audit_scores(timepoint);
    
    dist_matrix = read_distance_matrix(metricname, dataset, timepoint);
    parfor k = 1:length(dim_seq)
        dim = dim_seq(k);
        cv_outputs = [];
        for t = 1:repetitions
            
            preds = [];
            values = [];
            c = cvpartition(length(good_audit), 'KFold', folds);
            for i = 1:folds
                
                preds = cat(1, preds, get_preds(training(c,i), test(c,i), dist_matrix, good_audit, dim));
                values = cat(1, values, good_audit(test(c, i)));
            end
            [RHO,~] = corr(preds,values,'Type','Spearman');
            cv_outputs = [cv_outputs RHO];
        end

        rhos(k) = mean(cv_outputs);
        disp(["index" k "dim" dim "rho" rhos(k)]);
    end

    plot(dim_seq, rhos);
end

function ret = get_seq(length, lmin, lmax)
    min_multiplier = 1;
    max_multiplier = 2;

    while(true)
        list = get_seq_from_multiplier(length, lmin, max_multiplier);
        if max(list) >= lmax
            break
        else
            max_multiplier = max_multiplier * 2;
            min_multiplier = min_multiplier * 2;
        end
    end

    while(true)
        test_multiplier = (min_multiplier + max_multiplier)/2;
        list = get_seq_from_multiplier(length, lmin, test_multiplier);
        if max(list) == lmax
            ret = list;
            break
        elseif max(list) > lmax
            max_multiplier = test_multiplier;
        else
            min_multiplier = test_multiplier;
        end
    end
end

function ret = get_seq_from_multiplier(length, min, multiplier)
    value = min;
    list = zeros(length,1);
    for i = 1:length
        if i ~= 1
            value = max(value + 1, value * multiplier);
        end
        list(i) = round(value);
    end
    ret = list';
end

function main_cv()
    metricname = "Wasserstein";
    dataset = "IMAGEN";
    timepoint = "t2";
    min_dim = 1;
    max_dim = 720;
    n_dims = 720; %best if it's a multiple of the number of threads available
    folds = 10;
    repetitions = 10;

    dim_seq = get_seq(n_dims, min_dim, max_dim);

    graph_cv(metricname, dataset, timepoint, dim_seq, repetitions, folds)
end