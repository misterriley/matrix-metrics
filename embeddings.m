main()

function filepath = get_distance_filepath(metricname, dataset, timepoint)
    filename = strcat(timepoint, ".txt");
    filepath = get_generic_filepath("distances", metricname, dataset, filename);
end

function filepath = get_embedding_filepath(metricname, dataset, timepoint, dimension)
    filename = strcat(timepoint, "_dim_", int2str(dimension), ".txt");
    filepath = get_generic_filepath("embeddings", metricname, dataset, filename);
end

function filepath = get_generic_filepath(folder, metricname, dataset, filename)
    sep = "/";
    filepath = strcat(".", sep, folder, ...
        sep, metricname, ...
        sep, dataset, ...
        sep, filename);
end

function ret = read_distance_matrix(metricname, dataset, timepoint)
    filepath = get_distance_filepath(metricname, dataset, timepoint);
    ret = readmatrix(filepath);
end

function build_and_write_embedding(dm, metricname, dataset, timepoint, dimension)
    filepath = get_embedding_filepath(metricname, dataset, timepoint, dimension);
    if not(isfile(filepath))
        options = statset('MaxIter', 10000);
        [X,~] = mdscale(dm, dimension, 'Criterion','metricstress', 'Options', options);
        writematrix(X,filepath);
    end
end

function ret = read_embedding(metricname, dataset, timepoint, dimension)
    filepath = get_embedding_filepath(metricname, dataset, timepoint, dimension);
    if isfile(filepath)
        ret = readmatrix(filepath);
    else
        ret = false;
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

function graph_embedding_regressions(metricname, dataset, timepoint, max_dim)
    good_audit = get_good_audit_scores(timepoint);
    rs = zeros(max_dim);
    rs_adj = zeros(max_dim);
    dims = 1:max_dim;
    for dim = dims
        X = read_embedding(metricname, dataset, timepoint, dim);
        if X == false
            continue;
        end
        b = fitlm(X, good_audit);
        rs(dim) = sqrt(b.Rsquared.Ordinary);
        rs_adj(dim) = sqrt(b.Rsquared.Adjusted);
        disp(dim);
        disp(b);
    end

    rs_adj = rs_adj(rs ~= 0);
    dims = dims(rs ~= 0);
    rs = rs(rs ~= 0);

    display_graph(dims, rs, rs_adj);
end

function write_embeddings(metricname, dataset, max_dim)
    gcp;
    %parpool(20);
    for timepoint = ["t2", "t1", "t1_t2"]
        dm = read_distance_matrix(metricname, dataset, timepoint);
    
        parfor dim = 1:max_dim
            disp(strcat("starting ", timepoint, " ", int2str(dim)));
            %drawnow('update');
            build_and_write_embedding(dm, metricname, dataset, timepoint, dim);
            disp(strcat("finished ", timepoint, " ", int2str(dim)));
            drawnow('update');
        end
    end
end

function ret = get_master_data()
    addpath('G:\.shortcut-targets-by-id\1Y42MQjJzdev5CtNSh2pJh51BAqOrZiVX\IMAGEN\behav_variables')
    ret = readtable("IMAGEN_master_data_sheet_2022-01-03.xlsx");
end

function display_graph(dims, rs, rs_adj)
    h1 = plot(dims, rs);
    hold on;
    h2 = plot(dims, rs_adj);
    h3 = plot(dims, repelem(max(rs_adj, [], 'all'), length(dims)));
    legend([h1(1), h2(1), h3(1)], "|r|", "|r (adj)|", "max |r (adj)|", "Location", "southeast");
    
    hold off;
end

function main()
    metricname = "Wasserstein";
    dataset = "IMAGEN";
    timepoint = "t2";
    max_dim = 300;
    %write_embeddings(metricname, dataset, max_dim);
    graph_embedding_regressions(metricname, dataset, timepoint, max_dim);
end