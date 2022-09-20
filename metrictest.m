main()

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

function run_all_tests(mat_t1, mat_t2, origin)
    fid = fopen("output.txt", "a+");
    fprintf(fid, "\n\t---\nstarting tests for %s", origin);
    fclose(fid);
    Metrics.Angular_On_Z.run_test(mat_t1, mat_t2, origin)    
    Metrics.Angular.run_test(mat_t1, mat_t2, origin)
    Metrics.Euclidean_On_Z.run_test(mat_t1, mat_t2, origin)
    Metrics.Euclidean.run_test(mat_t1, mat_t2, origin)
    %Metrics.Hellinger.run_test(mat_t1, mat_t2, origin)
    %Metrics.KL.run_test(mat_t1, mat_t2, origin)
    Metrics.Rao.run_test(mat_t1, mat_t2, origin)
    Metrics.Wasserstein.run_test(mat_t1, mat_t2, origin)
end

function test_IMAGEN()
    [mat_t1, mat_t2] = load_IMAGEN("mid");
    %run_all_tests(mat_t1, mat_t2, "IMAGEN mid t1 predict t2")
    run_all_tests(mat_t2, mat_t1, "IMAGEN mid t2 predict t1")
    
    [mat_t1, mat_t2] = load_IMAGEN("sst");
    %run_all_tests(mat_t1, mat_t2, "IMAGEN sst t1 predict t2")
    run_all_tests(mat_t2, mat_t1, "IMAGEN sst t2 predict t1")
end

function ret = load_mls_by_time(full_folder, time_index)
    file_pattern = fullfile(full_folder, ...
        "*_matrix_" + time_index + "_matrix.txt");
    data_files = dir(file_pattern);
    ret = containers.Map();
    for f_index = 1:length(data_files)
        short_name = data_files(f_index).name;
        if ~endsWith(short_name, '.txt')
            continue
        end
        file_name = fullfile(full_folder, short_name);
        subj_id = extractBefore(short_name, "_");

        mat = tanh(readmatrix(file_name));
        for i = 1:size(mat, 1)
            mat(i,i) = 1;
        end

        ret(subj_id) = mat;
    end
end

function ret = load_mls(data_folder)
    directory = 'E:\dev\Repos\matrix-metrics\' + data_folder;
    time_index = 1;
    time_maps = {};
    while true
        fprintf("loading mls data time %d from %s folder\n", ...
            time_index, data_folder)
        time_map = load_mls_by_time(directory, time_index);
        if isempty(time_map)
            fprintf("no data found for time %d, breaking\n", time_index)
            break
        end
        time_maps{time_index} = time_map;
        time_index = time_index + 1;
    end
    
    all_keys = [];
    for i = 1:length(time_maps)
        all_keys = [all_keys keys(time_maps{i})];
    end

    all_keys = unique(all_keys);
    ret = NaN(length(time_maps),268,268,length(all_keys));
    for time_index = 1:length(time_maps)
        time_map = time_maps{time_index};
        for key_index = 1:length(all_keys)
            key = all_keys{key_index};
            if time_map.isKey(key)
                ret(time_index,:,:,key_index) = time_map(key);
            end
        end
    end
end

function test_mls_directory(dir)
    ret = load_mls(dir);
    n_time_indices = size(ret,1);
    for i = 1:n_time_indices
        for j = 1:n_time_indices
            if i == j
                continue
            end
            mat_t1 = squeeze(ret(i,:,:,:));
            mat_t2 = squeeze(ret(j,:,:,:));
            src = sprintf("mls %s t%d predict t%d", dir, i, j);
            run_all_tests(mat_t1, mat_t2, src)
        end
    end
end

function test_mls()
    for dir = [%"matrices" ...
            "matrices_gng" ...
            "matrices_gng_motion_checked" %...
            %"matrices_reward" ...
            %"matrices_reward_motion_checked"
            ]
        test_mls_directory(dir)
    end
end

function main()
    test_IMAGEN()
    %test_mls()
end
