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
    Metrics.Angular_On_Z.run_test(mat_t1, mat_t2, origin)    
    Metrics.Angular.run_test(mat_t1, mat_t2, origin)
    Metrics.Euclidean_On_Z.run_test(mat_t1, mat_t2, origin)
    Metrics.Euclidean.run_test(mat_t1, mat_t2, origin)
    Metrics.Hellinger.run_test(mat_t1, mat_t2, origin)
    Metrics.Rao.run_test(mat_t1, mat_t2, origin)
    Metrics.Wasserstein.run_test(mat_t1, mat_t2, origin)
end

function test_IMAGEN()
    [mat_t1, mat_t2] = load_IMAGEN("mid");
    run_all_tests(mat_t1, mat_t2, "IMAGEN mid")
    
    [mat_t1, mat_t2] = load_IMAGEN("sst");
    run_all_tests(mat_t1, mat_t2, "IMAGEN sst")
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
         file_name = fullfile(directory, short_name);
    end

end

function ret = load_mls(data_folder)
    directory = 'E:\dev\Repos\matrix-metrics\' + data_folder;
    time_index = 1;
    time_maps = []
    while true
        time_map = load_mls_by_time(directory, time_index)
        if size(time_map) == 0
            break
        end
    end

    
end

function test_mls()
    ret = load_mls("matrices");
    %run_all_test(mat_t1, mat_t2, "mls")
end

function main()
    %test_IMAGEN()
    test_mls()
end
