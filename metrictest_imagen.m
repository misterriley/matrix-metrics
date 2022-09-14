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

function load_Child_Risk()
    addpath("C:/Users/miste/developer/Repos/scratch/Child_Risk/Data")
end

function test_Child_Risk()
    [mat_t1, mat_t2] = load_Child_Risk();
    run_all_test(mat_t1, mat_t2, "Child Risk")
end

function main()
    test_IMAGEN()
    test_Child_Risk()
end
