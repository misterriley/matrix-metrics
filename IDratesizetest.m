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
    mat_t1 = load_IMAGEN_file(mid_or_sst, "bsl");

    disp("loading fu2")
    mat_t2 = load_IMAGEN_file(mid_or_sst, "fu2"); 
end

function main()
    [mt1, mt2] = load_IMAGEN("mid");

    good_i = [];
    for i = 1:size(mt1,3)
        if ~(any(any(isnan(mt1(:,:,i)))) || ...
                any(any(isnan(mt2(:,:,i)))))
            good_i = [good_i, i];
        end
    end
    n_good = length(good_i);

    utmask = triu(ones(size(mt1, 1)), 1) == 1;
    corr_matrix = zeros(n_good, n_good);
    for i = 1:size(corr_matrix, 1)
        m1index = good_i(i);
        m1 = mt1(:,:,m1index);
        ut1 = m1(utmask);
        for j = 1:size(corr_matrix, 2)
            m2index = good_i(j);
            m2 = mt2(:,:,m2index);
            ut2 = m2(utmask);

            corr_matrix(i,j) = corr(ut1,ut2);
        end
        fprintf("%d/%d\n", i, n_good)
    end

    n_repeats = 1000;
    sizes = [2, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900];
    for sz = sizes
        total_ids = 0;
        for r = 1:n_repeats
            sample_i = datasample(1:n_good,sz,'Replace',false);
            sample_matrix = corr_matrix(sample_i, sample_i);
            max_per_row = max(sample_matrix, [], 2);
            diagonals = diag(sample_matrix);
            n_ids = sum(max_per_row == diagonals);
            total_ids = total_ids + n_ids;
        end
        disp([sz, total_ids/(sz * n_repeats)])
    end
end