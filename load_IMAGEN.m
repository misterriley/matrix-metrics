function [mat_t1, mat_t2] = load_IMAGEN(mid_or_sst, bool_time_1, bool_time_2)
    addpath('./mats')
    
    mat_t1 = [];
    mat_t2 = [];

    if bool_time_1
        disp("loading bsl")
        mat = load_IMAGEN_file(mid_or_sst, "bsl");
        mat_t1 = fill_diags_with_ones(tanh(mat));
    end

    if bool_time_2
        disp("loading fu2")
        mat = load_IMAGEN_file(mid_or_sst, "fu2");
        mat_t2 = fill_diags_with_ones(tanh(mat));
    end
end

function ret = fill_diags_with_ones(mats)
    ret = mats;
    for i = 1:size(mats, 3)
        for j = 1:size(mats, 1)
            mats(j,j,i) = 1;
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