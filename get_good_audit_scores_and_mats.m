function [good_i, scores, mats] = get_good_audit_scores_and_mats(timepoint, mat_t1, mat_t2)
    if length(size(mat_t1)) <3 && length(size(mat_t2)) < 3
        [mat_t1, mat_t2] = load_IMAGEN("mid", contains(timepoint, "t1"), contains(timepoint, "t2"));
    end
    m_data = get_master_data();

    good_i = 1:size(m_data, 1);
    if length(size(mat_t2)) == 3
        good_i = remove_bad_indices(good_i, mat_t2, m_data, "audit_c_valid_fu2");
    end

    if length(size(mat_t1)) == 3
        good_i = remove_bad_indices(good_i, mat_t1, m_data, "audit_c_valid_bsl");
    end    

    good_t1 = m_data{good_i, "audit_c_total_score_bsl"};
    good_t2 = m_data{good_i, "audit_c_total_score_fu2"};
    if timepoint == "t1"
        scores = good_t1;
        mats = mat_t1(:,:,good_i);
    elseif timepoint == "t2"
        scores = good_t2;
        mats = mat_t2(:,:,good_i);
    else
        scores = cat(1, good_t1, good_t2);
        mats = cat(3, mat_t1(:,:,good_i), mat_t2(:,:,good_i));
    end
end

function ret = remove_bad_indices(index_list, mats, audit_scores, label)
    ret = index_list;
    for j = 1:length(index_list)
        i = index_list(j);
        if (any(any(isnan(mats(:,:,i))))) || ...
            audit_scores{i,label} ~= 1
            ret = ret(ret ~= i);
        end
    end
end

function ret = get_master_data()
    addpath('./mats')
    ret = readtable("IMAGEN_master_data_sheet_2022-01-03.xlsx");
end