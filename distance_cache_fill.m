speedtest()

main()

function speedtest()
    [~, landmarks_1] = load_IMAGEN("mid", true, true);
    indices_t1 = get_good_indices(landmarks_1);

    while true
        tic
        for i = 1:length(indices_t1)
            index_1 = indices_t1(i);
            mat1 = landmarks_1(:,:,index_1);
            index_2 = indices_t1(uint32(mod(7*i^2-1, length(indices_t1)))+1);
            mat2 = landmarks_1(:,:,index_2);
            d_squared = trace(mat1) + trace(mat2) - 2*trace(sqrtm(mat1*mat2));
            dist = real(sqrt(d_squared));
        end
        toc
    end
end

function main()
    metricname = "Wasserstein";
    
    [landmarks_1, landmarks_2] = load_IMAGEN("mid", true, true);
    indices_t1 = get_good_indices(landmarks_1);
    indices_t2 = get_good_indices(landmarks_2);

    dc = IMAGENDistanceCache();
    dc.load(metricname, "t1_t2_combined");

    pairs_to_calculate = [];
    disp("determining necessary calculations");

    for i_time = 1:2
        for j_time = 1:2
            pairs_to_calculate = cat(1, pairs_to_calculate, find_needed_calcs(i_time, j_time, indices_t1, indices_t2, dc));
        end
    end

    results = zeros(size(pairs_to_calculate,1),0);

    disp("calculating distances");
    parfor i = 1:length(pairs_to_calculate)
        x = pairs_to_calculate(i,:);
        
        results(i) = get_w_dist(x(1), x(2), x(3), x(4), landmarks_1, landmarks_2);
    end

    disp("finalizing and saving");
    for i = 1:length(pairs_to_calculate)
        disp(pairs_to_calculate(i,:));
        x = pairs_to_calculate(i,:);
        dc.setDistance(x(1), x(2), x(3), x(4), results(i));
    end

    dc.save();
end

function ret = find_needed_calcs(time_i, time_j, indices_1, indices_2, dc)
    disp([time_i, time_j]);
    ret = zeros(length(indices_1)*length(indices_2)/2,4);
    if time_i == 1
        indices_i = indices_1;
    else
        indices_i = indices_2;
    end

    if time_j == 1
        indices_j = indices_1;
    else
        indices_j = indices_2;
    end

    ret_index = 0;
    for i = indices_i
        for j = indices_j
            if i <= j
                continue
            end
            if dc.getDistance(i,time_i,j,time_j) == -1
                ret_index = ret_index + 1;
                ret(ret_index,:) = [i, time_i, j, time_j];
            end
        end
    end
    ret = ret(1:ret_index,:);
end

function indices = get_good_indices(mats)
    indices = [];
    for i = 1:size(mats, 3)
        if ~any(any(isnan(mats(:,:,i))))
            indices = [indices i];
        end
    end
end

function dist = get_w_dist(i, i_time, j, j_time, mats_t1, mats_t2)
    if i_time == 1
        mat1 = mats_t1(:,:,i);
    else
        mat1 = mats_t2(:,:,i);
    end

    if j_time == 1
        mat2 = mats_t1(:,:,j);
    else
        mat2 = mats_t2(:,:,j);
    end
    
    d_squared = trace(mat1) + trace(mat2) - 2*trace(sqrtm(mat1*mat2));
    dist = real(sqrt(d_squared));
end