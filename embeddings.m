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
    sep = "\\";
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

function main()
    gcp;
    metricname = "Wasserstein";
    dataset = "IMAGEN";
    timepoint = "t1";

    dm = read_distance_matrix(metricname, dataset, timepoint);

    parfor dim = 1:200
        disp(strcat("starting  ", int2str(dim)));
        drawnow('update');
        build_and_write_embedding(dm, metricname, dataset, timepoint, dim);
        disp(strcat("finishing ", int2str(dim)));
        drawnow('update');
    end
end