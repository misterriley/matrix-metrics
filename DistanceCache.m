classdef DistanceCache < handle
    properties
        distmatrix = [];
        metricname = "";
        dataset = "";
        timepoint = "";
    end
    methods
        function distance = getDistance(self, i, j)
            if i > size(self.distmatrix,1) || j > size(self.distmatrix,2)
                distance = -1;
            else
                distance = self.distmatrix(i,j);
            end
        end
        function setMatrixSize(self, n)
            self.distmatrix = -1*ones(n) + eye(n);
        end
        function setDistance(self, i, j, distance)
            self.distmatrix(i,j) = distance;
            self.distmatrix(j,i) = distance;
        end
        function save(self, metricname, dataset, timepoint)
            file = get_distance_filepath(metricname, dataset, timepoint);
            writematrix(self.distmatrix, file);
        end
        function load(self, metricname, dataset, timepoint)
            self.distmatrix = read_distance_matrix(metricname, dataset, timepoint);
        end
        function sc = subCache(self, indices)
            sc = DistanceCache();
            sc.distmatrix = self.distmatrix(indices, indices);
        end
    end
end

function ret = read_distance_matrix(metricname, dataset, timepoint)
    filepath = get_distance_filepath(metricname, dataset, timepoint);
    ret = readmatrix(filepath);
end

function filepath = get_distance_filepath(metricname, dataset, timepoint)
    filename = strcat(timepoint, ".txt");
    filepath = get_generic_filepath("distances", metricname, dataset, filename);
end

function filepath = get_generic_filepath(folder, metricname, dataset, filename)
    sep = "/";
    filepath = strcat(".", sep, folder, ...
        sep, metricname, ...
        sep, dataset, ...
        sep, filename);
end
