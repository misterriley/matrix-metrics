classdef IMAGENDistanceCache < handle
    properties
        dc = DistanceCache();
        n = 0;
        timepoint = "";
        metricname = "";
    end
    methods
        function distance = getDistance(self, i, i_time_index, j, j_time_index)
           i_idx = IMAGENDistanceCache.get_private_index(i, i_time_index, self.n);
           j_idx = IMAGENDistanceCache.get_private_index(j, j_time_index, self.n);
           distance = self.dc.getDistance(i_idx, j_idx);
        end
        function setDistance(self, i, i_time_index, j, j_time_index, distance)
            i_idx = IMAGENDistanceCache.get_private_index(i, i_time_index, self.n);
            j_idx = IMAGENDistanceCache.get_private_index(j, j_time_index, self.n);
            self.dc.setDistance(i_idx, j_idx,distance);
        end
        function save(self)
            self.dc.save(self.metricname, "IMAGEN", self.timepoint);
        end
        function load(self, metricname, timepoint)
            self.timepoint = timepoint;
            self.metricname = metricname;
            self.dc.load(metricname, "IMAGEN", timepoint);
            self.n = size(self.dc.distmatrix,1) / 2;
        end
        function set_n(self, n)
            self.n = n;
            self.dc.setMatrixSize(2*n);
        end
    end
    methods(Static)
        function idx = get_private_index(index, time_index, n)
            idx = index + (time_index - 1) * n;
        end
    end
end