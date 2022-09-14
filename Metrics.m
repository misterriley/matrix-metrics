classdef Metrics
    properties
        sim_count = 10000;
        print_every = 50000;
    end
    
    enumeration
        Angular,
        Angular_On_Z,
        Euclidean,
        Euclidean_On_Z,
        Hellinger,
        Wasserstein, 
        Rao
    end

    methods
        function run_test(obj, mat_t1, mat_t2, origin)
            good_i = [];
            for i = 1:size(mat_t1,3)
                if ~(any(any(isnan(mat_t1(:,:,i)))) || ...
                        any(any(isnan(mat_t2(:,:,i)))))
                    good_i = [good_i, i];
                end
            end
            
            n_closer = 0;
            n_farther = 0;
            while n_closer + n_farther < obj.sim_count
                is = datasample(good_i, 2);
            
                if is(1) == is(2)
                    continue
                end
            
                n = n_closer + n_farther;
                if mod(n + 1, obj.print_every) == 0
                    disp(n + 1)
                end
            
                x1 = mat_t1(:,:,is(1));
                x2 = mat_t2(:,:,is(1));
                y2 = mat_t2(:,:,is(2));
            
                dii = distance(obj, x1, x2);
                dij = distance(obj, x1, y2);
            
                if dii < dij
                    n_closer = n_closer + 1;
                elseif dii > dij
                    n_farther = n_farther + 1;
                else
                    n_farther = n_farther + .5;
                    n_closer = n_closer + .5;
                end
            end
            p = n_closer / obj.sim_count;
            q = 1 - p;
            se = sqrt(p * q / obj.sim_count);
            fprintf('%s %s\n', [string(obj) string(origin)])
            fprintf('p = %f +/- %f\n', [p se]);
        end

        function d = distance(obj, x1, x2)
            switch obj
                case Metrics.Euclidean
                    [ut1, ut2] = obj.triu_as_vector(x1, x2);
                    d = norm(ut1 - ut2);
                case Metrics.Euclidean_On_Z
                    [ut1, ut2] = obj.triu_as_vector(x1, x2);
                    [z1, z2] = obj.r_to_z(ut1, ut2);
                    d = norm(z1 - z2);
                case Metrics.Angular
                    [ut1, ut2] = obj.triu_as_vector(x1, x2);
                    d = acos(corr(ut1, ut2))/pi; % corr is not a proper metric
                case Metrics.Angular_On_Z
                    [ut1, ut2] = obj.triu_as_vector(x1, x2);
                    [z1, z2] = obj.r_to_z(ut1, ut2);
                    d = acos(corr(z1, z2))/pi;
                case Metrics.Hellinger
                    det1 = det(x1);
                    det2 = det(x2);
                    det12avg = det((x1 + x2)/2);
                    dsquared = 1 - ((det1*det2)^.25)/(det12avg^.5);
                    d = dsquared^.5;
                case Metrics.Wasserstein
                    sqrtx2 = sqrtm(x2);
                    d = trace(x1 + x2 - 2*sqrtm(sqrtx2 * x1 * sqrtx2));
                case Metrics.Rao
                    invsqrtx2 = (x2)^(-.5);
                    s = svd(invsqrtx2 * x1 * invsqrtx2);
                    logs = log(s);
                    d = sum(logs.^2);
            end
        end

        function [z1, z2] = r_to_z(~, x1, x2)
            z1 = atanh(x1);
            z2 = atanh(x2);
        end

        function [ut1, ut2] = triu_as_vector(~, x1,x2)
            utmask = triu(ones(size(x1, 1)), 1) == 1;
            ut1 = x1(utmask);
            ut2 = x2(utmask);
        end
    end
end
