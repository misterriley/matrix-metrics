cv %include functions from cv and other previous bits of code

main()

function main()
    
    timepoint = "t1";
    metricname = "Wasserstein"

    good_audit = get_good_audit_scores(timepoint);
    dist_matrix = read_distance_matrix(metricname, dataset, timepoint);
end