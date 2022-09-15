main()

function load_Child_Risk()
    addpath("C:/Users/miste/developer/Repos/scratch/Child_Risk/Data")
end

function test_Child_Risk()
    [mat_t1, mat_t2] = load_Child_Risk();
    run_all_test(mat_t1, mat_t2, "Child Risk")
end

function main()
    test_Child_Risk()
end