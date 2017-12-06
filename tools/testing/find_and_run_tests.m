function find_and_run_tests(this_fdr)
% Find all tests in all subdirectories and runs them

% check if there are tests in this directory
tests = dir([this_fdr,'/test_*.m']);

num_tests = size(tests,1);

for i = 1:num_tests
    results = runtests(tests(i).name);
    % display result in nicely arranged table
    num_subtests = size(results,2);
    for j = 1:num_subtests
        if (results(j).Passed ~= 1)
            warning(['Attention: ', results(j).Name,' did not pass!']);
        end
    end

    table(results)
end


allfiles = dir(this_fdr);
sub_dirs_flags = [allfiles.isdir];
sub_dirs = allfiles(sub_dirs_flags);

num_subdirs = size(sub_dirs,1);

% we start by 3 because 1 and 2 are reserved to the special dirs . and ..
for i = 3:num_subdirs
    find_and_run_tests([this_fdr,'/',sub_dirs(i).name]);
end


end

