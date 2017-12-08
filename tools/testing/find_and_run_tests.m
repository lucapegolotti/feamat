function [passed, failed, incomplete] = find_and_run_tests(this_fdr)
% Find all tests in all subdirectories and runs them

% check if there are tests in this directory
tests = dir([this_fdr,'/test_*.m']);

num_tests = size(tests,1);

passed = 0;
failed = 0;
incomplete = 0;

for i = 1:num_tests
    results = runtests(tests(i).name);
    % display result in nicely arranged table
    num_subtests = size(results,2);
    
    for j = 1:num_subtests
        if (results(j).Passed == 1)
            passed = passed + 1;
        else
            warning(['Attention: ', results(j).Name,' did not pass!']);
        end
        
        if (results(j).Failed == 1)
            failed = failed + 1;
        end
        
        if (results(j).Incomplete == 1)
            incomplete = incomplete + 1;
        end
        
    end
    
    table(results)
end

% find all subdirectories
allfiles = dir(this_fdr);
sub_dirs_flags = [allfiles.isdir];
sub_dirs = allfiles(sub_dirs_flags);

num_subdirs = size(sub_dirs,1);

% call recursively this function on all the subdirectories
% we start by 3 because 1 and 2 are reserved to the special dirs . and ..
for i = 3:num_subdirs
    [pa,fa,in] = find_and_run_tests([this_fdr,'/',sub_dirs(i).name]);
    passed = passed + pa;
    failed = failed + fa;
    incomplete = incomplete + in;
    
end

end

