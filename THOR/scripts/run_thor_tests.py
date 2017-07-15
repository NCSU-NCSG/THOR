import sys, os
path = os.getcwd() + "/../../scripts"
sys.path.append(path)
import thor_test_harness as th

# find all tests in sub folder
counter = 0
failures = 0
successes = 0
print "-"*80
all_test_files = th.find_all_tests(os.getcwd() + "/..")
for f in all_test_files:
    test_objects = th.parse_test_file(f)
    for t in test_objects:
        counter += 1
        name, result = t.execute('../')
        if result == "success":
            successes += 1
            print "Test ", counter, name, "success"
        else:
            failures += 1
            print "Test ", counter, name, "failure"
            print result

print "-"*80
print 'Successes: ', successes, '          ', 'Failures: ', failures
