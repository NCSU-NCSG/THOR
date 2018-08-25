import sys, os
path = os.getcwd() + "/../../scripts"
sys.path.append(path)
import thor_test_harness as th

# check command line args; currently hacky version of -j nproc
n_global_proc = None
argv = sys.argv
argv.pop(0)
for j in range(len(argv)):
    if argv[j] == "-j":
        if len(argv) <= j + 1: sys.exit("The command -j must be followed by a number")
        n_global_proc = int(argv[j + 1])

# find all tests in sub folder
counter = 0
failures = 0
successes = 0
print "-"*80
all_test_files = th.find_all_tests(os.getcwd() + "/..")
for f in all_test_files:
    test_objects = th.parse_test_file(f)
    for t in test_objects:
        if n_global_proc != None:
            t.set_nproc(n_global_proc)
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
