import sys, os, os.path
import numpy as np

""" Generic convert to type
:input: the variable to be converted
:type: type to convert to
:returns: the input converted to type
"""
def convert_to_type(input, type):
    try:
        if type == "string":
            return str(input)
        elif type == "int":
            return int(input)
        elif type == "float":
            return float(input)
        else:
            sys.exit("Unknown type " + str(type))
    except:
        sys.exit(str(input) + " cannot be converted to " + str(type))

""" Compares two csv files
:file: filename of test file
:gold_file: filename of gold file
:rel_tol: permissible relative error before asserting diff
:abs_tol: permissible absolute error before asserting diff
:returns: string, either "success" or error message
"""
def csv_diff(file, gold_file, rel_tol = 1.0e-6, abs_tol = 1.0e-6):
    # check if files exist
    if not os.path.isfile(file):
        return "file " + file + " not found"
    if not os.path.isfile(gold_file):
        return "gold file " + file + " not found"

    # check the data
    data = np.loadtxt(file, delimiter = ",", skiprows = 1)
    gold_data = np.loadtxt(gold_file, delimiter = ",", skiprows = 1)
    if len(gold_data.shape) == 2:
        if data.shape[0] != gold_data.shape[0]:
            return "Number of rows in file " + file + " (" + str(data.shape[0]) + " vs " + str(gold_data.shape[0]) + ") differs from gold file."
        if data.shape[1] != gold_data.shape[1]:
            return "Number of columns in file " + file + " (" + str(data.shape[1]) + " vs " + str(gold_data.shape[1]) + ") differs from gold file."
        for r in range(data.shape[0]):
            for c in range(data.shape[1]):
                if abs(data[r, c] - gold_data[r, c]) > abs_tol:
                    return "Entry (" + str(r) + "," + str(c) + ") in file " + file + " differs from gold file "\
                            + str(data[r, c]) + " != " + str(gold_data[r, c]) + " absolute diff: %8.4e" % (abs(data[r, c] - gold_data[r, c]))
                # relative check only if value is large enough
                if abs(gold_data[r, c]) > 1.0e-12:
                    if abs(1.0 - data[r, c] / gold_data[r, c]) > rel_tol:
                        return "Entry (" + str(r) + "," + str(c) + ") in file " + file + " differs from gold file "\
                                + str(data[r, c]) + " != " + str(gold_data[r, c]) + " relative diff: %8.4e" % (abs(1.0 - data[r, c] / gold_data[r, c]))
    else:
        if data.shape[0] != gold_data.shape[0]:
            return "Number of columns in file " + file + " (" + str(data.shape[0]) + " vs " + str(gold_data.shape[0]) + ") differs from gold file."
        for c in range(data.shape[0]):
            if abs(data[c] - gold_data[c]) > abs_tol:
                return "Entry (" + str(c) + ") in file " + file + " differs from gold file "\
                        + str(data[c]) + " != " + str(gold_data[c]) + " absolute diff: %8.4e" % (abs(data[c] - gold_data[c]))
            # relative check only if value is large enough
            if abs(gold_data[c]) > 1.0e-12:
                if abs(1.0 - data[c] / gold_data[c]) > rel_tol:
                    return "Entry (" + str(c) + ") in file " + file + " differs from gold file "\
                            + str(data[c]) + " != " + str(gold_data[c]) + " relative diff: %8.4e" % (abs(1.0 - data[c] / gold_data[c]))

    return "success"


""" A class to represent THOR cases. So far only CSV diffs are supported
:input_parameters: a dictionary storing parameter name -> value pairs
"""
class thor_test:

    parameters = {'directory' : 'string', 'name' : 'string', 'input' : "string",\
                  'gold_file_name' : "string", 'rel_tol' : "float", 'abs_tol' : "float", 'nproc' : "int"}
    defaults = {'directory' : '__required__', 'name' : '__required__', 'input' : "__required__",\
                'gold_file_name' : "__required__", 'rel_tol' : 1.0e-6, 'abs_tol' : 1.0e-6, 'nproc' : 1}

    def __init__(self, input_parameters):
        self._params = {}
        for p in input_parameters.keys():
            if p in self.parameters.keys():
                self._params[p] = convert_to_type(input_parameters[p], self.parameters[p])
            else:
                sys.exit("Unknown parameter named " + str(p))
        # work through defaults
        for p in self.defaults.keys():
            if not p in self._params.keys():
                if self.defaults[p] == "__required__":
                    sys.exit("Required parameter " + str(p) + " is missing")
                else:
                    self._params[p] = convert_to_type(self.defaults[p], self.parameters[p])

    """ Method for setting number of processors for this test
    :returns: None
    """
    def set_nproc(self, nproc):
        self._params['nproc'] = nproc

    """ Method for executing this THOR test
    :returns: csv diff results
    """
    def execute(self, exec_rel_path = ''):
        current_dir = os.getcwd()
        workdir = self._params['directory']
        nproc = self._params['nproc']
        if nproc == 1:
            command = os.getcwd() + "/" + exec_rel_path + "thor-1.0.exe" + " " + self._params['input'] + " > dummy"
        else:
            command = "mpiexec -n " + str(nproc) + " " + os.getcwd() + "/" + exec_rel_path + "thor-1.0.exe" + " " + self._params['input'] + " > dummy"
        test_file = workdir + self._params['gold_file_name']
        gold_file = workdir + "gold/" + self._params['gold_file_name']
        l = len(os.getcwd() + "/" + exec_rel_path)
        test_name = self._params['directory'][l:-1] + ":" + self._params['name']
        # execute command
        os.chdir(workdir)
        os.system(command)
        os.chdir(current_dir)
        # check csv diff and return
        res = csv_diff(test_file, gold_file, rel_tol = self._params['rel_tol'], abs_tol = self._params['abs_tol'])
        return test_name, res

""" Finds all files named test_input in subfolders
"""
def find_all_tests(root_dir = "."):
    all_test_files = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in [f for f in filenames if f == "test_input"]:
            all_test_files.append(os.path.join(dirpath, filename))
    return all_test_files

""" Parses a single test_input file and returns an array of tests to run
"""
def parse_test_file(file_name):
    directory = ""
    for j in range(len(file_name.split("/")) - 1):
        directory += file_name.split("/")[j] + "/"
    test_list = []
    input_parameters = {}
    lines = [line.rstrip('\n') for line in open(file_name)]
    reading_test = False
    for line in lines:
        if len(line.split()) > 0:
            if line.split()[0] == "start":
                if reading_test: sys.exit("Error when parsing " + file_name + " line " + line + " . Mismatch of start/end keywords?")
                reading_test = True
                input_parameters = {'directory' : directory}
                input_parameters['name'] = line.split()[-1]
            elif line.split()[0] == "end":
                if not reading_test: sys.exit("Error when parsing " + file_name + " line " + line + " . Mismatch of start/end keywords?")
                reading_test = False
                test_list.append(thor_test(input_parameters))
            else:
                if reading_test:
                    try:
                        if len(line.split("=")) != 2:
                            sys.exit("Something went wrong when parsing file " + file_name + " on line " + line +\
                                     ". Too many or too few entries.")
                        key = line.split("=")[0].strip()
                        val = line.split("=")[1].strip()
                        input_parameters[key] = val
                    except:
                        sys.exit("Something went wrong when parsing file " + file_name + " on line " + line)
                else:
                    # do nothing
                    pass
    return test_list
