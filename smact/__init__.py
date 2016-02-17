import os # get correct path for datafiles when called from another directory
module_directory = os.path.abspath( os.path.dirname(__file__))
data_directory = os.path.join(module_directory, 'data')
