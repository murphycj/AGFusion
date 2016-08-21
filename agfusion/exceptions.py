
class DataBaseError(Exception):
    def __init__(self,e):
        Exception.__init__(self,e)
        self.e = e

class GeneIDException(Exception):
    def __init__(self,dErrArguments):
        Exception.__init__(self,"my exception was raised with arguments {0}".format(dErrArguments))
        self.dErrorArguments = dErrorArguements

class JunctionException(Exception):
    def __init__(self):
        Exception.__init__(self,"Junction not within gene boundaries")
