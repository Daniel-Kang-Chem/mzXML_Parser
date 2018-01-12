from pyopenms import *

mzXML = MzXMLFile()
exp = MSExperiment()
mzXML.load("../data/B2_IDA.mzXML",exp)
mzXML.store("../data/B2_IDA.mzXML",exp)
