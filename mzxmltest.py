from pyteomics import mzxml, pepxml
import matplotlib.pyplot as plt

mzxmlFile = mzxml.read("B2_IDA.mzXML")
#pepxmlFile = pepxml.read("pepXML/B2_IDA_comet_tandem_irt_scored_ipro.pep.xml")
experiments = []
rt = 30.0

def filter_by_property(xmlfile, property, value):
  curScan = xmlfile.next()
  matched_scans = []
  if property == 'retention_time':
    filter_property = 'retention_time_sec'
  if property == 'charge': 
    filter_property = 'assumed_charge'
  while curScan:
    if (int(curScan[filter_property]) == value):
      matched_scans.append(curScan['index'])
    else:
        print('nope' + str(int(curScan[filter_property])))
    if (int(curScan[filter_property]) > value):
        break
    curScan = xmlfile.next()
  return matched_scans

def mzxml_filter_by_property(mzxml, property, value):
  curScan = mzxml.next()
  matched_scans = []
  if property == 'retention_time':
    filter_property = 'retentionTime'
  while curScan:
    if (int(curScan[filter_property]) == value):
      matched_scans.append(curScan['id'])
      print('we here')
      print(curScan['m/z array'])
      plt.plot(curScan['m/z array'], curScan['intensity array'])
      plt.show()
      print("found!")
    if (int(curScan[filter_property]) > value):
        break
    curScan = mzxml.next()
  return matched_scans

print(mzxml_filter_by_property(mzxmlFile, 'retention_time', rt))
#print(filter_by_property(pepxmlFile, 'retention_time', rt))
