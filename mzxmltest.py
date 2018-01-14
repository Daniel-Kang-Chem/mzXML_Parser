from pyteomics import mzxml,pepxml
mzxmlFile = mzxml.read("B2_IDA.mzXML")
pepxmlFile = pepxml.read("pepXML/B2_IDA_comet_tandem_irt_scored_ipro.pep.xml")
experiments = []
rt = 30.0
b = pepxmlFile.next()
while b:
 if (b['retention_time_sec'] == rt):
  print(b['index'])
  experiments.append(b['index'])
 else:
  print('sad bois' + str(b['retention_time_sec']))
 b = pepxmlFile.next()
print(experiments)

def filter(xmlfile, property, value):
 item = xmlfile.next()

 switch (property):
  case 'retention_time':
   filter_property = ('retention_time_sec',value)
   break
  case 'charge':
   filter_property = ('assumed_charge', value)
   break

 while item:
  if (item[filter_property[0]] == value):
   
