import csv
import requests #sudo easy_install pip & sudo pip install requests
import json

#List of diagnosis:
lyst=[]
with open('diagnosis.csv','r') as f:
     text = f.read().splitlines()
     for line in text:
         lyst.append(line)
#print lyst

#PREDICT AN ONTOLOGY FOR A TEXT VALUE e.g. diagnose
#QUERING API from ZOOMA; request from http://www.ebi.ac.uk/spot/zooma/v2/api
nomatch = []
matchlist = []

for entry in lyst: #Include the different diagnose in the API
    match = {}
    diagnose = {'propertyValue': entry}
    r = requests.get('http://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate', params=diagnose)
    x = r.json()

    if len(x) == 0:
        nomatch.append(entry) #200
        #print 'Empty result'
    else:
        match['originalEntry'] = entry #256
        match['zoomaReplay'] = x[0]['semanticTags']
        match['confidence'] = x[0]['confidence'] #HIGH=142, GOOD=92, MEDIUM=22, LOW=0
        matchlist.append(match)
#print matchlist

#QUERING API from OxO:
mappinglist = []

for replay in matchlist:
    mapping = {}
    mapping['originalEntry'] = replay['originalEntry']

    for i in replay['zoomaReplay']:
        identifier = {'ids': [(i.rsplit('/',1)[1]).replace("_",":")], 'mappingTarget':['NCIt'], 'distance':2}
        a = requests.post('http://www.ebi.ac.uk/spot/oxo/api/search', data=identifier)
        y = a.json()
        #print json.dumps(y, indent=2) #JSON visualization
        #if y['_embedded']['searchResults']:
        if y['_embedded']['searchResults'][0]['mappingResponseList']:
            mapping['NCItcode'] = y['_embedded']['searchResults'][0]['mappingResponseList'][0]['curie']
            mapping['NCItlabel'] = y['_embedded']['searchResults'][0]['label']
            mappinglist.append(mapping)

print(mappinglist)

#Export results #218 mappings out of 449
a=[]
b=[]
c=[]

for i in range (0, len(mappinglist)):
    a.append(mappinglist[i]['NCItlabel'])
    b.append(mappinglist[i]['NCItcode'])
    c.append(mappinglist[i]['originalEntry'])

lis=(a,b,c)
final=zip(*lis)
with open('results.csv','wb') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
    writer.writerows(final)
