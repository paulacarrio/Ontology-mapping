import csv
import requests #sudo easy_install pip & sudo pip install requests
import json
import urllib

#MAPPING ICDO MORPHOLOGY AND TOPOGRAPHY TO NCIt

# QUERING API OF ZOOMA. General request: http://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate
def getAnnotationFromZooma (value):
    diagnose = {'propertyValue': value}
    r = requests.get('http://snarf.ebi.ac.uk:8580/spot/zooma/v2/api/services/annotate', params=diagnose)
    x = r.json()
    if len(x) == 0:
        return None
    else:
        return x[0];

# QUERING API OF OXO
def getMappingFromOxO (uri):
    mapping = {}
    identifier = {'ids': [(uri.rsplit('/',1)[1]).replace("_",":")], 'mappingTarget':['NCIt'], 'distance':2}
    a = requests.post('http://www.ebi.ac.uk/spot/oxo/api/search', data=identifier)
    y = a.json()
    #print json.dumps(y, indent=2) #JSON visualization
    if y['_embedded']['searchResults'][0]['mappingResponseList']:
        mapping['NCItlabel'] = y['_embedded']['searchResults'][0]['label']
        print "mapping found! "+y['_embedded']['searchResults'][0]['mappingResponseList'][0]['curie']

        for i, val in enumerate(y['_embedded']['searchResults'][0]['mappingResponseList']):
            if val['curie'].replace('NCIt:C','').isnumeric():
                mapping['NCItcode'] = val['curie']
                print "returning mapping to "+val['curie']
                return mapping

    return None

# QUERING API OF OLS. http://snarf.ebi.ac.uk:8980/ols-beta/api/search?q=Squamous+cell+carcinoma&ontology=ncit&exact=true
def getMappingFromOLS (key):
    mapping = {}
    identifier = {'q': [key], 'ontology':['ncit'], 'exact':['True']}
    a = requests.get('http://snarf.ebi.ac.uk:8980/ols-beta/api/search', params=identifier)
    y = a.json()

    if y['response']['docs']:
        #print json.dumps(y, indent=2)
        mapping['NCItlabel'] = y['response']['docs'][0]['label']
        print "OLS mapping found " + y['response']['docs'][0]['label']
        mapping['NCItcode'] = y['response']['docs'][0]['short_form']
        print "OLS mapping found " + y['response']['docs'][0]['short_form']
        return mapping
    else:
        return None


# List of ICDO codes:
ICDO_M_T = {}
with open('input.csv','r') as csvfile:
    for line in csvfile.readlines():
        array = line.split(',')
        entry = array[0]
        nos_out = [' nos', ' NOS']
        for nos in nos_out:
            morph = array[1].replace(nos,'')
            top = array[3].replace(nos,'')
        icdmcode = array[2]
        icdtcode = array[4]
        Ncode = array[5]
        Nlabel = array[6]

        ICDO_M_T[entry+morph+top+Ncode+Nlabel] = { 'morph' : morph, 'top': top, 'icdmcode' : icdmcode, 'icdtcode' : icdtcode, 'Ncode' : Ncode, 'Nlabel' : Nlabel }

seen = {}
idMapping = {}

for entry in ICDO_M_T:
    ICDOM = ICDO_M_T[entry]['morph']
    ICDOT = ICDO_M_T[entry]['top']

#PIPELINE FOR MAPPING ICDOMORPHOLOGY TO NCIt:
    if ICDOM and ICDOM not in seen:
        print "Not seen before!"+ICDOM
        zoomaAnnotationICDOM = getAnnotationFromZooma(ICDOM)

        if not zoomaAnnotationICDOM:
            print "no match in Zooma for "+ICDOM

            print "searching ols for "+ICDOM+"..."
            ncitCode = getMappingFromOLS(ICDOM)
            if ncitCode:
                print "hit found in ols for "+ICDOM
                ICDO_M_T[entry]['morph_mapping'] = ncitCode

            seen[ICDOM] = None
        else:
            print "match in Zooma for "+ICDOM
            ICDO_M_T[entry]['morph_anno'] = zoomaAnnotationICDOM

            # check if NCIt term
            if zoomaAnnotationICDOM['derivedFrom']['provenance']['source']['uri'] == "http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl":
                tmpCode = zoomaAnnotationICDOM['semanticTags'][0]
                ICDO_M_T[entry]['morph_mapping'] = {'NCItcode' : "NCIt:"+tmpCode.rsplit('#',1)[1], 'NCItlabel' : ''}
            else:
                for uri in zoomaAnnotationICDOM['semanticTags']:
                    if uri not in idMapping:
                        ncitCode = getMappingFromOxO(uri)
                        if ncitCode:
                            ICDO_M_T[entry]['morph_mapping'] = ncitCode

            if 'morph_mapping' in ICDO_M_T[entry]:
                seen[ICDOM] = ICDO_M_T[entry]['morph_mapping']
            else:
                seen[ICDOM] = None
    elif ICDOM and seen[ICDOM]:
        ICDO_M_T[entry]['morph_mapping'] = seen[ICDOM]


#PIPELINE FOR MAPPING ICDOTOPOGRAPHY TO NCIt:
    if ICDOT not in seen:
        zoomaAnnotationICDOT = getAnnotationFromZooma(ICDOT)

        if not zoomaAnnotationICDOT:
            print "no match in Zooma for "+ICDOT

            print "searching ols for "+ICDOT+"..."
            ncitCode = getMappingFromOLS(ICDOT)
            if ncitCode:
                print "match found in OLS for "+ICDOT
                ICDO_M_T[entry]['morph_mapping'] = ncitCode

            seen[ICDOT] = None
        else:
            print "match for "+ICDOT
            ICDO_M_T[entry]['top_anno'] = zoomaAnnotationICDOT

            if zoomaAnnotationICDOT['derivedFrom']['provenance']['source']['uri'] == "http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl":
                tmpCode = zoomaAnnotationICDOT['semanticTags'][0]
                ICDO_M_T[entry]['top_mapping'] = {'NCItcode' : "NCIt:"+tmpCode.rsplit('#',1)[1], 'NCItlabel' : ''}
            else:
                for uri in zoomaAnnotationICDOT['semanticTags']:
                    if uri not in idMapping:
                        ncitCode = getMappingFromOxO(uri)
                        if ncitCode:
                             ICDO_M_T[entry]['top_mapping'] = ncitCode
            if 'top_mapping' in ICDO_M_T[entry]:
                seen[ICDOT] = ICDO_M_T[entry]['top_mapping']
            else:
                seen[ICDOT] = None

    elif seen[ICDOT]:
        ICDO_M_T[entry]['top_mapping'] = seen[ICDOT]


#CREATING A NEW ID ONTOLOGY FOR PAIRS OF ICD MORHPHOLOGY & TOPOGRAPHY
def getClassInMos(morphId, topId, morphNcit, topNcit):
    expression = "\nClass: {0} \n\
        Annotations: \n\
            rdfs:label \"{1}\" \n\
            ICD0TM:icdo_morphology_code \"{2}\"^^xsd:string,\n\
            ICD0TM:icdo_topography_code \"{3}\"^^xsd:string,\n\
                        \n\
        EquivalentTo: {4} and Thesaurus:R101 some {5}\n"

    id = urllib.quote(morphId+":"+topId)
    label = morphId+"/"+topId
    mnci = "Thesaurus:"+morphNcit
    tnci = "Thesaurus:"+topNcit
    x = [id, label, morphNcit, topNcit, mnci, tnci]
    return expression.format(*x)

# Export results
with open('output.csv', 'w') as f:
    csvw = csv.writer(f, delimiter=',', lineterminator='\n')

    for entry in ICDO_M_T:
        morph = ICDO_M_T[entry]['morph']
        morphCode = ICDO_M_T[entry]['morph_mapping']['NCItcode'].replace('NCIt:','') if 'morph_mapping' in ICDO_M_T[entry]  else ""
        morphCodeLabel = ICDO_M_T[entry]['morph_mapping']['NCItlabel'] if 'morph_mapping' in ICDO_M_T[entry] else ""
        top = ICDO_M_T[entry]['top']
        topCode = ICDO_M_T[entry]['top_mapping']['NCItcode'].replace('NCIt:','') if 'top_mapping' in ICDO_M_T[entry] else ""
        topCodeLabel = ICDO_M_T[entry]['top_mapping']['NCItlabel'] if 'top_mapping' in ICDO_M_T[entry]  else ""
        Ncode= ICDO_M_T[entry]['Ncode'].replace('"','')
        Nlabel = ICDO_M_T[entry]['Nlabel'].replace('"','')
        icdmcode = ICDO_M_T[entry]['icdmcode']
        icdtcode = ICDO_M_T[entry]['icdtcode']

        csvw.writerow([icdmcode, morph, morphCode, morphCodeLabel, icdtcode, top, topCode, topCodeLabel, Ncode, Nlabel])


seenNcitM = {}
seenNcitT = {}
classexpressions = []

for entry in ICDO_M_T:
    morph = ICDO_M_T[entry]['morph']
    morphCode = ICDO_M_T[entry]['morph_mapping']['NCItcode'].replace('NCIt:','') if 'morph_mapping' in ICDO_M_T[entry]  else ""
    morphCodeLabel = ICDO_M_T[entry]['morph_mapping']['NCItlabel'] if 'morph_mapping' in ICDO_M_T[entry] else ""
    top = ICDO_M_T[entry]['top']
    topCode = ICDO_M_T[entry]['top_mapping']['NCItcode'].replace('NCIt:','') if 'top_mapping' in ICDO_M_T[entry] else ""
    topCodeLabel = ICDO_M_T[entry]['top_mapping']['NCItlabel'] if 'top_mapping' in ICDO_M_T[entry]  else ""
    Ncode= ICDO_M_T[entry]['Ncode'].replace('"','')
    Nlabel = ICDO_M_T[entry]['Nlabel'].replace('"','')
    icdmcode = ICDO_M_T[entry]['icdmcode']
    icdtcode = ICDO_M_T[entry]['icdtcode']
    if morphCode and topCode:
        if morphCode not in seenNcitM:
            #print class decleration
            seenNcitM[morphCode] = 1
        if topCode not in seenNcitT:
            #print class decleration
            seenNcitT[topCode] = 1

            classexpressions.append(getClassInMos(icdmcode, icdtcode, morphCode,topCode))

    # Create an OWL file to read with Protege
    ontologyfile = open("icdo-ncit.owl","w")

    # Include prefix declerations and ontology and object property
    ontologyfile.write("Prefix: Thesaurus: <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#>")
    ontologyfile.write("Prefix: ICD0TM: <http://arraymap.org/ontology/icdo-tm-ncit.owl#>")
    ontologyfile.write("Prefix: rdfs: <http://www.w3.org/2000/01/rdf-schema#>")

    ontologyfile.write("Ontology: <http://arraymap.org/ontology/icdo-tm-ncit.owl#>")
    ontologyfile.write("AnnotationProperty: ICD0TM:icdo_morphology_code\n\
        Annotations: \n\
           rdfs:label \"ICD-O-3 Code morphology\"^^xsd:string")

    ontologyfile.write("AnnotationProperty: ICD0TM:icdo_topography_code\n\
       Annotations: \n\
          rdfs:label \"ICD-O-3 Code topography\"^^xsd:string")

    ontologyfile.write("ObjectProperty: Thesaurus:R101")

    # Include variant part "annotations":
    for i in classexpressions:
        ontologyfile.write(i)


ontologyfile.close()
