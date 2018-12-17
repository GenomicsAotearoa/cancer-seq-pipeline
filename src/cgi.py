from zipfile import ZipFile
import requests
import time

user="b.curran@auckland.ac.nz"
key = "67ad7874be54f64960a4";
url = "https://www.cancergenomeinterpreter.org/api/v1/"
uk= user+" "+key

#r = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
#jobs=r.json()
#print(jobs)
#for job in jobs:
#    print(job)
#    r = requests.delete('https://www.cancergenomeinterpreter.org/api/v1/'+jobId, headers=headers)
#    print(r.json())



#print("submitting job")
headers = {'Authorization': user+" "+key}
payload = {'cancer_type': 'CANCER', 'title': 'gis it a go'}

r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                headers=headers,
                files={
                        'mutations': open('mutations.tsv', 'rb')
                        },
                data=payload)
jobId = r.json()

#print("job id returned: "+jobId)


payload={'action':'logs'}
status='new'
while(status != "Done"):
    r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/'+jobId,
                headers=headers,
                params=payload
    ) 
    if r.json()['status'] != "Done":
        status=r.json()['status']
        #print(status)
        time.sleep(30)
    else:
        status = "Done"

#print(status)


payload={'action':'download'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/'+jobId,
                headers=headers,
                params=payload
                )

with open("file-"+jobId+".zip", 'wb') as fd:
    fd.write(r._content)

with ZipFile("file-"+jobId+".zip", 'r') as zip:
#    zip.printdir()
    zip.extract('mutation_analysis.tsv')
    zip.extract('drug_prescription.tsv')




#with open('file.zip', 'wb') as fd:
#    fd.write(r._content)
#httr::add_headers(
#
#    "Authorization" = paste(user, key, sep=" ",
#    mutations = "/projects/uoa00571/scriptes/reportGeneration/src/mutations")
#  )
#);
#POST(url, add_headers(Authorization =  paste(user, key, sep=" "), mutations = "/projects/uoa00571/scriptes/reportGeneration/src/mutations") , title = "Wibble" )
#
#POST(req)
r = requests.delete('https://www.cancergenomeinterpreter.org/api/v1/'+jobId, headers=headers)
#print(r.json())


