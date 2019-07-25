from zipfile import ZipFile
import requests
import time

#user="b.curran@auckland.ac.nz"
#key = "67ad7874be54f64960a4";
#url = "https://www.cancergenomeinterpreter.org/api/v1/"
#uk= user+" "+key

#r = requests.get('https://www.cancergenomeinterpreter.org/api/v1', headers=headers)
#jobs=r.json()
#print(jobs)
#for job in jobs:
#    print(job)
#    r = requests.delete('https://www.cancergenomeinterpreter.org/api/v1/'+jobId, headers=headers)
#    print(r.json())


def getDrugData(patientId):
#print("submitting job")
	print(patientId)
	user="b.curran@auckland.ac.nz"
	key = "67ad7874be54f64960a4";
	url = "https://www.cancergenomeinterpreter.org/api/v1/"
	uk= user+" "+key

	headers = {'Authorization': user+" "+key}
	payload = {'cancer_type': 'CANCER', 'title': 'gis it a go'}

	r = requests.post('https://www.cancergenomeinterpreter.org/api/v1', headers=headers, files={'mutations': open('mutations.tsv', 'rb')}, data=payload)
	jobId = r.json()

	payload={'action':'logs'}
	status='new'
	while(status != "Done"):
		r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/'+jobId, headers=headers, params=payload) 
   		if r.json()['status'] != "Done":
			status=r.json()['status']
       		 	time.sleep(30)
	    	else:
        		status = "Done"

	payload={'action':'download'}
	r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/'+jobId, headers=headers, params=payload )

	with open(patientId+".zip", 'wb') as fd:
		fd.write(r._content)

	with ZipFile( patientId+".zip", 'r') as zip:
		zip.extract('mutation_analysis.tsv',"../data/intermediate/"+patientId)
		zip.extract('drug_prescription.tsv',"../data/intermediate/"+patientId)

	r = requests.delete('https://www.cancergenomeinterpreter.org/api/v1/'+jobId, headers=headers)


