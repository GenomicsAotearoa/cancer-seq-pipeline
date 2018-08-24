import vcf


vcf_reader = vcf.Reader(open('../data/varscan-P1003C.indel.vcf', 'r'))
#print(vcf_reader.filters)
functions=["exonic","splicing","UTR5","upstream"]
for record in vcf_reader:
    #print(type(record.INFO['Func.refGene']))
    #print(record.INFO['Func.refGene'][0])
    if record.INFO['SS']=="1" and record.INFO['Func.refGene'][0] in functions:# and record.INFO['wgEncodeDacMapabilityConsensusExcludable'] == "." :
        print (record.INFO['wgEncodeDacMapabilityConsensusExcludable'][0])
    

  #filter(Func.refGene %in% c("exonic", "splicing", "UTR5", "upstream")) %>% 
  #filter(wgEncodeDacMapabilityConsensusExcludable == ".")
