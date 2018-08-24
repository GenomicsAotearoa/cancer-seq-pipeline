import sys
import re

'''
usage:
bgzip -d -c test.vcf.gz | python processVCF.py 
'''
sampleName = sys.argv[1]
#print('SAMPLE\tCHR\tPOS\tREF\tALT\tSS\tGENE\tGENE_REGION\tGENE_EXONFUC\tNRD\tNAD\tTRD\tTAD\tFILTER')

for line in sys.stdin:
    if line.startswith('#')==False:
        CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,NORMAL,TUMOR = line.split('\t')
        
        # Extract gene function information
        SS_search = re.search('SS=([0-9]);', INFO)
        GENE_search = re.search('Gene.refGene=(.*);GeneDetail.refGene', INFO)
        GENE_REGION_search = re.search('Func.refGene=(.*);Gene.refGene', INFO)
        GENE_EXONFUC_seach = re.search('ExonicFunc.refGene=(.*);AAChange.refGene=', INFO)
        
        SUPP_search = re.search('SUPPLEMENTARY', INFO)
        
        SUPP = '0'
        SS = SS_search.group(1)
        GENE = GENE_search.group(1)
        GENE_REGION = GENE_REGION_search.group(1)
        GENE_EXONFUC = GENE_EXONFUC_seach.group(1)
        if SUPP_search:
            SUPP = '1'
        
        NGT,NGQ,NDP,NRD,NAD,NFREQ,NDP4 = NORMAL.split(':')
        TGT,TGQ,TDP,TRD,TAD,TFREQ,TDP4 = TUMOR.split(':')
        
        if GENE_REGION == 'exonic':
            print(sampleName + '\t' + CHROM + '\t' + POS + '\t' + REF + '\t' + ALT + '\t' + SS + '\t' + SUPP + '\t' + GENE + '\t' + GENE_REGION + '\t' + GENE_EXONFUC + '\t' + NRD + '\t' + NAD + '\t' + TRD + '\t' + TAD + '\t' + FILTER)