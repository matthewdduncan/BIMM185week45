'''
Created on May 2, 2017

@author: Matthew Duncan
'''
from Bio import SeqIO
import gzip
import re

outgenes = open('files\\genestableatu.txt', 'w')
outfunctions = open('files\\functionsstableatu.txt', 'w')
outxref = open('files\\externaltableatu.txt', 'w')
outsynonyms = open('files\\synonymsstableatu.txt', 'w')
outreps = open('files\\repliconstableatu.txt', 'w')

taxid = ''
geneid = 4320
repid = 1;
genomid = 2;

#genbank files canhandle multiple datasets, this iterates through them (though there is only one here, for e coli k12)
for record in SeqIO.parse('atume.gbff', 'gb'):

	#each dataset has multiple features
	for feature in record.features:
	
		#gets the taxid for e coli
		if feature.type == 'source':
			taxid = feature.qualifiers['db_xref'][0]
			chrtype = feature.qualifiers.get('chromosome', [''])
			plasmtype = feature.qualifiers.get('plasmid', [''])
			if chrtype == '':
				type = 'plasmid'
			else:
				type = 'chromosome'
			repid += 1
			
			
			outreps.write(str(repid) + '\t' + str(genomid) + '\t' + type + '\t' + chrtype[0] + plasmtype[0] + '\n')
		#get information about each gene
		if feature.type == 'CDS':
			line= {}  #line is a dict where all information is stored.  this is necessary to beable to reorder everything.
			line['gene_id'] = geneid
			line['genome_id'] = genomid
			line['replicon_id'] = repid
			geneid += 1
			
			#can make a loop here instead of rewriting.  Extracts all relevant features from cds feature.qualifiers
			for x in ['protein_id', 'gene', 'locus_tag', 'gene_synonym', 'product', 'EC_number', 'db_xref']:
				first = True
				line[x] = ''
				
				#this loop adds multiple entries in cases where a feature has them.
				for ft in feature.qualifiers.get(x, '-'):
					if first:
						first = False
					else:
						line[x] += ','
					line[x] += ft
					
			#get features not in qualifiers
			line['coordinates'] = re.search('\[.*\]', str(feature.location)).group(0)
			line['strand'] = re.search('\(.*\)', str(feature.location)).group(0)[1]
			line['taxid'] = taxid
			line['chrtype'] = chrtype
			
			loc = re.findall('[0-9]+', line['coordinates'])
			length = int(loc[1]) - int(loc[0])
			line['length'] = length
			
			#print to output in the correct order, with tab spacing
			for i in ['gene_id', 'genome_id', 'replicon_id', 'locus_tag', 'gene', 'strand', 'length', 'product']:
				outgenes.write(str(line[i]) + '\t')
			outgenes.write('\n')
			
			for x in line['gene_synonym'].split('; '):
				if x == '-':
					continue
				outsynonyms.write(str(line['gene_id']) + '\t' + x + '\n')
				
			outfunctions.write(str(line['gene_id']) + '\t' + line['product'] + '\n')
			
			for x in line['db_xref'].split(','):
				if x == '-':
					continue
				outxref.write(str(line['gene_id']) + '\t' + x.split(':')[0] + '\t' + x.split(':')[1] + '\n')
				
			