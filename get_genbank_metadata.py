# Script to scrape genbank for geographical & clinical metadata of CMV sequences
# current form gets data from each sequence in a fasta alignment, could provide an id_list using Entrez db search etc.
import Bio
import os
import re
import numpy as npP
from Bio import SeqIO
from Bio import Entrez

Entrez.email = "oscarjamescharles@gmail.com"
# fasta file of minimum set
fasta_file = 'min-genbank.fas'
regex= '^([^.]+)' # everything before fullstop

start = 152
i = 0 
for record in SeqIO.parse(fasta_file, "fasta"):
    i = i+1
    if(i < start):
        continue
    # get genbank id
    #id = "KP745672"
    id = record.id
    id = re.search(regex, id)
    id = id.group()

    # get loc
    handle = Entrez.efetch(db="nucleotide", id=id, rettype = "gb")
    seq_record = SeqIO.read(handle, "gb")
    keys = seq_record.features[0].qualifiers
    country = keys.get('country')
    if country is None:
        country = "NA"
    else:
        country = ''.join(country)

    # get date
    date = keys.get('collection_date')
    if date is None:
        date = "NA"
    else:
        date = ''.join(date)


    # isolate, if exists then is clinical not passage
    isolate = keys.get('isolation_source')
    if isolate is None:
        isolate = "NA"
    else:
        isolate = ''.join(isolate)

    #strain - lab strains have this value
    strain = keys.get('strain')
    if strain is None:
        strain = "NA"
    else:
        strain = ''.join(strain)

    print(id+", "+country+", "+date+", "+isolate+", "+strain)
