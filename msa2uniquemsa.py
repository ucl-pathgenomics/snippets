from Bio import SeqIO
from collections import defaultdict

dedup_records = defaultdict(list)
for record in SeqIO.parse("u38_AA_msa2_from.fasta", "fasta"):
    # Use the sequence as the key and then have a list of id's as the value
    dedup_records[str(record.seq)].append(record.id)

with open("Output.fasta", 'w') as output:
    for seq, ids in dedup_records.items():
        # Join the ids and write them out as the fasta
        output.write(">{}\n".format('|'.join(ids)))
        output.write(seq + "\n")
