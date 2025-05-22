#Find potential fruit fly / yeast orthologs
#   Download FASTA files drosoph-ribosome.fasta.gz and yeast-ribosome.fasta.gz
#   from the course data-directory.
#   Uncompress and format each FASTA file for BLAST
#   Search fruit fly ribosomal proteins against yeast ribosomal proteins
#   For each fruit fly query, output the best yeast protein if it has a
#   significant E-value.
#       What ribosomal protein is most highly conserved between fruit fly
#       and yeast?

from Bio.Blast import NCBIXML

result_handle = open("blast_results.xml")
best_hits = {}                                      #store best yeast matches for each fruit fly query
for blast_result in NCBIXML.parse(result_handle):
    query_id = blast_result.query_id
    best_e_value = float('inf')                     #infinity place holder for tracking smallest E-value
    best_hit = None                                 #stores best yeast match for query
    for alignment in blast_result.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 1e-5:
                if hsp.expect < best_e_value:
                    best_e_value = hsp.expect
                    best_hit = alignment.title
    if best_hit:
        best_hits[query_id] = {                     #if best his is found for query, best yeast match and its E-value are stored
            'best_hit':best_hit,
            'e_value':best_e_value
        }

result_handle.close()

#prints the best yeast match for each fruit fly query
for query_id, hit_info in best_hits.items():
    print('Fruit Fly Query:',query_id)
    print('Best yeast protein match:',hit_info['best_hit'])
    print('E-value:',hit_info['e_value'])
    print()

#finds the overall most conserved ribosomal protein out of all of the queries
most_conserved_query = min(best_hits, key=lambda query: best_hits[query]['e_value'])
most_conserved_hit = best_hits[most_conserved_query]

print('Overall most highly conserved ribosomal protein:')
print('Fruit Fly Query:', most_conserved_query)
print('Best yeast match:', most_conserved_hit['best_hit'])
print('E-value:', most_conserved_hit['e_value'])
