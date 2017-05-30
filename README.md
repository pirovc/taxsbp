# TaxSBP

Vitor C. Piro (vitorpiro@gmail.com)


Implementation of the approximation algorithm for the hierarchically structured bin packing problem [1] based on the NCBI Taxonomy database [2].

Dependencies:
-------------
 - binpacking=1.3 (https://pypi.python.org/pypi/binpacking)

Input: 
------
 * nodes.dmp and a file with sequence information (identifier, length and taxonomic assignment)

nodes.dmp:

	wget -qO- ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz | tar xfz - nodes.dmp

sequence information (from NCBI refseq/genbank):

    # Bacteria - RefSeq - Complete Genomes (Download assembly reports adding taxonomic assignment at the end) 
    wget -qO- ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt | tail -n+3 |
    awk -F "\t" '$12=="Complete Genome" && $11=="latest"{url_count=split($20,url,"/"); print $6"\t"$20"/"url[url_count] "_assembly_report.txt"}' |
    parallel -j 24 --colsep "\t" 'wget -qO- {2} | grep "^[^#]" | tr -d "\r" | sed -e s/$/\\t{1}\\t`basename {2}`/' > refseq_bac_cg_ar.txt 2> refseq_bac_cg_ar.err

The assembly accession starts with a three letter prefix, GCA for GenBank assemblies and GCF for RefSeq assemblies. This is followed by an underscore and 9 digits. A version is then added to the accession. For example, the assembly accession for the GenBank version of the current public human reference assembly ( GRCh37.p2 ) is GCA_000001405.3 [https://www.ncbi.nlm.nih.gov/assembly/model/]

from new:

assembly_summary_file="refseq_bacteria_cg/2017-05-30_12-14-43_assembly_summary.txt"
assembly_reports_folder="refseq_bacteria_cg/files/"
cut -f 6,20 ${assembly_summary_file} | awk -F "\t" -v assembly_reports_folder="${assembly_reports_folder}" '{url_count=split($2,url,"/"); print $1 "\t" assembly_reports_folder url[url_count] "_assembly_report.txt"}' | parallel --colsep "\t" 'grep "^[^#]" {2} | tr -d "\r" | cut -f 7,9 | sed s/$/\\t{1}/' |

from _added.txt

added_accession_file="refseq_bacteria_cg/2017-05-30_13-01-54_added.txt"
assembly_summary_file="refseq_bacteria_cg/2017-05-30_13-01-54_assembly_summary.txt"
assembly_reports_folder="refseq_bacteria_cg/files/"
join <(sort -k 1,1 ${added_accession_file}) <(sort -k 1,1 ${assembly_summary_file}) -t$'\t' -o "2.6,2.20"  | awk -F "\t" -v assembly_reports_folder="${assembly_reports_folder}" '{url_count=split($2,url,"/"); print $1 "\t" assembly_reports_folder url[url_count] "_assembly_report.txt"}' | parallel --colsep "\t" 'grep "^[^#]" {2} | tr -d "\r" | cut -f 7,9 | sed s/$/\\t{1}/'

	
Output:
-------
 * A tab-separated file with sequence identifier and bin id

Running:
--------
	python3 TaxSBP.py -a refseq_bac_cg_ar.txt -n nodes.dmp -s 2 -b 50

References:
-----------

[1] Codenotti, B., De Marco, G., Leoncini, M., Montangero, M., & Santini, M. (2004). Approximation algorithms for a hierarchically structured bin packing problem. Information Processing Letters, 89(5), 215–221. http://doi.org/10.1016/j.ipl.2003.12.001

[2] Federhen, S. (2012). The NCBI Taxonomy database. Nucleic Acids Research, 40(D1), D136–D143. http://doi.org/10.1093/nar/gkr1178