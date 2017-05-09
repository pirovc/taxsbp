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