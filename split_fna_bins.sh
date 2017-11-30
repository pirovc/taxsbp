#!/bin/bash
# Usage: file.fasta.gz taxsbpbins.txt outputfolder/

# No prints/no gz
#awk -v taxsbpfile="../arc_bac_all_ar/2017-09-25_18-28-39_b59.txt" '
#BEGIN{while (getline < taxsbpfile){clust[$1]=$4;}; close(taxsbpfile);}
#/^>/ {filen=clust[substr($1,2)]}
#{print $0 > filen".fna"}' ../arc_bac_all_20170925.fna

zcat $1 | 
awk -v taxsbpfile=$2 -v outputfolder=$3 '
BEGIN{
	print "Parsing "taxsbpfile"...";
	while (getline < taxsbpfile){clust[$1]=$4;}; 
	close(taxsbpfile);
	l=length(clust)
	print "Done! "l" entries";
	print "Writing files ... ";
}
/^>/ {
filen=clust[substr($1,2)];
cont+=1;
if(cont % int(l/1000) == 0) {print cont" sequences out of "l" (" (cont/l)*100 "%)"};
}
{print $0 | " gzip >> "outputfolder"/"filen".fna.gz"}
END{print "Done!"}' 

