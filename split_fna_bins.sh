#!/bin/bash
# Usage: "file.fasta[.gz] [file2.fasta[.gz] ...]" bins.txt outputfolder/ [gzip]

zcat --force $1 | 
awk -v taxsbpfile=$2 -v outputfolder=$3 -v gzipcat=$(if [ $4 == "gzip" ]; then echo "gzip"; else echo "cat"; fi) -v ext=$(if [ $4 == "gzip" ]; then echo ".gz"; fi) '
BEGIN{
	print "Parsing "taxsbpfile"...";
	while (getline < taxsbpfile){clust[$1]=$4;}; 
	close(taxsbpfile);
	l=length(clust)
	print "Done! "l" entries";
	print "Writing files ... ";
}
/^>/ {
acc=substr($1,2);
filen=(acc in clust)?clust[acc]:"no_cluster";
cont+=1;
if(cont % int((l/1000)+1) == 0) {print cont" sequences out of "l" (" (cont/l)*100 "%)"};
}
{print $0 | gzipcat " >> "outputfolder"/"filen".fna"ext}
END{print "Done!"}' 

