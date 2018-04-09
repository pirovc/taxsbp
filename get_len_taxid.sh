#!/bin/bash
# Number of attempts to request data from e-utils
att=10

retrieve_nucleotide_fasta_xml()
{
	echo "$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=${1}")"
}

for ACC in $1
do
	# Try to retrieve information
	for i in $(seq 1 ${att});
	do
		xml_out="$(retrieve_nucleotide_fasta_xml "${ACC}")"
		taxid="$(echo "$xml_out" | grep -m 1 -oP '(?<=Name="TaxId" Type="Integer">)[^<]+')"
		# If taxid was found, break
		if [[ ! -z "${taxid}" ]]; then break; fi;
	done
	# If taxid was not found, add to the error list and continue
	if [[ -z "${taxid}" ]]; 
	then 
		nucl_error="${nucl_error}${ACC}"
		continue
	fi
	# Extract sequence length 
	len="$(echo "$xml_out" | grep -m 1 -oP '(?<=Name="Length" Type="Integer">)[^<]+')"
	
	# Print output to STDOUT
	echo ${ACC}$'\t'${len}$'\t'${taxid}
done

# Print errors to STDERR
if [ ! -z "${nucl_error}" ]
then
	(>&2 echo "Problems retrieving nucleotide information for the following entries: "${nucl_error}) 
fi