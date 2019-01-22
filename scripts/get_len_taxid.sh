#!/bin/bash
# Number of attempts to request data from e-utils
att=10
if [ ! -z "${NCBI_API_KEY}" ]
then
	api_key="&api_key=${NCBI_API_KEY}"
fi

retrieve_nucleotide_fasta_xml()
{
	echo "$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=${1}${2}")"
}

for ACC in $1
do
	xml_out=""
	taxid=""
	# Try to retrieve information
	for i in $(seq 1 ${att});
	do
		xml_out="$(retrieve_nucleotide_fasta_xml "${ACC}" "${api_key}")"
		taxid="$(echo "$xml_out" | grep -m 1 -oP '(?<=Name="TaxId" Type="Integer">)[^<]+')"
		# If taxid was found, break
		if [[ ! -z "${taxid}" ]]; then break; fi;
	done
	# If taxid was not found, add to the error list and continue
	if [[ -z "${taxid}" ]]; 
	then 
		error="${error} ${ACC}"
		continue
	fi
	# Extract sequence length 
	len="$(echo "$xml_out" | grep -m 1 -oP '(?<=Name="Length" Type="Integer">)[^<]+')"
	
	# Print output to STDOUT
	echo ${ACC}$'\t'${len}$'\t'${taxid}
done

# Print errors to STDERR
if [ ! -z "${error}" ]
then
	(>&2 echo "Failed to retrieve information: "${error})
	exit 1
fi
exit 0
