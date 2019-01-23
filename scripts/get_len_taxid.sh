#!/bin/bash
# Number of attempts to request data from e-utils
att=3
if [ ! -z "${NCBI_API_KEY}" ]
then
	api_key="&api_key=${NCBI_API_KEY}"
fi

retrieve_summary_xml()
{
	echo "$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=${1}${2}")"
}

retrieve_nucleotide_fasta_xml()
{
	echo "$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=xml&id=${1}${2}")"
}

for ACC in $1
do
	xml_out=""
	taxid=""
	len=""
	# Try to retrieve information
	for i in $(seq 1 ${att});
	do
		# First try to get from summary, lighter resource
		xml_out="$(retrieve_summary_xml "${ACC}" "${api_key}")"
		taxid="$(echo "${xml_out}" | grep -m 1 -oP '(?<=Name="TaxId" Type="Integer">)[^<]+')"
		# If taxid was found, get len and break the attempts
		if [[ ! -z "${taxid}" ]]; then 
			# Extract sequence length 
			len="$(echo "${xml_out}" | grep -m 1 -oP '(?<=Name="Length" Type="Integer">)[^<]+')"
			break
		else 
			# Look for error msg on xml_out
			error_xml_out="$(echo "${xml_out}" | grep -m 1 -oP '(?<=<ERROR>)[^<]+')"
			if [[ ! -z "${error_xml_out}" ]]; then # Found error msg
				# Try different resource
				xml_out="$(retrieve_nucleotide_fasta_xml "${ACC}" "${api_key}")"
				taxid="$(echo "${xml_out}" | grep -m 1 -oP '(?<=<TSeq_taxid>)[^<]+')"
				if [[ ! -z "${taxid}" ]]; then 
					# Extract sequence length 
					len="$(echo "${xml_out}" | grep -m 1 -oP '(?<=<TSeq_length>)[^<]+')"
					break
				else
					# more attempts
					continue 
				fi
			else
				# more attempts
				continue 
			fi
		fi
	done

	# If taxid or len was not found, add to the error list and continue
	if [[ -z "${taxid}" ]]; 
	then 
		error="${error} ${ACC}"
		continue
	fi

	# Print output to STDOUT
	echo ${ACC}$'\t'${len}$'\t'${taxid}
done

# Print errors to STDERR
if [ ! -z "${error}" ]
then
	(>&2 echo "Failed to get taxid and sequence length: "${error})
	exit 1
fi
exit 0
