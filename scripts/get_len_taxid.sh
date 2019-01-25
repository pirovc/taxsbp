#!/bin/bash
# Number of attempts to request data from e-utils
att=10
batch=200

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

get_lines()
{
	echo "$(sed -n "${2},$((${2}+${3}-1))p" ${1})"
}

batch_count=1
acc="$(get_lines ${1} 1 ${batch})"
while [[ ! -z "${acc}" ]];
do

	for i in $(seq 1 ${att});
	do
		# First try to get from summary, lighter resource 
		# replacing \n for , 
		xml_summary="$(retrieve_summary_xml "${acc//$'\n'/,}" "${api_key}")"
		acc_summary="$(echo "${xml_summary}" | grep -oP '(?<=Name="AccessionVersion" Type="String">)[^<]+')"

		# if no accession returned, try again
		if [[ -z "${acc_summary}" ]]; then 
			continue
		else
			len_summary="$(echo "${xml_summary}" | grep -oP '(?<=Name="Length" Type="Integer">)[^<]+')"
			taxid_summary="$(echo "${xml_summary}" | grep -oP '(?<=Name="TaxId" Type="Integer">)[^<]+')"
			out_summary="$(paste <(echo "$acc_summary") <(echo "$len_summary") <(echo "$taxid_summary") --delimiters '\t')"
			echo "${out_summary}"
			break
		fi
	done

	# if there are less output lines than input accessions, get accessions difference (or return is empty)
	if [[ "$(echo "${acc_summary}" | wc -l)" -lt "$(echo "${acc}" | wc -l)" ||  -z "${acc_summary}"  ]];
	then
		acc_diff="$(diff --changed-group-format='%>' --unchanged-group-format='' <(echo "${acc_summary}") <(echo "$acc"))"
	else
		acc_diff=""
	fi

	# If there are accessions left
	if [[ ! -z "${acc_diff}" ]]; then 
		for i in $(seq 1 ${att});
		do
			# try another method
			xml_fetch="$(retrieve_nucleotide_fasta_xml "${acc_diff//$'\n'/,}" "${api_key}")"
			acc_fetch="$(echo "${xml_fetch}" | grep -oP '(?<=<TSeq_accver>)[^<]+')"
			# if no accession returned, try again
			if [[ -z "${acc_fetch}" ]]; then 
				continue
			else
				len_fetch="$(echo "${xml_fetch}" | grep -oP '(?<=<TSeq_length>)[^<]+')"
				taxid_fetch="$(echo "${xml_fetch}" | grep -oP '(?<=<TSeq_taxid>)[^<]+')"
				out_fetch="$(paste <(echo "$acc_fetch") <(echo "$len_fetch") <(echo "$taxid_fetch") --delimiters '\t')"
				echo "${out_fetch}"
				break
			fi
		done
		# if there are less output lines than input accessions, get accessions missing (or return is empty)
		if [[ "$(echo "${acc_fetch}" | wc -l)" -lt "$(echo "${acc_diff}" | wc -l)" ||  -z "${acc_fetch}" ]];
		then
			acc_missing="$(diff --changed-group-format='%>' --unchanged-group-format='' <(echo "${acc_fetch}") <(echo "$acc_diff"))"
			failed="${failed}${acc_missing}\n"
		fi
	fi

	# read new batch
	batch_count=$((batch_count+batch))
	acc="$(get_lines ${1} ${batch_count} ${batch})"
done

# Print errors to STDERR
if [[ ! -z "${failed}" ]];
then
	(>&2 printf "Failed to get taxid and sequence length:\n${failed}")
	exit 1
fi
exit 0

# for ACC in $1
# do
# 	xml_out=""
# 	taxid=""
# 	len=""
# 	# Try to retrieve information
# 	for i in $(seq 1 ${att});
# 	do
# 		# First try to get from summary, lighter resource
# 		xml_out="$(retrieve_summary_xml "${ACC}" "${api_key}")"
# 		taxid="$(echo "${xml_out}" | grep -m 1 -oP '(?<=Name="TaxId" Type="Integer">)[^<]+')"
# 		# If taxid was found, get len and break the attempts
# 		if [[ ! -z "${taxid}" ]]; then 
# 			# Extract sequence length 
# 			len="$(echo "${xml_out}" | grep -m 1 -oP '(?<=Name="Length" Type="Integer">)[^<]+')"
# 			break
# 		else 
# 			# Look for error msg on xml_out
# 			error_xml_out="$(echo "${xml_out}" | grep -m 1 -oP '(?<=<ERROR>)[^<]+')"
# 			if [[ ! -z "${error_xml_out}" ]]; then # Found error msg
# 				# Try different resource
# 				xml_out="$(retrieve_nucleotide_fasta_xml "${ACC}" "${api_key}")"
# 				taxid="$(echo "${xml_out}" | grep -m 1 -oP '(?<=<TSeq_taxid>)[^<]+')"
# 				if [[ ! -z "${taxid}" ]]; then 
# 					# Extract sequence length 
# 					len="$(echo "${xml_out}" | grep -m 1 -oP '(?<=<TSeq_length>)[^<]+')"
# 					break
# 				else
# 					# more attempts
# 					continue 
# 				fi
# 			else
# 				# more attempts
# 				continue 
# 			fi
# 		fi
# 	done

# 	# If taxid or len was not found, add to the error list and continue
# 	if [[ -z "${taxid}" ]]; 
# 	then 
# 		error="${error} ${ACC}"
# 		continue
# 	fi

# 	# Print output to STDOUT
# 	echo ${ACC}$'\t'${len}$'\t'${taxid}
# done

# # Print errors to STDERR
# if [ ! -z "${error}" ]
# then
# 	(>&2 echo "Failed to get taxid and sequence length: "${error})
# 	exit 1
# fi
# exit 0
