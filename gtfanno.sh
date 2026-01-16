#!/bin/bash


info(){
	printf "\e[44m%s\e[0m\n" "$*"
}
ok(){
	printf "\e[42m%s\e[0m\n" "$*"
}
warn(){
	printf "\e[33m%s\e[0m\n" "$*"
}

usage() {
	cat <<\END_HELP
   _____ _______ ______                      
  / ____|__   __|  ____|                     
 | |  __   | |  | |__ __ _ _ __  _ __   ___  
 | | |_ |  | |  |  __/ _` | '_ \| '_ \ / _ \ 
 | |__| |  | |  | | | (_| | | | | | | | (_) |
  \_____|  |_|  |_|  \__,_|_| |_|_| |_|\___/  
                                             
This script is for generating TSS, exon, intron, intergenic region annotations from GTF file.
The input GTF file requires chromosome name with "chr" prefix.

Usage:
	bash gtfanno.sh -f <GTF file>			# Output file to current working directory
	bash gtfanno.sh -f <GTF file> -o <output directory> -k   # Includes scaffolds
	
Parameters:
	-f	Required. Path to GTF file. Either gzipped or plain file is accepted.
	-p	Optional. Prefix to output files. Default the same as the prefix of GTF file name.
	-o	Optional. Output directory. Default the current working directory.
	-k	Optional. Also include scaffolds.
	-u	Optional. Upstream distance from TSS (bp). Default: 2000.
	-d	Optional. Downstgream distance from TSS (bp). Default: 500.
	-s	Optional. The local chromosome size file for calculating intergenic area. If left empty, it will be automatically downloaded from https://github.com/igvteam/igv/tree/maaster/genomes/sizes
	-h	Print this help message.

END_HELP
	exit "$1";
}

check_cmd() {
	local uninstalled=()
	for c in bedtools; do
		if ! command -v $c > /dev/null; then
			uninstalled+=("$c")
		fi
	done
	if [[ ${#uninstalled[@]} -gt 0 ]]; then
		warn "Required command[s] not found:" "${uninstalled[@]}"
		exit 1
	fi
}
check_cmd

include_scaffold=false
outdir=$(pwd)
tss_radius=300

while getopts f:p:o:ku:d:s:h opt; do
case ${opt} in
	f) gtf_file=${OPTARG};;
	p) prefix=${OPTARG};;
	o) outdir=${OPTARG};;
	k) include_scaffold=true;;
	u) up_dist=${OPTARG};;
	d) down_dist=${OPTARG};;
	s) size_file=${OPTARG};;
	h) usage 0;;
	*) usage 1;;
esac
done

if [[ -z $gtf_file ]];then
	usage 1
fi

regiondir=$outdir/region
extradir=$outdir/extra
tmpdir=$outdir/.tmp

mkdir -p "$regiondir"
mkdir -p "$extradir"
mkdir -p "$tmpdir"

if [[ -z $prefix ]];then
	prefix=$(basename "$gtf_file" | awk -F '.gtf' '{print $1}')
fi

load_gtf() {
	if [[ $(file -b "$gtf_file") == *gzip* ]];then
		if [ "$(uname)" = "Linux" ]; then
			zcat "$gtf_file"
		elif [ "$(uname)" = "Darwin" ]; then
			gzcat "$gtf_file"
		fi
	else
		cat "$gtf_file"
	fi
}

if load_gtf | head | grep -E "GRCh38|hg38" > /dev/null; then
	genome=hg38
elif load_gtf | head | grep -E "GRCh37|hg19" > /dev/null; then
	genome=hg19
elif load_gtf | head | grep -E "GRCm39|mm39" > /dev/null; then
	genome=mm39
elif load_gtf | head | grep -E "GRCh38|mm10" > /dev/null; then
	genome=mm10
fi

genome_size_url="https://raw.githubusercontent.com/igvteam/igv/refs/heads/main/genomes/sizes/$genome.chrom.sizes"

chr_file=$outdir/$prefix.chr.bed
tss_file=$regiondir/$prefix.tss_u${up_dist}_d${down_dist}.bed
exon_file=$regiondir/$prefix.exon_no_tss.bed
intron_file=$regiondir/$prefix.intron.bed
intergenic_file=$regiondir/$prefix.intergenic.bed
three_utr_file=$extradir/$prefix.3utr.bed
three_utr_notss_file=$extradir/$prefix.3utr_no_tss.bed
five_utr_file=$extradir/$prefix.5utr.bed
five_utr_notss_file=$extradir/$prefix.5utr_no_tss.bed
cds_file=$extradir/$prefix.cds.bed
cds_notss_file=$extradir/$prefix.cds_no_tss.bed

#####################################################
#													#
#			1) GTF to BED conversion				#
#													#
#####################################################

info "Start converting GTF to BED: $(warn $chr_file)"

# row example
# 1	2	3	4	5	6	7	8	9
# chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.1";
gtf_to_bed() {
	 awk -F '\t' -v OFS='\t' '{
		split($9, attr, ";")
		for (i = 1; i <= length(attr); i++) {
			if (attr[i] ~ /gene_name/) {
				split(attr[i], gene_name, " ")
				gsub("\"", "", gene_name[2])
				print $1, $4-1, $5, $1":"$4-1"-"$5":"$3":"gene_name[2]":"$7, "0", $7
				break
			}
		}
	}'
}

filter_scaffold() {
	if $include_scaffold; then
		grep ""
	else
		grep -E '^(chr[0-9]+|[0-9]+|chr[XYM])'
	fi
}

sort_bed() {
	awk -F '\t' -v OFS='\t' '{
		if ($1 == "chrM")
			print $0, "chrZ"
		else
			print $0, $1
	}' | sort -S1G -k7,7V -k2,2n -k3,3n | awk -F '\t' -v OFS='\t' '{
		$NF=""
		sub(/\t$/, "")
		print $0
	}'
}

load_gtf "$gtf_file" | gtf_to_bed | filter_scaffold | sort_bed | uniq > $chr_file
ok "Job finished"

#####################################################
#													#
#				2) TSS annotation					#
#													#
#####################################################

info "Start generate TSS $tss_radius annotation: $(warn $tss_file)"

get_tss_from_chr() {
	awk -F '\t' -v u="$up_dist" -v d="$down_dist" -v OFS='\t' '{
		if ($6 == "+") {
			s = $2 - u
			e = $2 + d
		} else {
			s = $3 - d
			e = $3 + u
		}
		if (s < 0) s = 0;
		print $1, s, e, $4, $5, $6
	}'
}

bedtools_merge_with_genename() {
    sort_bed | bedtools merge -s -c 4,6 -o distinct,distinct -i stdin | \
	awk -F '\t' -v OFS='\t' -v region="$1" '{
		split($4, name, ":")
        print $1, $2, $3, $1":"$2"-"$3":"region":"name[4]":"$5, "0", $5
    }' | sort_bed
}

bedtools_merge_simple() {
    sort_bed | bedtools merge -s -c 6 -o distinct -i stdin | \
	awk -F '\t' -v OFS='\t' -v region="$1" '{
        print $1, $2, $3, $1":"$2"-"$3":"region":"$4, "0", $4
    }' | sort_bed
}


grep transcript "$chr_file" | get_tss_from_chr | bedtools_merge_with_genename tss > $tss_file

ok "Job finished."

#####################################################
#													#
#				3) exon annotation					#
#													#
#####################################################

info "Start generate exon (excluding TSS) annotation: $(warn $exon_file)"

exon_tmp=$tmpdir/exon.tmp
grep "exon" $chr_file > $exon_tmp
bedtools subtract -s -a $exon_tmp -b $tss_file | bedtools_merge_with_genename exon > $exon_file

ok "Job finished."

#####################################################
#													#
#				4) intron annotation				#
#													#
#####################################################

info "Start generate intron (excluding TSS) annotation: $(warn $intron_file)"

exon_tss_tmp=$tmpdir/exon_tss.tmp
gene_tmp=$tmpdir/gene.tmp

grep "gene" $chr_file > $gene_tmp
cat $exon_tmp $tss_file | bedtools_merge_simple exontss > $exon_tss_tmp

bedtools subtract -s -a $gene_tmp -b $exon_tss_tmp | bedtools_merge_with_genename intron > $intron_file

ok "Job finished."

#####################################################
#													#
#			5) intergenic annotation				#
#													#
#####################################################

info "Start generate intergenic annotation: $(warn $intergenic_file)"

gene_tss_fwd_tmp=$tmpdir/gene_tss_fwd.tmp
gene_tss_rev_tmp=$tmpdir/gene_tss_rev.tmp
intergenic_fwd_tmp=$tmpdir/intergenic_fwd.tmp
intergenic_rev_tmp=$tmpdir/intergenic_rev.tmp
genome_size_tmp=$tmpdir/genome_size.tmp

if [[ -z $size_file ]] || [[ ! -f $size_file ]]; then
	warn "Local genome size file not found. Start downloading from $genome_size_url"
	curl $genome_size_url > $genome_size_tmp
	ok "$genome size file is downloaded to $genome_size_tmp"
else
	cat $size_file > genome_size_tmp
fi


cat $gene_tmp $tss_file | awk -F '\t' -v OFS='\t' '{if ($6 == "+") {print $0}}' | \
	bedtools_merge_simple > $gene_tss_fwd_tmp
cat $gene_tmp $tss_file | awk -F '\t' -v OFS='\t' '{if ($6 == "-") {print $0}}' | \
	bedtools_merge_simple > $gene_tss_rev_tmp

bedtools complement -i $gene_tss_fwd_tmp -g $genome_size_tmp | \
	awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1":"$2"-"$3":+", "0", "+"}' > $intergenic_fwd_tmp
bedtools complement -i $gene_tss_rev_tmp -g $genome_size_tmp | \
	awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1":"$2"-"$3":-", "0", "-"}' > $intergenic_rev_tmp

cat $intergenic_fwd_tmp $intergenic_rev_tmp | bedtools_merge_simple intergenic > $intergenic_file

ok "Job finished."

#####################################################
#													#
#				6) UTRs						#
#													#
#####################################################

info "Start generate 3'UTR and 5'UTR annotation: $(warn $three_utr_file $three_utr_notss_file $five_utr_file $five_utr_notss_file)"
utr_tmp=$tmpdir/utr.tmp

load_gtf $gtf_file | grep chr | awk '$3 == "UTR" || $3 == "CDS"' | python3 -c '''
import sys

gtf = dict(CDS={}, UTR={})

for line in sys.stdin:
    content = line.strip().split("\t")
    feature_type = content[2]
    transcript_id = ""
    gene_name = ""
    for field in content[8].strip(";").split(";"):
        k, v = field.strip().split(" ")
        if k == "transcript_id":
            transcript_id = v.replace("\"", "")
        elif k == "gene_name":
            gene_name = v.replace("\"", "")
            break

    if gtf[feature_type].get(transcript_id) is None:
        gtf[feature_type][transcript_id] = []
    gtf[feature_type][transcript_id].append([content[i] for i in [0, 2, 3, 4, 6]] + [gene_name])

for k in gtf["UTR"].keys():
    for u in gtf["UTR"][k]:
        if k not in gtf["CDS"]:
            continue
        if gtf["CDS"][k][0][4] == "+":
            utr_type = "three_prime_utr" if int(u[2]) > max([int(c[3]) for c in gtf["CDS"][k]]) else "five_prime_utr"
        else:
            utr_type = "three_prime_utr" if int(u[3]) < min([int(c[2]) for c in gtf["CDS"][k]]) else "five_prime_utr"

        print("{}\t{}\t{}\t{}\t{}\t{}".format(
            u[0],
            int(u[2]) - 1,
            u[3],
            "{}:{}-{}:{}:{}:{}:{}".format(u[0], str(int(u[2]) - 1), u[3], utr_type, k, gtf["UTR"][k][0][5], gtf["UTR"][k][0][4]),
            0,
            gtf["UTR"][k][0][4]
        ))
''' > $utr_tmp

grep five_prime_utr $utr_tmp > $five_utr_file
bedtools subtract -s -a $five_utr_file -b $tss_file | sort_bed |\
	bedtools merge -s -i stdin -c 4,6 -o distinct,distinct | awk -F '\t' -v OFS='\t' '{
		split($4, name, ":"); 
		print $1, $2, $3, $1":"$2"-"$3":5utr:"name[5]":"$5, "0", $5
	}' > $five_utr_notss_file

grep three_prime_utr $utr_tmp > $three_utr_file
bedtools subtract -s -a $three_utr_file -b $tss_file | sort_bed |\
	bedtools merge -s -i stdin -c 4,6 -o distinct,distinct | awk -F '\t' -v OFS='\t' '{
		split($4, name, ":"); 
		print $1, $2, $3, $1":"$2"-"$3":3utr:"name[5]":"$5, "0", $5
	}' > $three_utr_notss_file


ok "Job finished."


#####################################################
#													#
#					7) CDS							#
#													#
#####################################################

info "Start generate CDS annotation: $(warn $cds_file $cds_notss_file)"

load_gtf $gtf_file | filter_scaffold | awk '$3 == "CDS"' | gtf_to_bed | sort_bed | uniq > $cds_file
bedtools subtract -s -a $cds_file  -b $tss_file | sort_bed |\
	bedtools merge -s -i stdin -c 4,6 -o distinct,distinct | awk -F '\t' -v OFS='\t' '{
		split($4, name, ":"); 
		print $1, $2, $3, $1":"$2"-"$3":cds:"name[4]":"$5, "0", $5
	}' > $cds_notss_file

ok "Job finished."

#####################################################
#                                                   #
#                   8) Merge regions                #
#                                                   #
#####################################################

merged_file=$outdir/$prefix.region.bed
info "Start merging all region files: $(warn $merged_file)"

cat $tss_file $exon_file $intron_file $intergenic_file | sort_bed > $merged_file

ok "Job finished."

#####################################################
#													#
#					9) clean						#
#													#
#####################################################

rm -rf "$tmpdir"
