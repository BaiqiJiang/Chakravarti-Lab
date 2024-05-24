[[[Introduction]]]

We proudly introduce the calculate_concor tool! Developed in Python and based on cyvcf2, this tool is designed to compare two specified samples in two VCF files, calculate Genotype Concordance and Non-Reference Genotype Concordance, and provide corresponding statistical data.

Compared to GATK 4.2.3.0's Picard Genotype Concordance, our advantages include:

1. Users can specify INFO intervals (e.g., AF or RAF) to calculate concordance for each interval.
2. Thanks to the Cython framework, it remains faster than scripts based on GATK even on less powerful computing resources. In tests with three pairs of similar-sized datasets, tools combining bcftools with GATK took an average of 5 minutes per task on a single-core 64RAM CPU in an HPC environment. In contrast, our tool, running on a single-core 10RAM CPU in a local environment, completed the same tasks in an average of 9.03 seconds.
3. Built on the Python architecture for ease of maintenance and understanding.
4. Current tools from GATK have unresolved bugs that inaccurately reflect Genotype Concordance and TN/TP/FN/FP data. We have fixed this issue.
5. When using bcftools to divide intervals, if a variant's frequency is exactly on the boundary of an interval, it could be counted twice. We solved this issue by assigning a unique identifier to each variant.

[[[Dependencies]]]

1. Python [3.8.19]
2. htslib [1.19.1]
3. cyvcf2 [0.30.28]
4. matplotlib [3.7.5]

[[[Computing Environment Required]]]

Currently, this tool requires you to have a relative large RAM, depends on your call VCF's file size. For example, if your call VCF's file size is 32GB, you may want to use a 40GB RAM to complete your task. We strongly recommend you to use a High-Performance Computing Core (HPCC) instead of a personal computer. 

If you don't have an applicable computing environment that meets the requirements, we are currently working on a major update for significantly saving users' RAM. Once completed, for most of the common VCF files, a 16GB RAM might be enough. This update is expected to complete by August, 2024.

[[[Usage Example]]]

python calculate_concor.py --truth_vcf HG00099ori.chr21.vcf.gz --call_vcf HG00099cp.chr21.vcf.gz --sample_truth HG00099ori --sample_call HG00099cp --info_field AF --intervals "0-0.01" "0.01-0.02" "0.02-0.05" "0.05-0.1" "0.1-0.2" "0.2-0.5" "0.5-1.0" --output_path output.tsv

[[[Required Parameters]]]

For brief versions of this section, refer to the -h (help) function.

--truth_vcf: Path to the VCF file containing the truth set or as defined by the user, in .vcf.gz format with a .tbi index file in the same directory.
--call_vcf: Path to the VCF file for comparison. Note the distinction from the truth set. This file must also be in .vcf.gz format with a .tbi index.
--sample_truth: The sample name within the --truth_vcf.
--sample_call: The sample name within the --call_vcf.
--info_field: The INFO field for interval filtering, such as AF or RAF. This parameter must exist under INFO in both VCF files.
--intervals: The filtering intervals, automatically filtered according to --info_field.
--output_path: The path to the output file, it must be set as a TSV file.

[[[Optional Parameters]]]

For brief versions of this section, refer to the -h (help) function.

--stat_counts: Default is False. When True, includes counts of TP, TN, FP, FN in the output file.
--generate_figs: Default is False. When True, requires --figs_path to be specified for generating concordance scatter plots based on variant types and info intervals.
--figs_path: Default is empty. Must be filled when --generate_figs is True. The path where charts will be stored.
--get_info_from_call: Default is False. When True, get INFOs from call VCF instead of truth VCF. We strongly recommend you to ignore this parameter, unless you DO NOT have INFOs for your truth VCF.

[[[Output File Explanations]]]

1. Variant_Type: Corresponds to the variant type.
2. INFO_Interval: Corresponds to the INFO interval.
3. G_Concordance: Genotype concordance, calculated as TP + TN / TP + TN + FP + FN. If the denominator is 0, the output is 0.
4. Non_Ref_G_Concordance: Non-reference genotype concordance, calculated as TP / TP + FP + FN. If the denominator is 0, the output is 0.
5. TP_Counts: The number of TP. At the same locus, if both the truth and call samples are ALT, it counts as TP. Phase is not considered; e.g., 0|1 and 1|0 are also viewed as TP.
6. TN_Counts: The number of TN. At the same locus, if both the truth and call samples are REF, it counts as TN.
7. FP_Counts: The number of FP. At the same locus, if the truth sample is REF while the call sample is ALT, it counts as FP.
8. FN_Counts: The number of FN. At the same locus, if the truth sample is ALT while the call sample is REF, it counts as FN.
9. Not_Matched: The number of loci detected in the truth set but not in the call set. Here, loci refer to the matching of CHROM, POS, REF, and ALT parameters. For example, for a specific data point, if the CHROM, POS, REF, and ALT are exactly the same for two samples, they belong to the same locus.

[[[Acknowledgements]]]

I sincerely appreciate the guidance and contributions of Dr. Or Yaacov in the conceptualization, development, and debugging of this tool. Additionally, I am grateful for the opportunity to work under Dr. Aravinda Chakravarti and his lab. Lastly, I thank the High Performance Computing Core (HPCC) facility at NYU Langone Health and NYU Grossman School of Medicine for their support with computational resources.

[[[References]]]

1. Pedersen BS, Quinlan AR. cyvcf2: fast, flexible variant analysis with Python. Bioinformatics. 2017;33(12):1867-1869. doi:10.1093/bioinformatics/btx057

[[[Contact Information]]]

If you have questions about using the tool, request for functionality updates, bug fixes, or have any other inquiries, please contact:

Baiqi Jiang
Email: baiqijiang0108@gmail.com (preferred)
Cell/Text: 8573487988
