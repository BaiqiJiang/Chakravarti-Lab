import cyvcf2
import argparse

### Current Version: [1.0.0]
### Last Update: 04/03/2024
### Important
### For future developers in our lab, please refer to the readme.txt to understand these codes quickly without reading each line
### Req libs versions info:
### cyvcf2 [0.30.28]
### htslib [1.19.1]
### python [3.8.19]
### matplotlib [3.7.5]
### Copyright: Baiqi Jiang & Chakravarti Lab @2024

def calculate_concordance(truth_vcf_path, call_vcf_path, sample_truth, sample_call, info_field, intervals, output_path):
    """A function, or, method that calculates the genotype (GT) concordance and non-reference genotype concordance between two VCF files."""

    # Set counters
    counters = {interval: {'SNP': {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0, 'not_matched': 0}, 'INDEL': {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0, 'not_matched': 0}} for interval in intervals}

    # Load truth_vcf and call_vcf files
    truth_vcf = cyvcf2.VCF(truth_vcf_path)
    call_vcf = cyvcf2.VCF(call_vcf_path)

    # Create a dictionary to store all call variants into a specific format
    call_variants = {}

    for var in call_vcf:
        variant_key = f'{var.CHROM}:{var.POS}-{var.REF}-{",".join(var.ALT)}'
        call_variants[variant_key] = var

    # Locate samples
    truth_sample_index = truth_vcf.samples.index(sample_truth)
    call_sample_index = call_vcf.samples.index(sample_call)

    # Go through each line in truth vcf
    # Locate interval
    # Find the location that matches with the call variants dict
    # Define variant type
    # Calculations
    for truth_var in truth_vcf:
        variant_freq = truth_var.INFO.get(info_field)
        interval_key = None

        for interval in intervals:
            start, end = map(float, interval.split('-'))

            if start <= variant_freq <= end:
                interval_key = interval
                break

        if not interval_key:
            continue

        variant_key = f'{truth_var.CHROM}:{truth_var.POS}-{truth_var.REF}-{",".join(truth_var.ALT)}'
        variant_type = 'SNP' if truth_var.is_snp else 'INDEL'

        if variant_key in call_variants:
            matched_var = call_variants[variant_key]
            truth_gt = truth_var.genotypes[truth_sample_index][:2]
            call_gt = matched_var.genotypes[call_sample_index][:2]

            truth_simplified = [1 if gt > 0 else 0 for gt in truth_gt]
            call_simplified = [1 if gt > 0 else 0 for gt in call_gt]

            if truth_simplified == call_simplified:
                if truth_simplified == [0, 0]:
                    counters[interval_key][variant_type]['tn'] += 1
                else:
                    counters[interval_key][variant_type]['tp'] += 1
            else:
                if truth_simplified == [0, 0]:
                    counters[interval_key][variant_type]['fp'] += 1
                elif call_simplified == [0, 0]:
                    counters[interval_key][variant_type]['fn'] += 1
                elif sum(truth_simplified) == 1 and sum(call_simplified) == 1:
                    counters[interval_key][variant_type]['tp'] += 1
        else:
            counters[interval_key][variant_type]['not_matched'] += 1

    # Write the results into a file
    with open(output_path, 'w') as f_out:
        f_out.write(f"Variant_Type\t{info_field}_Interval\tG_Concordance\tNon_Ref_G_Concordance\tNot_Matched\n")

        for interval, data in counters.items():
            for var_type, counts in data.items():
                total = sum(counts.values()) - counts['not_matched']
                g_concordance = (counts['tp'] + counts['tn']) / total if total > 0 else 0
                non_ref_g_concordance = counts['tp'] / (counts['tp'] + counts['fp'] + counts['fn']) if (counts['tp'] + counts['fp'] + counts['fn']) > 0 else 0
                f_out.write(f"{var_type}\t{interval}\t{g_concordance:.10f}\t{non_ref_g_concordance:.10f}\t{counts['not_matched']}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Calculate genotype concordance, non reference genotype concordance and unmatched variant counts based on selected info intervals between two VCF files.")
    parser.add_argument("--truth_vcf", required = True, help = "Path to the truth VCF file.")
    parser.add_argument("--call_vcf", required = True, help = "Path to the call VCF file.")
    parser.add_argument("--sample_truth", required = True, help = "The sample name in the truth VCF.")
    parser.add_argument("--sample_call", required = True, help = "The sample name in the call VCF.")
    parser.add_argument("--info_field", required = True, help = "The INFO field to use for filtering intervals (e.g., AF or RAF).")
    parser.add_argument("--intervals", required = True, nargs = '+', help = "List of --info_field intervals.")
    parser.add_argument("--output_path", required = True, help = "Path to the output TSV file.")

    args = parser.parse_args()

    calculate_concordance(args.truth_vcf, args.call_vcf, args.sample_truth, args.sample_call, args.info_field, args.intervals, args.output_path)
