import cyvcf2
import argparse
import matplotlib.pyplot as plt
import os

### Current Version: [1.2.0]
### Last Update: 04/04/2024
### Important
### For future developers in our lab, please refer to the readme.txt to understand these codes quickly without reading each line
### Req libs versions info:
### cyvcf2 [0.30.28]
### htslib [1.19.1]
### python [3.8.19]
### matplotlib [3.7.5]
### Copyright: Baiqi Jiang & Chakravarti Lab @2024

def plot_concordance(figs_path, intervals, concordance_data, variant_type, concordance_type, info_field):
    """A function that generate stat figures based on the calculated concordances, variant types and selected intervals."""

    # Check path if not exist then create
    if not os.path.exists(figs_path):
        os.makedirs(figs_path)
    
    # Plot
    plt.figure()
    plt.plot(intervals, concordance_data, marker = 'o', linestyle = '-')
    plt.title(f"{variant_type} {concordance_type} Concordance by {info_field} Intervals")
    plt.xlabel(f"{info_field} Intervals")
    plt.ylabel(f"{concordance_type} Concordance (%)")
    plt.grid(True)
    plt.xticks(rotation = 45)
    plt.tight_layout()
    
    # Save files to path
    plt.savefig(f'{figs_path}/{variant_type}_{concordance_type}_Concordance_{info_field}.tiff', format = 'tiff', dpi = 600)
    plt.close()

def calculate_concordance(truth_vcf_path, call_vcf_path, sample_truth, sample_call, info_field, intervals, output_path, stat_counts, generate_figs, figs_path):
    """A function, or, method that calculates the genotype (GT) concordance and non-reference genotype concordance between two VCF files."""

    # Set counters
    counters = {interval: {'SNP': {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0, 'not_matched': 0},
                           'INDEL': {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0, 'not_matched': 0}}
                for interval in intervals}

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
        headers = ["Variant_Type", f"{info_field}_Interval", "G_Concordance", "Non_Ref_G_Concordance"]

        if stat_counts:
            headers.extend(["TP_Counts", "TN_Counts", "FP_Counts", "FN_Counts"])

        headers.append("Not_Matched")
        f_out.write("\t".join(headers) + "\n")

        for interval, data in counters.items():
            for var_type, counts in data.items():
                total = sum(counts.values()) - counts['not_matched']
                g_concordance = (counts['tp'] + counts['tn']) / total if total > 0 else 0
                non_ref_g_concordance = counts['tp'] / (counts['tp'] + counts['fp'] + counts['fn']) if (counts['tp'] + counts['fp'] + counts['fn']) > 0 else 0
                line = [var_type, interval, f"{g_concordance:.10f}", f"{non_ref_g_concordance:.10f}"]

                if stat_counts:
                    line.extend([str(counts['tp']), str(counts['tn']), str(counts['fp']), str(counts['fn'])])
                    
                line.append(str(counts['not_matched']))
                f_out.write("\t".join(line) + "\n")

    if generate_figs:
        for var_type in ['SNP', 'INDEL']:
            g_concordance_data = []
            non_ref_g_concordance_data = []

            for interval in intervals:
                data = counters[interval][var_type]
                total = sum(data.values()) - data['not_matched']
                g_concordance = (data['tp'] + data['tn']) / total if total > 0 else 0
                non_ref_g_concordance = data['tp'] / (data['tp'] + data['fp'] + data['fn']) if (data['tp'] + data['fp'] + data['fn']) > 0 else 0
                g_concordance_data.append(g_concordance * 100)
                non_ref_g_concordance_data.append(non_ref_g_concordance * 100)
        
            plot_concordance(figs_path, intervals, g_concordance_data, var_type, 'G', info_field)
            plot_concordance(figs_path, intervals, non_ref_g_concordance_data, var_type, 'Non_Ref_G', info_field)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Calculate genotype concordance, non reference genotype concordance and unmatched variant counts based on selected info intervals between two VCF files.")
    parser.add_argument("--truth_vcf", required = True, help = "Path to the truth VCF file.")
    parser.add_argument("--call_vcf", required = True, help = "Path to the call VCF file.")
    parser.add_argument("--sample_truth", required = True, help = "The sample name in the truth VCF.")
    parser.add_argument("--sample_call", required = True, help = "The sample name in the call VCF.")
    parser.add_argument("--info_field", required = True, help = "The INFO field to use for filtering intervals (e.g., AF or RAF).")
    parser.add_argument("--intervals", required = True, nargs = '+', help = "List of --info_field intervals.")
    parser.add_argument("--output_path", required = True, help = "Path to the output TSV file.")
    parser.add_argument("--stat_counts", type = bool, default = False, help = "Default is False, when True, include detailed TP, TN, FP, FN counts in output.")
    parser.add_argument("--generate_figs", type = bool, default = False, help = "Default is False, when True, generate concordance figures based on concordance data, variant types and selected info intervals.")
    parser.add_argument("--figs_path", type = str, help = "Path to save the generated figures.")
    
    args = parser.parse_args()

    if args.generate_figs and not args.figs_path:
        parser.error("--figs_path is required when --generate_figs is set to True.")

    calculate_concordance(args.truth_vcf, args.call_vcf, args.sample_truth, args.sample_call, args.info_field, args.intervals, args.output_path, args.stat_counts, args.generate_figs, args.figs_path)