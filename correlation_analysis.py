#!/usr/bin/env python3

import csv
import math
import statistics

def pearson_correlation(x, y):
    """Calculate Pearson correlation coefficient between two lists"""
    if len(x) != len(y):
        raise ValueError("Lists must have same length")

    n = len(x)
    if n == 0:
        return 0

    # Calculate means
    mean_x = sum(x) / n
    mean_y = sum(y) / n

    # Calculate correlation
    numerator = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
    sum_sq_x = sum((x[i] - mean_x) ** 2 for i in range(n))
    sum_sq_y = sum((y[i] - mean_y) ** 2 for i in range(n))

    denominator = math.sqrt(sum_sq_x * sum_sq_y)

    if denominator == 0:
        return 0

    return numerator / denominator

def main():
    print("Reading data from 1_np_cutoff_results.csv...")

    # Read CSV data
    data = []
    with open("1_np_cutoff_results.csv", 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            data.append({
                'ENSG': row['ENSG'],
                'gene': row['gene'],
                'Trait': row['Trait'],
                'page_rank': float(row['page.rank']),
                'degree': int(row['degree']),
                'emp_pval': float(row['emp_pval']),
                'qval': float(row['qval'])
            })

    print(f"Data dimensions: {len(data)} rows")

    # Calculate -log10 transformations
    for row in data:
        row['neg_log10_emp_pval'] = -math.log10(row['emp_pval'] + 1e-10)
        row['neg_log10_qval'] = -math.log10(row['qval'] + 1e-10)

    # Group by trait
    traits = {}
    for row in data:
        trait = row['Trait']
        if trait not in traits:
            traits[trait] = []
        traits[trait].append(row)

    print(f"Number of unique traits/diseases: {len(traits)}")

    # Calculate correlations for each disease
    results = []

    print("\nCalculating correlations for each disease...")
    for trait, trait_data in traits.items():
        # Extract data for correlation
        page_ranks = [row['page_rank'] for row in trait_data]
        neg_log10_emp_pvals = [row['neg_log10_emp_pval'] for row in trait_data]
        neg_log10_qvals = [row['neg_log10_qval'] for row in trait_data]

        # Calculate correlations
        corr_emp_pval = pearson_correlation(page_ranks, neg_log10_emp_pvals)
        corr_qval = pearson_correlation(page_ranks, neg_log10_qvals)

        result = {
            'Trait': trait,
            'corr_prop_emp_pval': corr_emp_pval,
            'corr_prop_qval': corr_qval,
            'n_genes': len(trait_data)
        }
        results.append(result)

        print(f"  {trait} - Genes: {len(trait_data)}, "
              f"Corr(prop_score, -log10_emp_pval): {corr_emp_pval:.4f}, "
              f"Corr(prop_score, -log10_qval): {corr_qval:.4f}")

    # Save results to CSV
    with open("disease_correlations.csv", 'w', newline='') as file:
        fieldnames = ['Trait', 'corr_prop_emp_pval', 'corr_prop_qval', 'n_genes']
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result)

    print("\nResults saved to disease_correlations.csv")

    # Summary statistics
    emp_pval_corrs = [r['corr_prop_emp_pval'] for r in results if not math.isnan(r['corr_prop_emp_pval'])]
    qval_corrs = [r['corr_prop_qval'] for r in results if not math.isnan(r['corr_prop_qval'])]

    print("\n=== SUMMARY STATISTICS ===")
    print("Correlation with -log10(empirical P-value):")
    print(f"  Mean: {statistics.mean(emp_pval_corrs):.4f}")
    print(f"  Median: {statistics.median(emp_pval_corrs):.4f}")
    print(f"  Range: {min(emp_pval_corrs):.4f} to {max(emp_pval_corrs):.4f}")

    print("\nCorrelation with -log10(q-value):")
    print(f"  Mean: {statistics.mean(qval_corrs):.4f}")
    print(f"  Median: {statistics.median(qval_corrs):.4f}")
    print(f"  Range: {min(qval_corrs):.4f} to {max(qval_corrs):.4f}")

    # Top and bottom correlations
    results_sorted_emp = sorted(results, key=lambda x: x['corr_prop_emp_pval'], reverse=True)

    print("\n=TOP 10 highest corr. (Empirical P-value)=")
    for i, result in enumerate(results_sorted_emp[:10], 1):
        print(f"{i:2d}. {result['Trait'][:50]:<50} {result['corr_prop_emp_pval']:7.4f} ({result['n_genes']} genes)")

    print("\n=TOP 10 lowest corr. (Empirical P-value) ===")
    for i, result in enumerate(results_sorted_emp[-10:], 1):
        print(f"{i:2d}. {result['Trait'][:50]:<50} {result['corr_prop_emp_pval']:7.4f} ({result['n_genes']} genes)")

    print("\nAnalysis complete!")
    print("Note: For visualizations, consider installing matplotlib/seaborn or use R for plotting.")

if __name__ == "__main__":
    main()