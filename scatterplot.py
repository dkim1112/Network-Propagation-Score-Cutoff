#!/usr/bin/env python3

import csv
import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

def create_scatter_plot(page_ranks, neg_log10_vals, trait_name, val_type):
    """Create scatter plot and save as PNG"""

    plt.figure(figsize=(10, 6))

    # Create scatter plot
    plt.scatter(page_ranks, neg_log10_vals, alpha=0.6, s=1)

    # Add trend line
    import numpy as np
    z = np.polyfit(page_ranks, neg_log10_vals, 1)
    p = np.poly1d(z)
    plt.plot(page_ranks, p(page_ranks), "r--", alpha=0.8, linewidth=1)

    # Calculate correlation
    def pearson_correlation(x, y):
        n = len(x)
        mean_x = sum(x) / n
        mean_y = sum(y) / n
        numerator = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
        sum_sq_x = sum((x[i] - mean_x) ** 2 for i in range(n))
        sum_sq_y = sum((y[i] - mean_y) ** 2 for i in range(n))
        denominator = math.sqrt(sum_sq_x * sum_sq_y)
        return numerator / denominator if denominator != 0 else 0

    corr = pearson_correlation(page_ranks, neg_log10_vals)

    plt.xlabel('Propagation Score', fontsize=12)
    plt.ylabel(f'-log10({val_type})', fontsize=12)
    plt.title(f'{trait_name}\nCorrelation: {corr:.4f}', fontsize=14)
    plt.grid(True, alpha=0.3)

    # Save plot
    filename = f"scatter_{trait_name}_{val_type.replace(' ', '_').replace('-', '_')}.png"
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved: {filename}")
    return filename

def find_extremes(page_ranks, neg_log10_vals, trait_name, val_type):
    """Find and display extreme values"""

    combined = list(zip(page_ranks, neg_log10_vals))

    # Sort by propagation score (highest to lowest)
    combined_by_prop = sorted(combined, key=lambda x: x[0], reverse=True)

    # Sort by -log10 values (highest to lowest)
    combined_by_sig = sorted(combined, key=lambda x: x[1], reverse=True)

    print(f"\n=== EXTREME VALUES for {trait_name} ({val_type}) ===")

    print(f"\nTop 10 by Propagation Score:")
    for i, (prop, sig) in enumerate(combined_by_prop[:10], 1):
        print(f"  {i:2d}. Prop: {prop:.6f}, -log10: {sig:.4f}")

    print(f"\nTop 10 by -log10 {val_type}:")
    for i, (prop, sig) in enumerate(combined_by_sig[:10], 1):
        print(f"  {i:2d}. Prop: {prop:.6f}, -log10: {sig:.4f}")

    return combined_by_prop[:10], combined_by_sig[:10]

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

    print(f"Data dimensions: {len(data)} rows")
    print(f"Number of unique traits/diseases: {len(traits)}")

    # Create combined subplot figure
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    fig.suptitle('Propagation Score vs -log10(P-values) for All Diseases', fontsize=16, y=0.98)

    trait_names = list(traits.keys())
    all_extremes = {}

    # Plot for each trait
    for i, trait in enumerate(trait_names):
        print(f"\n{'='*60}")
        print(f"ANALYZING: {trait}")
        print(f"{'='*60}")

        trait_data = traits[trait]
        page_ranks = [row['page_rank'] for row in trait_data]
        neg_log10_emp_pvals = [row['neg_log10_emp_pval'] for row in trait_data]
        neg_log10_qvals = [row['neg_log10_qval'] for row in trait_data]

        # Calculate correlations
        def pearson_correlation(x, y):
            n = len(x)
            mean_x = sum(x) / n
            mean_y = sum(y) / n
            numerator = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
            sum_sq_x = sum((x[i] - mean_x) ** 2 for i in range(n))
            sum_sq_y = sum((y[i] - mean_y) ** 2 for i in range(n))
            denominator = math.sqrt(sum_sq_x * sum_sq_y)
            return numerator / denominator if denominator != 0 else 0

        corr_emp = pearson_correlation(page_ranks, neg_log10_emp_pvals)
        corr_q = pearson_correlation(page_ranks, neg_log10_qvals)

        # Calculate top 5% threshold
        import numpy as np
        prop_threshold = np.percentile(page_ranks, 95)  # Top 5%

        # Separate data into top 5% and rest
        top5_indices = [j for j, p in enumerate(page_ranks) if p >= prop_threshold]
        rest_indices = [j for j, p in enumerate(page_ranks) if p < prop_threshold]

        # Top 5% data
        top5_props_emp = [page_ranks[j] for j in top5_indices]
        top5_vals_emp = [neg_log10_emp_pvals[j] for j in top5_indices]
        top5_props_q = [page_ranks[j] for j in top5_indices]
        top5_vals_q = [neg_log10_qvals[j] for j in top5_indices]

        # Calculate top 5% correlations
        top5_corr_emp = pearson_correlation(top5_props_emp, top5_vals_emp) if len(top5_props_emp) > 1 else 0
        top5_corr_q = pearson_correlation(top5_props_q, top5_vals_q) if len(top5_props_q) > 1 else 0

        # Empirical P-value plot (top row)
        ax_emp = axes[0, i]

        # Plot rest of data in gray
        rest_props_emp = [page_ranks[j] for j in rest_indices]
        rest_vals_emp = [neg_log10_emp_pvals[j] for j in rest_indices]
        ax_emp.scatter(rest_props_emp, rest_vals_emp, alpha=0.3, s=0.3, c='gray')

        # Plot top 5% in red
        ax_emp.scatter(top5_props_emp, top5_vals_emp, alpha=0.8, s=1.5, c='red')

        # Add trend line for ALL data (but make it more subtle)
        if len(page_ranks) > 1:
            z_emp = np.polyfit(page_ranks, neg_log10_emp_pvals, 1)
            p_emp = np.poly1d(z_emp)
            x_range = np.linspace(min(page_ranks), max(page_ranks), 100)
            ax_emp.plot(x_range, p_emp(x_range), "b-", alpha=0.5, linewidth=1, label='All data')

        # Add trend line for top 5% if enough points
        if len(top5_props_emp) > 1:
            z_top = np.polyfit(top5_props_emp, top5_vals_emp, 1)
            p_top = np.poly1d(z_top)
            x_range_top = np.linspace(min(top5_props_emp), max(top5_props_emp), 50)
            ax_emp.plot(x_range_top, p_top(x_range_top), "r--", alpha=0.8, linewidth=2, label='Top 5%')

        # Clean up trait name for title
        clean_trait = trait.replace('finngen_R12_', '').replace('_', ' ')
        ax_emp.set_title(f'{clean_trait}\nAll r={corr_emp:.3f}, Top5% r={top5_corr_emp:.3f}', fontsize=9)
        ax_emp.set_xlabel('Propagation Score', fontsize=8)
        ax_emp.set_ylabel('-log10(Emp P-val)', fontsize=8)
        ax_emp.grid(True, alpha=0.3)

        # Q-value plot (bottom row)
        ax_q = axes[1, i]

        # Plot rest of data in gray
        rest_props_q = [page_ranks[j] for j in rest_indices]
        rest_vals_q = [neg_log10_qvals[j] for j in rest_indices]
        ax_q.scatter(rest_props_q, rest_vals_q, alpha=0.3, s=0.3, c='gray')

        # Plot top 5% in red
        ax_q.scatter(top5_props_q, top5_vals_q, alpha=0.8, s=1.5, c='red')

        # Add trend line for ALL data
        if len(page_ranks) > 1:
            z_q = np.polyfit(page_ranks, neg_log10_qvals, 1)
            p_q = np.poly1d(z_q)
            x_range = np.linspace(min(page_ranks), max(page_ranks), 100)
            ax_q.plot(x_range, p_q(x_range), "b-", alpha=0.5, linewidth=1)

        # Add trend line for top 5%
        if len(top5_props_q) > 1:
            z_top_q = np.polyfit(top5_props_q, top5_vals_q, 1)
            p_top_q = np.poly1d(z_top_q)
            x_range_top = np.linspace(min(top5_props_q), max(top5_props_q), 50)
            ax_q.plot(x_range_top, p_top_q(x_range_top), "r--", alpha=0.8, linewidth=2)

        ax_q.set_title(f'All r={corr_q:.3f}, Top5% r={top5_corr_q:.3f}', fontsize=9)
        ax_q.set_xlabel('Propagation Score', fontsize=8)
        ax_q.set_ylabel('-log10(Q-val)', fontsize=8)
        ax_q.grid(True, alpha=0.3)

        print(f"  Correlations - All data: Emp P-val r={corr_emp:.3f}, Q-val r={corr_q:.3f}")
        print(f"  Correlations - Top 5%: Emp P-val r={top5_corr_emp:.3f}, Q-val r={top5_corr_q:.3f}")
        print(f"  Top 5% threshold: {prop_threshold:.6f} ({len(top5_props_emp)} points)")

        # Find extremes
        print("\n--- EMPIRICAL P-VALUES ---")
        ext_prop_emp, ext_sig_emp = find_extremes(page_ranks, neg_log10_emp_pvals, trait, "Empirical P-value")

        print("\n--- Q-VALUES ---")
        ext_prop_qval, ext_sig_qval = find_extremes(page_ranks, neg_log10_qvals, trait, "Q-value")

        all_extremes[trait] = {
            'emp_by_prop': ext_prop_emp,
            'emp_by_sig': ext_sig_emp,
            'qval_by_prop': ext_prop_qval,
            'qval_by_sig': ext_sig_qval
        }

    # Adjust layout and save combined plot
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)

    combined_filename = "combined_scatter_plots.png"
    plt.savefig(combined_filename, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n{'='*80}")
    print(f"COMBINED SCATTER PLOT SAVED: {combined_filename}")
    print(f"{'='*80}")

    # Summary
    print(f"\nSUMMARY: EXTREME VALUES ACROSS ALL DISEASES")
    print(f"{'='*80}")

    for trait, extremes in all_extremes.items():
        print(f"\n{trait}:")
        print(f"  Highest Prop Score (Emp): {extremes['emp_by_prop'][0][0]:.6f}")
        print(f"  Most Significant (Emp): -log10 = {extremes['emp_by_sig'][0][1]:.4f}")
        print(f"  Highest Prop Score (Qval): {extremes['qval_by_prop'][0][0]:.6f}")
        print(f"  Most Significant (Qval): -log10 = {extremes['qval_by_sig'][0][1]:.4f}")

    print(f"\n{'='*80}")
    print(f"ANALYSIS COMPLETE! Check: {combined_filename}")
    print(f"{'='*80}")

if __name__ == "__main__":
    main()