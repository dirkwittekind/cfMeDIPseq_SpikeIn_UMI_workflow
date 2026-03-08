#!/usr/bin/env python3
"""
MEDIPIPE ML Analysis Report Generator

Generates comprehensive ML discrimination report with:
- Model performance summary (L2 + Elastic-Net comparison)
- Exact permutation p-values
- Robustness check results
- Feature importance / driver genes
- Brief interpretation of results

Author: MEDIPIPE
"""

import argparse
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import numpy as np


def load_json_safe(filepath: str) -> dict:
    """Safely load JSON file."""
    if not os.path.exists(filepath):
        return {}
    try:
        with open(filepath) as f:
            return json.load(f)
    except:
        return {}


def load_tsv_safe(filepath: str) -> pd.DataFrame:
    """Safely load TSV file."""
    if not os.path.exists(filepath):
        return pd.DataFrame()
    try:
        return pd.read_csv(filepath, sep='\t')
    except:
        return pd.DataFrame()


def collect_ml_results(ml_output_dir: str) -> Dict:
    """Collect all ML analysis results."""
    results = {
        'report_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'output_dir': ml_output_dir
    }
    
    # Find all exact permutation results
    exact_perm_dir = os.path.join(ml_output_dir, 'exact_permutation')
    if os.path.exists(exact_perm_dir):
        results['exact_permutation'] = {}
        
        # Load all summary files
        for f in os.listdir(exact_perm_dir):
            if f.endswith('.summary.tsv'):
                feature_space = f.replace('.summary.tsv', '')
                df = load_tsv_safe(os.path.join(exact_perm_dir, f))
                if len(df) > 0:
                    results['exact_permutation'][feature_space] = df.iloc[0].to_dict()
        
        # Load JSON p-values for more detail
        for f in os.listdir(exact_perm_dir):
            if f.endswith('.exact_p_values.json'):
                feature_space = f.replace('.exact_p_values.json', '')
                pvals = load_json_safe(os.path.join(exact_perm_dir, f))
                if pvals and feature_space in results['exact_permutation']:
                    results['exact_permutation'][feature_space].update(pvals)
    
    # Load robustness checks
    robustness_dir = os.path.join(ml_output_dir, 'robustness')
    if os.path.exists(robustness_dir):
        robustness_file = os.path.join(robustness_dir, 'robustness_checks.json')
        results['robustness'] = load_json_safe(robustness_file)
        
        # Also try TSV
        if not results['robustness']:
            robustness_tsv = os.path.join(robustness_dir, 'robustness_checks.tsv')
            df = load_tsv_safe(robustness_tsv)
            if len(df) > 0:
                results['robustness'] = df.to_dict('records')
    
    # Load feature drivers
    drivers_dir = os.path.join(ml_output_dir, 'feature_drivers')
    if os.path.exists(drivers_dir):
        drivers_file = os.path.join(drivers_dir, 'univariate_top_windows.tsv')
        df = load_tsv_safe(drivers_file)
        if len(df) > 0:
            # Get top 20 drivers
            results['top_drivers'] = df.head(20).to_dict('records')
    
    # Load DMR results
    dmr_dirs = [d for d in os.listdir(ml_output_dir) if d.startswith('DMR')]
    if dmr_dirs:
        dmr_dir = os.path.join(ml_output_dir, dmr_dirs[0])
        dmr_summary = os.path.join(dmr_dir, 'dmr.summary.tsv')
        df = load_tsv_safe(dmr_summary)
        if len(df) > 0:
            results['dmr_summary'] = df.iloc[0].to_dict() if len(df) == 1 else df.to_dict('records')
    
    # Load CV results
    cv_dirs = [d for d in os.listdir(ml_output_dir) if d.startswith('windows') and '__' in d]
    if cv_dirs:
        cv_dir = os.path.join(ml_output_dir, cv_dirs[0])
        cv_summary = os.path.join(cv_dir, 'ml_cv.summary.tsv')
        df = load_tsv_safe(cv_summary)
        if len(df) > 0:
            results['cv_summary'] = df.iloc[0].to_dict() if len(df) == 1 else df.to_dict('records')
    
    return results


def generate_interpretation(results: Dict) -> str:
    """Generate brief interpretation of ML results."""
    interp = []
    
    # Check exact permutation results
    if 'exact_permutation' in results:
        perm_results = results['exact_permutation']
        
        # Find best performing feature space
        best_space = None
        best_auroc = 0
        significant_spaces = []
        
        for space, data in perm_results.items():
            auroc = data.get('observed_auroc', 0)
            p_val = data.get('exact_p_auroc', 1)
            
            if auroc > best_auroc:
                best_auroc = auroc
                best_space = space
            
            if p_val < 0.05:
                significant_spaces.append((space, auroc, p_val))
        
        if significant_spaces:
            interp.append(f"**Statistical Significance**: {len(significant_spaces)} feature space(s) achieved statistically significant discrimination (p < 0.05).")
            for space, auroc, pval in significant_spaces:
                interp.append(f"  - {space}: AUROC = {auroc:.3f}, exact p = {pval:.4f}")
        else:
            interp.append("**No Significant Signal**: None of the tested feature spaces achieved statistically significant discrimination (all p > 0.05).")
        
        if best_space:
            interp.append(f"\n**Best Performer**: {best_space} with AUROC = {best_auroc:.3f}")
    
    # Robustness interpretation
    if 'robustness' in results and results['robustness']:
        rob = results['robustness']
        if isinstance(rob, dict):
            # Check depth confounder
            depth_p = rob.get('depth_correlation_p', None)
            if depth_p is not None:
                if depth_p > 0.05:
                    interp.append(f"\n**Depth Confounder**: ✓ No significant correlation between classifier score and sequencing depth (p = {depth_p:.3f})")
                else:
                    interp.append(f"\n**Depth Confounder**: ⚠️ WARNING - Significant correlation with sequencing depth (p = {depth_p:.3f}). Results may be confounded.")
    
    # Feature driver interpretation
    if 'top_drivers' in results and results['top_drivers']:
        drivers = results['top_drivers']
        hypo_count = sum(1 for d in drivers if d.get('direction', '') == 'hypomethylation' or d.get('cohens_d', 0) < 0)
        hyper_count = len(drivers) - hypo_count
        
        interp.append(f"\n**Methylation Direction**: Among top 20 driver features, {hypo_count} show hypomethylation and {hyper_count} show hypermethylation in the case group.")
        
        # Top gene
        top = drivers[0]
        gene = top.get('nearest_gene', top.get('gene', 'Unknown'))
        cohens_d = top.get('cohens_d', 0)
        direction = "hypermethylation" if cohens_d > 0 else "hypomethylation"
        interp.append(f"  - Top driver: {gene} ({direction}, Cohen's d = {cohens_d:+.2f})")
    
    # Model comparison if both L2 and ElasticNet results
    if 'exact_permutation' in results:
        l2_results = {k: v for k, v in results['exact_permutation'].items() if '.L2' in k or k.endswith('_windows')}
        enet_results = {k: v for k, v in results['exact_permutation'].items() if '.ElasticNet' in k}
        
        if l2_results and enet_results:
            l2_key = list(l2_results.keys())[0]
            enet_key = list(enet_results.keys())[0]
            l2_auroc = l2_results[l2_key].get('observed_auroc', 0)
            enet_auroc = enet_results[enet_key].get('observed_auroc', 0)
            
            interp.append(f"\n**Model Comparison**:")
            interp.append(f"  - L2 Logistic Regression: AUROC = {l2_auroc:.3f}")
            interp.append(f"  - Elastic-Net: AUROC = {enet_auroc:.3f}")
            
            if abs(l2_auroc - enet_auroc) < 0.02:
                interp.append("  - Models perform similarly, suggesting robust signal")
            elif l2_auroc > enet_auroc:
                interp.append("  - L2 regularization performs better, suggesting distributed signal across many features")
            else:
                interp.append("  - Elastic-Net performs better, suggesting sparse signal concentrated in fewer features")
    
    # Clinical interpretation
    if 'exact_permutation' in results:
        has_significant = any(
            data.get('exact_p_auroc', 1) < 0.05 
            for data in results['exact_permutation'].values()
        )
        
        interp.append("\n**Clinical Implications**:")
        if has_significant:
            interp.append("  - The significant discrimination suggests cfMeDIP-seq may have diagnostic utility for this comparison")
            interp.append("  - Validation in an independent cohort is required before clinical application")
        else:
            interp.append("  - No significant discrimination was achieved in this cohort")
            interp.append("  - Larger sample sizes or alternative feature spaces may be needed")
    
    return "\n".join(interp)


def generate_html_report(results: Dict, output_file: str, case_group: str, ctrl_group: str):
    """Generate comprehensive HTML ML report."""
    
    interpretation = generate_interpretation(results)
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>MEDIPIPE ML Discrimination Report - {case_group} vs {ctrl_group}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; line-height: 1.6; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 15px; }}
        h2 {{ color: #34495e; margin-top: 40px; border-left: 4px solid #3498db; padding-left: 15px; }}
        h3 {{ color: #7f8c8d; }}
        .summary-box {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 8px; margin: 20px 0; }}
        .summary-box h2 {{ color: white; border-left-color: white; }}
        .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }}
        .metric-card {{ background: #f8f9fa; padding: 15px; border-radius: 5px; text-align: center; border-left: 4px solid #3498db; }}
        .metric-value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
        .metric-label {{ font-size: 12px; color: #7f8c8d; text-transform: uppercase; }}
        .significant {{ border-left-color: #27ae60; }}
        .not-significant {{ border-left-color: #e74c3c; }}
        table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 10px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .interpretation {{ background: #fff3cd; border: 1px solid #ffc107; border-radius: 5px; padding: 20px; margin: 20px 0; }}
        .interpretation h3 {{ margin-top: 0; color: #856404; }}
        .badge {{ display: inline-block; padding: 3px 8px; border-radius: 3px; font-size: 11px; font-weight: bold; }}
        .badge-success {{ background: #d4edda; color: #155724; }}
        .badge-danger {{ background: #f8d7da; color: #721c24; }}
        .badge-warning {{ background: #fff3cd; color: #856404; }}
        pre {{ background: #f4f4f4; padding: 15px; border-radius: 5px; overflow-x: auto; white-space: pre-wrap; }}
        .footer {{ margin-top: 40px; padding-top: 20px; border-top: 2px solid #eee; color: #7f8c8d; font-size: 12px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 MEDIPIPE ML Discrimination Report</h1>
        <p><strong>Comparison:</strong> {case_group} vs {ctrl_group}</p>
        <p><strong>Report Generated:</strong> {results['report_date']}</p>
        
        <div class="summary-box">
            <h2 style="margin-top: 0;">Executive Summary</h2>
"""
    
    # Add executive summary metrics
    if 'exact_permutation' in results:
        for space, data in results['exact_permutation'].items():
            auroc = data.get('observed_auroc', 0)
            pval = data.get('exact_p_auroc', 1)
            sig_class = "significant" if pval < 0.05 else "not-significant"
            sig_badge = '<span class="badge badge-success">SIGNIFICANT</span>' if pval < 0.05 else '<span class="badge badge-danger">NOT SIGNIFICANT</span>'
            
            html += f"""
            <div class="metric-grid">
                <div class="metric-card {sig_class}">
                    <div class="metric-value">{auroc:.3f}</div>
                    <div class="metric-label">AUROC ({space})</div>
                </div>
                <div class="metric-card {sig_class}">
                    <div class="metric-value">{pval:.4f}</div>
                    <div class="metric-label">Exact p-value</div>
                </div>
                <div class="metric-card">
                    {sig_badge}
                    <div class="metric-label" style="margin-top: 5px;">at α=0.05</div>
                </div>
            </div>
"""
    
    html += """        </div>
        
        <div class="interpretation">
            <h3>📊 Interpretation</h3>
            <pre>""" + interpretation + """</pre>
        </div>
        
        <h2>📈 Exact Permutation Test Results</h2>
"""
    
    # Add detailed permutation results table
    if 'exact_permutation' in results:
        html += """
        <table>
            <tr>
                <th>Feature Space</th>
                <th>Model</th>
                <th>AUROC</th>
                <th>Exact p-value</th>
                <th>Null Mean ± SD</th>
                <th>Significant?</th>
            </tr>
"""
        for space, data in results['exact_permutation'].items():
            auroc = data.get('observed_auroc', 0)
            pval = data.get('exact_p_auroc', 1)
            null_mean = data.get('null_auroc_mean', 0.5)
            null_std = data.get('null_auroc_std', 0)
            model = data.get('model_type', 'L2')
            sig = "✓ Yes" if pval < 0.05 else "✗ No"
            
            html += f"""
            <tr>
                <td>{space}</td>
                <td>{model}</td>
                <td><strong>{auroc:.4f}</strong></td>
                <td>{pval:.4f}</td>
                <td>{null_mean:.3f} ± {null_std:.3f}</td>
                <td>{sig}</td>
            </tr>
"""
        html += """        </table>
"""
    
    # Add top drivers table
    if 'top_drivers' in results:
        html += """
        <h2>🧬 Top Discriminative Features</h2>
        <table>
            <tr>
                <th>Rank</th>
                <th>Feature/Gene</th>
                <th>Cohen's d</th>
                <th>Direction</th>
                <th>p-value</th>
            </tr>
"""
        for i, driver in enumerate(results['top_drivers'][:10], 1):
            gene = driver.get('nearest_gene', driver.get('gene', driver.get('feature_id', 'Unknown')))
            cohens_d = driver.get('cohens_d', 0)
            direction = "↑ Hyper" if cohens_d > 0 else "↓ Hypo"
            pval = driver.get('p_value', driver.get('pval', 'N/A'))
            
            html += f"""
            <tr>
                <td>{i}</td>
                <td><strong>{gene}</strong></td>
                <td>{cohens_d:+.3f}</td>
                <td>{direction}</td>
                <td>{pval}</td>
            </tr>
"""
        html += """        </table>
"""
    
    # Add robustness section
    if 'robustness' in results and results['robustness']:
        html += """
        <h2>✅ Robustness Checks</h2>
        <div class="metric-grid">
"""
        rob = results['robustness']
        if isinstance(rob, dict):
            for key, value in rob.items():
                if isinstance(value, (int, float)):
                    html += f"""
            <div class="metric-card">
                <div class="metric-value">{value:.3f}</div>
                <div class="metric-label">{key.replace('_', ' ').title()}</div>
            </div>
"""
        html += """        </div>
"""
    
    # Footer
    html += f"""
        <div class="footer">
            <p>Generated by MEDIPIPE cfMeDIP-seq Analysis Pipeline</p>
            <p>Analysis performed using L2 Logistic Regression and Elastic-Net with exact permutation testing</p>
            <p>Report version: 1.0 | {datetime.now().strftime('%Y-%m-%d')}</p>
        </div>
    </div>
</body>
</html>
"""
    
    with open(output_file, 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(
        description="Generate MEDIPIPE ML analysis report"
    )
    parser.add_argument("--ml-output-dir", required=True, 
                        help="MEDIPIPE ML output directory (e.g., output/ml_discrimination/Outputs/AEG_vs_CTRL)")
    parser.add_argument("--report-dir", required=True, help="Output directory for report")
    parser.add_argument("--case-group", default="Case", help="Case group name")
    parser.add_argument("--ctrl-group", default="Control", help="Control group name")
    parser.add_argument("--format", choices=['html', 'tsv', 'json', 'all'], default='all',
                        help="Output format")
    args = parser.parse_args()
    
    print(f"Collecting ML results from: {args.ml_output_dir}")
    results = collect_ml_results(args.ml_output_dir)
    
    # Create output directory
    os.makedirs(args.report_dir, exist_ok=True)
    
    # Generate reports
    if args.format in ['html', 'all']:
        html_file = os.path.join(args.report_dir, f"ml_report_{args.case_group}_vs_{args.ctrl_group}.html")
        generate_html_report(results, html_file, args.case_group, args.ctrl_group)
        print(f"  HTML report: {html_file}")
    
    if args.format in ['json', 'all']:
        json_file = os.path.join(args.report_dir, f"ml_report_{args.case_group}_vs_{args.ctrl_group}.json")
        # Convert non-serializable items
        def convert_to_serializable(obj):
            if isinstance(obj, (np.integer, np.floating)):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, pd.DataFrame):
                return obj.to_dict('records')
            return obj
        
        serializable_results = json.loads(
            json.dumps(results, default=convert_to_serializable)
        )
        with open(json_file, 'w') as f:
            json.dump(serializable_results, f, indent=2)
        print(f"  JSON report: {json_file}")
    
    if args.format in ['tsv', 'all']:
        # Create summary TSV
        summary_rows = []
        if 'exact_permutation' in results:
            for space, data in results['exact_permutation'].items():
                summary_rows.append({
                    'feature_space': space,
                    'model_type': data.get('model_type', 'L2'),
                    'observed_auroc': data.get('observed_auroc'),
                    'exact_p_auroc': data.get('exact_p_auroc'),
                    'significant': data.get('exact_p_auroc', 1) < 0.05
                })
        
        if summary_rows:
            tsv_file = os.path.join(args.report_dir, f"ml_summary_{args.case_group}_vs_{args.ctrl_group}.tsv")
            pd.DataFrame(summary_rows).to_csv(tsv_file, sep='\t', index=False)
            print(f"  TSV summary: {tsv_file}")
    
    # Print interpretation to console
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    print(generate_interpretation(results))
    print("="*70)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
