#!/usr/bin/env python3
"""
MEDIPIPE QC Report Generator
Generates comprehensive QC reports including:
- Deduplicated read counts
- Individual molecules (UMI)
- Methylation QC metrics (relH, GoGe)
- Overamplification signs
- Spike-in DNA recovery

Author: MEDIPIPE Team
"""

import argparse
import os
import re
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False
    print("Warning: pandas not available, some features limited")


class QCMetrics:
    """Container for sample QC metrics"""
    def __init__(self, sample_id):
        self.sample_id = sample_id
        self.raw_reads = None
        self.trimmed_reads = None
        self.aligned_reads = None
        self.dedup_reads = None
        self.unique_molecules = None  # UMI-based
        self.duplication_rate = None
        self.relH = None  # Relative CpG enrichment
        self.GoGe = None  # Global methylation
        self.saturation = None
        self.enrichment_score = None
        self.spikein_recovery = None
        self.spikein_expected = None
        # New spike-in specific metrics
        self.spikein_relH = None  # relH from spike-in reference
        self.spikein_GoGe = None  # GoGe from spike-in reference
        self.spikein_meth_reads = None  # Methylated spike-in read count
        self.spikein_unmeth_reads = None  # Unmethylated spike-in read count
        self.spikein_meth_ratio = None  # Methylated / total spike-in ratio
        self.spikein_total_reads = None  # Total spike-in reads
        self.spikein_status = None  # PASS/FAIL status
        # Additional alignment metrics
        self.insert_size_avg = None
        self.insert_size_std = None
        self.avg_quality = None
        self.error_rate = None
        self.overamplification_flag = False
        self.overamplification_details = ""
        # Tissue of origin deconvolution
        self.too_top_tissues = []  # List of (tissue, proportion) tuples
        self.too_brain_total = None  # Combined brain fraction
        self.too_blood_total = None  # Combined blood fraction
        self.too_liver = None
        self.too_other = None
        self.too_status = None  # Normal/Elevated/etc
        self.qc_pass = True
        self.qc_warnings = []
        
    def to_dict(self):
        return {
            "sample_id": self.sample_id,
            "raw_reads": self.raw_reads,
            "trimmed_reads": self.trimmed_reads,
            "aligned_reads": self.aligned_reads,
            "dedup_reads": self.dedup_reads,
            "unique_molecules": self.unique_molecules,
            "duplication_rate": self.duplication_rate,
            "relH": self.relH,
            "GoGe": self.GoGe,
            "saturation": self.saturation,
            "enrichment_score": self.enrichment_score,
            "spikein_recovery": self.spikein_recovery,
            "spikein_expected": self.spikein_expected,
            "spikein_relH": self.spikein_relH,
            "spikein_GoGe": self.spikein_GoGe,
            "spikein_meth_reads": self.spikein_meth_reads,
            "spikein_unmeth_reads": self.spikein_unmeth_reads,
            "spikein_meth_ratio": self.spikein_meth_ratio,
            "spikein_total_reads": self.spikein_total_reads,
            "spikein_status": self.spikein_status,
            "insert_size_avg": self.insert_size_avg,
            "insert_size_std": self.insert_size_std,
            "avg_quality": self.avg_quality,
            "error_rate": self.error_rate,
            "overamplification_flag": self.overamplification_flag,
            "overamplification_details": self.overamplification_details,
            "too_top_tissues": "; ".join([f"{t}:{p:.1%}" for t, p in self.too_top_tissues]) if self.too_top_tissues else None,
            "too_brain_total": self.too_brain_total,
            "too_blood_total": self.too_blood_total,
            "too_liver": self.too_liver,
            "too_status": self.too_status,
            "qc_pass": self.qc_pass,
            "qc_warnings": "; ".join(self.qc_warnings)
        }


class QCReportGenerator:
    """Generate comprehensive QC reports from MEDIPIPE output"""
    
    # QC thresholds
    THRESHOLDS = {
        "min_dedup_reads": 5_000_000,
        "max_duplication_rate": 0.9,
        "min_relH": 1.5,
        "max_relH": 10.0,
        "min_GoGe": 1.0,
        "min_spikein_recovery": 0.5,
        "min_spikein_relH": 2.5,  # Spike-in relH threshold
        "min_spikein_GoGe": 1.5,  # Spike-in GoGe threshold
        "min_spikein_meth_ratio": 0.7,  # Meth ratio threshold
        "overamp_umi_ratio": 10  # reads per UMI
    }
    
    def __init__(self, output_dir, samples=None, spikein_dir=None):
        self.output_dir = Path(output_dir)
        # Spike-in output may be in a separate directory
        if spikein_dir:
            self.spikein_dir = Path(spikein_dir)
        else:
            # Auto-detect spike-in directory
            parent = self.output_dir.parent
            if (parent / "output_spikein").exists():
                self.spikein_dir = parent / "output_spikein"
            else:
                self.spikein_dir = self.output_dir
        self.samples = samples or self._detect_samples()
        self.metrics = {}
        # Load tissue of origin data once for all samples
        self.too_data = self._load_tissue_of_origin()
        
    def _detect_samples(self):
        """Auto-detect samples from output directory"""
        samples = set()
        
        # Check various output folders
        for pattern in ["meth_qc_quant/*_meth_qc.txt", "dedup_bam_pe/*_dedup.bam"]:
            for f in self.output_dir.glob(pattern):
                sample = re.sub(r'_(meth_qc|dedup).*$', '', f.stem)
                samples.add(sample)
        
        return sorted(samples)
    
    def _load_tissue_of_origin(self):
        """Load tissue of origin deconvolution data from NNLS output"""
        too_data = {}  # sample_id -> {tissue: proportion}
        
        # Look for NNLS proportions file
        too_files = [
            self.output_dir / "tissue_deconv/nnls_proportions.tsv",
            self.output_dir / "tissue_of_origin/nnls_proportions.tsv",
            self.output_dir / "tissue_of_origin/too_fractions.tsv",
        ]
        
        for too_file in too_files:
            if too_file.exists():
                try:
                    with open(too_file) as f:
                        header = next(f).strip().split('\t')
                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) >= 3:
                                # Format: sample\ttissue\tproportion
                                sample_raw = parts[0].strip("'\"")
                                tissue = parts[1]
                                proportion = float(parts[2])
                                
                                # Clean sample name (remove .bw extension etc)
                                sample_id = self._clean_sample_name(sample_raw)
                                
                                if sample_id not in too_data:
                                    too_data[sample_id] = {}
                                too_data[sample_id][tissue] = proportion
                    break
                except Exception as e:
                    print(f"Warning: Error loading TOO data: {e}")
        
        # Also try wide format
        too_wide_files = [
            self.output_dir / "tissue_deconv/nnls_proportions_wide.tsv",
            self.output_dir / "tissue_of_origin/nnls_proportions_wide.tsv",
        ]
        
        if not too_data:  # Only if long format not found
            for too_file in too_wide_files:
                if too_file.exists():
                    try:
                        with open(too_file) as f:
                            header = next(f).strip().split('\t')
                            tissues = header[1:]  # First column is sample_id
                            for line in f:
                                parts = line.strip().split('\t')
                                if len(parts) >= 2:
                                    sample_raw = parts[0].strip("'\"")
                                    sample_id = self._clean_sample_name(sample_raw)
                                    
                                    too_data[sample_id] = {}
                                    for i, tissue in enumerate(tissues):
                                        if i + 1 < len(parts):
                                            try:
                                                too_data[sample_id][tissue] = float(parts[i + 1])
                                            except ValueError:
                                                pass
                        break
                    except Exception as e:
                        print(f"Warning: Error loading TOO wide data: {e}")
        
        return too_data
    
    def _clean_sample_name(self, raw_name):
        """Clean sample name from bigwig filename to match sample ID"""
        # Remove common suffixes
        name = raw_name
        for suffix in ['.bw', '.bigwig', '.bam', '_dedup', '_sorted']:
            name = name.replace(suffix, '')
        # Match lane-tag patterns like '014-CSF_S1_L003' -> '014-CSF'
        name = re.sub(r'_S\d+_L\d{3}$', '', name)
        # Also handle patterns like '149_S7_L004' -> '149'
        name = re.sub(r'_S\d+_L\d{3}.*$', '', name)
        return name
    
    def _parse_tissue_of_origin(self, m):
        """Parse tissue of origin data for a sample"""
        if not self.too_data:
            return
        
        # Try to find matching sample in TOO data
        sample_data = None
        for too_sample, data in self.too_data.items():
            # Try exact match first
            if too_sample == m.sample_id:
                sample_data = data
                break
            # Try partial match (sample ID contained in TOO sample name)
            if m.sample_id in too_sample or too_sample in m.sample_id:
                sample_data = data
                break
        
        if not sample_data:
            return
        
        # Calculate aggregate tissue fractions
        brain_tissues = ['brain_neuron', 'brain_oligodendrocyte']
        blood_tissues = ['blood_monocyte', 'blood_granulocyte', 'blood_bcell', 
                        'blood_bcell_memory', 'blood_tcell_cd4', 'blood_tcell_cd4_central_memory',
                        'blood_tcell_cd8', 'blood_tcell_cd8_effector', 'blood_nkcell']
        
        m.too_brain_total = sum(sample_data.get(t, 0) for t in brain_tissues)
        m.too_blood_total = sum(sample_data.get(t, 0) for t in blood_tissues)
        m.too_liver = sample_data.get('liver', 0)
        
        # Calculate other (non-brain, non-blood, non-liver)
        m.too_other = 1.0 - m.too_brain_total - m.too_blood_total - m.too_liver
        if m.too_other < 0:
            m.too_other = 0  # Rounding errors
        
        # Get top 5 tissues
        sorted_tissues = sorted(sample_data.items(), key=lambda x: x[1], reverse=True)
        m.too_top_tissues = [(t, p) for t, p in sorted_tissues[:5] if p > 0.01]  # >1% only
        
        # Determine status based on expected cfDNA composition
        # Normal plasma cfDNA is predominantly from blood cells
        if m.too_brain_total > 0.20:  # >20% brain-derived
            m.too_status = "ELEVATED_BRAIN"
        elif m.too_liver > 0.30:  # >30% liver
            m.too_status = "ELEVATED_LIVER"
        elif m.too_blood_total < 0.30:  # <30% blood (unusual for plasma)
            m.too_status = "LOW_BLOOD"
        else:
            m.too_status = "NORMAL"
    
    def collect_metrics(self, sample_id):
        """Collect all QC metrics for a sample"""
        m = QCMetrics(sample_id)
        
        # 1. Methylation QC (relH, GoGe)
        self._parse_meth_qc(m)
        
        # 2. Deduplication stats
        self._parse_dedup_stats(m)
        
        # 3. UMI stats
        self._parse_umi_stats(m)
        
        # 4. Spike-in recovery
        self._parse_spikein_stats(m)
        
        # 5. Tissue of origin
        self._parse_tissue_of_origin(m)
        
        # 6. Check for overamplification
        self._check_overamplification(m)
        
        # 7. Apply QC thresholds
        self._apply_thresholds(m)
        
        self.metrics[sample_id] = m
        return m
    
    def _parse_meth_qc(self, m):
        """Parse methylation QC file"""
        qc_file = self.output_dir / f"meth_qc_quant/{m.sample_id}_meth_qc.txt"
        
        if not qc_file.exists():
            m.qc_warnings.append("Missing meth_qc file")
            return
            
        try:
            with open(qc_file) as f:
                content = f.read()
            
            # Parse relH (relative CpG enrichment)
            match = re.search(r'relH[:\s]+([0-9.]+)', content, re.IGNORECASE)
            if match:
                m.relH = float(match.group(1))
            
            # Parse GoGe (global methylation)
            match = re.search(r'GoGe[:\s]+([0-9.]+)', content, re.IGNORECASE)
            if match:
                m.GoGe = float(match.group(1))
            
            # Parse saturation
            match = re.search(r'saturation[:\s]+([0-9.]+)', content, re.IGNORECASE)
            if match:
                m.saturation = float(match.group(1))
            
            # Parse enrichment score
            match = re.search(r'enrichment[:\s]+([0-9.]+)', content, re.IGNORECASE)
            if match:
                m.enrichment_score = float(match.group(1))
                
        except Exception as e:
            m.qc_warnings.append(f"Error parsing meth_qc: {e}")
    
    def _parse_dedup_stats(self, m):
        """Parse deduplication statistics from samtools stats format"""
        stats_file = self.output_dir / f"dedup_bam_pe/{m.sample_id}_dedup.bam.stats.txt"
        
        if not stats_file.exists():
            # Try alternative locations
            for alt in ["stats/{}_dedup_stats.txt", "bam_stats/{}.stats", "dedup_bam_pe/{}_dedup.bam.stats.txt"]:
                alt_file = self.output_dir / alt.format(m.sample_id)
                if alt_file.exists():
                    stats_file = alt_file
                    break
        
        if not stats_file.exists():
            m.qc_warnings.append("Missing dedup stats file")
            return
            
        try:
            with open(stats_file) as f:
                content = f.read()
            
            # Parse samtools stats format (SN lines)
            for line in content.split('\n'):
                if line.startswith('SN\t'):
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        key = parts[1].rstrip(':')
                        value = parts[2].split('#')[0].strip()  # Remove comments
                        
                        if key == 'raw total sequences':
                            m.aligned_reads = int(value)
                        elif key == 'sequences':
                            m.dedup_reads = int(value)  # After dedup, this is the count
                        elif key == 'reads mapped':
                            if not m.dedup_reads:
                                m.dedup_reads = int(value)
                        elif key == 'insert size average':
                            m.insert_size_avg = float(value)
                        elif key == 'insert size standard deviation':
                            m.insert_size_std = float(value)
                        elif key == 'average quality':
                            m.avg_quality = float(value)
                        elif key == 'error rate':
                            m.error_rate = float(value)
            
            # For deduplicated BAMs, aligned_reads = dedup_reads (no duplication info in stats)
            # The file is already deduplicated, so we report the final count
            if m.dedup_reads and not m.aligned_reads:
                m.aligned_reads = m.dedup_reads
                
        except Exception as e:
            m.qc_warnings.append(f"Error parsing dedup stats: {e}")
    
    def _parse_umi_stats(self, m):
        """Parse UMI-based statistics"""
        # Check for UMI tools output
        umi_log = self.output_dir / f"umi_extract/{m.sample_id}_umi_extract.log"
        
        if not umi_log.exists():
            # Try dedup log which may have UMI info
            umi_log = self.output_dir / f"dedup_bam_pe/{m.sample_id}_dedup.log"
        
        if not umi_log.exists():
            return  # UMI stats optional
            
        try:
            with open(umi_log) as f:
                content = f.read()
            
            # Parse unique UMI count
            match = re.search(r'unique UMIs?[:\s]+(\d+)', content, re.IGNORECASE)
            if match:
                m.unique_molecules = int(match.group(1))
            
            # Alternative: reads per UMI
            match = re.search(r'mean reads per UMI[:\s]+([0-9.]+)', content, re.IGNORECASE)
            if match and m.dedup_reads:
                reads_per_umi = float(match.group(1))
                m.unique_molecules = int(m.dedup_reads / reads_per_umi) if reads_per_umi > 0 else None
                
        except Exception as e:
            m.qc_warnings.append(f"Error parsing UMI stats: {e}")
    
    def _parse_spikein_stats(self, m):
        """Parse spike-in recovery statistics from dedicated spike-in output"""
        
        # Look for spike-in QC file (meth_qc.txt)
        spikein_qc_files = [
            self.spikein_dir / f"meth_qc_quant_spikein/{m.sample_id}_meth_qc.txt",
            self.output_dir / f"spikein_qc/{m.sample_id}_spikein_recovery.txt",
            self.output_dir / f"meth_qc_quant/{m.sample_id}_spikein.txt",
        ]
        
        # Parse spike-in QC metrics (relH, GoGe)
        for qc_file in spikein_qc_files:
            if qc_file.exists():
                try:
                    with open(qc_file) as f:
                        lines = f.readlines()
                    
                    if len(lines) >= 2:
                        # TSV format: header + data row
                        header = lines[0].strip().split('\t')
                        values = lines[1].strip().split('\t')
                        
                        data = dict(zip(header, values))
                        
                        # Extract relH and GoGe from spike-in analysis
                        if 'enrichment.relH' in data:
                            m.spikein_relH = float(data['enrichment.relH'])
                        if 'enrichment.GoGe' in data:
                            m.spikein_GoGe = float(data['enrichment.GoGe'])
                        if 'saturation.numberReads' in data:
                            m.saturation = int(data['saturation.numberReads'])
                    break
                except Exception as e:
                    m.qc_warnings.append(f"Error parsing spike-in QC: {e}")
        
        # Look for spike-in quantification file (meth_quant.tsv)
        spikein_quant_files = [
            self.spikein_dir / f"meth_qc_quant_spikein/{m.sample_id}_meth_quant.tsv",
            self.output_dir / f"spikein_quant/{m.sample_id}_meth_quant.tsv",
        ]
        
        for quant_file in spikein_quant_files:
            if quant_file.exists():
                try:
                    meth_count = 0
                    unmeth_count = 0
                    
                    with open(quant_file) as f:
                        next(f)  # Skip header
                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                region = parts[0]
                                count = int(parts[3])
                                
                                if region.endswith('_meth'):
                                    meth_count += count
                                elif region.endswith('_unmeth'):
                                    unmeth_count += count
                    
                    m.spikein_meth_reads = meth_count
                    m.spikein_unmeth_reads = unmeth_count
                    m.spikein_total_reads = meth_count + unmeth_count
                    
                    if m.spikein_total_reads > 0:
                        m.spikein_meth_ratio = meth_count / m.spikein_total_reads
                    
                    # Determine spike-in status
                    m.spikein_status = self._evaluate_spikein_status(m)
                    
                    break
                except Exception as e:
                    m.qc_warnings.append(f"Error parsing spike-in quant: {e}")
    
    def _evaluate_spikein_status(self, m):
        """Evaluate spike-in QC status based on thresholds"""
        T = self.THRESHOLDS
        issues = []
        
        if m.spikein_relH and m.spikein_relH < T["min_spikein_relH"]:
            issues.append(f"low relH ({m.spikein_relH:.2f})")
        
        if m.spikein_GoGe and m.spikein_GoGe < T["min_spikein_GoGe"]:
            issues.append(f"low GoGe ({m.spikein_GoGe:.2f})")
        
        if m.spikein_meth_ratio and m.spikein_meth_ratio < T["min_spikein_meth_ratio"]:
            issues.append(f"low meth ratio ({m.spikein_meth_ratio:.2f})")
        
        if m.spikein_total_reads and m.spikein_total_reads < 10:
            issues.append(f"low reads ({m.spikein_total_reads})")
        
        if not issues:
            return "PASS"
        else:
            return f"WARNING: {', '.join(issues)}"
    
    def _check_overamplification(self, m):
        """Check for signs of overamplification"""
        reasons = []
        
        # High duplication rate
        if m.duplication_rate and m.duplication_rate > 0.85:
            reasons.append(f"High duplication rate: {m.duplication_rate:.1%}")
        
        # Low unique molecules relative to reads
        if m.unique_molecules and m.dedup_reads:
            reads_per_molecule = m.dedup_reads / m.unique_molecules
            if reads_per_molecule > self.THRESHOLDS["overamp_umi_ratio"]:
                reasons.append(f"High reads/UMI ratio: {reads_per_molecule:.1f}")
        
        # Very high relH can indicate overamplification
        if m.relH and m.relH > self.THRESHOLDS["max_relH"]:
            reasons.append(f"Unusually high relH: {m.relH:.2f}")
        
        if reasons:
            m.overamplification_flag = True
            m.overamplification_details = "; ".join(reasons)
    
    def _apply_thresholds(self, m):
        """Apply QC pass/fail thresholds"""
        T = self.THRESHOLDS
        
        if m.dedup_reads and m.dedup_reads < T["min_dedup_reads"]:
            m.qc_pass = False
            m.qc_warnings.append(f"Low read count: {m.dedup_reads:,}")
        
        if m.duplication_rate and m.duplication_rate > T["max_duplication_rate"]:
            m.qc_pass = False
            m.qc_warnings.append(f"Excessive duplication: {m.duplication_rate:.1%}")
        
        if m.relH and m.relH < T["min_relH"]:
            m.qc_warnings.append(f"Low CpG enrichment (relH): {m.relH:.2f}")
        
        if m.GoGe and m.GoGe < T["min_GoGe"]:
            m.qc_warnings.append(f"Low GoGe: {m.GoGe:.2f}")
        
        if m.spikein_recovery and m.spikein_recovery < T["min_spikein_recovery"]:
            m.qc_warnings.append(f"Low spike-in recovery: {m.spikein_recovery:.1%}")
        
        if m.overamplification_flag:
            m.qc_warnings.append("Signs of overamplification detected")
    
    def generate_report(self, output_path=None):
        """Generate the full QC report"""
        # Collect metrics for all samples
        for sample in self.samples:
            self.collect_metrics(sample)
        
        # Generate outputs
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Text report
        report = self._generate_text_report(timestamp)
        
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Save text report
            with open(output_path, 'w') as f:
                f.write(report)
            
            # Save TSV summary if pandas available
            if HAS_PANDAS:
                tsv_path = output_path.with_suffix('.tsv')
                df = pd.DataFrame([m.to_dict() for m in self.metrics.values()])
                df.to_csv(tsv_path, sep='\t', index=False)
            
            # Save JSON for programmatic access
            json_path = output_path.with_suffix('.json')
            with open(json_path, 'w') as f:
                json.dump({
                    "timestamp": timestamp,
                    "samples": [m.to_dict() for m in self.metrics.values()],
                    "thresholds": self.THRESHOLDS
                }, f, indent=2)
        
        return report
    
    def _generate_text_report(self, timestamp):
        """Generate human-readable text report"""
        lines = [
            "=" * 80,
            "MEDIPIPE COMPREHENSIVE QC REPORT",
            "=" * 80,
            f"Generated: {timestamp}",
            f"Output Directory: {self.output_dir}",
            f"Samples: {len(self.samples)}",
            "",
            "-" * 80,
            "SUMMARY",
            "-" * 80,
        ]
        
        # Summary statistics
        passed = sum(1 for m in self.metrics.values() if m.qc_pass)
        overamp = sum(1 for m in self.metrics.values() if m.overamplification_flag)
        
        lines.extend([
            f"Samples Passing QC: {passed}/{len(self.metrics)}",
            f"Samples with Overamplification Signs: {overamp}",
            ""
        ])
        
        # Per-sample details
        lines.extend([
            "-" * 80,
            "PER-SAMPLE METRICS",
            "-" * 80,
            ""
        ])
        
        for sample_id, m in sorted(self.metrics.items()):
            status = "✅ PASS" if m.qc_pass else "❌ FAIL"
            overamp_status = "⚠️ OVERAMP" if m.overamplification_flag else ""
            
            lines.append(f"## {sample_id} [{status}] {overamp_status}")
            lines.append("")
            
            # Read counts
            lines.append("  Read Counts:")
            if m.dedup_reads:
                lines.append(f"    Deduplicated Reads: {m.dedup_reads:,}")
            if m.aligned_reads and m.aligned_reads != m.dedup_reads:
                lines.append(f"    Pre-dedup Reads: {m.aligned_reads:,}")
            if m.unique_molecules:
                lines.append(f"    Unique Molecules (UMI): {m.unique_molecules:,}")
            if m.duplication_rate:
                lines.append(f"    Duplication Rate: {m.duplication_rate:.1%}")
            
            # Alignment quality metrics
            if m.insert_size_avg or m.avg_quality:
                lines.append("  Alignment Quality:")
                if m.insert_size_avg:
                    std_str = f"{m.insert_size_std:.1f}" if m.insert_size_std else "0.0"
                    lines.append(f"    Insert Size: {m.insert_size_avg:.1f} ± {std_str} bp")
                if m.avg_quality:
                    lines.append(f"    Avg Base Quality: {m.avg_quality:.1f}")
                if m.error_rate:
                    lines.append(f"    Error Rate: {m.error_rate:.4f}")
            
            # Methylation QC
            lines.append("  Methylation QC:")
            if m.relH:
                lines.append(f"    relH (CpG enrichment): {m.relH:.2f}")
            if m.GoGe:
                lines.append(f"    GoGe (global methylation): {m.GoGe:.2f}")
            if m.saturation:
                lines.append(f"    Saturation: {m.saturation:.2f}")
            if m.enrichment_score:
                lines.append(f"    Enrichment Score: {m.enrichment_score:.2f}")
            
            # Spike-in Analysis
            if m.spikein_relH or m.spikein_meth_ratio:
                lines.append("  Spike-in Analysis:")
                if m.spikein_relH:
                    lines.append(f"    relH: {m.spikein_relH:.2f} (threshold >2.5)")
                if m.spikein_GoGe:
                    lines.append(f"    GoGe: {m.spikein_GoGe:.2f} (threshold >1.5)")
                if m.spikein_meth_reads is not None:
                    lines.append(f"    Methylated Reads: {m.spikein_meth_reads}")
                if m.spikein_unmeth_reads is not None:
                    lines.append(f"    Unmethylated Reads: {m.spikein_unmeth_reads}")
                if m.spikein_total_reads:
                    lines.append(f"    Total Spike-in Reads: {m.spikein_total_reads}")
                if m.spikein_meth_ratio:
                    lines.append(f"    Meth Ratio: {m.spikein_meth_ratio:.2f} (threshold >0.7)")
                if m.spikein_status:
                    status_icon = "✅" if m.spikein_status == "PASS" else "⚠️"
                    lines.append(f"    Status: {status_icon} {m.spikein_status}")
            elif m.spikein_recovery:
                lines.append("  Spike-in:")
                lines.append(f"    Recovery: {m.spikein_recovery:.1%}")
                if m.spikein_expected:
                    lines.append(f"    Expected: {m.spikein_expected}")
            
            # Tissue of Origin Analysis
            if m.too_top_tissues:
                lines.append("  Tissue of Origin (NNLS Deconvolution):")
                # Status indicator
                if m.too_status:
                    status_icon = "✅" if m.too_status == "NORMAL" else "⚠️"
                    lines.append(f"    Status: {status_icon} {m.too_status}")
                # Aggregate fractions
                lines.append(f"    Brain-derived: {m.too_brain_total:.1%}")
                lines.append(f"    Blood-derived: {m.too_blood_total:.1%}")
                lines.append(f"    Liver: {m.too_liver:.1%}")
                if m.too_other and m.too_other > 0.01:
                    lines.append(f"    Other tissues: {m.too_other:.1%}")
                # Top tissues breakdown
                lines.append("    Top Contributing Tissues:")
                for tissue, prop in m.too_top_tissues:
                    tissue_display = tissue.replace('_', ' ').title()
                    lines.append(f"      - {tissue_display}: {prop:.1%}")
            
            # Overamplification
            if m.overamplification_flag:
                lines.append("  ⚠️ Overamplification Details:")
                lines.append(f"    {m.overamplification_details}")
            
            # Warnings
            if m.qc_warnings:
                lines.append("  Warnings:")
                for w in m.qc_warnings:
                    lines.append(f"    - {w}")
            
            lines.append("")
        
        # Thresholds used
        lines.extend([
            "-" * 80,
            "QC THRESHOLDS APPLIED",
            "-" * 80,
        ])
        for k, v in self.THRESHOLDS.items():
            lines.append(f"  {k}: {v}")
        
        lines.append("=" * 80)
        
        return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Generate comprehensive QC reports for MEDIPIPE samples"
    )
    parser.add_argument("--output-dir", "-d", required=True,
                       help="MEDIPIPE output directory")
    parser.add_argument("--spikein-dir", "-k",
                       help="Spike-in output directory (default: auto-detect)")
    parser.add_argument("--samples", "-s", nargs="+",
                       help="Specific samples to include (default: auto-detect)")
    parser.add_argument("--report", "-r", default="qc_report.txt",
                       help="Output report filename")
    parser.add_argument("--format", choices=["text", "all"], default="all",
                       help="Output format (text or all formats)")
    
    args = parser.parse_args()
    
    # Determine output path
    output_path = Path(args.output_dir) / args.report
    
    # Generate report
    generator = QCReportGenerator(args.output_dir, args.samples, args.spikein_dir)
    report = generator.generate_report(output_path)
    
    print(report)
    print(f"\n✅ Report saved to: {output_path}")
    if HAS_PANDAS:
        print(f"   TSV: {output_path.with_suffix('.tsv')}")
    print(f"   JSON: {output_path.with_suffix('.json')}")


if __name__ == "__main__":
    main()
