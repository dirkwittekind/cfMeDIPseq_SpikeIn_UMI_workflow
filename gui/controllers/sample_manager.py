"""
MEDIPIPE Sample Manager
Handles sample discovery, validation, and file management.
"""

import os
import re
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass


@dataclass
class Sample:
    """Represents a single sample"""
    sample_id: str
    group: str = "default"
    r1_path: str = ""
    r2_path: str = ""
    is_valid: bool = False
    validation_message: str = ""


class SampleManager:
    """Manages sample discovery and validation"""
    
    # Patterns to match FASTQ filenames and extract sample ID + read number
    # Order matters - more specific patterns should come first
    FASTQ_PATTERNS = [
        # Illumina with .merge suffix: 148_S6_L004_R1_001.merge.fastq.gz -> sample_id=148, read=1
        r"(.+?)_S\d+_L\d+_R([12])_\d+\.merge\.(?:fastq|fq)(?:\.gz)?$",
        # Standard Illumina: sample_S1_L001_R1_001.fastq.gz -> sample_id=sample, read=1
        r"(.+?)_S\d+_L\d+_R([12])_\d+\.(?:fastq|fq)(?:\.gz)?$",
        # Illumina without lane: sample_S1_R1_001.fastq.gz
        r"(.+?)_S\d+_R([12])_\d+\.(?:fastq|fq)(?:\.gz)?$",
        # Simple with .merge: sample_R1.merge.fastq.gz
        r"(.+?)_R([12])\.merge\.(?:fastq|fq)(?:\.gz)?$",
        # Simple underscore: sample_1.fastq.gz or sample_R1.fastq.gz
        r"(.+?)_R?([12])\.(?:fastq|fq)(?:\.gz)?$",
        # Dot separator: sample.R1.fastq.gz
        r"(.+?)\.R([12])\.(?:fastq|fq)(?:\.gz)?$",
    ]
    
    def __init__(self):
        self.samples: List[Sample] = []
        self.samples_file: Optional[str] = None
    
    def discover_samples(self, input_dir: str, paired_end: bool = True) -> List[Sample]:
        """
        Discover FASTQ files in a directory and pair them.
        
        Args:
            input_dir: Directory to search for FASTQ files
            paired_end: Whether to look for paired-end files
        
        Returns:
            List of discovered Sample objects
        """
        if not os.path.exists(input_dir):
            return []
        
        # Find all FASTQ files
        fastq_files = []
        for ext in ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]:
            fastq_files.extend(Path(input_dir).rglob(ext))
        
        # Group files by sample ID
        sample_files: Dict[str, Dict[str, str]] = {}
        
        for fq_path in fastq_files:
            filename = fq_path.name
            sample_id, read_num = self._parse_fastq_name(filename)
            
            if sample_id:
                if sample_id not in sample_files:
                    sample_files[sample_id] = {"R1": "", "R2": ""}
                
                key = f"R{read_num}" if read_num else "R1"
                sample_files[sample_id][key] = str(fq_path)
        
        # Create Sample objects
        self.samples = []
        for sample_id, files in sorted(sample_files.items()):
            sample = Sample(
                sample_id=sample_id,
                r1_path=files.get("R1", ""),
                r2_path=files.get("R2", "") if paired_end else ""
            )
            self._validate_sample(sample, paired_end)
            self.samples.append(sample)
        
        return self.samples
    
    def _parse_fastq_name(self, filename: str) -> Tuple[Optional[str], Optional[str]]:
        """Parse FASTQ filename to extract sample ID and read number"""
        for pattern in self.FASTQ_PATTERNS:
            match = re.match(pattern, filename, re.IGNORECASE)
            if match:
                return match.group(1), match.group(2)
        
        # Fallback: use filename without extension as sample ID
        # Handle .merge.fastq.gz, .fastq.gz, .fq.gz, etc.
        base = filename
        for ext in [".merge.fastq.gz", ".merge.fq.gz", ".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
            if base.lower().endswith(ext):
                base = base[:-len(ext)]
                break
        return base, None
    
    def _validate_sample(self, sample: Sample, paired_end: bool) -> None:
        """Validate a sample's files"""
        errors = []
        
        if not sample.r1_path:
            errors.append("R1 file not found")
        elif not os.path.exists(sample.r1_path):
            errors.append(f"R1 file missing: {sample.r1_path}")
        
        if paired_end:
            if not sample.r2_path:
                errors.append("R2 file not found")
            elif not os.path.exists(sample.r2_path):
                errors.append(f"R2 file missing: {sample.r2_path}")
        
        sample.is_valid = len(errors) == 0
        sample.validation_message = "; ".join(errors) if errors else "Valid"
    
    def add_sample(self, sample: Sample) -> None:
        """Add a sample to the list"""
        self.samples.append(sample)
    
    def add_files(self, file_paths: List[str], paired_end: bool = True) -> List[Sample]:
        """
        Add specific FASTQ files as samples (only the selected files, not the whole directory).
        
        Args:
            file_paths: List of paths to FASTQ files
            paired_end: Whether to pair R1/R2 files
        
        Returns:
            List of newly added Sample objects
        """
        # Group files by sample ID
        sample_files: Dict[str, Dict[str, str]] = {}
        
        for fq_path in file_paths:
            if not os.path.isfile(fq_path):
                continue
            
            filename = os.path.basename(fq_path)
            sample_id, read_num = self._parse_fastq_name(filename)
            
            if sample_id:
                if sample_id not in sample_files:
                    sample_files[sample_id] = {"R1": "", "R2": ""}
                
                key = f"R{read_num}" if read_num else "R1"
                sample_files[sample_id][key] = str(fq_path)
        
        # Create Sample objects for new samples only
        new_samples = []
        existing_ids = {s.sample_id for s in self.samples}
        
        for sample_id, files in sorted(sample_files.items()):
            if sample_id in existing_ids:
                # Update existing sample if it has missing paths
                for existing in self.samples:
                    if existing.sample_id == sample_id:
                        if not existing.r1_path and files.get("R1"):
                            existing.r1_path = files["R1"]
                        if not existing.r2_path and files.get("R2") and paired_end:
                            existing.r2_path = files["R2"]
                        self._validate_sample(existing, paired_end)
                        break
            else:
                sample = Sample(
                    sample_id=sample_id,
                    r1_path=files.get("R1", ""),
                    r2_path=files.get("R2", "") if paired_end else ""
                )
                self._validate_sample(sample, paired_end)
                self.samples.append(sample)
                new_samples.append(sample)
        
        return new_samples
    
    def remove_sample(self, sample_id: str) -> bool:
        """Remove a sample by ID"""
        for i, sample in enumerate(self.samples):
            if sample.sample_id == sample_id:
                del self.samples[i]
                return True
        return False
    
    def update_sample_group(self, sample_id: str, group: str) -> bool:
        """Update a sample's group"""
        for sample in self.samples:
            if sample.sample_id == sample_id:
                sample.group = group
                return True
        return False
    
    def load_samples_tsv(self, tsv_path: str, paired_end: bool = True) -> List[Sample]:
        """Load samples from a TSV file"""
        if not os.path.exists(tsv_path):
            raise FileNotFoundError(f"Samples file not found: {tsv_path}")
        
        df = pd.read_csv(tsv_path, sep="\t")
        
        self.samples = []
        for _, row in df.iterrows():
            sample = Sample(
                sample_id=str(row.get("sample_id", row.get("sample", ""))),
                group=str(row.get("group", "default")),
                r1_path=str(row.get("R1", row.get("fastq1", ""))),
                r2_path=str(row.get("R2", row.get("fastq2", ""))) if paired_end else ""
            )
            self._validate_sample(sample, paired_end)
            self.samples.append(sample)
        
        self.samples_file = tsv_path
        return self.samples
    
    def save_samples_tsv(self, tsv_path: str) -> str:
        """Save samples to a TSV file"""
        if not self.samples:
            raise ValueError("No samples to save")
        
        data = []
        for sample in self.samples:
            data.append({
                "sample_id": sample.sample_id,
                "group": sample.group,
                "R1": sample.r1_path,
                "R2": sample.r2_path
            })
        
        df = pd.DataFrame(data)
        
        # Ensure directory exists
        os.makedirs(os.path.dirname(tsv_path), exist_ok=True)
        
        df.to_csv(tsv_path, sep="\t", index=False)
        self.samples_file = tsv_path
        return tsv_path
    
    def get_valid_samples(self) -> List[Sample]:
        """Get only valid samples"""
        return [s for s in self.samples if s.is_valid]
    
    def get_invalid_samples(self) -> List[Sample]:
        """Get only invalid samples"""
        return [s for s in self.samples if not s.is_valid]
    
    def get_sample_ids(self) -> List[str]:
        """Get list of all sample IDs"""
        return [s.sample_id for s in self.samples]
    
    def get_groups(self) -> List[str]:
        """Get unique groups"""
        return list(set(s.group for s in self.samples))
    
    def validate_all(self, paired_end: bool = True) -> Tuple[int, int]:
        """
        Validate all samples and return (valid_count, invalid_count)
        """
        valid = 0
        invalid = 0
        
        for sample in self.samples:
            self._validate_sample(sample, paired_end)
            if sample.is_valid:
                valid += 1
            else:
                invalid += 1
        
        return valid, invalid
    
    def discover_processed_samples(
        self,
        bam_dir: str = None,
        bigwig_dir: str = None,
        annotation_file: str = None
    ) -> List[Sample]:
        """
        Discover samples from processed files (BAM/BigWig) and annotation.
        
        This is used for post-analysis workflows like ML discrimination,
        where we start from processed files rather than raw FASTQ.
        
        Args:
            bam_dir: Directory containing deduplicated BAM files
            bigwig_dir: Directory containing BigWig files
            annotation_file: Sample annotation TSV with sample_id and group columns
        
        Returns:
            List of discovered Sample objects with groups assigned
        """
        self.samples = []
        sample_groups = {}
        
        # Load groups from annotation file if provided
        if annotation_file and os.path.exists(annotation_file):
            sample_groups = self._load_groups_from_annotation(annotation_file)
        
        # Discover from BigWig directory
        if bigwig_dir and os.path.exists(bigwig_dir):
            for bw_file in Path(bigwig_dir).glob("*.bw"):
                sample_id = bw_file.stem
                group = sample_groups.get(sample_id, self._infer_group_from_name(sample_id))
                
                sample = Sample(
                    sample_id=sample_id,
                    group=group,
                    r1_path=str(bw_file),  # Store bigwig path in r1_path
                    is_valid=True,
                    validation_message="BigWig file"
                )
                self.samples.append(sample)
        
        # Discover from BAM directory (if no bigwigs found)
        if not self.samples and bam_dir and os.path.exists(bam_dir):
            for bam_file in Path(bam_dir).glob("*.bam"):
                sample_id = bam_file.stem
                group = sample_groups.get(sample_id, self._infer_group_from_name(sample_id))
                
                sample = Sample(
                    sample_id=sample_id,
                    group=group,
                    r1_path=str(bam_file),  # Store bam path in r1_path
                    is_valid=os.path.exists(str(bam_file) + ".bai"),
                    validation_message="BAM file" if os.path.exists(str(bam_file) + ".bai") else "BAM index missing"
                )
                self.samples.append(sample)
        
        return self.samples
    
    def _load_groups_from_annotation(self, annotation_file: str) -> Dict[str, str]:
        """
        Load sample-to-group mapping from annotation file.
        
        Supports various column naming conventions:
        - sample_id, sample, Sample, SAMPLE
        - group, Group, GROUP, condition, Condition
        """
        try:
            df = pd.read_csv(annotation_file, sep="\t")
            df.columns = [c.strip() for c in df.columns]
            
            # Find sample ID column
            sample_col = None
            for col in ["sample_id", "sample", "Sample", "SAMPLE", "SampleID"]:
                if col in df.columns:
                    sample_col = col
                    break
            
            # Find group column
            group_col = None
            for col in ["group", "Group", "GROUP", "condition", "Condition", "class", "Class"]:
                if col in df.columns:
                    group_col = col
                    break
            
            if sample_col is None or group_col is None:
                return {}
            
            # Build mapping
            sample_groups = {}
            for _, row in df.iterrows():
                sample_id = str(row[sample_col]).strip()
                group = str(row[group_col]).strip()
                # Normalize control group names
                if group.lower() in ["kontrolle", "control", "ctrl", "healthy", "normal"]:
                    group = "CTRL"
                sample_groups[sample_id] = group
            
            return sample_groups
        
        except Exception as e:
            print(f"Warning: Could not load annotation file: {e}")
            return {}
    
    def _infer_group_from_name(self, sample_id: str) -> str:
        """
        Try to infer group from sample ID naming convention.
        
        Common patterns:
        - CTRL_001, Control_1 -> CTRL
        - AEG_001, Case_1 -> AEG (or other case name)
        - Sample names containing disease/condition keywords
        """
        sample_lower = sample_id.lower()
        
        # Control patterns
        control_patterns = ["ctrl", "control", "healthy", "normal", "hd", "nc"]
        for pattern in control_patterns:
            if pattern in sample_lower:
                return "CTRL"
        
        # Common disease/case patterns - extract the prefix
        # E.g., AEG_001 -> AEG, CRC_S1 -> CRC
        parts = sample_id.split("_")
        if len(parts) >= 2:
            prefix = parts[0]
            # Check if prefix looks like a group (2-4 uppercase letters)
            if prefix.isupper() and 2 <= len(prefix) <= 6:
                return prefix
        
        return "unknown"
    
    def detect_available_groups(self) -> List[str]:
        """
        Detect all unique groups from current samples.
        Automatically merges similar group names.
        """
        groups = set()
        for sample in self.samples:
            if sample.group and sample.group != "unknown":
                groups.add(sample.group)
        
        # Sort with CTRL/Control first if present
        group_list = sorted(list(groups))
        if "CTRL" in group_list:
            group_list.remove("CTRL")
            group_list = ["CTRL"] + group_list
        
        return group_list
    
    def get_samples_by_group(self, group: str) -> List[Sample]:
        """Get all samples belonging to a specific group."""
        return [s for s in self.samples if s.group == group]
    
    def get_group_counts(self) -> Dict[str, int]:
        """Get count of samples per group."""
        counts = {}
        for sample in self.samples:
            group = sample.group or "unknown"
            counts[group] = counts.get(group, 0) + 1
        return counts
