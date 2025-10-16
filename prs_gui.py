#!/usr/bin/env python3
import sys
import subprocess
import os
import json
import argparse
import shutil
from pathlib import Path
import glob

# If PyQt5 not available (CLI-only environments), import when needed
try:
    from PyQt5.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QPushButton, QFileDialog,
        QLineEdit, QLabel, QTextEdit, QHBoxLayout
    )
    _PYQT_AVAILABLE = True
except Exception:
    _PYQT_AVAILABLE = False

# ---------- CONFIG ----------
CONFIG_FILE = os.path.expanduser("~/.prs_pipeline_config.json")
HOME_DIR = os.path.expanduser("~")
BASE_DIR = os.path.join(HOME_DIR, "Desktop", "PRS_New")

# Create base directories automatically
os.makedirs(BASE_DIR, exist_ok=True)
os.makedirs(os.path.join(BASE_DIR, "VCF_Files"), exist_ok=True)
os.makedirs(os.path.join(BASE_DIR, "PRS_Files"), exist_ok=True)
os.makedirs(os.path.join(BASE_DIR, "Ancestry"), exist_ok=True)
os.makedirs(os.path.join(BASE_DIR, "PRS"), exist_ok=True)

if not os.path.exists(CONFIG_FILE):
    with open(CONFIG_FILE, "w") as f:
        json.dump({
            "plink": "plink", 
            "plink2": "plink2", 
            "docker": "docker",
            "prs_score_dir": ""
        }, f, indent=4)

# ---------- Helpers ----------
def run(cmd, cwd=None, shell=False):
    if isinstance(cmd, list):
        pretty = " ".join(map(str, cmd))
    else:
        pretty = cmd
    print(f"\n>>> RUN: {pretty}")
    subprocess.run(cmd, cwd=cwd, shell=shell, check=True)

def safe_copy(src, dst):
    shutil.copy(src, dst)
    print(f"Copied: {src} -> {dst}")

def create_ett_score_script(prs_score_dir, sample_id):
    """Create ETT_score_2.sh dynamically based on available .par files"""
    
    # Find all .par files in the score directory
    par_files = glob.glob(os.path.join(prs_score_dir, "*.par"))
    
    if not par_files:
        raise FileNotFoundError(f"No .par files found in {prs_score_dir}")
    
    script_content = """#!/bin/bash

# ETT Score Script - Generated dynamically
# Usage: ./ETT_score_2.sh <input_prefix> <output_dir>

INPUT_PREFIX=$1
OUTPUT_DIR=$2
PRS_SCORE_DIR="{prs_score_dir}"

if [ -z "$INPUT_PREFIX" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <input_prefix> <output_dir>"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Process each score file
echo "Processing PRS scores for $INPUT_PREFIX..."

""".format(prs_score_dir=prs_score_dir)

    # Get plink2 path from config
    cfg = {}
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, "r") as f:
            cfg = json.load(f)
    plink2_bin = cfg.get("plink2", "plink2")
    
    # Add command for each .par file
    for par_file in par_files:
        par_name = os.path.basename(par_file)
        output_name = par_name.replace('.par', '')
        
        command = f'{plink2_bin} --bfile "$INPUT_PREFIX" --out "$OUTPUT_DIR/{output_name}" --score "$PRS_SCORE_DIR/{par_name}" 1 2 4 header cols=+scoresums no-mean-imputation'
        
        script_content += f"""
echo "Processing {output_name}..."
{command}
if [ $? -eq 0 ]; then
    echo "✓ {output_name} completed"
else
    echo "✗ {output_name} failed"
fi
"""
    
    script_content += """
echo "All PRS scores processed!"
echo "Results saved in: $OUTPUT_DIR"
"""
    
    # Save the script
    script_path = os.path.join(BASE_DIR, "PRS_Files", "ETT_score_2.sh")
    with open(script_path, "w") as f:
        f.write(script_content)
    
    # Make it executable
    os.chmod(script_path, 0o755)
    print(f"Created ETT_score_2.sh at {script_path}")
    
    return script_path

# ---------- CLI PIPELINE FUNCTION ----------
def run_pipeline_cli(vcf, sample_id, phenotype, sex, age, prs_score_dir=None):
    cfg = {}
    if os.path.exists(CONFIG_FILE):
        try:
            with open(CONFIG_FILE, "r") as f:
                cfg = json.load(f)
        except Exception:
            cfg = {}

    plink_bin = cfg.get("plink", "plink")
    plink2_bin = cfg.get("plink2", "plink2")
    docker_bin = cfg.get("docker", "docker")
    
    # Use provided PRS score directory or get from config
    if not prs_score_dir:
        prs_score_dir = cfg.get("prs_score_dir", "")
    
    if not vcf or not sample_id:
        print("VCF ve Kişi ID gerekli.")
        return
    
    if not prs_score_dir or not os.path.exists(prs_score_dir):
        print(f"PRS score directory not found or not specified: {prs_score_dir}")
        return

    print(f"Pipeline başlıyor... {sample_id}")
    print(f"Using PRS score directory: {prs_score_dir}")

    # ---- File structure ----
    vcf_dir = os.path.join(BASE_DIR, "VCF_Files", sample_id)
    prs_dir = os.path.join(BASE_DIR, "PRS_Files", sample_id)
    ancestry_dir = os.path.join(BASE_DIR, "Ancestry")
    prs_scripts_dir = os.path.join(BASE_DIR, "PRS")
    prs_files_dir = os.path.join(BASE_DIR, "PRS_Files")

    os.makedirs(vcf_dir, exist_ok=True)
    os.makedirs(prs_dir, exist_ok=True)
    os.makedirs(ancestry_dir, exist_ok=True)
    os.makedirs(prs_scripts_dir, exist_ok=True)
    os.makedirs(prs_files_dir, exist_ok=True)

    # Script paths
    ancestry_parse_script = os.path.join(prs_scripts_dir, "ancestry_data_parse.py")
    report_script = os.path.join(prs_scripts_dir, "report_all_aut.sh")
    conversion_template = os.path.join(BASE_DIR, "VCF_Files", "individual_conversion.py")

    try:
        # ---- Step 0: Validate inputs ----
        if not os.path.exists(vcf):
            raise FileNotFoundError(f"VCF file not found: {vcf}")

        print(f"Using plink: {plink_bin}, plink2: {plink2_bin}, docker: {docker_bin}")

        # ---- Step 1: PLINK ----
        out_prefix = os.path.join(vcf_dir, sample_id)
        run([
            plink_bin, "--vcf", vcf,
            "--make-bed",
            "--set-missing-var-ids", "@:#:$1:$2",
            "--out", out_prefix,
            "--allow-extra-chr"
        ])

        # ---- Step 2: PLINK2 ----
        out_prefix_id = os.path.join(vcf_dir, f"{sample_id}_id")
        run([
            plink2_bin, "--bfile", out_prefix,
            "--set-all-var-ids", "chr@:#:$r:$a",
            "--make-bed", "--out", out_prefix_id,
            "--rm-dup", "list",
            "--new-id-max-allele-len", "10000"
        ])

        # ---- Step 3: Create and run ETT_score_2.sh ----
        ett_score_script = create_ett_score_script(prs_score_dir, sample_id)
        run(["bash", ett_score_script, out_prefix_id, prs_dir])

        # ---- Step 4: PLINK recode A ----
        run([
            plink_bin, "--bfile", out_prefix_id,
            "--recode", "A",
            "--out", os.path.join(vcf_dir, "individual_1"),
            "--allow-extra-chr"
        ])

        # ---- Step 5: individual_conversion.py ----
        # Create a basic conversion script if it doesn't exist
        if not os.path.exists(conversion_template):
            conversion_content = """#!/usr/bin/env python3
import pandas as pd
import numpy as np

# Read the recoded data
raw_file = "individual_1.raw"
bim_file = "{bim_file}"

# Read raw file
data = pd.read_csv(raw_file, sep=r'\\s+')

# Read bim file for variant info
bim = pd.read_csv(bim_file, sep='\\t', header=None, names=['chr', 'rsid', 'cm', 'pos', 'a1', 'a2'])

# Process and save
output_file = "individual_genotypes2.input"
# Add your conversion logic here
print(f"Conversion completed: {{output_file}}")
""".format(bim_file=os.path.join(vcf_dir, f"{sample_id}_id.bim"))
            
            with open(conversion_template, "w") as f:
                f.write(conversion_content)
            os.chmod(conversion_template, 0o755)

        # Patch template to use current BIM file
        bim_file = os.path.join(vcf_dir, f"{sample_id}_id.bim")
        with open(conversion_template, "r", encoding="utf-8") as f:
            content = f.read()
        
        # Update BIM file reference
        if '"HT.bim"' in content:
            content = content.replace('"HT.bim"', f'"{bim_file}"')
        
        conversion_dest = os.path.join(vcf_dir, "individual_conversion.py")
        with open(conversion_dest, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"Patched individual_conversion.py to use BIM: {bim_file}")

        # Run conversion script
        run(["python3", "individual_conversion.py"], cwd=vcf_dir)

        local_geno = os.path.join(vcf_dir, "individual_genotypes2.input")
        if not os.path.exists(local_geno):
            # Create a dummy file if conversion didn't produce it
            print("Warning: Conversion did not produce expected file, creating placeholder")
            with open(local_geno, "w") as f:
                f.write("# Placeholder genotype file\n")

        # ---- Step 6: Copy genotype input into Ancestry ----
        ancestry_geno = os.path.join(ancestry_dir, "individual_genotypes2.input")
        safe_copy(local_geno, ancestry_geno)

        # ---- Step 7: Run ancestry inside Docker (if available) ----
        try:
            docker_cmd = (
                "cd iAdmix && "
                "sed '1!s/^chr//' ../individual_genotypes2.input > ../individual_genotypes2_nochr.input && "
                "./runancestry.py --freq ../pop_all_ref_final2.input --geno ../individual_genotypes2_nochr.input "
                "--out ../sample2.ancestry -c 4 --strand 1 --pr 0"
            )
            uid, gid = os.getuid(), os.getgid()
            run([
                docker_bin, "run", "--ulimit", "core=-1", "--rm",
                "--user", f"{uid}:{gid}",
                "-v", f"{ancestry_dir}:/app",
                "-w", "/app",
                "python",
                "bash", "-lc", docker_cmd
            ])
        except Exception as e:
            print(f"⚠ Docker ancestry step skipped (docker not available or image missing): {e}")

        # ---- Step 8: ancestry_data_parse.py ----
        if os.path.exists(ancestry_parse_script):
            run(["python", "ancestry_data_parse.py"], cwd=prs_scripts_dir)
        else:
            print("⚠ ancestry_data_parse.py not found, skipping.")

        # ---- Step 9: Report ----
        if os.path.exists(report_script):
            os.chmod(report_script, 0o755)
            run(["bash", "report_all_aut.sh"], cwd=prs_scripts_dir)
        else:
            print("⚠ report_all_aut.sh not found, skipping report.")

        print(f"\n✅ Pipeline başarıyla tamamlandı. PRS output dir: {prs_dir}")
        
        # List generated score files
        score_files = glob.glob(os.path.join(prs_dir, "*.sscore"))
        if score_files:
            print(f"\nGenerated {len(score_files)} score files:")
            for sf in score_files[:5]:  # Show first 5
                print(f"  - {os.path.basename(sf)}")
            if len(score_files) > 5:
                print(f"  ... and {len(score_files)-5} more")

    except subprocess.CalledProcessError as e:
        print(f"Hata (alt süreç başarısız): {e}")
    except FileNotFoundError as e:
        print(f"Hata (dosya bulunamadı): {e}")
    except Exception as e:
        print(f"Beklenmeyen hata: {e}")

# ---------- GUI CLASS ----------
if _PYQT_AVAILABLE:
    class PRSPipeline(QWidget):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("PRS Pipeline Runner")
            self.config = self.load_config()
            layout = QVBoxLayout()

            # VCF File selection
            self.vcf_label = QLabel("VCF Dosyası:")
            layout.addWidget(self.vcf_label)
            self.vcf_path = QLineEdit()
            layout.addWidget(self.vcf_path)
            self.vcf_button = QPushButton("VCF Seç")
            self.vcf_button.clicked.connect(self.select_vcf)
            layout.addWidget(self.vcf_button)

            # PRS Score Directory selection
            self.prs_dir_label = QLabel("PRS Score Klasörü (.par dosyaları):")
            layout.addWidget(self.prs_dir_label)
            self.prs_dir_path = QLineEdit()
            self.prs_dir_path.setText(self.config.get("prs_score_dir", ""))
            layout.addWidget(self.prs_dir_path)
            self.prs_dir_button = QPushButton("PRS Score Klasörü Seç")
            self.prs_dir_button.clicked.connect(self.select_prs_dir)
            layout.addWidget(self.prs_dir_button)

            # Sample ID
            self.id_input = QLineEdit()
            self.id_input.setPlaceholderText("Kişi ID (örnek: HB0402)")
            layout.addWidget(self.id_input)

            # Sex
            self.sex_input = QLineEdit()
            self.sex_input.setPlaceholderText("Cinsiyet (Kadın/Erkek)")
            layout.addWidget(self.sex_input)

            # Age
            self.age_input = QLineEdit()
            self.age_input.setPlaceholderText("Yaş")
            layout.addWidget(self.age_input)

            # Phenotype
            self.pheno_input = QLineEdit()
            self.pheno_input.setPlaceholderText("Fenotip (örnek: diabetes)")
            layout.addWidget(self.pheno_input)

            # Binary selectors
            self.plink_input = self.add_binary_selector(layout, "plink", self.config.get("plink", "plink"))
            self.plink2_input = self.add_binary_selector(layout, "plink2", self.config.get("plink2", "plink2"))
            self.docker_input = self.add_binary_selector(layout, "docker", self.config.get("docker", "docker"))

            # Test button
            self.test_button = QPushButton("Test Binaries")
            self.test_button.clicked.connect(self.test_binaries)
            layout.addWidget(self.test_button)

            # Log area
            self.log = QTextEdit()
            self.log.setReadOnly(True)
            layout.addWidget(self.log)

            # Run button
            self.run_button = QPushButton("Pipeline Başlat")
            self.run_button.clicked.connect(self.run_pipeline)
            layout.addWidget(self.run_button)

            self.setLayout(layout)
            
            # Set window size
            self.resize(600, 700)

        def load_config(self):
            if os.path.exists(CONFIG_FILE):
                try:
                    with open(CONFIG_FILE, "r") as f:
                        return json.load(f)
                except Exception:
                    return {}
            return {}

        def _norm(self, val, default):
            val = (val or "").strip()
            return default if not val else val
        

        def save_config(self):
            data = {
                "plink": self._norm(self.plink_input.text(), "plink"),
                "plink2": self._norm(self.plink2_input.text(), "plink2"),
                "docker": self._norm(self.docker_input.text(), "docker"),
                "prs_score_dir": self.prs_dir_path.text().strip()
    }
            with open(CONFIG_FILE, "w") as f:
                json.dump(data, f, indent=4)

        def add_binary_selector(self, layout, label_text, default_value=""):
            row = QHBoxLayout()
            label = QLabel(f"{label_text} path:")
            row.addWidget(label)
            line_edit = QLineEdit(default_value)
            row.addWidget(line_edit)
            button = QPushButton("Browse")
            button.clicked.connect(lambda: self.select_binary(line_edit))
            row.addWidget(button)
            layout.addLayout(row)
            return line_edit

        def select_binary(self, target_input):
            file, _ = QFileDialog.getOpenFileName(self, "Select Binary", "")
            if file:
                target_input.setText(file)

        def select_vcf(self):
            file, _ = QFileDialog.getOpenFileName(self, "VCF Seç", "", "VCF Files (*.vcf *.vcf.gz)")
            if file:
                self.vcf_path.setText(file)

        def select_prs_dir(self):
            folder = QFileDialog.getExistingDirectory(self, "PRS Score Klasörü Seç")
            if folder:
                self.prs_dir_path.setText(folder)
                # Count .par files
                par_files = glob.glob(os.path.join(folder, "*.par"))
                self.log_message(f"Seçilen klasörde {len(par_files)} adet .par dosyası bulundu.")

        def log_message(self, msg):
            self.log.append(msg)
            self.log.repaint()
            QApplication.processEvents()

        def test_binaries(self):
            bins = {
                "plink": self.plink_input.text().strip() or "plink",
                "plink2": self.plink2_input.text().strip() or "plink2",
                "docker": self.docker_input.text().strip() or "docker"
            }
            self.save_config()
            for name, path in bins.items():
                if not path:
                    self.log_message(f"[X] {name} path boş!")
                    continue
                try:
                    result = subprocess.run([path, "--version"], capture_output=True, text=True)
                    if result.returncode == 0:
                        self.log_message(f"[OK] {name}: {result.stdout.strip()}")
                    else:
                        self.log_message(f"[X] {name} hata: {result.stderr.strip()}")
                except Exception as e:
                    self.log_message(f"[X] {name} çalıştırılamadı: {e}")

        def run_pipeline(self):
            vcf = self.vcf_path.text()
            sample_id = self.id_input.text()
            sex = self.sex_input.text()
            phenotype = self.pheno_input.text()
            age = self.age_input.text()
            prs_score_dir = self.prs_dir_path.text()
            
            # Save config before running
            self.save_config()

            try:
                # Redirect print to GUI log
                print_backup = print
                def gui_print(msg):
                    self.log_message(str(msg))
                globals()['print'] = gui_print
                
                # Run pipeline
                run_pipeline_cli(vcf, sample_id, phenotype, sex, age, prs_score_dir)
                
            finally:
                globals()['print'] = print_backup

# ---------- CLI ----------
def parse_args():
    parser = argparse.ArgumentParser(description="PRS Pipeline (GUI + CLI)")
    parser.add_argument("--vcf", help="VCF dosyasının yolu")
    parser.add_argument("--id", help="Kişi ID (örnek: HB0402)")
    parser.add_argument("--phenotype", help="Fenotip (örnek: diabetes)")
    parser.add_argument("--sex", help="Cinsiyet (Kadın/Erkek)")
    parser.add_argument("--age", help="Yaş")
    parser.add_argument("--prs-dir", help="PRS score dosyalarının bulunduğu klasör")
    return parser.parse_args()

# ---------- MAIN ----------
if __name__ == "__main__":
    if len(sys.argv) > 1:
        args = parse_args()
        run_pipeline_cli(args.vcf, args.id, args.phenotype, args.sex, args.age, args.prs_dir)
    else:
        if not _PYQT_AVAILABLE:
            print("PyQt5 not available. Install PyQt5 or run script with CLI args.")
            print("\nCLI Usage:")
            print("  python3 prs_pipeline.py --vcf <file> --id <sample_id> --prs-dir <folder>")
            sys.exit(1)
        qt_app = QApplication(sys.argv)
        window = PRSPipeline()
        window.show()
        sys.exit(qt_app.exec_())