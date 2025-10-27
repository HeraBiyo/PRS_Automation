#!/usr/bin/env python3
import sys
import subprocess
import os
import json
import argparse
import shutil
from pathlib import Path

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

if not os.path.exists(CONFIG_FILE):
    with open(CONFIG_FILE, "w") as f:
        json.dump({"plink": "plink", "plink2": "plink2", "docker": "docker"}, f, indent=4)

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

# ---------- CLI PIPELINE FUNCTION ----------
def run_pipeline_cli(vcf, sample_id, phenotype, sex, age):
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

    if not vcf or not sample_id:
        print("VCF ve Kişi ID gerekli.")
        return

    print(f"Pipeline başlıyor... {sample_id}")

    # ---- File structure ----
    vcf_dir = os.path.join(BASE_DIR, "VCF_Files", sample_id)
    prs_dir = os.path.join(BASE_DIR, "PRS_Files", sample_id)
    ancestry_dir = os.path.join(BASE_DIR, "Ancestry")
    prs_scripts_dir = os.path.join(BASE_DIR, "PRS")
    prs_files_dir = os.path.join(BASE_DIR, "PRS_Files")

    os.makedirs(vcf_dir, exist_ok=True)
    os.makedirs(prs_dir, exist_ok=True)
    os.makedirs(ancestry_dir, exist_ok=True)

    # Resolve script paths
    ett_score_script = os.path.join(prs_files_dir, "ETT_score_2.sh")
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

        # ---- Step 3: PRS scoring ----
        if not os.path.exists(ett_score_script):
            raise FileNotFoundError(f"ETT_score_2.sh not found at {ett_score_script}")
        os.chmod(ett_score_script, 0o755)
        run(["bash", ett_score_script, out_prefix_id, prs_dir])

        # ---- Step 4: PLINK recode A ----
        run([
            plink_bin, "--bfile", out_prefix_id,
            "--recode", "A",
            "--out", os.path.join(vcf_dir, "individual_1"),
            "--allow-extra-chr"
        ])

        # ---- Step 5: individual_conversion.py ----
        if not os.path.exists(conversion_template):
            raise FileNotFoundError(f"individual_conversion.py not found at {conversion_template}")

        # Patch template to use current BIM file
        bim_file = os.path.join(vcf_dir, f"{sample_id}_id.bim")
        with open(conversion_template, "r", encoding="utf-8") as f:
            content = f.read()
        content = content.replace('"HT.bim"', f'"{bim_file}"')  # Replace placeholder with actual BIM

        conversion_dest = os.path.join(vcf_dir, "individual_conversion.py")
        with open(conversion_dest, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"Patched individual_conversion.py to use BIM: {bim_file}")

        # Run conversion script
        run(["python3", "individual_conversion.py"], cwd=vcf_dir)

        local_geno = os.path.join(vcf_dir, "individual_genotypes2.input")
        if not os.path.exists(local_geno):
            raise FileNotFoundError(f"Conversion did not produce {local_geno}")

        # ---- Step 6: Copy genotype input into Ancestry ----
        ancestry_geno = os.path.join(ancestry_dir, "individual_genotypes2.input")
        safe_copy(local_geno, ancestry_geno)

        # ---- Step 7: Run ancestry inside Docker ----
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
            "ancestry-py27",
            "bash", "-lc", docker_cmd
        ])

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

            self.vcf_label = QLabel("VCF Dosyası:")
            layout.addWidget(self.vcf_label)
            self.vcf_path = QLineEdit()
            layout.addWidget(self.vcf_path)
            self.vcf_button = QPushButton("VCF Seç")
            self.vcf_button.clicked.connect(self.select_vcf)
            layout.addWidget(self.vcf_button)

            self.id_input = QLineEdit()
            self.id_input.setPlaceholderText("Kişi ID (örnek: HB0402)")
            layout.addWidget(self.id_input)

            self.sex_input = QLineEdit()
            self.sex_input.setPlaceholderText("Cinsiyet (Kadın/Erkek)")
            layout.addWidget(self.sex_input)

            self.age_input = QLineEdit()
            self.age_input.setPlaceholderText("Yaş")
            layout.addWidget(self.age_input)

            self.pheno_input = QLineEdit()
            self.pheno_input.setPlaceholderText("Fenotip (örnek: diabetes)")
            layout.addWidget(self.pheno_input)

            self.plink_input = self.add_binary_selector(layout, "plink", self.config.get("plink", "plink"))
            self.plink2_input = self.add_binary_selector(layout, "plink2", self.config.get("plink2", "plink2"))
            self.docker_input = self.add_binary_selector(layout, "docker", self.config.get("docker", "docker"))

            self.test_button = QPushButton("Test Binaries")
            self.test_button.clicked.connect(self.test_binaries)
            layout.addWidget(self.test_button)

            self.log = QTextEdit()
            self.log.setReadOnly(True)
            layout.addWidget(self.log)

            self.run_button = QPushButton("Pipeline Başlat")
            self.run_button.clicked.connect(self.run_pipeline)
            layout.addWidget(self.run_button)

            self.setLayout(layout)

        def load_config(self):
            if os.path.exists(CONFIG_FILE):
                try:
                    with open(CONFIG_FILE, "r") as f:
                        return json.load(f)
                except Exception:
                    return {}
            return {}

        def save_config(self):
            data = {
                "plink": self.plink_input.text(),
                "plink2": self.plink2_input.text(),
                "docker": self.docker_input.text()
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

        def log_message(self, msg):
            self.log.append(msg)
            self.log.repaint()

        def test_binaries(self):
            bins = {
                "plink": self.plink_input.text(),
                "plink2": self.plink2_input.text(),
                "docker": self.docker_input.text()
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

            try:
                print_backup = print
                def gui_print(msg):
                    self.log_message(str(msg))
                globals()['print'] = gui_print
                run_pipeline_cli(vcf, sample_id, phenotype, sex, age)
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
    return parser.parse_args()

# ---------- MAIN ----------
if __name__ == "__main__":
    if len(sys.argv) > 1:
        args = parse_args()
        run_pipeline_cli(args.vcf, args.id, args.phenotype, args.sex, args.age)
    else:
        if not _PYQT_AVAILABLE:
            print("PyQt5 not available. Install PyQt5 or run script with CLI args.")
            sys.exit(1)
        qt_app = QApplication(sys.argv)
        window = PRSPipeline()
        window.show()
        sys.exit(qt_app.exec_())