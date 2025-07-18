# -*- coding: utf-8 -*-
import sys
import re
import tempfile
import webbrowser
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QTextEdit, QVBoxLayout,
    QPushButton, QSpinBox, QTableWidget, QTableWidgetItem,
    QMessageBox, QComboBox, QDialog
)
from PyQt5.QtGui import QFont
from Bio.Blast import NCBIWWW, NCBIXML

# 제한효소 사전
RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "XhoI": "CTCGAG",
    "NotI": "GCGGCCGC"
}

def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def calc_gc_content(seq):
    return round((sum(1 for base in seq.upper() if base in "GC") / len(seq)) * 100, 2)

def calc_tm(seq):
    seq = seq.upper()
    N = len(seq)
    if N < 14:
        return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))
    else:
        return round(64.9 + 41 * (seq.count("G") + seq.count("C") - 16.4) / N, 2)

def dimer_score(seq):
    return sum(1 for i in range(len(seq)-3) if reverse_complement(seq[i:i+4]) in seq)

def hairpin_score(seq):
    return sum(1 for i in range(len(seq)-3) if seq[i:i+4] == reverse_complement(seq[i:i+4]))

def insert_re_site(primer, enzyme, direction):
    site = RESTRICTION_ENZYMES.get(enzyme, "")
    return (site + primer) if direction == "5'" else (primer + site)

def design_primers(full_seq, target_seq, f_len, r_len, f_enzyme, r_enzyme):
    full_seq = re.sub(r"\s+", "", full_seq.upper())
    target_seq = re.sub(r"\s+", "", target_seq.upper())

    start = full_seq.find(target_seq)
    if start == -1:
        raise ValueError("❌ 표적 서열이 전체 유전자 서열에 존재하지 않습니다.")
    end = start + len(target_seq)

    if start < f_len or end + r_len >= len(full_seq):
        raise ValueError("❌ 프라이머 길이에 맞게 앞뒤 여유 공간이 부족합니다.")

    f_core = full_seq[start - f_len:start]
    r_core = reverse_complement(full_seq[end:end + r_len])

    if RESTRICTION_ENZYMES[f_enzyme] in f_core or RESTRICTION_ENZYMES[r_enzyme] in r_core:
        raise ValueError("❌ 제한효소 서열이 프라이머에 중복됩니다.")

    f_final = insert_re_site(f_core, f_enzyme, "5'")
    r_final = insert_re_site(r_core, r_enzyme, "5'")

    return [
        ("Forward", f_final, len(f_final), calc_gc_content(f_final), calc_tm(f_final), dimer_score(f_final), hairpin_score(f_final)),
        ("Reverse", r_final, len(r_final), calc_gc_content(r_final), calc_tm(r_final), dimer_score(r_final), hairpin_score(r_final))
    ]

def show_fasta_popup(forward, reverse):
    fasta = f">Forward_Primer\n{forward}\n>Reverse_Primer\n{reverse}"
    dlg = QDialog()
    dlg.setWindowTitle("FASTA Format")
    layout = QVBoxLayout()
    edit = QTextEdit()
    edit.setText(fasta)
    layout.addWidget(edit)
    dlg.setLayout(layout)
    dlg.resize(400, 200)
    dlg.exec_()

def run_blast(seq):
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", seq)
        blast_record = NCBIXML.read(result_handle)
        html = "<html><body><h2>BLAST 결과</h2><ul>"
        for alignment in blast_record.alignments[:5]:
            ident = alignment.hsps[0].identities / alignment.hsps[0].align_length * 100
            color = "green" if ident > 95 else "orange" if ident > 80 else "red"
            html += f"<li style='color:{color}'>{alignment.title} ({ident:.1f}% match)</li>"
        html += "</ul></body></html>"
        path = tempfile.mktemp(suffix=".html")
        with open(path, "w") as f:
            f.write(html)
        webbrowser.open(path)
    except Exception as e:
        print("BLAST 실패:", e)

class PrimerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 PCR 프라이머 설계기")
        self.resize(800, 700)
        layout = QVBoxLayout()

        font = QFont("Courier New", 10)

        layout.addWidget(QLabel("전체 유전자 서열 (5'→3')"))
        self.seq_input = QTextEdit()
        self.seq_input.setFont(font)
        layout.addWidget(self.seq_input)

        layout.addWidget(QLabel("표적 유전자 서열"))
        self.target_input = QTextEdit()
        self.target_input.setFont(font)
        layout.addWidget(self.target_input)

        self.f_len_input = QSpinBox()
        self.f_len_input.setRange(10, 40)
        self.f_len_input.setValue(20)
        self.r_len_input = QSpinBox()
        self.r_len_input.setRange(10, 40)
        self.r_len_input.setValue(20)

        layout.addWidget(QLabel("Forward Primer 길이"))
        layout.addWidget(self.f_len_input)
        layout.addWidget(QLabel("Reverse Primer 길이"))
        layout.addWidget(self.r_len_input)

        self.f_enzyme_box = QComboBox()
        self.f_enzyme_box.addItems(RESTRICTION_ENZYMES.keys())
        layout.addWidget(QLabel("Forward Primer 제한효소"))
        layout.addWidget(self.f_enzyme_box)

        self.r_enzyme_box = QComboBox()
        self.r_enzyme_box.addItems(RESTRICTION_ENZYMES.keys())
        layout.addWidget(QLabel("Reverse Primer 제한효소"))
        layout.addWidget(self.r_enzyme_box)

        self.btn_generate = QPushButton("프라이머 생성")
        self.btn_generate.clicked.connect(self.generate_primers)
        layout.addWidget(self.btn_generate)

        self.btn_fasta = QPushButton("FASTA 보기")
        self.btn_fasta.clicked.connect(self.show_fasta)
        layout.addWidget(self.btn_fasta)

        self.btn_blast = QPushButton("BLAST 실행")
        self.btn_blast.clicked.connect(self.run_blast_align)
        layout.addWidget(self.btn_blast)

        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(["Primer", "Sequence", "Length", "GC%", "Tm", "DimerScore", "HairpinScore"])
        layout.addWidget(self.table)

        self.setLayout(layout)

    def generate_primers(self):
        try:
            full_seq = self.seq_input.toPlainText()
            target_seq = self.target_input.toPlainText()
            f_len = self.f_len_input.value()
            r_len = self.r_len_input.value()
            f_enzyme = self.f_enzyme_box.currentText()
            r_enzyme = self.r_enzyme_box.currentText()

            results = design_primers(full_seq, target_seq, f_len, r_len, f_enzyme, r_enzyme)
            self.results = results  # 저장

            self.table.setRowCount(len(results))
            for i, row in enumerate(results):
                for j, val in enumerate(row):
                    self.table.setItem(i, j, QTableWidgetItem(str(val)))

        except Exception as e:
            QMessageBox.warning(self, "에러", str(e))

    def show_fasta(self):
        try:
            f = self.results[0][1]
            r = self.results[1][1]
            show_fasta_popup(f, r)
        except:
            QMessageBox.information(self, "오류", "먼저 프라이머를 생성하세요")

    def run_blast_align(self):
        try:
            f = self.results[0][1]
            run_blast(f)
        except:
            QMessageBox.information(self, "오류", "먼저 프라이머를 생성하세요")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PrimerApp()
    window.show()
    sys.exit(app.exec_())
