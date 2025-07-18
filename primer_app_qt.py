# ===== Biopython Import 시도 =====
try:
    from Bio.Blast import NCBIWWW, NCBIXML
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

import sys
import re
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QTextEdit, QVBoxLayout,
    QPushButton, QSpinBox, QTableWidget, QTableWidgetItem,
    QMessageBox, QComboBox, QHBoxLayout
)
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt

# ===== 제한효소 데이터 =====
RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "XhoI": "CTCGAG",
    "NotI": "GCGGCCGC",
    "HindIII": "AAGCTT"
}

# ===== 유틸리티 함수 =====
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def calc_gc_content(seq):
    gc_count = sum(1 for base in seq.upper() if base in "GC")
    return round((gc_count / len(seq)) * 100, 2)

def calc_tm(seq):
    seq = seq.upper()
    N = len(seq)
    if N < 14:
        return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))
    else:
        return round(64.9 + 41 * (seq.count("G") + seq.count("C") - 16.4) / N, 2)

def check_gc_clamp(seq):
    return 1 <= sum(1 for base in seq[-5:] if base in "GC") <= 2

def check_repeat(seq):
    return bool(re.search(r"(AT){3,}|(TA){3,}|A{5,}|T{5,}|G{5,}|C{5,}", seq))

def has_self_complementarity(primer):
    rev_comp = reverse_complement(primer)
    return any(rev_comp[i:i+4] in primer for i in range(len(primer) - 3))

def score_self_complementarity(primer):
    count = 0
    rev_comp = reverse_complement(primer)
    for i in range(len(primer) - 3):
        window = rev_comp[i:i+4]
        if window in primer:
            count += 1
    return count

def is_primer_unique(seq, primer):
    return seq.count(primer) == 1

# ===== 프라이머 설계 함수 =====
def design_primers(full_seq, target_seq, f_len, r_len, enzyme, insert_side):
    full_seq = re.sub(r"\s+", "", full_seq.upper())
    target_seq = re.sub(r"\s+", "", target_seq.upper())

    start = full_seq.find(target_seq)
    if start == -1:
        raise ValueError("❌ 표적 서열이 전체 유전자 서열에 없습니다.")

    end = start + len(target_seq)
    if start < f_len or end + r_len > len(full_seq):
        raise ValueError("❌ 프라이머 설계에 필요한 여유 염기 수 부족")

    f_primer = full_seq[start - f_len:start]
    r_primer_raw = full_seq[end:end + r_len]
    r_primer = reverse_complement(r_primer_raw)

    enzyme_seq = RESTRICTION_ENZYMES.get(enzyme, "")

    if insert_side == "5'":
        f_primer = enzyme_seq + f_primer
        r_primer = enzyme_seq + r_primer
    else:
        f_primer = f_primer + enzyme_seq
        r_primer = r_primer + enzyme_seq

    if enzyme_seq in f_primer[f_len:] or enzyme_seq in r_primer[r_len:]:
        raise ValueError("⚠️ 제한효소 서열이 프라이머 내부에 중복 포함되어 있습니다.")

    return [
        ("Forward", f_primer, len(f_primer), calc_gc_content(f_primer), calc_tm(f_primer),
         check_gc_clamp(f_primer), check_repeat(f_primer), score_self_complementarity(f_primer)),

        ("Reverse", r_primer, len(r_primer), calc_gc_content(r_primer), calc_tm(r_primer),
         check_gc_clamp(r_primer), check_repeat(r_primer), score_self_complementarity(r_primer)),
    ]

# ===== UI 클래스 =====
class PrimerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 전문가용 PCR 프라이머 설계기")
        self.resize(800, 700)

        layout = QVBoxLayout()

        self.seq_input = QTextEdit()
        self.seq_input.setPlaceholderText("전체 유전자 서열 (5'→3')")
        layout.addWidget(self.seq_input)

        self.target_input = QTextEdit()
        self.target_input.setPlaceholderText("표적 DNA 서열")
        layout.addWidget(self.target_input)

        f_layout = QHBoxLayout()
        self.f_len = QSpinBox(); self.f_len.setRange(10, 40); self.f_len.setValue(20)
        self.r_len = QSpinBox(); self.r_len.setRange(10, 40); self.r_len.setValue(20)
        f_layout.addWidget(QLabel("Forward 길이:")); f_layout.addWidget(self.f_len)
        f_layout.addWidget(QLabel("Reverse 길이:")); f_layout.addWidget(self.r_len)
        layout.addLayout(f_layout)

        self.enzyme_box = QComboBox()
        self.enzyme_box.addItems(RESTRICTION_ENZYMES.keys())
        layout.addWidget(QLabel("제한효소 선택:"))
        layout.addWidget(self.enzyme_box)

        self.side_box = QComboBox()
        self.side_box.addItems(["5'", "3'"])
        layout.addWidget(QLabel("제한효소 삽입 위치:"))
        layout.addWidget(self.side_box)

        self.button = QPushButton("프라이머 생성")
        self.button.clicked.connect(self.generate_primers)
        layout.addWidget(self.button)

        self.blast_button = QPushButton("NCBI BLAST 실행")
        self.blast_button.clicked.connect(self.run_blast)
        if not BIOPYTHON_AVAILABLE:
            self.blast_button.setEnabled(False)
            self.blast_button.setToolTip("Biopython 미설치됨 – BLAST 사용 불가")
        layout.addWidget(self.blast_button)

        self.table = QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels([
            "Primer", "Sequence", "Length", "GC_Content", "Tm", "GC_Clamp", "Repeats", "DimerScore"
        ])
        layout.addWidget(self.table)

        self.setLayout(layout)

    def generate_primers(self):
        try:
            full_seq = self.seq_input.toPlainText()
            target_seq = self.target_input.toPlainText()
            f_len = self.f_len.value()
            r_len = self.r_len.value()
            enzyme = self.enzyme_box.currentText()
            if enzyme != "없음" and enzyme in enzyme_dict:
                insert_side = self.side_box.currentText()

            result = design_primers(full_seq, target_seq, f_len, r_len, enzyme, insert_side)
            self.table.setRowCount(len(result))
            for i, row in enumerate(result):
                for j, val in enumerate(row):
                    self.table.setItem(i, j, QTableWidgetItem(str(val)))
        except Exception as e:
            QMessageBox.warning(self, "에러", str(e))

    def run_blast(self):
        if not BIOPYTHON_AVAILABLE:
            QMessageBox.warning(self, "BLAST 불가", "Biopython이 설치되지 않아 BLAST 실행이 불가능합니다.")
            return

        try:
            primer_seq = self.table.item(0, 1).text()
            result_handle = NCBIWWW.qblast("blastn", "nt", primer_seq)
            blast_record = NCBIXML.read(result_handle)

            top_hits = []
            for alignment in blast_record.alignments[:5]:
                top_hits.append(f"{alignment.title}\nScore: {alignment.hsps[0].score}")

            QMessageBox.information(self, "BLAST 결과 (Top 5)", "\n\n".join(top_hits))
        except Exception as e:
            QMessageBox.critical(self, "BLAST 오류", f"오류 발생: {e}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PrimerApp()
    window.show()
    sys.exit(app.exec_())
