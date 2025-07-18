import sys
import re
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QTextEdit, QVBoxLayout,
    QPushButton, QSpinBox, QComboBox, QTableWidget, QTableWidgetItem, QMessageBox, QLineEdit
)
from PyQt5.QtGui import QFont
import math

# ===== 유틸리티 함수 =====
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def calc_gc_content(seq):
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return round((gc_count / len(seq)) * 100, 2)

def calc_tm(seq):
    seq = seq.upper()
    n = len(seq)
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    if n <= 14:
        return 4 * gc + 2 * at
    else:
        return round(64.9 + 41 * (gc - 16.4) / n, 2)

def check_gc_clamp(seq):
    return seq[-1] in "GC"

def check_repeat(seq):
    return bool(re.search(r"(.)\1{4,}", seq))

def has_self_complementarity(seq):
    rev = reverse_complement(seq)
    for i in range(4, 7):
        if rev[:i] in seq or rev[-i:] in seq:
            return True
    return False

def is_valid_gc_content(seq):
    gc = calc_gc_content(seq)
    return 40 <= gc <= 60

# ===== 프라이머 설계 함수 =====
def design_strict_primers(full_seq, target_seq, primer_length=20):
    full_seq = re.sub(r"\s+", "", full_seq.upper())
    target_seq = re.sub(r"\s+", "", target_seq.upper())

    start = full_seq.find(target_seq)
    if start == -1:
        raise ValueError("❌ 전체 서열에서 타겟 서열을 찾을 수 없습니다.")
    end = start + len(target_seq)

    if start < 0 or end + primer_length > len(full_seq):
        raise ValueError("❌ 프라이머를 설계하기에 충분한 영역이 없습니다.")

    f_primer = full_seq[start : start + primer_length]
    r_raw = full_seq[end - primer_length : end]
    r_primer = reverse_complement(r_raw)

    f_tm = calc_tm(f_primer)
    r_tm = calc_tm(r_primer)

    f_gc = calc_gc_content(f_primer)
    r_gc = calc_gc_content(r_primer)

    f_warn, r_warn = [], []
    if not is_valid_gc_content(f_primer): f_warn.append("⚠ GC% 비정상")
    if not is_valid_gc_content(r_primer): r_warn.append("⚠ GC% 비정상")
    if not check_gc_clamp(f_primer): f_warn.append("⚠ 3' 끝에 G/C 없음")
    if not check_gc_clamp(r_primer): r_warn.append("⚠ 3' 끝에 G/C 없음")
    if check_repeat(f_primer): f_warn.append("⚠ 반복염기")
    if check_repeat(r_primer): r_warn.append("⚠ 반복염기")
    if has_self_complementarity(f_primer): f_warn.append("⚠ 자기 상보성")
    if has_self_complementarity(r_primer): r_warn.append("⚠ 자기 상보성")

    return [
        ("Forward", f_primer, len(f_primer), f_gc, f_tm, ", ".join(f_warn) or "OK"),
        ("Reverse", r_primer, len(r_primer), r_gc, r_tm, ", ".join(r_warn) or "OK")
    ]

# ===== UI 앱 =====
class PrimerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 고정밀 PCR 프라이머 설계기")
        self.resize(750, 600)
        layout = QVBoxLayout()
        font = QFont("Courier New", 10)

        layout.addWidget(QLabel("전체 유전자 서열 (5'→3')"))
        self.seq_input = QTextEdit()
        self.seq_input.setFont(font)
        layout.addWidget(self.seq_input)

        layout.addWidget(QLabel("표적 DNA 서열"))
        self.target_input = QTextEdit()
        self.target_input.setFont(font)
        layout.addWidget(self.target_input)

        self.len_input = QSpinBox()
        self.len_input.setRange(15, 40)
        self.len_input.setValue(20)
        layout.addWidget(QLabel("프라이머 길이 (15~40):"))
        layout.addWidget(self.len_input)

        self.button = QPushButton("프라이머 생성")
        self.button.clicked.connect(self.generate_primers)
        layout.addWidget(self.button)

        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["Primer", "Sequence", "Length", "GC%", "Tm", "Warnings"])
        layout.addWidget(self.table)

        self.setLayout(layout)

    def generate_primers(self):
        seq = self.seq_input.toPlainText()
        target = self.target_input.toPlainText()
        length = self.len_input.value()

        try:
            results = design_strict_primers(seq, target, length)
            self.table.setRowCount(len(results))
            for i, row in enumerate(results):
                for j, val in enumerate(row):
                    self.table.setItem(i, j, QTableWidgetItem(str(val)))
        except Exception as e:
            QMessageBox.critical(self, "에러", str(e))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PrimerApp()
    window.show()
    sys.exit(app.exec_())
