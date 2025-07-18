import sys
import re
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QTextEdit, QVBoxLayout,
    QPushButton, QSpinBox, QComboBox, QTableWidget, QTableWidgetItem, QMessageBox
)
from PyQt5.QtGui import QFont

# ===== 유틸리티 함수 =====
def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def calc_gc_content(seq):
    gc_count = sum(1 for base in seq.upper() if base in "GC")
    return round((gc_count / len(seq)) * 100, 2)

def calc_tm(seq):
    seq = seq.upper()
    return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C")) - 5

def check_gc_clamp(seq):
    return 1 <= sum(1 for base in seq[-5:] if base in "GC") <= 2

def check_repeat(seq):
    return bool(re.search(r"(AT){3,}|(TA){3,}|A{5,}|T{5,}|G{5,}|C{5,}", seq))

def is_primer_unique(seq, primer):
    return seq.count(primer) == 1

def has_self_complementarity(primer):
    rev_comp = reverse_complement(primer)
    return any(rev_comp[i:i+4] in primer for i in range(len(primer) - 3))

# 제한효소 사전
RE_SITES = {
    "(없음)": "",
    "EcoRI (GAATTC)": "GAATTC",
    "XhoI (CTCGAG)": "CTCGAG",
    "BamHI (GGATCC)": "GGATCC",
    "HindIII (AAGCTT)": "AAGCTT"
}

# ===== 프라이머 설계 함수 =====
def design_primers(full_seq, target_seq, primer_length=20, enzyme_f="", enzyme_r=""):
    full_seq = re.sub(r"\s+", "", full_seq.upper())
    target_seq = re.sub(r"\s+", "", target_seq.upper())

    start_pos = full_seq.find(target_seq)
    if start_pos == -1:
        raise ValueError("❌ 표적 서열이 전체 유전자 서열에 존재하지 않습니다.")

    end_pos = start_pos + len(target_seq) - 1

    if start_pos < primer_length or end_pos + primer_length >= len(full_seq):
        raise ValueError("❌ 프라이머 길이에 맞게 타겟 앞뒤 여유 공간이 부족합니다.")

    f_core = full_seq[start_pos - primer_length : start_pos]
    r_core_raw = full_seq[end_pos + 1 : end_pos + 1 + primer_length]
    r_core = reverse_complement(r_core_raw)

    # 제한효소 서열 추가
    f_primer = enzyme_f + f_core if enzyme_f else f_core
    r_primer = enzyme_r + r_core if enzyme_r else r_core

    tm_f = calc_tm(f_core)
    tm_r = calc_tm(r_core)

    warning_messages = []

    if abs(tm_f - tm_r) > 5:
        warning_messages.append("⚠️ Tm 차이가 큽니다 (권장 ≤ 5℃)")

    if not is_primer_unique(full_seq, f_core):
        warning_messages.append("⚠️ Forward 프라이머가 유전자 내 여러 위치에 존재합니다.")

    if not is_primer_unique(full_seq, reverse_complement(r_core)):
        warning_messages.append("⚠️ Reverse 프라이머가 유전자 내 여러 위치에 존재합니다.")

    if has_self_complementarity(f_core):
        warning_messages.append("⚠️ Forward 프라이머에 self-complementary 구조 가능성")

    if has_self_complementarity(r_core):
        warning_messages.append("⚠️ Reverse 프라이머에 self-complementary 구조 가능성")

    amplicon_length = end_pos + primer_length - (start_pos - primer_length) + 1

    result = [
        ("Forward", f_primer, len(f_primer), calc_gc_content(f_core),
         tm_f, check_gc_clamp(f_core), check_repeat(f_core), "Yes" if is_primer_unique(full_seq, f_core) else "No"),

        ("Reverse", r_primer, len(r_primer), calc_gc_content(r_core),
         tm_r, check_gc_clamp(r_core), check_repeat(r_core), "Yes" if is_primer_unique(full_seq, reverse_complement(r_core)) else "No"),

        ("Amplicon Length", amplicon_length, "", "", "", "", "", "")
    ]

    return result, warning_messages

# ===== UI 앱 클래스 =====
class PrimerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 전문가용 PCR 프라이머 설계기")
        self.resize(750, 700)

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
        self.len_input.setRange(10, 40)
        self.len_input.setValue(20)
        layout.addWidget(QLabel("프라이머 길이 (10~40):"))
        layout.addWidget(self.len_input)

        self.enzyme_f_input = QComboBox()
        self.enzyme_f_input.addItems(RE_SITES.keys())
        layout.addWidget(QLabel("Forward 제한효소 (선택):"))
        layout.addWidget(self.enzyme_f_input)

        self.enzyme_r_input = QComboBox()
        self.enzyme_r_input.addItems(RE_SITES.keys())
        layout.addWidget(QLabel("Reverse 제한효소 (선택):"))
        layout.addWidget(self.enzyme_r_input)

        self.button = QPushButton("프라이머 생성")
        self.button.clicked.connect(self.generate_primers)
        layout.addWidget(self.button)

        self.table = QTableWidget(0, 8)
        self.table.setHorizontalHeaderLabels(
            ["Primer", "Sequence", "Length", "GC_Content", "Tm", "GC_Clamp", "Repeats", "Unique"]
        )
        layout.addWidget(self.table)

        self.setLayout(layout)

    def generate_primers(self):
        full_seq = self.seq_input.toPlainText()
        target_seq = self.target_input.toPlainText()
        primer_len = self.len_input.value()
        enzyme_f = RE_SITES[self.enzyme_f_input.currentText()]
        enzyme_r = RE_SITES[self.enzyme_r_input.currentText()]

        try:
            results, warnings = design_primers(full_seq, target_seq, primer_len, enzyme_f, enzyme_r)
            self.table.setRowCount(len(results))
            for i, row in enumerate(results):
                for j, val in enumerate(row):
                    self.table.setItem(i, j, QTableWidgetItem(str(val)))
            if warnings:
                QMessageBox.warning(self, "주의 사항", "\n".join(warnings))
        except Exception as e:
            QMessageBox.critical(self, "에러", str(e))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PrimerApp()
    window.show()
    sys.exit(app.exec_())
