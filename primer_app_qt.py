import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QTextEdit, QVBoxLayout,
    QPushButton, QSpinBox, QTableWidget, QTableWidgetItem
)
import re

def reverse_complement(seq):
    comp_table = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp_table)[::-1]

def calc_gc_content(seq):
    seq = seq.upper()
    gc_count = sum(1 for base in seq if base in "GC")
    return round((gc_count / len(seq)) * 100, 2)

def calc_tm(seq):
    seq = seq.upper()
    A = seq.count("A")
    T = seq.count("T")
    G = seq.count("G")
    C = seq.count("C")
    return 2 * (A + T) + 4 * (G + C) - 5

def check_gc_clamp(seq):
    tail = seq[-5:]
    gc_count = sum(1 for base in tail if base in "GC")
    return 1 <= gc_count <= 2

def check_repeat(seq):
    return bool(re.search(r"(AT){3,}|(TA){3,}|A{5,}|T{5,}|G{5,}|C{5,}", seq))

def design_primers(full_seq, target_seq, primer_length=20):
    full_seq = re.sub(r"\s+", "", full_seq.upper())
    target_seq = re.sub(r"\s+", "", target_seq.upper())

    start_pos = full_seq.find(target_seq)
    if start_pos == -1:
        raise ValueError("❌ 표적 서열을 전체 유전자에서 찾을 수 없습니다.")
    
    end_pos = start_pos + len(target_seq) - 1

    if start_pos < primer_length or end_pos + primer_length > len(full_seq):
        raise ValueError("❌ 프라이머 길이에 맞는 위치가 충분하지 않습니다.")

    f_primer = full_seq[start_pos - primer_length:start_pos]
    r_primer_raw = full_seq[end_pos + 1:end_pos + 1 + primer_length]
    r_primer = reverse_complement(r_primer_raw)

    return [
        ("Forward", f_primer, len(f_primer), calc_gc_content(f_primer), calc_tm(f_primer), check_gc_clamp(f_primer), check_repeat(f_primer)),
        ("Reverse", r_primer, len(r_primer), calc_gc_content(r_primer), calc_tm(r_primer), check_gc_clamp(r_primer), check_repeat(r_primer)),
    ]

class PrimerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 PCR 프라이머 만들기...?")
        self.resize(700, 600)
        layout = QVBoxLayout()

        self.seq_label = QLabel("전체 유전자 서열 (5'→3')")
        self.seq_input = QTextEdit()
        layout.addWidget(self.seq_label)
        layout.addWidget(self.seq_input)

        self.target_label = QLabel("표적 DNA 서열")
        self.target_input = QTextEdit()
        layout.addWidget(self.target_label)
        layout.addWidget(self.target_input)

        self.len_label = QLabel("프라이머 길이 (10~40):")
        self.len_input = QSpinBox()
        self.len_input.setMinimum(10)
        self.len_input.setMaximum(40)
        self.len_input.setValue(20)
        layout.addWidget(self.len_label)
        layout.addWidget(self.len_input)

        self.button = QPushButton("프라이머 생성")
        self.button.clicked.connect(self.generate_primers)
        layout.addWidget(self.button)

        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(["Primer", "Sequence", "Length", "GC_Content", "Tm", "GC_Clamp", "Repeats"])
        layout.addWidget(self.table)

        self.setLayout(layout)

    def generate_primers(self):
        full_seq = self.seq_input.toPlainText()
        target_seq = self.target_input.toPlainText()
        primer_len = self.len_input.value()

        try:
            result = design_primers(full_seq, target_seq, primer_len)
            self.table.setRowCount(len(result))
            for i, row in enumerate(result):
                for j, val in enumerate(row):
                    self.table.setItem(i, j, QTableWidgetItem(str(val)))
        except Exception as e:
            self.table.setRowCount(1)
            self.table.setItem(0, 0, QTableWidgetItem(str(e)))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PrimerApp()
    window.show()
    sys.exit(app.exec_())
