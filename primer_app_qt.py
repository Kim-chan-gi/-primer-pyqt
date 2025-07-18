import sys
import re
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QTextEdit, QVBoxLayout,
    QPushButton, QSpinBox, QTableWidget, QTableWidgetItem
)

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

# ===== 프라이머 설계 함수 =====
def design_primers(full_seq, target_seq, primer_length=20):
    full_seq = re.sub(r"\s+", "", full_seq.upper())
    target_seq = re.sub(r"\s+", "", target_seq.upper())

    start_pos = full_seq.find(target_seq)
    if start_pos == -1:
        raise ValueError("❌ 표적 서열이 전체 유전자 서열에 존재하지 않습니다.")

    end_pos = start_pos + len(target_seq) - 1

    # 양쪽 프라이머 길이가 확보 가능한지 검사
    if start_pos + primer_length > len(full_seq) or end_pos - primer_length + 1 < 0:
        raise ValueError("❌ 프라이머 길이에 맞게 충분한 여유가 없습니다.")

    # Forward Primer: 타겟의 시작 위치부터
    f_primer = full_seq[start_pos : start_pos + primer_length]

    # Reverse Primer: 타겟의 끝 포함해 그 앞에서 primer_length 추출 → 역상보
    r_primer_raw = full_seq[end_pos - primer_length + 1 : end_pos + 1]
    r_primer = reverse_complement(r_primer_raw)

    return [
        ("Forward", f_primer, len(f_primer), calc_gc_content(f_primer),
         calc_tm(f_primer), check_gc_clamp(f_primer), check_repeat(f_primer)),
        ("Reverse", r_primer, len(r_primer), calc_gc_content(r_primer),
         calc_tm(r_primer), check_gc_clamp(r_primer), check_repeat(r_primer)),
    ]

# ===== UI 앱 클래스 =====
class PrimerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬 PCR 프라이머 설계기")
        self.resize(700, 600)

        layout = QVBoxLayout()

        # 전체 유전자 서열 입력
        layout.addWidget(QLabel("전체 유전자 서열 (5'→3')"))
        self.seq_input = QTextEdit()
        layout.addWidget(self.seq_input)

        # 타겟 서열 입력
        layout.addWidget(QLabel("표적 DNA 서열"))
        self.target_input = QTextEdit()
        layout.addWidget(self.target_input)

        # 프라이머 길이 입력
        self.len_input = QSpinBox()
        self.len_input.setRange(10, 40)
        self.len_input.setValue(20)
        layout.addWidget(QLabel("프라이머 길이 (10~40):"))
        layout.addWidget(self.len_input)

        # 버튼
        self.button = QPushButton("프라이머 생성")
        self.button.clicked.connect(self.generate_primers)
        layout.addWidget(self.button)

        # 출력 테이블
        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(
            ["Primer", "Sequence", "Length", "GC_Content", "Tm", "GC_Clamp", "Repeats"]
        )
        layout.addWidget(self.table)

        self.setLayout(layout)

    def generate_primers(self):
        full_seq = self.seq_input.toPlainText()
        target_seq = self.target_input.toPlainText()
        primer_len = self.len_input.value()

        try:
            results = design_primers(full_seq, target_seq, primer_len)
            self.table.setRowCount(len(results))
            for i, row in enumerate(results):
                for j, val in enumerate(row):
                    self.table.setItem(i, j, QTableWidgetItem(str(val)))
        except Exception as e:
            self.table.setRowCount(1)
            self.table.setItem(0, 0, QTableWidgetItem("⚠️ " + str(e)))

# ===== 실행 진입점 =====
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PrimerApp()
    window.show()
    sys.exit(app.exec_())
