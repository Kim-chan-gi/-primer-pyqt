import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QTextEdit, QVBoxLayout, QPushButton, QSpinBox, QTableWidget, QTableWidgetItem
)
import re

def reverse_complement(seq):
    comp_table = str.maketrans("ACGT", "TGCA")
    comp_seq = seq.translate(comp_table)
    return comp_seq[::-1]

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

    data = {
        "Primer": ["Forward", "Reverse"],
        "Sequence": [f_primer, r_primer],
        "Length": [len(f_primer), len(r_primer)],
        "GC_Content": [calc_gc_content(f_primer), calc_gc_content(r_primer)],
        "Tm": [calc_tm(f_primer), calc_tm(r_primer)],
        "GC_Clamp": [check_gc_clamp(f_primer), check_gc_clamp(r_primer)],
        "Repeats": [check_repeat(f_primer), check_repeat(r_primer)]
    }
    return pd.DataFrame(data)

st.title("🧬 PCR 프라이머 설계기")

full_seq = st.text_area("전체 유전자 서열 (5'→3')", height=100)
target_seq = st.text_area("표적 DNA 서열", height=100)
primer_length = st.slider("프라이머 길이", min_value=10, max_value=40, value=20)
run_button = st.button("프라이머 생성")

if run_button:
    try:
        df = design_primers(full_seq, target_seq, primer_length)
        st.dataframe(df)

        df_filtered = df[
            (df["GC_Content"] >= 40) & (df["GC_Content"] <= 60) &
            (df["Tm"] >= 52) & (df["Tm"] <= 65) &
            (df["GC_Clamp"]) & (~df["Repeats"])
        ]
        if len(df_filtered) < 2:
            st.warning("⚠️ 프라이머 조건이 적합하지 않거나 1개 이상 조건을 충족하지 않습니다.")
        else:
            st.success("✅ 모든 생물학적 조건을 충족하는 프라이머입니다.")
    except Exception as e:
        st.error(f"❌ 오류 발생: {e}")
