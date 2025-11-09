import io
from typing import List, Dict, Tuple
import pandas as pd
import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.pyplot as plt

# ------------------------------
# 유틸 함수
# ------------------------------
def clean_sequence(seq: str) -> str:
    """대문자화 + ATGC만 남김"""
    seq = (seq or "").upper()
    return "".join([b for b in seq if b in ["A", "T", "G", "C"]])

def count_bases(seq: str) -> Dict[str, int]:
    return {b: seq.count(b) for b in ["A", "T", "G", "C"]}

def base_ratios(counts: Dict[str, int]) -> Dict[str, float]:
    total = sum(counts.values())
    return {b: (counts[b] / total) if total else 0.0 for b in counts}

def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    g = seq.count("G")
    c = seq.count("C")
    return (g + c) / len(seq) * 100

def read_fasta_bytes(file_bytes) -> str:
    """업로더에서 받은 바이너리 FASTA를 첫 레코드로 읽어 문자열 반환"""
    try:
        handle = io.StringIO(file_bytes.getvalue().decode("utf-8"))
        record = next(SeqIO.parse(handle, "fasta"))
        return str(record.seq)
    except Exception:
        return ""

def find_mutations(query: str, ref: str) -> List[Tuple[int, str, str]]:
    """위치, ref, query를 담은 변이 리스트 반환(0-index)"""
    muts = []
    L = min(len(query), len(ref))
    for i in range(L):
        if query[i] != ref[i]:
            muts.append((i, ref[i], query[i]))
    return muts

def kmers_count(seq: str, k: int) -> Dict[str, int]:
    counts = {}
    if k <= 0 or len(seq) < k:
        return counts
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        counts[kmer] = counts.get(kmer, 0) + 1
    return counts

def find_orfs(seq: str, min_nt: int = 90) -> List[Dict]:
    """ATG에서 시작해 TAA/TAG/TGA에서 끝나는 ORF 탐색(3 프레임, + strand)"""
    stops = {"TAA", "TAG", "TGA"}
    results = []
    for frame in range(3):
        i = frame
        while i + 3 <= len(seq):
            codon = seq[i:i+3]
            if codon == "ATG":
                j = i + 3
                while j + 3 <= len(seq):
                    stop = seq[j:j+3]
                    if stop in stops:
                        length = j + 3 - i
                        if length >= min_nt:
                            prot = Seq(seq[i:j+3]).translate(to_stop=True)
                            results.append({
                                "frame": f"+{frame}",
                                "start": i,
                                "end": j+3,
                                "length_nt": length,
                                "length_aa": len(prot),
                                "protein_preview": str(prot)[:20] + ("..." if len(prot) > 20 else "")
                            })
                        i = j + 3
                        break
                    j += 3
                else:
                    i += 3
            else:
                i += 3
    return results

# ------------------------------
# UI
# ------------------------------
st.set_page_config(page_title="DNA Analyzer (BioPython)", layout="wide")
st.title("🧬 DNA Analyzer (BioPython) — 심화 웹앱")

with st.sidebar:
    st.header("설정")
    min_orf_len = st.number_input("최소 ORF 길이 (nt)", min_value=60, max_value=3000, value=90, step=30)
    k_for_kmer = st.slider("k-mer 크기", min_value=2, max_value=3, value=2)
    show_plots = st.checkbox("그래프 표시", value=True)
    st.caption("참고: 업로드가 없으면 텍스트 입력을 사용함.")

# 입력 섹션
st.subheader("1) 서열 입력/업로드")
col1, col2 = st.columns(2)

with col1:
    seq_text = st.text_area(
        "분석할 DNA 서열(텍스트 붙여넣기 또는 FASTA 업로드 사용)",
        height=150,
        placeholder="ATGCGTACGTT... (A/T/G/C 외 문자는 자동 제거됨)"
    )
    upload_seq = st.file_uploader("또는 FASTA 업로드", type=["fa", "fasta"], key="query")

with col2:
    ref_text = st.text_area(
        "기준(reference) 서열(선택)",
        height=150,
        placeholder="비교 기준 서열이 있으면 입력 또는 FASTA 업로드"
    )
    upload_ref = st.file_uploader("또는 기준 FASTA 업로드", type=["fa", "fasta"], key="ref")

# 데이터 준비
raw_seq = read_fasta_bytes(upload_seq) if upload_seq else seq_text
raw_ref = read_fasta_bytes(upload_ref) if upload_ref else ref_text

seq = clean_sequence(raw_seq)
ref = clean_sequence(raw_ref)

if not seq:
    st.info("분석할 서열을 입력하거나 업로드해 주세요.")
    st.stop()

st.success(f"입력 서열 길이: {len(seq)} nt" + (f" · 기준 서열 길이: {len(ref)} nt" if ref else ""))

# ------------------------------
# 분석 1: 빈도/비율/GC
# ------------------------------
st.subheader("2) 염기 조성 · 비율 · GC%")
counts = count_bases(seq)
ratios = base_ratios(counts)
gc = gc_content(seq)

left, right = st.columns([1, 1])
with left:
    df_comp = pd.DataFrame({
        "base": ["A", "T", "G", "C"],
        "count": [counts["A"], counts["T"], counts["G"], counts["C"]],
        "ratio(%)": [round(ratios["A"]*100, 2),
                     round(ratios["T"]*100, 2),
                     round(ratios["G"]*100, 2),
                     round(ratios["C"]*100, 2)]
    })
    st.dataframe(df_comp, use_container_width=True)

with right:
    st.metric(label="GC Content", value=f"{gc:.2f}%")

if show_plots:
    fig = plt.figure()
    plt.bar(df_comp["base"], df_comp["count"])
    plt.title("Base Counts")
    plt.xlabel("Base")
    plt.ylabel("Count")
    st.pyplot(fig)

# ------------------------------
# 분석 2: 변이 검출
# ------------------------------
st.subheader("3) 변이 검출(선택)")
if ref:
    muts = find_mutations(seq, ref)
    st.write(f"총 변이 수: **{len(muts)}**")
    if muts:
        df_muts = pd.DataFrame(muts, columns=["pos(0-based)", "ref", "query"])
        st.dataframe(df_muts, use_container_width=True, height=240)
        csv = df_muts.to_csv(index=False).encode("utf-8")
        st.download_button("변이 목록 CSV 다운로드", csv, file_name="mutations.csv", mime="text/csv")
else:
    st.caption("기준 서열이 없으므로 변이 검출을 건너뜀.")

# ------------------------------
# 분석 3: ORF 탐색
# ------------------------------
st.subheader("4) ORF 탐색(+ strand)")
orfs = find_orfs(seq, min_nt=int(min_orf_len))
if orfs:
    df_orf = pd.DataFrame(orfs)
    st.dataframe(df_orf, use_container_width=True, height=260)
    csv_orf = df_orf.to_csv(index=False).encode("utf-8")
    st.download_button("ORF 목록 CSV 다운로드", csv_orf, file_name="orfs.csv", mime="text/csv")
else:
    st.write("조건을 만족하는 ORF가 없음.")

# ------------------------------
# 분석 4: k-mer 빈도
# ------------------------------
st.subheader("5) k-mer 빈도 분석")
kcounts = kmers_count(seq, k_for_kmer)
if kcounts:
    df_k = pd.DataFrame(sorted(kcounts.items(), key=lambda x: (-x[1], x[0])),
                        columns=[f"{k_for_kmer}-mer", "count"])
    st.dataframe(df_k, use_container_width=True, height=260)
    if show_plots:
        fig2 = plt.figure()
        plt.bar(df_k.iloc[:20, 0], df_k.iloc[:20, 1])
        plt.title(f"Top {min(20, len(df_k))} {k_for_kmer}-mers")
        plt.xlabel("k-mer")
        plt.ylabel("Count")
        plt.xticks(rotation=90)
        st.pyplot(fig2)
else:
    st.caption("서열 길이가 k보다 짧아 k-mer 분석 불가.")

# ------------------------------
# 시퀀스 보기 & 리버스컴플리먼트
# ------------------------------
st.subheader("6) 시퀀스 도구")
with st.expander("리버스-컴플리먼트 보기"):
    rc = str(Seq(seq).reverse_complement())
    st.code(rc, language="text")
    st.download_button("RC 서열 다운로드", rc, file_name="reverse_complement.txt", mime="text/plain")

st.success("분석 완료.")
