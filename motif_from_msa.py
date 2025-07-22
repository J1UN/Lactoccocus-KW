import argparse
from collections import OrderedDict

# Step 1: parse MSA
def parse_clustal_msa(path):
    seqs = OrderedDict()
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("CLUSTAL") or line[0] in "*: .":
                continue
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            seq_id, seq = parts
            if seq_id not in seqs:
                seqs[seq_id] = ""
            seqs[seq_id] += seq
    return seqs

# Step 2: calculate conserved positions
def find_conserved_positions(seqs, threshold=0.8, min_block=10):
    seq_ids = list(seqs.keys())
    n = len(seq_ids)
    aln_len = len(next(iter(seqs.values())))
    conserved = []

    for i in range(aln_len):
        chars = [seqs[seq_id][i] for seq_id in seq_ids]
        if '-' in chars:
            continue
        major = max(set(chars), key=chars.count)
        count = chars.count(major)
        if count / n >= threshold:
            conserved.append(i)

    # group into blocks
    blocks = []
    temp = []
    for i in conserved:
        if not temp or i == temp[-1] + 1:
            temp.append(i)
        else:
            if len(temp) >= min_block:
                blocks.append((temp[0], temp[-1] + 1))
            temp = [i]
    if len(temp) >= min_block:
        blocks.append((temp[0], temp[-1] + 1))
    return blocks

# Step 3: extract
def extract_conserved(seqs, blocks):
    trimmed = OrderedDict()
    for sid, seq in seqs.items():
        trimmed_seq = ''.join(seq[start:end] for start, end in blocks)
        trimmed[sid] = trimmed_seq
    return trimmed

# Step 4: write FASTA
def write_fasta(seqs, out_path):
    with open(out_path, 'w') as f:
        for sid, seq in seqs.items():
            f.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

# Entry point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract conserved blocks from Clustal-format MSA file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--msa_path", required=True, help="Input Clustal MSA file path")
    parser.add_argument("-o", "--output_path", required=True, help="Output FASTA file path for conserved regions")
    parser.add_argument("-t", "--threshold", type=float, default=0.5,
                        help="Threshold for conservation (e.g., 0.5 means ≥50% of sequences)")
    parser.add_argument("-b", "--min_block_length", type=int, default=10,
                        help="Minimum length of continuous conserved block to include")

    args = parser.parse_args()

    try:
        seqs = parse_clustal_msa(args.msa_path)
        blocks = find_conserved_positions(seqs, threshold=args.threshold, min_block=args.min_block_length)
        conserved_seqs = extract_conserved(seqs, blocks)
        write_fasta(conserved_seqs, args.output_path)

        print(f"Conserved regions (≥{args.threshold*100:.0f}%) written to: {args.output_path}")
        print(f"Included {len(blocks)} blocks, total length: {sum(end-start for start, end in blocks)}")
    except Exception as e:
        print("An error occurred during processing:")
        print(str(e))
        parser.print_help()
