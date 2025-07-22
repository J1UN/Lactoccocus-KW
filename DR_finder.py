import os
import re
import argparse
from Bio import SeqIO

def find_direct_repeats(seq1_id, seq1, seq2_id, seq2, seq_starts):
    """Find direct repeats (DRs) between two sequences."""
    window_sizes = range(40, 9, -1)  # Search window sizes from 40bp to 10bp
    direct_repeats = []
    detected_positions = set()
    
    start1, start2 = seq_starts[seq1_id], seq_starts[seq2_id]

    for size in window_sizes:
        for i in range(len(seq1) - size + 1):
            if (seq1_id, i) in detected_positions:
                continue
            window = seq1[i:i+size]
            for j in range(len(seq2) - size + 1):
                if (seq2_id, j) in detected_positions:
                    continue
                window2 = seq2[j:j+size]
                identity = sum(1 for x, y in zip(window, window2) if x == y) / size
                if identity >= 0.9:
                    # Save only if the repeat contains TTTTTG or CAAAAA
                    if 'TTTTTG' in window or 'CAAAAA' in window:
                        pos_1, pos_2 = start1 + i, start2 + j
                        dr_length = len(window)
                        direct_repeats.append(
                            (seq1_id, window, pos_1, seq2_id, window2, pos_2, dr_length, identity)
                        )
                        # Mark these positions to avoid overlapping detection
                        for k in range(size):
                            detected_positions.add((seq1_id, i + k))
                            detected_positions.add((seq2_id, j + k))
    return direct_repeats

def find_and_save_direct_repeats(input_file, output_file):
    """Find direct repeats from a FASTA file and save them to a TSV output file."""
    sequences, seq_starts = {}, {}

    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            header = record.id
            match = re.search(r':(\d+)-\d+', header)
            start_position = int(match.group(1)) if match else 1
            sequences[header] = str(record.seq)
            seq_starts[header] = start_position

    # Ensure exactly two sequences are provided
    if len(sequences) != 2:
        raise ValueError(f"Input file '{input_file}' must contain exactly two sequences.")

    direct_repeat_records = []
    seq_ids = list(sequences.keys())
    seq1_id, seq2_id = seq_ids[0], seq_ids[1]
    seq1, seq2 = sequences[seq1_id], sequences[seq2_id]

    direct_repeat_records = find_direct_repeats(seq1_id, seq1, seq2_id, seq2, seq_starts)

    with open(output_file, "w") as f:
        f.write("SeqID1\tRepeat1\tPos1\tSeqID2\tRepeat2\tPos2\tDR_len\tper_match\n")
        for record in direct_repeat_records:
            seqid1, dr1, pos1, seqid2, dr2, pos2, dr_length, identity = record
            line = f"{seqid1}\t{dr1}\t{pos1}\t{seqid2}\t{dr2}\t{pos2}\t{dr_length}\t{identity:.2f}\n"
            f.write(line)

def main():
    parser = argparse.ArgumentParser(description="""
    Find direct repeats (DRs) between two regions in FASTA files.
    Input FASTA must contain exactly two sequences (subjective regions).
    Only DRs containing 'TTTTTG' or 'CAAAAA' are saved.
    """)
    parser.add_argument("-i", "--input_dir", required=True, help="Path to input directory containing FASTA (.txt) files.")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to output directory for TSV result files.")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        if file_name.endswith(".txt"):
            input_path = os.path.join(input_dir, file_name)
            output_path = os.path.join(output_dir, file_name.replace(".txt", ".tsv"))
            try:
                find_and_save_direct_repeats(input_path, output_path)
                print(f"Processed: {file_name}")
            except ValueError as e:
                print(f"Error in {file_name}: {e}")

if __name__ == "__main__":
    main()

# EOF
