#!/usr/bin/env python
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os
import tempfile
import subprocess

current_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(current_folder)

def get_hmm_lengths(hmm_database_path):
    """
    return a dictionary of hmm model's lengths
    key is hmm NAME, value is models' length
    """
    hmm_len = {}
    with open(hmm_database_path, "r") as hf:
        cur_name = ""
        for line in hf.readlines():
            if line.startswith("NAME"):
                cur_name = line.split(" ")[-1].rstrip()
            if line.startswith("LENG") and len(cur_name) > 0:
                hmm_len[cur_name] = int(line.split(" ")[-1].rstrip())
    return hmm_len

def hmmscan(sequence, hmm_file, hmm_len, cutoff):
    """
    given an AA sequence and a .hmm filepath, perform hmmscan
    and parse its result into JSON-formatted list
    example JSON:
    [
        {
            "model_start": 5,
            "model_end": 100,
            "seq_start": 10,
            "seq_end": 100,
            "bitscore": 230,
            "sequence": "AMAMAMAMAMMAMAMAA"
        }
    ]
    """
    result = []

    with tempfile.TemporaryDirectory() as temp_dir:
        fasta_path = os.path.join(temp_dir, "temp.fa")
        fasta_file = open(fasta_path, "w")
        fasta_file.write(">temp\n{}\n".format(sequence))

        fasta_file.close()
        result_path = os.path.join(temp_dir, "temp.output")
        command = "hmmscan -o {} -T {} --incT {} {} {}".format(result_path, cutoff, cutoff, hmm_file, fasta_path)
        subprocess.check_output(command, shell=True)
        for runresult in SearchIO.parse(result_path, 'hmmer3-text'):
            for hsp in runresult.hsps:
                padding_left = ""
                padding_right = ""
                for i in range(0, hsp.hit_start):
                    padding_left += "-"
                for i in range(hsp.hit_end, hmm_len[hsp.hit_id] + 1):
                    padding_right += "-"
                if hsp.hit_strand != hsp.query_strand:
                    padding_left, padding_right = padding_right, padding_left
                seq = ""
                hmmseq = str(hsp.hit.seq)
                queryseq = str(hsp.query.seq)
                for i in range(0, len(hmmseq)):
                    if hmmseq[i] != '.':
                        seq += queryseq[i]
                seq = "".join([padding_left, seq, padding_right])
                result.append({
                    "name": hsp.hit_id,
                    "model_start": hsp.hit_start,
                    "model_end": hsp.hit_end,
                    "seq_start": hsp.query_start,
                    "seq_end": hsp.query_end,
                    "bitscore": hsp.bitscore,
                    "seq": seq
                })

    # filter overlapping hits
    temp_result = []
    result = sorted(result, key=lambda obj:obj["seq_start"])
    for hit in result:
        if len(temp_result) > 0:
            if (hit["seq_start"] - temp_result[-1]["seq_end"]) < 0: # is overlapping
                if hit["bitscore"] > temp_result[-1]["bitscore"]:
                    temp_result.pop()
                else:
                    continue
        temp_result.append(hit)
    result = temp_result
                
    return result

if __name__ == "__main__":
    sequence = "%s/examples/lipstatin.fasta" % parent_folder
    hmm_file = "%s/data/AMP-binding_3.2.hmm" % parent_folder
    # hmmlen = get_hmm_lengths("%s/data/AMP-binding.hmm" % parent_folder)
    hmmlen = 418
    cutoff = 300
    # print(hmmlen)
    hmmscan(sequence, hmm_file, hmmlen, cutoff)
