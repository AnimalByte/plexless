#!/usr/bin/env python3
"""Compare equivalent transformations on synthetic FASTQ; no speed claims.

Usage: python3 scripts/compare_fastp.py --fastp /path/to/fastp \
    --plexless target/debug/plexless [--save-golden tests/fixtures/fastp_v1_4_0.json]
Only sequences and qualities are compared: merge/header annotations differ.
"""
import argparse
import json
import random
import subprocess
import tempfile
from pathlib import Path

def dna(n):
    rng = random.Random(1874)
    return ''.join(rng.choice('ACGT') for _ in range(n))
def rc(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
def record(sequence,quality=None,mate=1):
    return f'@read/{mate}\n{sequence}\n+\n{quality or "I"*len(sequence)}\n'
def records(path):
    lines=path.read_text().splitlines()
    return [[lines[i+1],lines[i+3]] for i in range(0,len(lines),4)]
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastp',required=True,type=Path)
    parser.add_argument('--plexless',required=True,type=Path)
    parser.add_argument('--save-golden',type=Path)
    args=parser.parse_args()
    fastp=str(args.fastp.resolve()); plexless=str(args.plexless.resolve())
    version=subprocess.run([fastp,'--version'],capture_output=True,text=True,check=True)
    sequence=dna(80)
    cases=[
        ('global',record(sequence),None,['--trim_front1','3','--trim_tail1','5'],['--trim-front1','3','--trim-tail1','5']),
        ('front',record(sequence,'!'*8+'I'*72),None,['--cut_front'],['--cut-front']),
        ('tail',record(sequence,'I'*72+'!'*8),None,['--cut_tail'],['--cut-tail']),
        ('right',record(sequence,'I'*35+'!'*8+'I'*37),None,['--cut_right'],['--cut-right']),
        ('polyg',record(dna(48)+'CC'+'G'*15),None,['--trim_poly_g'],['--trim-poly-g']),
        ('polyx',record('ATTTTAAAAAAAAAATAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAT'),None,['--trim_poly_x'],['--trim-poly-x']),
        ('umi',record('ACGT'+sequence),None,['--umi','--umi_loc','read1','--umi_len','4'],['--umi','--umi-loc','read1','--umi-len','4']),
        ('merge',record(dna(220)[:150]),record(rc(dna(220)[70:]),mate=2),['--merge'],['--merge']),
        ('overlap_trim',record(dna(100)+'CGATTGACGATCCGAATGCA'),record(rc(dna(100))+'TGCAGGATCGTACCAGTACT',mate=2),[],[]),
    ]
    golden={'reference_version':(version.stdout+version.stderr).strip(),'source':'https://github.com/OpenGene/fastp/releases/tag/v1.4.0','cases':[]}
    with tempfile.TemporaryDirectory(prefix='plexless-fastp-compare-') as tmp:
        for name,first,second,fast_flags,rust_flags in cases:
            root=Path(tmp)/name;root.mkdir();(root/'in1.fq').write_text(first)
            base=[fastp,'-w','1','-Q','-L','--dont_eval_duplication','-i',str(root/'in1.fq'),'-o',str(root/'fastp1.fq'),'-j',str(root/'fastp.json'),'-h',str(root/'fastp.html')]
            # Sequence adapters/automatic discovery stay off for these paired
            # overlap comparisons: fastp's default overlap remains enabled.
            if second is None:base+=['-A']
            else: (root/'in2.fq').write_text(second);base+=['-I',str(root/'in2.fq'),'-O',str(root/'fastp2.fq')]
            cmd=[plexless,'trim','--output',str(root/'rust'),'--out1',str(root/'rust1.fq')]
            cmd+=['--r1',str(root/'in1.fq'),'--r2',str(root/'in2.fq'),'--out2',str(root/'rust2.fq')] if second else ['--reads',str(root/'in1.fq')]
            if name=='merge':base+=['--merged_out',str(root/'fastp_merged.fq')];cmd+=['--merged-out',str(root/'rust_merged.fq')]
            for command in [base+fast_flags,cmd+rust_flags]:
                run=subprocess.run(command,capture_output=True,text=True)
                if run.returncode:raise RuntimeError(f'{name}: {command}\n{run.stderr}')
            outputs=[]
            for suffix in (['merged'] if name=='merge' else ['1','2'] if second else ['1']):
                fast_path=root/(f'fastp_{suffix}.fq' if suffix=='merged' else f'fastp{suffix}.fq')
                rust_path=root/(f'rust_{suffix}.fq' if suffix=='merged' else f'rust{suffix}.fq')
                expected=records(fast_path);actual=records(rust_path)
                if actual!=expected:raise AssertionError(f'{name} {suffix}: fastp {expected!r}; Plexless {actual!r}')
                outputs.append(expected)
            golden['cases'].append({'name':name,'r1':first,'r2':second,'args':rust_flags,'outputs':outputs})
            print(f'{name}: sequences and qualities agree')
    if args.save_golden:args.save_golden.write_text(json.dumps(golden,indent=2)+'\n')
if __name__=='__main__':main()
