#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
원본 PDB를 Rosetta surface_docking 입력 형식으로 정리:
- 표면(resname 리스트) 라인에 chain ID 부여(기본 Z)
- (옵션) 표면 resSeq 재번호(1부터)
- 단백질(기본 A 등) -> 표면(Z) 순서로 정렬
- 섹션 사이 TER/END 정리
- HETATM/ATOM 고정폭 포맷 유지
"""

import argparse
import re
from typing import List, Tuple

AA3 = {
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
}

def parse_args():
    p = argparse.ArgumentParser(description="Fix PDB for Rosetta surface_docking.")
    p.add_argument("-i","--inp", required=True, help="input PDB")
    p.add_argument("-o","--out", required=True, help="output PDB")
    p.add_argument("--surface-res", nargs="+", required=True,
                   help="surface residue names (e.g. TI4 PSF CAL)")
    p.add_argument("--surface-chain", default="Z", help="chain ID to assign to surface (default Z)")
    p.add_argument("--protein-chains", nargs="+", default=["A"],
                   help="protein chain IDs to keep as protein section (default A). \
                         비어 있으면 자동 분류(AA3면 단백질로).")
    p.add_argument("--renumber-surface", action="store_true",
                   help="surface residues를 1부터 연속 번호로 재부여")
    p.add_argument("--start-resseq", type=int, default=1,
                   help="renumber-surface 시작 번호 (default 1)")
    return p.parse_args()

def is_atom_line(line:str)->bool:
    rec = line[:6]
    return rec.startswith("ATOM  ") or rec.startswith("HETATM")

def get_fields(line:str)->Tuple[str,int,str,str,str,str]:
    """쉬운 파싱용: (recName, serial, atomName, resName, chainID, resSeq)"""
    rec = line[:6]
    serial = int(line[6:11])
    atom  = line[12:16]
    resn  = line[17:21].strip()
    chain = line[21].strip()
    resi  = line[22:26].strip()
    return rec, serial, atom, resn, chain, resi

def format_atom_line(line:str, *, resn=None, chain=None, resi=None)->str:
    """필요한 칼럼만 덮어써서 PDB 고정폭 포맷 유지"""
    # pad line to 80 cols
    s = list(line.rstrip("\n"))
    if len(s) < 80:
        s += [" "] * (80 - len(s))
    if resn is not None:
        s[17:21] = f"{resn:>3} " if len(resn)==3 else f"{resn:>4}"
    if chain is not None and len(chain)==1:
        s[21] = chain
    if resi is not None:
        s[22:26] = f"{int(resi):>4}"
    return "".join(s) + "\n"

def main():
    args = parse_args()
    surf_set = set([r.strip() for r in args.surface_res])
    protein_chains = None if not args.protein_chains else set(args.protein_chains)

    prot_lines: List[str] = []
    surf_lines: List[str] = []
    other_lines: List[str] = []  # REMARK, CRYST1 등

    with open(args.inp, "r") as f:
        for line in f:
            if not is_atom_line(line):
                # TER/END 등은 나중에 다시 찍을 것이므로 저장해뒀다가 제일 앞의 헤더만 살림
                other_lines.append(line)
                continue

            rec, serial, atom, resn, chain, resi = get_fields(line)

            # 표면 분류: resname이 surface 목록에 있으면 표면
            if resn in surf_set:
                # 체인ID 강제 부여
                line = format_atom_line(line, chain=args.surface_chain)
                surf_lines.append(line)
            else:
                # 단백질/그 외: 지정 체인에 속하면 단백질, 아니면 단백질로 추정(AA3)
                is_protein = False
                if protein_chains:
                    is_protein = (chain in protein_chains)
                else:
                    is_protein = (resn in AA3)
                (prot_lines if is_protein else prot_lines).append(line)  # 기본은 단백질로 취급

    # 표면 resSeq 재번호(옵션)
    if args.renumber_surface and surf_lines:
        new_lines = []
        resseq_map = {}
        next_resi = args.start_resseq
        prev_key = None
        # resSeq는 (resName, chain, resSeq) 조합으로 바뀌어야 하지만,
        # 표면은 일반적으로 각 원자가 자기 residue를 따르므로 resName+oldResSeq 기준 매핑
        for line in surf_lines:
            rec, serial, atom, resn, chain, resi = get_fields(line)
            key = (resn, chain, resi)
            if key not in resseq_map:
                resseq_map[key] = next_resi
                next_resi += 1
            line = format_atom_line(line, resi=resseq_map[key])
            new_lines.append(line)
        surf_lines = new_lines

    # 출력 조립: 헤더(비원자 라인 중 CRYST1/REMARK 등 일부) -> 단백질 -> TER -> 표면 -> TER -> END
    with open(args.out, "w") as w:
        # 헤더 중 CRYST1/REMARK/HETNAM이 있으면 유지
        for line in other_lines:
            if line.startswith(("REMARK", "CRYST1", "SCALE", "ORIGX", "HETNAM", "FORMUL", "HEADER", "TITLE")):
                w.write(line)
        if prot_lines:
            for l in prot_lines:
                w.write(l)
            w.write("TER\n")
        if surf_lines:
            for l in surf_lines:
                w.write(l)
            w.write("TER\n")
        w.write("END\n")

if __name__ == "__main__":
    main()