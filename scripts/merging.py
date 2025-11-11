#!/usr/bin/env python3
# merge_pdb.py
import os
import sys

def merge_pdb(slab_file, protein_file, output_file):
    with open(slab_file, "r") as f1, open(protein_file, "r") as f2:
        slab_lines = f1.readlines()
        prot_lines = f2.readlines()

    merged = []

    # slab 파트: END, TER 같은 건 빼고 atom/hetatm만 남김
    for line in slab_lines:
        if line.startswith(("ATOM", "HETATM")):
            merged.append(line)

    # protein 파트도 마찬가지
    for line in prot_lines:
        if line.startswith(("ATOM", "HETATM", "SSBOND")):
            merged.append(line)

    # 마지막에 END 붙이기
    merged.append("END\n")

    with open(output_file, "w") as f:
        f.writelines(merged)

    print(f"[OK] merged pdb saved to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python merging.py <slab_file> <protein_file>")
        sys.exit(1)

    slab_file = sys.argv[1]
    protein_file = sys.argv[2]

    # Extract metal name from slab_file path
    # 예: TiAlV_slab.pdb -> TiAlV
    # 예: fix_Ti.pdb -> Ti (기존 형식도 지원)
    slab_basename = os.path.basename(slab_file)
    
    if '_slab.pdb' in slab_basename:
        # TiAlV_slab.pdb 형식
        metal_name = slab_basename.replace('_slab.pdb', '')
    elif slab_basename.startswith('fix_'):
        # fix_Ti.pdb 형식 (기존 호환)
        metal_name = slab_basename.split('_')[1].split('.')[0]
    else:
        # 기타 형식: 첫 번째 언더스코어 전까지
        metal_name = slab_basename.split('_')[0].replace('.pdb', '')

    # Extract protein name from protein_file path
    # 예: albumin_clean.pdb -> albumin
    # 예: albumin_A_0001.pdb -> albumin
    protein_basename = os.path.basename(protein_file)
    protein_name = protein_basename.split('_')[0]

    # Build output file path (절대 경로로 변환)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    output_file = os.path.join(project_root, "input", "merge_pdb", f"merged_{protein_name}_{metal_name}.pdb")
    
    # 출력 디렉토리 생성
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    merge_pdb(slab_file, protein_file, output_file)