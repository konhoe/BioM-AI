#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
surface_vector_generator.py
- PDB 파일에서 표면 법선 벡터(surface normal vector) 계산
- 표면 원자들의 좌표로부터 평면을 피팅하여 법선 방향 계산

Usage:
  python surface_vector_generator.py --pdb cleaned_surface.pdb --output surface_vectors.txt
"""

import argparse
import numpy as np
from pathlib import Path

def parse_args():
    ap = argparse.ArgumentParser(description="PDB 파일로부터 surface vector 생성")
    ap.add_argument("--pdb", required=True, help="표면 PDB 파일")
    ap.add_argument("--output", default="surface_vectors.txt", help="출력 파일")
    ap.add_argument("--method", choices=["pca", "z-axis", "centroid"], default="pca",
                    help="계산 방법: pca(주성분분석), z-axis(단순 Z축), centroid(중심점 기반)")
    return ap.parse_args()

def read_surface_coords(pdb_file):
    """PDB에서 표면 원자 좌표 읽기"""
    coords = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            
            try:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
            except (ValueError, IndexError):
                continue
    
    return np.array(coords)

def calculate_surface_normal_pca(coords):
    """
    PCA(주성분분석)를 사용하여 표면 법선 계산
    - 좌표들이 이루는 평면에 수직인 벡터를 찾음
    """
    # 중심화
    centroid = coords.mean(axis=0)
    centered = coords - centroid
    
    # 공분산 행렬의 고유값 분해
    cov_matrix = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
    
    # 가장 작은 고유값에 해당하는 고유벡터가 법선 방향
    normal = eigenvectors[:, 0]
    
    # Z축 양의 방향을 향하도록 정규화
    if normal[2] < 0:
        normal = -normal
    
    # 단위 벡터로 정규화
    normal = normal / np.linalg.norm(normal)
    
    return normal, centroid

def calculate_surface_normal_z_axis(coords):
    """
    간단한 방법: Z축을 법선으로 가정
    """
    centroid = coords.mean(axis=0)
    normal = np.array([0.0, 0.0, 1.0])
    return normal, centroid

def calculate_surface_normal_centroid(coords):
    """
    중심점 기반: 가장 위/아래 원자들의 평균 위치로 방향 결정
    """
    centroid = coords.mean(axis=0)
    
    # Z 좌표 기준으로 정렬
    z_coords = coords[:, 2]
    z_sorted_indices = np.argsort(z_coords)
    
    # 상위 10%와 하위 10% 원자들
    n_sample = max(10, len(coords) // 10)
    top_atoms = coords[z_sorted_indices[-n_sample:]]
    bottom_atoms = coords[z_sorted_indices[:n_sample]]
    
    # 방향 계산
    top_center = top_atoms.mean(axis=0)
    bottom_center = bottom_atoms.mean(axis=0)
    
    normal = top_center - bottom_center
    normal = normal / np.linalg.norm(normal)
    
    return normal, centroid

def generate_surface_vectors_file(normal, centroid, output_file, coords):
    """
    Rosetta surface_docking에 필요한 형식으로 파일 생성
    
    파일 형식: 3개의 벡터 (각 줄에 x y z)
    - 첫 번째: 표면 중심점 좌표
    - 두 번째: 표면 상의 다른 점 (중심에서 X방향)
    - 세 번째: 표면 상의 또 다른 점 (중심에서 Y방향)
    """
    
    # 표면 평면 상에서 직교하는 두 벡터 찾기
    # 법선 벡터에 수직인 벡터들을 생성
    
    # 임의의 벡터 선택 (법선과 평행하지 않은)
    if abs(normal[0]) < 0.9:
        arbitrary = np.array([1.0, 0.0, 0.0])
    else:
        arbitrary = np.array([0.0, 1.0, 0.0])
    
    # 법선에 수직인 첫 번째 벡터
    tangent1 = np.cross(normal, arbitrary)
    tangent1 = tangent1 / np.linalg.norm(tangent1)
    
    # 법선과 tangent1에 수직인 두 번째 벡터
    tangent2 = np.cross(normal, tangent1)
    tangent2 = tangent2 / np.linalg.norm(tangent2)
    
    # 표면 상의 세 점 생성 (거리 10Å 정도)
    point1 = centroid
    point2 = centroid + tangent1 * 10.0
    point3 = centroid + tangent2 * 10.0
    
    with open(output_file, 'w') as f:
        # 3개의 점을 각 줄에 출력 (주석 없이)
        f.write(f"{point1[0]:8.3f} {point1[1]:8.3f} {point1[2]:8.3f}\n")
        f.write(f"{point2[0]:8.3f} {point2[1]:8.3f} {point2[2]:8.3f}\n")
        f.write(f"{point3[0]:8.3f} {point3[1]:8.3f} {point3[2]:8.3f}\n")
    
    print(f"Surface vector 파일 생성: {output_file}")
    print(f"  Point 1 (중심): ({point1[0]:.3f}, {point1[1]:.3f}, {point1[2]:.3f})")
    print(f"  Point 2 (X방향): ({point2[0]:.3f}, {point2[1]:.3f}, {point2[2]:.3f})")
    print(f"  Point 3 (Y방향): ({point3[0]:.3f}, {point3[1]:.3f}, {point3[2]:.3f})")
    print(f"  표면 법선: ({normal[0]:.6f}, {normal[1]:.6f}, {normal[2]:.6f})")

def analyze_surface_geometry(coords):
    """표면 기하학 분석 정보 출력"""
    print("\n=== 표면 기하학 분석 ===")
    print(f"총 원자 수: {len(coords)}")
    
    # 좌표 범위
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    ranges = maxs - mins
    
    print(f"X 범위: {mins[0]:.2f} ~ {maxs[0]:.2f} (폭: {ranges[0]:.2f} Å)")
    print(f"Y 범위: {mins[1]:.2f} ~ {maxs[1]:.2f} (폭: {ranges[1]:.2f} Å)")
    print(f"Z 범위: {mins[2]:.2f} ~ {maxs[2]:.2f} (높이: {ranges[2]:.2f} Å)")
    
    # 평면도 판단
    if ranges[2] < 5.0:
        print(f"→ 평평한 표면 (Z 변화 < 5Å)")
    else:
        print(f"→ 굴곡진 표면 (Z 변화 > 5Å)")
    
    print()

def main():
    args = parse_args()
    
    print(f"PDB 파일 읽기: {args.pdb}")
    coords = read_surface_coords(args.pdb)
    
    if len(coords) == 0:
        print("ERROR: PDB 파일에서 좌표를 읽을 수 없습니다.")
        return 1
    
    # 표면 분석
    analyze_surface_geometry(coords)
    
    # 법선 벡터 계산
    print(f"계산 방법: {args.method}")
    
    if args.method == "pca":
        normal, centroid = calculate_surface_normal_pca(coords)
        print("PCA를 사용하여 최적 평면 피팅")
    elif args.method == "z-axis":
        normal, centroid = calculate_surface_normal_z_axis(coords)
        print("Z축을 법선으로 가정")
    elif args.method == "centroid":
        normal, centroid = calculate_surface_normal_centroid(coords)
        print("상하 원자들의 중심점으로 방향 계산")
    
    # 파일 생성
    generate_surface_vectors_file(normal, centroid, args.output, coords)
    
    print("\n사용법:")
    print(f"  surface_docking ... -in:file:surface_vectors {args.output}")
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())