#!/usr/bin/env bash
set -euo pipefail

# Rosetta Score Extraction Script
# 사용법: ./extract_scores.sh <experiment_dir>

if [[ $# -lt 1 ]]; then
    echo "사용법: $0 <experiment_dir>"
    echo "예시: $0 ../experiments/surfdock_merged_complex_20241002_171500"
    exit 1
fi

EXP_DIR="$1"

# 환경 설정
ROSETTA_HOME="${ROSETTA_HOME:-/Users/junyoung/Desktop/BioMatAI/rosetta/rosetta.binary.m1.release-371/main}"
SCRIPT_DIR="$(cd "$(dirname "$0")"; pwd)"

# 실험 디렉토리 절대경로화
if [[ ! "$EXP_DIR" =~ ^/ ]]; then
    EXP_DIR="$SCRIPT_DIR/$EXP_DIR"
    EXP_DIR="$(cd "$EXP_DIR" && pwd)"
fi

echo "========================================"
echo "Rosetta Score Extraction"
echo "========================================"
echo "실험 디렉토리: $EXP_DIR"

# 입력 확인
if [[ ! -d "$EXP_DIR" ]]; then
    echo "[ERROR] 실험 디렉토리가 없습니다: $EXP_DIR"
    exit 1
fi

if [[ ! -d "$EXP_DIR/output" ]]; then
    echo "[ERROR] output 디렉토리가 없습니다: $EXP_DIR/output"
    exit 1
fi

# PDB 파일 확인
PDB_FILES=($EXP_DIR/output/complex_*.pdb)
if [[ ! -f "${PDB_FILES[0]}" ]]; then
    echo "[ERROR] complex_*.pdb 파일이 없습니다"
    exit 1
fi

PDB_COUNT=${#PDB_FILES[@]}
echo "[OK] PDB 파일 ${PDB_COUNT}개 발견"

# Params 파일 확인
PARAMS_FILES=($EXP_DIR/input/*.params)
if [[ ! -f "${PARAMS_FILES[0]}" ]]; then
    echo "[WARN] Params 파일이 없습니다. 표준 residue만 스코어링합니다."
    PARAMS_FLAGS=""
else
    echo "[OK] Params 파일들:"
    PARAMS_FLAGS=""
    for params in "${PARAMS_FILES[@]}"; do
        if [[ ! "$(basename "$params")" =~ _centroid\.params$ ]]; then
            echo "     - $(basename "$params")"
            PARAMS_FLAGS="$PARAMS_FLAGS -in:file:extra_res_fa $params"
        fi
    done
fi

# Score 바이너리 확인
SCORE_BIN=""
for bin in score_jd2.static.macosclangrelease score_jd2.macosclangrelease; do
    if [[ -f "$ROSETTA_HOME/source/bin/$bin" ]]; then
        SCORE_BIN="$ROSETTA_HOME/source/bin/$bin"
        break
    fi
done

if [[ -z "$SCORE_BIN" ]]; then
    echo "[ERROR] score_jd2 바이너리를 찾을 수 없습니다"
    exit 1
fi

echo "[OK] 스코어링 바이너리: $(basename "$SCORE_BIN")"

# 기존 score.sc 백업
if [[ -f "$EXP_DIR/output/score.sc" ]]; then
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    mv "$EXP_DIR/output/score.sc" "$EXP_DIR/output/score_${TIMESTAMP}.sc.bak"
    echo "[INFO] 기존 score.sc를 백업했습니다"
fi

# 스코어 추출 실행
echo ""
echo "========================================"
echo "스코어 추출 시작"
echo "========================================"
echo "시작 시간: $(date)"

cd "$EXP_DIR/output"

"$SCORE_BIN" \
    -in:file:s complex_*.pdb \
    -out:file:scorefile score.sc \
    $PARAMS_FLAGS \
    > ../logs/score_extraction.log 2>&1

EXIT_CODE=$?

echo "완료 시간: $(date)"
echo "종료 코드: $EXIT_CODE"

# 결과 확인
echo ""
echo "========================================"
echo "결과 요약"
echo "========================================"

if [[ -f "score.sc" ]]; then
    SCORE_COUNT=$(tail -n +2 score.sc | wc -l | tr -d ' ')
    echo "[SUCCESS] score.sc 생성 완료"
    echo "스코어링된 구조 수: ${SCORE_COUNT}개"
    echo ""
    echo "상위 5개 구조 (낮은 점수가 좋음):"
    echo "----------------------------------------"
    tail -n +2 score.sc | sort -k2 -n | head -5 | awk '{printf "%-25s %12.2f\n", $NF, $2}'
    echo ""
    
    # 최고 구조 정보
    BEST_LINE=$(tail -n +2 score.sc | sort -k2 -n | head -1)
    BEST_NAME=$(echo "$BEST_LINE" | awk '{print $NF}')
    BEST_SCORE=$(echo "$BEST_LINE" | awk '{print $2}')
    
    echo "최고 점수 구조: $BEST_NAME"
    echo "Total Score: $BEST_SCORE"
    echo ""
    echo "PyMOL로 확인: pymol output/${BEST_NAME%_*}.pdb"
else
    echo "[ERROR] score.sc 파일이 생성되지 않았습니다"
    echo "로그 확인: cat $EXP_DIR/logs/score_extraction.log"
    exit 1
fi

echo "========================================"
echo "Score 파일 위치: $EXP_DIR/output/score.sc"
echo "========================================"