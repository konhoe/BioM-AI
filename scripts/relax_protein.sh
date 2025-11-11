#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Rosetta Protein Clean + Relax 자동 실행 스크립트
# 사용법:
#   ./relax_protein.sh ../input/protein/albumin/albumin.pdb A [nstruct]
# 예시:
#   ./relax_protein.sh ../input/protein/albumin/albumin.pdb A 1
# ============================================================

if [[ $# -lt 2 ]]; then
    echo "사용법: $0 <protein_pdb> <chain_id> [nstruct]"
    exit 1
fi

PROTEIN_PDB="$1"
CHAIN_ID="$2"
NSTRUCT="${3:-1}"

# ------------------------------------------------------------
# 환경 변수
# ------------------------------------------------------------
ROSETTA_HOME="${ROSETTA_HOME:-/Users/junyoung/Desktop/BioMatAI/rosetta/rosetta.binary.m1.release-371/main}"
SCRIPT_DIR="$(cd "$(dirname "$0")"; pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.."; pwd)"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# ------------------------------------------------------------
# 파일 경로 처리
# ------------------------------------------------------------
if [[ ! "$PROTEIN_PDB" =~ ^/ ]]; then
    PROTEIN_PDB="$(cd "$(dirname "$PROTEIN_PDB")" && pwd)/$(basename "$PROTEIN_PDB")"
fi

BASENAME=$(basename "$PROTEIN_PDB" .pdb)
EXP_NAME="relax_${BASENAME}_${TIMESTAMP}"
WORK_DIR="$ROOT_DIR/experiments/$EXP_NAME"
mkdir -p "$WORK_DIR"/{input,output,logs}

echo "========================================"
echo "Rosetta Protein Clean + Relax"
echo "========================================"
echo "입력 파일: $PROTEIN_PDB"
echo "체인: $CHAIN_ID"
echo "구조 수: $NSTRUCT"
echo "작업 디렉토리: $WORK_DIR"
echo "========================================"

# ------------------------------------------------------------
# clean_pdb.py 수행
# ------------------------------------------------------------
CLEAN_BIN="$ROSETTA_HOME/tools/protein_tools/scripts/clean_pdb.py"
if [[ ! -f "$CLEAN_BIN" ]]; then
    echo "[ERROR] clean_pdb.py를 찾을 수 없습니다: $CLEAN_BIN"
    exit 1
fi

echo "[INFO] clean_pdb.py 실행 중..."
python "$CLEAN_BIN" "$PROTEIN_PDB" "$CHAIN_ID" 2>&1 | tee "$WORK_DIR/logs/clean_pdb.log"

CLEANED_PDB_DIR="$(dirname "$PROTEIN_PDB")"
CLEANED_PDB="$CLEANED_PDB_DIR/${BASENAME}_${CHAIN_ID}.pdb"

if [[ ! -f "$CLEANED_PDB" ]]; then
    echo "[ERROR] clean_pdb 결과 PDB를 찾을 수 없습니다: $CLEANED_PDB"
    echo "로그 확인: tail -50 $WORK_DIR/logs/clean_pdb.log"
    exit 1
fi

# 깨끗한 파일을 작업 디렉토리로 복사
cp "$CLEANED_PDB" "$WORK_DIR/input/"
echo "[OK] clean 완료 → $WORK_DIR/input/$(basename "$CLEANED_PDB")"

# ------------------------------------------------------------
# Relax 실행 설정
# ------------------------------------------------------------
RELAX_INPUT="$WORK_DIR/input/$(basename "$CLEANED_PDB")"
FLAGS_FILE="$WORK_DIR/input/flags_relax"

cat > "$FLAGS_FILE" <<EOF
# Database & input
-database $ROSETTA_HOME/database
-s $RELAX_INPUT
-nstruct $NSTRUCT

# Output
-out:pdb
-out:path:all output/
-out:file:scorefile output/relax_score.sc
-overwrite

# Scoring
-score:weights ref2015

# Relax options
-relax:fast
-relax:constrain_relax_to_start_coords
-relax:coord_constrain_sidechains

# General
-run:jran $(date +%s)
EOF

echo "[OK] flags 파일 생성 완료: $FLAGS_FILE"

# ------------------------------------------------------------
# Relax 바이너리 탐색
# ------------------------------------------------------------
POSSIBLE_BINS=(
  "$ROSETTA_HOME/source/bin/relax.static.macosclangrelease"
  "$ROSETTA_HOME/source/bin/relax.macosclangrelease"
)
BIN=""
for b in "${POSSIBLE_BINS[@]}"; do
  if [[ -f "$b" ]]; then BIN="$b"; break; fi
done
if [[ -z "$BIN" ]]; then
  echo "[ERROR] Relax 바이너리를 찾을 수 없습니다 (${POSSIBLE_BINS[*]})"
  exit 1
fi
echo "[OK] 바이너리 확인: $(basename "$BIN")"

# ------------------------------------------------------------
# 실행
# ------------------------------------------------------------
cd "$WORK_DIR"
echo ""
echo "========================================"
echo "Relax 실행 시작"
echo "========================================"
echo "시작 시간: $(date)"
"$BIN" @input/flags_relax 2>&1 | tee logs/relax.log
EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "========================================"
echo "완료 (종료 코드: $EXIT_CODE)"
echo "========================================"

# ------------------------------------------------------------
# 결과 확인
# ------------------------------------------------------------
if [[ $EXIT_CODE -ne 0 ]]; then
  echo "[ERROR] Relax 실패"
  echo "로그 확인: tail -100 $WORK_DIR/logs/relax.log"
  exit 1
fi

TOTAL_PDBS=$(ls output/*.pdb 2>/dev/null | wc -l | tr -d ' ')
if [[ $TOTAL_PDBS -le 0 ]]; then
  echo "[ERROR] PDB 출력이 없습니다"
  echo "로그 확인: tail -50 $WORK_DIR/logs/relax.log"
  exit 1
fi

echo "[SUCCESS] Relax 완료!"
echo "생성된 구조 수: $TOTAL_PDBS 개"
ls -lh output/*.pdb | awk '{print "  " $9 " (" $5 ")"}'

if [[ -f output/relax_score.sc ]]; then
  echo ""
  echo "상위 5개 스코어:"
  tail -n +2 output/relax_score.sc | sort -k2 -n | head -5
else
  echo "[WARN] score.sc 파일이 없습니다. 로그에서 SCORE 라인 검색 권장."
fi

echo ""
echo "결과 파일:"
echo "  $WORK_DIR/output/"
echo "========================================"