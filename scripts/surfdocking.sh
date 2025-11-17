#!/usr/bin/env bash
set -euo pipefail

# Rosetta Surface Docking
# 사용법:
#   1) DB 모드(로제타 DB에 params 설치해둔 경우):
#        ./surfdocking.sh <merged_pdb> DB [nstruct] [protein_name]
#   2) 로컬 라이브러리에서 자동 선택(플랫 폴더):
#        ./surfdocking.sh <merged_pdb> <params_dir> [nstruct] [protein_name]
#      -> PDB의 HETATM resname과 일치하는 *.params만 자동 선택하여 -extra_res_fa에 추가
#
#   예:
#        ./surfdocking.sh input/merge_pdb/merged_albumin_TiAl.pdb params/TiAl 1 albumin
#        ./surfdocking.sh input/merge_pdb/merged_1H3Y_TiAl.pdb    params/TiAl 5 1H3Y

if [[ $# -lt 2 ]]; then
    echo "사용법: $0 <merged_pdb> <params_dir|DB> [nstruct] [protein_name]"
    echo "예시:  $0 input/merge_pdb/merged_complex.pdb DB 1 albumin"
    echo "       $0 input/merge_pdb/merged_complex.pdb params/TiAl 5 1H3Y"
    exit 1
fi

MERGED_PDB="$1"
PARAMS_SRC="$2"           # 디렉토리 또는 DB
NSTRUCT="${3:-1}"         # 기본 1개 구조 생성
PROT_NAME="${4:-albumin}" # 기본 단백질 이름: albumin

# 환경 설정
ROSETTA_HOME="${ROSETTA_HOME:-/Users/junyoung/Desktop/BioMatAI/rosetta/rosetta.binary.m1.release-371/main}"
SCRIPT_DIR="$(cd "$(dirname "$0")"; pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.."; pwd)"

# --- merged pdb 절대경로 정규화 (현재 위치와 무관하게 동작) ---
if [[ "$MERGED_PDB" != /* ]]; then
  MERGED_PDB="$(cd "$(dirname "$MERGED_PDB")" 2>/dev/null && pwd)/$(basename "$MERGED_PDB")"
fi

# --- params dir (DB 모드가 아닐 때만) ---
PARAMS_DIR=""
if [[ "$PARAMS_SRC" != "DB" ]]; then
  PARAMS_DIR="$PARAMS_SRC"
  # ~ 확장 보정
  PARAMS_DIR="${PARAMS_DIR/#\~/$HOME}"

  # 상대경로면 ROOT_DIR 기준으로 절대경로화
  if [[ "$PARAMS_DIR" != /* ]]; then
    PARAMS_DIR="$ROOT_DIR/$PARAMS_DIR"
  fi
fi

echo "========================================"
echo "Rosetta Surface Docking"
echo "========================================"
echo "병합된 PDB  : $MERGED_PDB"
echo "Params 소스 : $PARAMS_SRC"
echo "구조 수     : $NSTRUCT"
echo "단백질 이름 : $PROT_NAME"
echo "========================================"

# 입력 확인
if [[ ! -f "$MERGED_PDB" ]]; then
  echo "[ERROR] 병합된 PDB 파일이 없습니다: $MERGED_PDB"
  exit 1
fi

if [[ "$PARAMS_SRC" != "DB" ]]; then
  if [[ ! -d "$PARAMS_DIR" ]]; then
    echo "[ERROR] Params 디렉토리가 없습니다: $PARAMS_DIR"
    exit 1
  fi
fi

# 실험 디렉토리
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
EXP_NAME="surfdock_$(basename "$MERGED_PDB" .pdb)_${TIMESTAMP}"
WORK_DIR="$ROOT_DIR/experiments/$EXP_NAME"
mkdir -p "$WORK_DIR"/{input,output,logs}
cd "$WORK_DIR"

# 입력 복사
cp "$MERGED_PDB" "input/complex.pdb"

# 체인 ID 추출
PROTEIN_CHAIN=$(grep '^ATOM'   input/complex.pdb | head -1 | awk '{print substr($0,22,1)}' || true)
SURFACE_CHAIN=$(grep '^HETATM' input/complex.pdb | head -1 | awk '{print substr($0,22,1)}' || true)
[[ -z "$PROTEIN_CHAIN" ]] && PROTEIN_CHAIN="A"
[[ -z "$SURFACE_CHAIN" ]] && SURFACE_CHAIN="Z"
echo "[INFO] Protein chain: $PROTEIN_CHAIN"
echo "[INFO] Surface/Metal chain: $SURFACE_CHAIN"

# surface_vectors (현재는 Ti 기본 사용, 나중에 필요하면 metal별로 분기)
SURFVECS_PATH="input/surface_vectors.txt"
USER_SURFACE_VECTORS="$ROOT_DIR/input/metal/Ti/Ti_surface_vectors.txt"
if [[ -f "$USER_SURFACE_VECTORS" ]]; then
    cp "$USER_SURFACE_VECTORS" "$SURFVECS_PATH"
    echo "[OK] surface_vectors 사용: $SURFVECS_PATH (사용자 제공)"
else
    echo "[INFO] surface_vectors 자동 생성"
    GEN_BIN="$ROSETTA_HOME/source/bin/surface_docking.static.macosclangrelease"
    if [[ ! -f "$GEN_BIN" ]]; then
        echo "[ERROR] surface_docking 바이너리 없음: $GEN_BIN"
        exit 1
    fi
    "$GEN_BIN" -in:file:s input/complex.pdb -in:file:generate_surface_vectors > logs/generate_surface_vectors.log 2>&1 || true
    if   [[ -f "complex.surfvec" ]];        then mv "complex.surfvec"        "$SURFVECS_PATH"
    elif [[ -f "input/complex.surfvec" ]];  then mv "input/complex.surfvec"  "$SURFVECS_PATH"; fi
    [[ -f "$SURFVECS_PATH" ]] || { echo "[ERROR] surface_vectors 생성 실패. logs/generate_surface_vectors.log 확인"; exit 1; }
    echo "[OK] surface_vectors 생성 완료: $SURFVECS_PATH"
fi

# Fragment(옵션) - input/protein/<PROT_NAME>/<PROT_NAME>_3mers, _9mers 사용
FRAG3_NAME="${PROT_NAME}_3mers"
FRAG9_NAME="${PROT_NAME}_9mers"
FRAG_DIR="$ROOT_DIR/input/protein/${PROT_NAME}"
FRAG3_SRC="$FRAG_DIR/${FRAG3_NAME}"
FRAG9_SRC="$FRAG_DIR/${FRAG9_NAME}"
USE_FRAGMENTS=false

if [[ -f "$FRAG3_SRC" && -f "$FRAG9_SRC" ]]; then
    cp "$FRAG3_SRC" "input/${FRAG3_NAME}"
    cp "$FRAG9_SRC" "input/${FRAG9_NAME}"
    USE_FRAGMENTS=true
    echo "[OK] Fragment 파일 복사 완료 (protein = $PROT_NAME)"
else
    echo "[INFO] ${PROT_NAME} fragment 파일 없음 - fragment-free 모드"
fi

# --- HETATM resname 기반으로 필요한 params만 자동 선택 ---
declare -a CHOSEN_PARAMS=()
if [[ "$PARAMS_SRC" != "DB" ]]; then
    echo "[INFO] PDB(HETATM)에서 필요한 params 자동 선택"
    # HETATM의 residue name(18-20열)만 수집. 물/중복 제거.
    RESNAMES=$(awk '/^HETATM/ {r=substr($0,18,3); if (r!="HOH" && r!="WAT" && r!="DOD") print r}' input/complex.pdb | sort | uniq)
    if [[ -z "${RESNAMES:-}" ]]; then
        echo "[WARN] HETATM resname이 없어 로컬 params 추가가 필요 없어 보임."
    else
        echo "[INFO] 발견된 HETATM resname: $RESNAMES"
        for r in $RESNAMES; do
            cand1="$PARAMS_DIR/${r}.params"
            cand2="$PARAMS_DIR/${r}0.params"   # 예: TI → TI0
            chosen=""
            if   [[ -f "$cand1" ]]; then chosen="$cand1"
            elif [[ -f "$cand2" ]]; then chosen="$cand2"; fi
            if [[ -n "$chosen" ]]; then
                cp "$chosen" "input/"
                CHOSEN_PARAMS+=("input/$(basename "$chosen")")
                echo "  ✓ $(basename "$chosen")"
            else
                echo "  ⚠️  매칭 params 없음: $r  (디렉토리: $PARAMS_DIR)"
            fi
        done
        if [[ ${#CHOSEN_PARAMS[@]} -eq 0 ]]; then
            echo "[WARN] 선택된 params가 없습니다. DB에 설치된 params만으로 진행합니다."
        fi
    fi
fi

# Flags 파일 생성
FLAGS_FILE="input/flags"
cat > "$FLAGS_FILE" << EOF
# Rosetta Database
-database $ROSETTA_HOME/database

# Input
-s input/complex.pdb
-nstruct $NSTRUCT

# Output
-out:pdb
-out:file:scorefile output/score.sc
-out:path:all output/
-overwrite

# Score Function
-score:weights ref2015

# Surface Docking
-docking:partners ${PROTEIN_CHAIN}_${SURFACE_CHAIN}
-in:file:surface_vectors $SURFVECS_PATH
-docking:low_res_protocol_only false
-docking:dock_mcm_first_cycles 0
-include_surfaces

# Packing
-packing:ex1
-packing:ex2aro
EOF

# Fragment 옵션
if [[ "$USE_FRAGMENTS" == "true" ]]; then
cat >> "$FLAGS_FILE" << EOF
# Fragment Files
-in:file:frag3 input/${FRAG3_NAME}
-in:file:frag9 input/${FRAG9_NAME}
EOF
fi

# 선택된 params만 -extra_res_fa로 추가
if [[ ${#CHOSEN_PARAMS[@]} -gt 0 ]]; then
    echo "# Metal/Ligand Params (Full Atom)" >> "$FLAGS_FILE"
    for pf in "${CHOSEN_PARAMS[@]}"; do
        echo "-extra_res_fa $pf" >> "$FLAGS_FILE"
    done
fi

# 일반 옵션
cat >> "$FLAGS_FILE" << EOF

# General Options
-run:preserve_header
-run:jran $(date +%s)
EOF

echo "[OK] 설정 파일 생성 완료: $FLAGS_FILE"
echo "----------------------------------------"
cat "$FLAGS_FILE"
echo "----------------------------------------"

# 실행
BIN="$ROSETTA_HOME/source/bin/surface_docking.static.macosclangrelease"
[[ -f "$BIN" ]] || { echo "[ERROR] surface_docking 바이너리 없음: $BIN"; exit 1; }

echo ""
echo "========================================"
echo "Surface Docking 실행 시작"
echo "========================================"
echo "시작 시간: $(date)"
echo "작업 디렉토리: $WORK_DIR"
echo "Fragment 모드: $USE_FRAGMENTS"
echo ""

"$BIN" @input/flags 2>&1 | tee logs/surfdock.log &
ROSETTA_PID=$!
echo "프로세스 ID: $ROSETTA_PID"
echo ""

wait $ROSETTA_PID
EXIT_CODE=$?

echo ""
echo "========================================"
echo "완료 (종료 코드: $EXIT_CODE)"
echo "========================================"
echo "종료 시간: $(date)"
echo ""

if [[ -f "output/score.sc" ]]; then
    TOTAL_PDBS=$(ls output/*.pdb 2>/dev/null | wc -l | tr -d ' ')
    echo "✓ 생성된 구조 수: ${TOTAL_PDBS}개"
    if [[ -s "output/score.sc" && $TOTAL_PDBS -gt 0 ]]; then
        echo ""
        echo "상위 5개 구조 (total_score 기준):"
        echo "----------------------------------------"
        tail -n +2 "output/score.sc" | sort -k2 -n | head -5 | awk '{printf "  %s: %.2f\n", $NF, $2}'
        echo "----------------------------------------"
    fi
    echo ""
    echo "결과 위치: $WORK_DIR/output/"
else
    echo "❌ ERROR - 결과 생성 실패"
    echo ""
    echo "마지막 50줄 로그:"
    echo "----------------------------------------"
    tail -50 logs/surfdock.log
    echo "----------------------------------------"
    echo ""
    echo "전체 로그: $WORK_DIR/logs/surfdock.log"
fi

echo "========================================"
