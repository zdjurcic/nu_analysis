#!/bin/bash
# Sync between workshop local_full/ and nu_analysis/2x2_mc/
# local_full is the dev sandbox; 2x2_mc is the publishable copy (no ROOTs)
# Usage: ./sync_local_full.sh [to_nu | to_local | status]
set -e
WS="/Users/zdjurcic/ZelimirDjurcic/DUNE/DUNE_related_Computing/service-enabled-science-jobs-data-ai-agents_Sep_30_2026/Service_Enabled_Science/05_Agentic_Workflows/caf_analysis_code/local_full"
NU="/Users/zdjurcic/ZelimirDjurcic/nu_analysis/2x2_mc"
MODE=${1:-status}
if [ "$MODE" = "to_nu" ]; then
  echo "→ Sync local_full → 2x2_mc (no ROOTs, keep plots)"
  rsync -av --exclude='caf_MC/*.root' --exclude='__pycache__' --exclude='*.pyc' --exclude='test_*.py' --exclude='fix_syntax2.py' --exclude='run_full_*.py' \
    "$WS/analyze_verified.py" "$WS/analyze_verified_minerva.py" "$WS/analyze_full.py" "$WS/nominal_mc.cc" "$WS/README.md" "$WS/analysis_summary.html" "$NU/" 2>&1 | tail -20
  rsync -av "$WS/plots/" "$NU/plots/" --include='*.png' --include='*.json' 2>&1 | tail -5
  rsync -av "$WS/plots_verified/" "$NU/plots_verified/" 2>&1 | tail -5
  rsync -av "$WS/deps/duneanaobj/LinkDef.h" "$WS/deps/duneanaobj/LinkDef_full.h" "$NU/deps/duneanaobj/" 2>&1 | tail -5
  rsync -av "$WS/deps/duneanaobj/libduneanaobj*.so" "$NU/deps/duneanaobj/" 2>&1 | tail -5
  echo "Done. Review with: cd $NU && git status"
elif [ "$MODE" = "to_local" ]; then
  echo "→ Sync 2x2_mc → local_full"
  rsync -av "$NU/analyze_verified.py" "$NU/analyze_verified_minerva.py" "$NU/analyze_full.py" "$NU/nominal_mc.cc" "$NU/README.md" "$NU/analysis_summary.html" "$WS/" 2>&1 | tail -20
  rsync -av "$NU/deps/duneanaobj/LinkDef.h" "$NU/deps/duneanaobj/LinkDef_full.h" "$WS/deps/duneanaobj/" 2>&1 | tail -5
  echo "Done. local_full updated."
elif [ "$MODE" = "status" ]; then
  echo "Diff local_full ↔ 2x2_mc (ignoring ROOTs):"
  diff -rq --exclude='*.root' --exclude='__pycache__' --exclude='*.pyc' --exclude='test_*.py' --exclude='fix_syntax2.py' --exclude='run_full_*.py' \
    "$WS" "$NU" 2>&1 | head -50 || true
  echo ""
  echo "Usage: $0 [to_nu|to_local|status]"
else
  echo "Unknown mode $MODE — use to_nu, to_local, or status"
  exit 1
fi
