"""Bulut dağıtımı giriş noktası: ana arayüzü (../app.py) küçültülmüş veriyle çalıştırır.

- Koşu sonuçları: deploy/runs/ (deploy/make_deploy_data.py ile üretilir, ~9 MB, depoda)
- Canlı mod: 2000 test parçacığı (bulut makinesinde akıcı kalsın diye)

Çalıştırma:  streamlit run deploy/app.py
Streamlit Community Cloud'da "Main file path" olarak deploy/app.py verilir.
"""

import os
import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

DEPLOY_LIVE_N = 2000

os.environ.setdefault("MW_RUNS_DIR", str(HERE / "runs"))
os.environ.setdefault("MW_LIVE_N", str(DEPLOY_LIVE_N))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

runpy.run_path(str(ROOT / "app.py"), run_name="__main__")
