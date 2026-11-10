from pathlib import Path
from datetime import datetime, timezone
import hashlib
import requests

ROOT = Path(__file__).resolve().parents[1]
URL = "https://bigg.ucsd.edu/static/models/e_coli_core.xml"
OUT = ROOT / "data" / "e_coli_core.xml"

def main():
    response = requests.get(URL, timeout=60)
    response.raise_for_status()
    if b"<sbml" not in response.content[:5000].lower():
        raise RuntimeError("Downloaded content does not look like an SBML XML file.")
    OUT.write_bytes(response.content)
    digest = hashlib.sha256(response.content).hexdigest()
    metadata = (
        f"URL: {URL}\n"
        f"Retrieved UTC: {datetime.now(timezone.utc).isoformat()}\n"
        f"SHA256: {digest}\n"
        f"Bytes: {len(response.content)}\n"
    )
    (ROOT / "data" / "e_coli_core_download_metadata.txt").write_text(metadata)
    print(f"Saved {OUT}")
    print(metadata)

if __name__ == "__main__":
    main()
