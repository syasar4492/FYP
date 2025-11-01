# Flow cytometry: Size vs fluorescence scatter from an FCS file

from pathlib import Path
import matplotlib.pyplot as plt
import fcsparser

# ---- 1) Set FCS file path (raw string avoids \0 issues on Windows) ----
path = Path(r"D:\ICL Module Notes and more\Year 4\FYP\Software Automation\Raw Flow Data\06.11.2025_TE varying NaCl concs_EFG\0.5M NaCl_F2.fcs")

# ---- 2) Load FCS ----
meta, data = fcsparser.parse(str(path), reformat_meta=True)
print("Channels:", list(data.columns))

# ---- 3) Robust channel picking helpers ----
def first_match(columns, patterns):
    """
    Return the first column whose uppercase name contains any of the given
    uppercase pattern substrings, in order of patterns.
    """
    cols = list(columns)
    up = [c.upper() for c in cols]
    for pat in patterns:
        pat = pat.upper()
        for i, cu in enumerate(up):
            if pat in cu:
                return cols[i]
    return None

# Prefer area signals for scatter
fsc = first_match(data.columns, ["FSC-A", "FSC"])
if fsc is None:
    raise ValueError("Could not find an FSC channel. Found: " + str(list(data.columns)))

# Common aliases for FITC/GFP on blue 488 nm laser (adjust/extend as needed)
fluor_priority = [
    "FITC", "GFP",
    "B525-A", "B525", "B530", "B515", "530/30", "525/50",
    "BL1-A", "BL1",
    "FL1-A", "FL1",
]

fluor = first_match(data.columns, fluor_priority)

# If still not found, fall back to any '-A' fluorescence channel (excluding FSC/SSC/Time/Width)
if fluor is None:
    area_candidates = [
        c for c in data.columns
        if c.upper().endswith("-A")
        and not any(key in c.upper() for key in ["FSC", "SSC", "TIME", "WIDTH"])
    ]
    if area_candidates:
        fluor = area_candidates[0]

if fluor is None:
    raise ValueError(
        "Could not find a fluorescence channel. "
        f"Available columns: {list(data.columns)}"
    )

print(f"Using X={fsc}, Y={fluor}")

# ---- 4) Plot ----
plt.scatter(data[fsc], data[fluor], s=2)
plt.xlabel(fsc)
plt.ylabel(fluor)
plt.title(f"{path.name}: {fsc} vs {fluor}")
plt.tight_layout()
plt.show()
