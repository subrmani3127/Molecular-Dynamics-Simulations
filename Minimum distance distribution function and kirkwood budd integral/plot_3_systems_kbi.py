import json
import matplotlib.pyplot as plt

# --- Define your JSON files and labels/colors ---
json_files = {

    "C1":       {"file": "PR_PB_41_results.json",           "color": "black"},
    "C2":       {"file": "PR_PB_130_results.json",          "color": "red"},
    "C3":       {"file": "PR_PB_434_results.json",          "color": "blue"},
}

# --- Load KBI data ---
kbi_data = {}

for label, info in json_files.items():
    with open(info["file"], 'r') as f:
        data = json.load(f)
        kbi_data[label] = {
            "distance": data["d"],
            "kbi": data["kb"],
            "color": info["color"]
        }

# --- Plot ---
plt.figure(figsize=(10, 6))

for label, dataset in kbi_data.items():
    plt.plot(dataset["distance"], dataset["kbi"], label=label,
             linewidth=2, color=dataset["color"])

# --- Style the plot ---
plt.xlabel("Distance (Å)", fontsize=16, fontweight='bold')
plt.ylabel("Kirkwood-Buff Integral (cm³/mol)", fontsize=14, fontweight='bold')
#plt.title("KBI vs Distance", fontsize=20, fontweight='bold')
plt.legend(fontsize=10)
plt.xticks(fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')

# Set border width
for spine in plt.gca().spines.values():
    spine.set_linewidth(2.0)

plt.grid(False)
plt.tight_layout()
plt.savefig("kbi1.svg",format="svg",bbox_inches="tight",pad_inches=0.05)

plt.show()
