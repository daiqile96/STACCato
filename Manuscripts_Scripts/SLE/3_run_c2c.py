import pandas as pd
import scanpy as sc
import plotnine as p9
import matplotlib.pyplot as plt

import liana as li
import cell2cell as c2c
import decoupler as dc  # needed for pathway enrichment

import warnings
warnings.filterwarnings('ignore')
from collections import defaultdict
import pickle
import os
import sys

# ============================================================================
# PARSE COMMAND-LINE ARGUMENTS
# ============================================================================
if len(sys.argv) != 3:
    print("Usage: python script.py <data_dir> <output_dir>")
    sys.exit(1)

data_dir = sys.argv[1]
output_dir = sys.argv[2]

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

print("="*80)
print("CONFIGURATION")
print("="*80)
print(f"Data directory: {data_dir}")
print(f"Output directory: {output_dir}")
print("="*80 + "\n")

# ============================================================================
# LOAD TENSOR AND METADATA
# ============================================================================
tensor_path = data_dir + "/Tensor.pkl"
meta_path = data_dir + "/Tensor-Metadata.pkl"

print("Loading tensor and metadata...")
with open(tensor_path, "rb") as f:
    tensor = pickle.load(f)

with open(meta_path, "rb") as f:
    meta_tensor = pickle.load(f)

print("Metadata:")
print(meta_tensor)

# ============================================================================
# RUN TENSOR DECOMPOSITION WITH AUTOMATIC RANK SELECTION
# ============================================================================
print("\n" + "="*80)
print("Running Tensor-Cell2Cell decomposition...")
print("="*80)

tensor = c2c.analysis.run_tensor_cell2cell_pipeline(
    tensor,
    meta_tensor,
    upper_rank=15,  # Maximum rank to test in elbow analysis
    rank=None,  # Automatically determine optimal rank
    tf_optimization='regular',  # 'regular' or 'robust'
    random_state=0,  # Random seed for reproducibility
    device='cpu',
    elbow_metric='error',  # Metric for elbow analysis
    smooth_elbow=False,  # Whether to smooth the elbow curve
)

# ============================================================================
# SAVE OPTIMIZED RANK AND ELBOW PLOT
# ============================================================================
print(f"\n{'='*80}")
print(f"OPTIMAL RANK SELECTED: {tensor.rank}")
print(f"{'='*80}\n")

# Save the optimal rank to a text file
with open(f"{output_dir}/optimal_rank.txt", "w") as f:
    f.write(f"Optimal Rank: {tensor.rank}\n")
    f.write(f"Upper Rank Tested: 15\n")
    f.write(f"Optimization Method: regular\n")
    f.write(f"Random State: 0\n")

print(f"✓ Optimal rank saved to: {output_dir}/optimal_rank.txt")

# Save elbow plot if available
if hasattr(tensor, 'elbow_metric_raw'):
    plt.figure(figsize=(10, 6))
    
    # Plot elbow curve
    if tensor.elbow_metric_raw.ndim == 2:
        # Multiple runs - plot mean and std
        mean_error = tensor.elbow_metric_raw.mean(axis=1)
        std_error = tensor.elbow_metric_raw.std(axis=1)
        ranks = range(1, len(mean_error) + 1)
        
        plt.plot(ranks, mean_error, 'b-', linewidth=2, label='Mean Error')
        plt.fill_between(ranks, mean_error - std_error, mean_error + std_error, 
                         alpha=0.2, color='blue', label='±1 SD')
    else:
        # Single run
        ranks = range(1, len(tensor.elbow_metric_raw) + 1)
        plt.plot(ranks, tensor.elbow_metric_raw, 'b-', linewidth=2)
    
    # Mark the selected rank
    plt.axvline(x=tensor.rank, color='red', linestyle='--', linewidth=2, 
                label=f'Selected Rank = {tensor.rank}')
    
    plt.xlabel('Rank (Number of Factors)', fontsize=12)
    plt.ylabel('Reconstruction Error', fontsize=12)
    plt.title('Elbow Analysis for Optimal Rank Selection', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    elbow_plot_path = f"{output_dir}/elbow_plot.png"
    plt.savefig(elbow_plot_path, dpi=300, bbox_inches='tight')
    print(f"✓ Elbow plot saved to: {elbow_plot_path}")
    plt.close()

# ============================================================================
# SAVE TENSOR DECOMPOSITION RESULTS
# ============================================================================
print("\nSaving tensor decomposition results...")

# Save the complete tensor object
tensor_result_path = f"{output_dir}/tensor_decomposed.pkl"
with open(tensor_result_path, "wb") as f:
    pickle.dump(tensor, f)
print(f"✓ Complete tensor object saved to: {tensor_result_path}")
