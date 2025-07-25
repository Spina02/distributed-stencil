import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob

# --- Color palette variables ---
ideal_color         = '#1C1F75'     # Dark Navy (for ideal lines)
strong_color        = '#3F6DB5'     # Medium Blue
omp_color           = '#7BB8E8'     # Light Blue
weak_color          = '#B8C8E8'     # Very Light Blue
# strong_eff_color    = '#1C1F75'     # Dark Navy (for efficiency)
# omp_eff_color       = '#2B2E77'     # Navy Purple
computation_color   = '#E8B8C8'     # Light Pink (lighter for computation)
communication_color = '#C8A8E8'     # Light Purple (lighter for communication)
total_color         = '#855988'     # Medium Purple (darker for total)

def read_data(data_folder):
    """Read all CSV files from the data folder"""
    csv_files = glob.glob(os.path.join(data_folder, "*.csv"))
    data = {}
    for file_path in csv_files:
        filename = os.path.basename(file_path)
        test_type = filename.replace('_results.csv', '')
        try:
            df = pd.read_csv(file_path)
            for col in ['Nodes', 'TotalTasks', 'ThreadsPerTask', 'TotalTime', 'ComputationTime', 'CommunicationTime']:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
            data[test_type] = df.dropna()
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}")
    return data

def plot_strong_speedup(data, output_folder):
    if 'strong' in data and not data['strong'].empty:
        strong_data = data['strong']
        grouped = strong_data.groupby('Nodes')['TotalTime'].mean().sort_index()
        nodes = grouped.index.values
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        ideal_speedup = nodes / nodes[0]
        plt.figure(figsize=(9,7))
        plt.plot(nodes, speedup, color=strong_color, marker='o', linewidth=2.5, markersize=8, label='Measured Speedup')  # Azzurro
        plt.plot(nodes, ideal_speedup, color=ideal_color, linestyle='--', linewidth=2, label='Ideal Speedup')  # Giallo limone
        plt.xlabel('Number of nodes', fontsize=18)
        plt.ylabel('Speedup', fontsize=18)
        plt.title('Strong Scaling Speedup (as a function of nodes)', fontsize=24)
        plt.legend(fontsize=18)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xticks(nodes)
        plt.ylim(bottom=0)
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(output_folder, 'speedup_strong.png'), dpi=300, bbox_inches='tight')
        plt.close()

def plot_omp_speedup(data, output_folder):
    if 'omp' in data and not data['omp'].empty:
        omp_data = data['omp']
        grouped = omp_data.groupby('ThreadsPerTask')['TotalTime'].mean().sort_index()
        threads = grouped.index.values
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        ideal_speedup = threads
        plt.figure(figsize=(9,7))
        plt.plot(threads, speedup, color=omp_color, marker='o', linewidth=2.5, markersize=8, label='Measured Speedup')  # Verde acqua
        plt.plot(threads, ideal_speedup, color=ideal_color, linestyle='--', linewidth=2, label='Ideal Speedup')  # Giallo limone
        plt.xlabel('Number of Threads', fontsize=18)
        plt.ylabel('Speedup', fontsize=18)
        plt.title('OpenMP Thread Scaling Speedup', fontsize=24)
        plt.legend(fontsize=18)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xticks(threads)
        plt.ylim(bottom=0)
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(output_folder, 'speedup_omp.png'), dpi=300, bbox_inches='tight')
        plt.close()

def plot_strong_efficiency(data, output_folder):
    if 'strong' in data and not data['strong'].empty:
        strong_data = data['strong']
        grouped = strong_data.groupby('Nodes')['TotalTime'].mean().sort_index()
        nodes = grouped.index.values
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        efficiency = speedup / (nodes / nodes[0])
        plt.figure(figsize=(9,7))
        plt.plot(nodes, efficiency, color=strong_color, marker='o', linewidth=2.5, markersize=8, label='Measured Efficiency')  # Turchese
        plt.axhline(y=1.0, color=ideal_color, linestyle='--', linewidth=2, label='Ideal Efficiency')  # Giallo limone
        plt.xlabel('Number of nodes', fontsize=18)
        plt.ylabel('Efficiency', fontsize=18)
        plt.title('Strong Scaling Efficiency (as a function of nodes)', fontsize=24)
        plt.legend(fontsize=18)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xticks(nodes)
        plt.ylim(0, 1.2)
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(output_folder, 'efficiency_strong.png'), dpi=300, bbox_inches='tight')
        plt.close()

def plot_omp_efficiency(data, output_folder):
    if 'omp' in data and not data['omp'].empty:
        omp_data = data['omp']
        grouped = omp_data.groupby('ThreadsPerTask')['TotalTime'].mean().sort_index()
        threads = grouped.index.values
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        efficiency = speedup / threads
        plt.figure(figsize=(9,7))
        plt.plot(threads, efficiency, color=omp_color, marker='o', linewidth=2.5, markersize=8, label='Measured Efficiency')  # Verde lime
        plt.axhline(y=1.0, color=ideal_color, linestyle='--', linewidth=2, label='Ideal Efficiency')  # Giallo limone
        plt.xlabel('Number of Threads', fontsize=18)
        plt.ylabel('Efficiency', fontsize=18)
        plt.title('OpenMP Thread Scaling Efficiency', fontsize=24)
        plt.legend(fontsize=18)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xticks(threads)
        plt.ylim(0, 1.2)
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(output_folder, 'efficiency_omp.png'), dpi=300, bbox_inches='tight')
        plt.close()

def plot_weak_efficiency(data, output_folder):
    if 'weak' in data and not data['weak'].empty:
        weak_data = data['weak']
        grouped = weak_data.groupby('Nodes')['TotalTime'].mean().sort_index()
        nodes = grouped.index.values
        baseline_time = grouped.iloc[0]
        efficiency = baseline_time / grouped.values
        plt.figure(figsize=(9,7))
        plt.plot(nodes, efficiency, color=weak_color, marker='o', linewidth=2.5, markersize=8, label='Measured Efficiency')  # Corallo
        plt.axhline(y=1.0, color=ideal_color, linestyle='--', linewidth=2, label='Ideal Efficiency')  # Giallo limone
        plt.xlabel('Number of Nodes', fontsize=18)
        plt.ylabel('Efficiency', fontsize=18)
        plt.title('Weak Scaling Efficiency', fontsize=24)
        plt.legend(fontsize=18)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xticks(nodes)
        plt.ylim(0, 1.2)
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(output_folder, 'efficiency_weak.png'), dpi=300, bbox_inches='tight')
        plt.close()

def plot_strong_time(data, output_folder):
    if 'strong' in data and not data['strong'].empty:
        strong_data = data['strong']
        grouped = strong_data.groupby('Nodes').mean(numeric_only=True).sort_index()
        nodes = grouped.index
        plt.figure(figsize=(9,7))
        plt.plot(nodes, grouped['ComputationTime'], color=computation_color, marker='s', linewidth=2, markersize=6, label='Computation Time')  # Azzurro
        plt.plot(nodes, grouped['CommunicationTime'], color=communication_color, marker='^', linewidth=2, markersize=6, label='Communication Time')  # Verde acqua
        plt.plot(nodes, grouped['TotalTime'], color=total_color, marker='o', linewidth=2, markersize=8, label='Total Time')  # Blu cielo
        plt.xlabel('Number of nodes', fontsize=18)
        plt.ylabel('Time (seconds)', fontsize=18)
        plt.title('Strong Scaling - Execution Time (as a function of nodes)', fontsize=24)
        plt.legend(fontsize=18)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xticks(nodes)
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(output_folder, 'time_strong.png'), dpi=300, bbox_inches='tight')
        plt.close()

def plot_omp_time(data, output_folder):
    if 'omp' in data and not data['omp'].empty:
        omp_data = data['omp']
        grouped = omp_data.groupby('ThreadsPerTask').mean(numeric_only=True).sort_index()
        threads = grouped.index
        plt.figure(figsize=(9,7))
        plt.plot(threads, grouped['ComputationTime'], color=computation_color, marker='s', linewidth=2, markersize=6, label='Computation Time')  # Azzurro
        plt.plot(threads, grouped['CommunicationTime'], color=communication_color, marker='^', linewidth=2, markersize=6, label='Communication Time')  # Verde acqua
        plt.plot(threads, grouped['TotalTime'], color=total_color, marker='o', linewidth=2, markersize=8, label='Total Time')  # Blu cielo
        plt.xlabel('Number of Threads', fontsize=18)
        plt.ylabel('Time (seconds)', fontsize=18)
        plt.title('OpenMP Scaling - Execution Time', fontsize=24)
        plt.legend(fontsize=18)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xticks(threads)
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(output_folder, 'time_omp.png'), dpi=300, bbox_inches='tight')
        plt.close()

def plot_weak_time(data, output_folder):
    if 'weak' in data and not data['weak'].empty:
        weak_data = data['weak']
        grouped = weak_data.groupby('Nodes').mean(numeric_only=True).sort_index()
        nodes = grouped.index
        plt.figure(figsize=(9,7))
        plt.plot(nodes, grouped['ComputationTime'], color=computation_color, marker='s', linewidth=2, markersize=6, label='Computation Time')  # Azzurro
        plt.plot(nodes, grouped['CommunicationTime'], color=communication_color, marker='^', linewidth=2, markersize=6, label='Communication Time')  # Verde acqua
        plt.plot(nodes, grouped['TotalTime'], color=total_color, marker='o', linewidth=2, markersize=8, label='Total Time')  # Blu cielo
        plt.xlabel('Number of Nodes', fontsize=18)
        plt.ylabel('Time (seconds)', fontsize=18)
        plt.title('Weak Scaling - Execution Time', fontsize=24)
        plt.legend(fontsize=18)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xticks(nodes)
        plt.tight_layout(pad=2.0)
        plt.savefig(os.path.join(output_folder, 'time_weak.png'), dpi=300, bbox_inches='tight')
        plt.close()

def plot_strong_overlap_comparison(overlap_path, no_overlap_path, output_folder):
    """Plot TotalTime vs Nodes for strong scaling with and without overlap."""
    import pandas as pd
    import os
    import matplotlib.pyplot as plt
    # Read both CSVs
    df_overlap = pd.read_csv(overlap_path)
    df_no_overlap = pd.read_csv(no_overlap_path)
    # Ensure numeric
    for df in [df_overlap, df_no_overlap]:
        for col in ['Nodes', 'TotalTime']:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
    # Group by Nodes and average if needed
    overlap_grouped = df_overlap.groupby('Nodes')['TotalTime'].mean().sort_index()
    no_overlap_grouped = df_no_overlap.groupby('Nodes')['TotalTime'].mean().sort_index()
    nodes = sorted(set(overlap_grouped.index).intersection(no_overlap_grouped.index))
    plt.figure(figsize=(8, 6))
    plt.plot(nodes, [overlap_grouped[n] for n in nodes], marker='o', linewidth=2.5, label='With Overlap')
    plt.plot(nodes, [no_overlap_grouped[n] for n in nodes], marker='o', linewidth=2.5, label='No Overlap')
    plt.xlabel('Number of Nodes', fontsize=18)
    plt.ylabel('Total Execution Time (s)', fontsize=18)
    plt.title('Strong Scaling: Overlap vs No Overlap', fontsize=24)
    plt.legend(fontsize=18)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xticks(nodes)
    plt.tight_layout(pad=2.0)
    plt.savefig(os.path.join(output_folder, 'strong_overlap_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot the results of the parallel stencil code')
    parser.add_argument('--data-folder', type=str, default='data', help='Path to the data folder')
    parser.add_argument('--output-folder', type=str, default='plots', help='Path to the output folder')
    parser.add_argument('--plot-type', type=str, nargs='+', default=['all'], help='Type(s) of plot to generate', choices=['speedup', 'efficiency', 'time', 'overlap', 'all'])
    parser.add_argument('--strong-overlap', type=str, default='data/strong_results', help='Path to strong_results.csv (with overlap)')
    parser.add_argument('--strong-no-overlap', type=str, default='data/strong_results_no_overlap.csv', help='Path to strong_results_no_overlap.csv (no overlap)')
    args = parser.parse_args()
    os.makedirs(args.output_folder, exist_ok=True)
    plot_types = args.plot_type
    data = read_data(args.data_folder)
    if 'all' in plot_types or 'speedup' in plot_types:
        print("Generating strong speedup plot...")
        plot_strong_speedup(data, args.output_folder)
        print("Generating OpenMP speedup plot...")
        plot_omp_speedup(data, args.output_folder)
    if 'all' in plot_types or 'efficiency' in plot_types:
        print("Generating strong efficiency plot...")
        plot_strong_efficiency(data, args.output_folder)
        print("Generating OpenMP efficiency plot...")
        plot_omp_efficiency(data, args.output_folder)
        print("Generating weak efficiency plot...")
        plot_weak_efficiency(data, args.output_folder)
    if 'all' in plot_types or 'time' in plot_types:
        print("Generating strong time plot...")
        plot_strong_time(data, args.output_folder)
        print("Generating OpenMP time plot...")
        plot_omp_time(data, args.output_folder)
        print("Generating weak time plot...")
        plot_weak_time(data, args.output_folder)
    if 'all' in plot_types or 'overlap' in plot_types:
        if args.strong_overlap and args.strong_no_overlap:
            print("Generating strong overlap comparison plot...")
            plot_strong_overlap_comparison(args.strong_overlap, args.strong_no_overlap, args.output_folder)
        else:
            print("Error: --strong-overlap and --strong-no-overlap must be provided for overlap plot.")

if __name__ == "__main__":
    main()