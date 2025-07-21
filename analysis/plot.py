import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob

def read_data(data_folder):
    """Read all CSV files from the data folder"""
    csv_files = glob.glob(os.path.join(data_folder, "*.csv"))
    data = {}
    
    for file_path in csv_files:
        filename = os.path.basename(file_path)
        test_type = filename.replace('_results.csv', '')
        try:
            df = pd.read_csv(file_path)
            data[test_type] = df
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}")
    
    return data

def plot_speedup(data_folder, output_folder):
    """Plot speedup for strong scaling and OpenMP scaling"""
    data = read_data(data_folder)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Strong scaling speedup
    if 'strong' in data:
        strong_data = data['strong']
        grouped = strong_data.groupby('TotalTasks')['TotalTime'].mean().sort_index()
        processors = grouped.index
        p_base = processors[0]
        nodes = processors / p_base
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        ideal_speedup = nodes  # Ideal speedup = nodi
        ax1.plot(nodes, speedup, 'bo-', linewidth=2, markersize=8, label='Measured Speedup')
        ax1.plot(nodes, ideal_speedup, 'r--', linewidth=2, label='Ideal Speedup')
        ax1.set_xlabel('Number of nodes')
        ax1.set_ylabel('Speedup')
        ax1.set_title('Strong Scaling Speedup (as a function of nodes)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xticks(nodes)
        ax1.set_yticks(np.arange(0, max(ideal_speedup) * 1.2, 2))
    
    # OpenMP scaling speedup
    if 'omp' in data:
        omp_data = data['omp']
        # CORREZIONE: Ordina SEMPRE i dati raggruppati
        grouped = omp_data.groupby('ThreadsPerTask')['TotalTime'].mean().sort_index()
        
        threads = grouped.index
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        
        ax2.plot(threads, speedup, 'go-', linewidth=2, markersize=8, label='Measured Speedup')
        ax2.plot(threads, threads, 'r--', linewidth=2, label='Ideal Speedup')
        ax2.set_xlabel('Number of Threads')
        ax2.set_ylabel('Speedup')
        ax2.set_title('OpenMP Thread Scaling Speedup')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xticks(threads)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, 'speedup_plot.png'), dpi=300, bbox_inches='tight')

def plot_efficiency(data_folder, output_folder):
    """Plot efficiency for strong, weak, and OpenMP scaling"""
    data = read_data(data_folder)
    
    # AGGIUNTA: Un subplot in più per l'efficienza di scaling debole
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22, 6))
    
    # Strong scaling efficiency
    if 'strong' in data:
        strong_data = data['strong']
        grouped = strong_data.groupby('TotalTasks')['TotalTime'].mean().sort_index()
        processors = grouped.index
        p_base = processors[0]
        nodes = processors / p_base
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        efficiency = speedup / nodes
        ax1.plot(nodes, efficiency, 'bo-', linewidth=2, markersize=8, label='Measured Efficiency')
        ax1.axhline(y=1.0, color='r', linestyle='--', linewidth=2, label='Ideal Efficiency')
        ax1.set_xlabel('Number of nodes')
        ax1.set_ylabel('Efficiency')
        ax1.set_title('Strong Scaling Efficiency (as a function of nodes)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xticks(nodes)
        ax1.set_ylim(0, 1.2)

    # OpenMP scaling efficiency
    if 'omp' in data:
        omp_data = data['omp']
        # CORREZIONE: Ordina SEMPRE i dati raggruppati
        grouped = omp_data.groupby('ThreadsPerTask')['TotalTime'].mean().sort_index()
        
        threads = grouped.index
        baseline_time = grouped.iloc[0]
        
        speedup = baseline_time / grouped.values
        efficiency = speedup / threads
        
        ax2.plot(threads, efficiency, 'go-', linewidth=2, markersize=8, label='Measured Efficiency')
        ax2.axhline(y=1.0, color='r', linestyle='--', linewidth=2, label='Ideal Efficiency')
        ax2.set_xlabel('Number of Threads')
        ax2.set_ylabel('Efficiency')
        ax2.set_title('OpenMP Thread Scaling Efficiency')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xticks(threads)
        ax2.set_ylim(0, 1.2)

    # NUOVO E CONSIGLIATO: Weak scaling efficiency
    if 'weak' in data:
        weak_data = data['weak']
        grouped = weak_data.groupby('TotalTasks')['TotalTime'].mean().sort_index()
        
        processors = grouped.index
        baseline_time = grouped.iloc[0]
        
        # L'efficienza nello scaling debole è T_base / T_p
        efficiency = baseline_time / grouped.values
        
        ax3.plot(processors, efficiency, 'mo-', linewidth=2, markersize=8, label='Measured Efficiency')
        ax3.axhline(y=1.0, color='r', linestyle='--', linewidth=2, label='Ideal Efficiency')
        ax3.set_xlabel('Number of Tasks (Processes)')
        ax3.set_ylabel('Efficiency')
        ax3.set_title('Weak Scaling Efficiency')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        ax3.set_xticks(processors)
        ax3.set_ylim(0, 1.2)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, 'efficiency_plot.png'), dpi=300, bbox_inches='tight')

def plot_time(data_folder, output_folder):
    """Plot execution time for all scaling tests"""
    data = read_data(data_folder)
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    # Strong scaling time
    if 'strong' in data:
        strong_data = data['strong']
        grouped = strong_data.groupby('TotalTasks').mean(numeric_only=True).sort_index()
        processors = grouped.index
        p_base = processors[0]
        nodes = processors / p_base
        ax1.plot(nodes, grouped['ComputationTime'], 'bs-', linewidth=2, markersize=6, label='Computation Time')
        ax1.plot(nodes, grouped['CommunicationTime'], 'r^-', linewidth=2, markersize=6, label='Communication Time')
        ax1.plot(nodes, grouped['TotalTime'], 'go-', linewidth=2, markersize=8, label='Total Time')
        ax1.set_xlabel('Number of nodes')
        ax1.set_ylabel('Time (seconds)')
        ax1.set_title('Strong Scaling - Execution Time (as a function of nodes)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xticks(nodes)
    
    # OpenMP scaling time
    if 'omp' in data:
        omp_data = data['omp']
        # CORREZIONE: Aggiunto numeric_only=True
        grouped = omp_data.groupby('ThreadsPerTask').mean(numeric_only=True).sort_index()
        
        ax2.plot(grouped.index, grouped['ComputationTime'], 'bs-', linewidth=2, markersize=6, label='Computation Time')
        ax2.plot(grouped.index, grouped['CommunicationTime'], 'r^-', linewidth=2, markersize=6, label='Communication Time')
        ax2.plot(grouped.index, grouped['TotalTime'], 'go-', linewidth=2, markersize=8, label='Total Time')
        ax2.set_xlabel('Number of Threads')
        ax2.set_ylabel('Time (seconds)')
        ax2.set_title('OpenMP Scaling - Execution Time')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xticks(grouped.index)
    
    # Weak scaling time
    if 'weak' in data:
        weak_data = data['weak']
        # CORREZIONE: Aggiunto numeric_only=True
        grouped = weak_data.groupby('TotalTasks').mean(numeric_only=True).sort_index()
        
        ax3.plot(grouped.index, grouped['ComputationTime'], 'bs-', linewidth=2, markersize=6, label='Computation Time')
        ax3.plot(grouped.index, grouped['CommunicationTime'], 'r^-', linewidth=2, markersize=6, label='Communication Time')
        ax3.plot(grouped.index, grouped['TotalTime'], 'go-', linewidth=2, markersize=8, label='Total Time')
        ax3.set_xlabel('Number of Tasks')
        ax3.set_ylabel('Time (seconds)')
        ax3.set_title('Weak Scaling - Execution Time')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        ax3.set_xticks(grouped.index)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, 'time_plot.png'), dpi=300, bbox_inches='tight')
    # plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot the results of the parallel stencil code')
    parser.add_argument('--data-folder', type=str, default='data', help='Path to the data folder')
    parser.add_argument('--output-folder', type=str, default='plots', help='Path to the output folder')
    parser.add_argument('--plot-type', type=str, nargs='+', default=['all'], help='Type(s) of plot to generate', choices=['speedup', 'efficiency', 'time', 'all'])
    
    args = parser.parse_args()
    
    os.makedirs(args.output_folder, exist_ok=True)
    
    types = args.plot_type
    
    if 'all' in types or 'speedup' in types:
        plot_speedup(args.data_folder, args.output_folder)
    if 'all' in types or 'efficiency' in types:
        plot_efficiency(args.data_folder, args.output_folder)
    if 'all' in types or 'time' in types:
        plot_time(args.data_folder, args.output_folder)

if __name__ == "__main__":
    main()