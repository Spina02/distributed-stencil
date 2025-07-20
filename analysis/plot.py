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
        # Group by TaskTotali and calculate average time
        grouped = strong_data.groupby('TaskTotali')['TempoTotale'].mean()
        
        # Calculate speedup (baseline is single task)
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        processors = grouped.index
        
        ax1.plot(processors, speedup, 'bo-', linewidth=2, markersize=8, label='Measured Speedup')
        ax1.plot(processors, processors, 'r--', linewidth=2, label='Ideal Speedup')
        ax1.set_xlabel('Number of Tasks')
        ax1.set_ylabel('Speedup')
        ax1.set_title('Strong Scaling Speedup')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
    
    # OpenMP scaling speedup
    if 'omp' in data:
        omp_data = data['omp']
        # Group by ThreadPerTask and calculate average time
        grouped = omp_data.groupby('ThreadPerTask')['TempoTotale'].mean()
        
        # Calculate speedup (baseline is single thread)
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        threads = grouped.index
        
        ax2.plot(threads, speedup, 'go-', linewidth=2, markersize=8, label='Measured Speedup')
        ax2.plot(threads, threads, 'r--', linewidth=2, label='Ideal Speedup')
        ax2.set_xlabel('Number of Threads')
        ax2.set_ylabel('Speedup')
        ax2.set_title('OpenMP Thread Scaling Speedup')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, 'speedup_plot.png'), dpi=300, bbox_inches='tight')
    plt.show()

def plot_efficiency(data_folder, output_folder):
    """Plot efficiency for strong scaling and OpenMP scaling"""
    data = read_data(data_folder)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Strong scaling efficiency
    if 'strong' in data:
        strong_data = data['strong']
        # Group by TaskTotali and calculate average time
        grouped = strong_data.groupby('TaskTotali')['TempoTotale'].mean()
        
        # Calculate efficiency
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        processors = grouped.index
        efficiency = speedup / processors
        
        ax1.plot(processors, efficiency, 'bo-', linewidth=2, markersize=8, label='Measured Efficiency')
        ax1.axhline(y=1.0, color='r', linestyle='--', linewidth=2, label='Ideal Efficiency')
        ax1.set_xlabel('Number of Tasks')
        ax1.set_ylabel('Efficiency')
        ax1.set_title('Strong Scaling Efficiency')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xscale('log')
    
    # OpenMP scaling efficiency
    if 'omp' in data:
        omp_data = data['omp']
        # Group by ThreadPerTask and calculate average time
        grouped = omp_data.groupby('ThreadPerTask')['TempoTotale'].mean()
        
        # Calculate efficiency
        baseline_time = grouped.iloc[0]
        speedup = baseline_time / grouped.values
        threads = grouped.index
        efficiency = speedup / threads
        
        ax2.plot(threads, efficiency, 'go-', linewidth=2, markersize=8, label='Measured Efficiency')
        ax2.axhline(y=1.0, color='r', linestyle='--', linewidth=2, label='Ideal Efficiency')
        ax2.set_xlabel('Number of Threads')
        ax2.set_ylabel('Efficiency')
        ax2.set_title('OpenMP Thread Scaling Efficiency')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, 'efficiency_plot.png'), dpi=300, bbox_inches='tight')
    plt.show()

def plot_time(data_folder, output_folder):
    """Plot execution time for all scaling tests"""
    data = read_data(data_folder)
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    # Strong scaling time
    if 'strong' in data:
        strong_data = data['strong']
        # Group by TaskTotali and calculate average time
        grouped = strong_data.groupby('TaskTotali')['TempoTotale'].mean()
        
        ax1.plot(grouped.index, grouped.values, 'bo-', linewidth=2, markersize=8, label='Total Time')
        ax1.set_xlabel('Number of Tasks')
        ax1.set_ylabel('Time (seconds)')
        ax1.set_title('Strong Scaling - Execution Time')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
    
    # OpenMP scaling time
    if 'omp' in data:
        omp_data = data['omp']
        # Group by ThreadPerTask and calculate average time
        grouped = omp_data.groupby('ThreadPerTask')['TempoTotale'].mean()
        
        ax2.plot(grouped.index, grouped.values, 'go-', linewidth=2, markersize=8, label='Total Time')
        ax2.set_xlabel('Number of Threads')
        ax2.set_ylabel('Time (seconds)')
        ax2.set_title('OpenMP Scaling - Execution Time')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
    
    # Weak scaling time
    if 'weak' in data:
        weak_data = data['weak']
        # Group by TaskTotali and calculate average time
        grouped = weak_data.groupby('TaskTotali')['TempoTotale'].mean()
        
        ax3.plot(grouped.index, grouped.values, 'ro-', linewidth=2, markersize=8, label='Total Time')
        ax3.set_xlabel('Number of Tasks')
        ax3.set_ylabel('Time (seconds)')
        ax3.set_title('Weak Scaling - Execution Time')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        ax3.set_xscale('log')
        ax3.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, 'time_plot.png'), dpi=300, bbox_inches='tight')
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot the results of the parallel stencil code')
    parser.add_argument('--data-folder', type=str, default='data/local_results', help='Path to the data folder')
    parser.add_argument('--output-folder', type=str, default='plots', help='Path to the output folder')
    parser.add_argument('--plot-type', type=str, nargs='+', default=['all'], help='Type(s) of plot to generate', choices=['speedup', 'efficiency', 'time', 'all'])
    
    args = parser.parse_args()
    
    # Create output folder if it doesn't exist
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