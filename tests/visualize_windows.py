#!/usr/bin/env python3
"""
Visualize window selection results from Fortran test program
"""
import numpy as np
import matplotlib.pyplot as plt
from obspy import read

def read_window_results(filename='window_results.txt'):
    """Read window selection results"""
    windows = []
    metadata = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Number of windows:'):
                metadata['n_windows'] = int(line.split(':')[1].strip())
            elif line.startswith('# tstart:'):
                metadata['tstart'] = float(line.split(':')[1].strip())
            elif line.startswith('# tend:'):
                metadata['tend'] = float(line.split(':')[1].strip())
            elif line.startswith('# noise_level:'):
                metadata['noise_level'] = float(line.split(':')[1].strip())
            elif not line.startswith('#') and line:
                parts = line.split()
                windows.append({
                    'number': int(parts[0]),
                    'start': float(parts[1]),
                    'end': float(parts[2]),
                    'duration': float(parts[3]),
                    'avg_cc': float(parts[4]),
                    'avg_shift': float(parts[5])
                })
    
    return windows, metadata

def read_cc_shift_data(filename='window_cc_shift.txt'):
    """Read cross-correlation and time shift data"""
    data = np.loadtxt(filename, comments='#')
    return {
        'time': data[:, 0],
        'cc': data[:, 1],
        'time_shift': data[:, 2]
    }

def plot_window_selection(obs_file='../example_data/obs.II.ABKT.LHZ.sac',
                          syn_file='../example_data/syn.II.ABKT.LHZ.sac',
                          window_file='window_results.txt',
                          cc_shift_file='window_cc_shift.txt'):
    """Plot window selection results"""
    
    # Read SAC files
    obs = read(obs_file)[0]
    syn = read(syn_file)[0]
    obs.filter('bandpass', freqmin=1/150, freqmax=1/50, corners=2, zerophase=True)
    syn.filter('bandpass', freqmin=1/150, freqmax=1/50, corners=2, zerophase=True)

    # Read window results
    windows, metadata = read_window_results(window_file)
    cc_data = read_cc_shift_data(cc_shift_file)
    
    # Create figure
    fig, axes = plt.subplots(4, 1, figsize=(14, 12))
    
    # Plot 1: Waveforms with windows
    ax = axes[0]
    times = obs.times()
    ax.plot(times, obs.data, 'k', alpha=0.7, linewidth=0.8, label='Observed')
    ax.plot(times, syn.data, 'r--', alpha=0.7, linewidth=0.8, label='Synthetic')
    
    # Plot time window bounds
    ax.axvline(x=metadata['tstart'], color='b', linestyle='--', linewidth=1.5, 
               alpha=0.5, label='tstart')
    ax.axvline(x=metadata['tend'], color='m', linestyle='--', linewidth=1.5, 
               alpha=0.5, label='tend')
    
    # Highlight selected windows
    for win in windows:
        ax.axvspan(win['start'], win['end'], alpha=0.3, color='green')
    
    ax.set_xlim(0, cc_data['time'][-1])
    ax.set_ylabel('Amplitude')
    ax.legend(loc='upper right')
    ax.set_title('Waveforms with Selected Windows (green shading)')
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Cross-correlation
    ax = axes[1]
    ax.plot(cc_data['time'], cc_data['cc'], 'k', linewidth=1, alpha=0.7, 
            label='Cross-correlation')
    ax.axhline(y=0.8, color='r', linestyle='--', linewidth=2, alpha=0.5, 
               label='CC Threshold = 0.8')
    
    # Highlight windows
    for win in windows:
        win_mask = (cc_data['time'] >= win['start']) & (cc_data['time'] <= win['end'])
        ax.scatter(cc_data['time'][win_mask], cc_data['cc'][win_mask], 
                  c='green', s=10, alpha=0.6, zorder=5)
    
    ax.set_xlim(0, cc_data['time'][-1])
    ax.set_ylabel('Cross-correlation')
    ax.legend(loc='upper right')
    ax.set_title('Cross-correlation Coefficient')
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Time shift
    ax = axes[2]
    ax.plot(cc_data['time'], cc_data['time_shift'], 'b', linewidth=1, alpha=0.7, 
            label='Time shift')
    ax.axhline(y=15, color='r', linestyle='--', linewidth=2, alpha=0.5)
    ax.axhline(y=-15, color='r', linestyle='--', linewidth=2, alpha=0.5, 
               label='Shift Threshold = ±15s')
    
    # Highlight windows
    for win in windows:
        win_mask = (cc_data['time'] >= win['start']) & (cc_data['time'] <= win['end'])
        ax.scatter(cc_data['time'][win_mask], cc_data['time_shift'][win_mask], 
                  c='green', s=10, alpha=0.6, zorder=5)
    
    ax.set_xlim(0, cc_data['time'][-1])
    ax.set_ylabel('Time shift (s)')
    ax.legend(loc='upper right')
    ax.set_title('Time Shift')
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Window statistics
    ax = axes[3]
    if windows:
        win_numbers = [w['number'] for w in windows]
        durations = [w['duration'] for w in windows]
        avg_ccs = [w['avg_cc'] for w in windows]
        
        ax2 = ax.twinx()
        bars = ax.bar(win_numbers, durations, alpha=0.6, color='green', label='Duration')
        line = ax2.plot(win_numbers, avg_ccs, 'ro-', linewidth=2, markersize=8, 
                       label='Avg CC')
        
        ax.set_xlabel('Window Number')
        ax.set_ylabel('Duration (s)', color='green')
        ax2.set_ylabel('Average CC', color='red')
        ax.tick_params(axis='y', labelcolor='green')
        ax2.tick_params(axis='y', labelcolor='red')
        ax.set_title('Window Statistics')
        ax.grid(True, alpha=0.3)
        
        # Combined legend
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
    
    plt.tight_layout()
    plt.savefig('window_selection_results.png', dpi=150, bbox_inches='tight')
    print("Figure saved as: window_selection_results.png")
    plt.show()
    
    # Print summary
    print("\n" + "="*70)
    print("Window Selection Summary")
    print("="*70)
    print(f"Total windows found: {len(windows)}")
    print(f"Time window: {metadata['tstart']:.2f} - {metadata['tend']:.2f} s")
    print(f"Noise level: {metadata['noise_level']:.4f}")
    
    if windows:
        total_duration = sum(w['duration'] for w in windows)
        coverage = total_duration / (metadata['tend'] - metadata['tstart']) * 100
        print(f"Total selected duration: {total_duration:.2f} s")
        print(f"Coverage: {coverage:.1f}%")
        print(f"\nAverage CC: {np.mean([w['avg_cc'] for w in windows]):.3f}")
        print(f"Average time shift: {np.mean([w['avg_shift'] for w in windows]):.3f} s")
    print("="*70)

if __name__ == '__main__':
    plot_window_selection()
