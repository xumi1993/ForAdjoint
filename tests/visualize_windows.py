import numpy as np
import obspy
import pyflex
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import subprocess

def run_fortran_window_selection(obs_file, syn_file, min_period, max_period, is_split_phases):
    # Prepare command
    fstr_is_split_phases = ".true." if is_split_phases else ".false."
    cmd = [
        "../bin/xex_win_sel",
        obs_file,
        syn_file,
        str(min_period),
        str(max_period), 
        fstr_is_split_phases
    ]
    # Run Fortran executable
    subprocess.run(cmd, check=True)

def read_fortran_windows(filename):
    windows = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip().startswith('#') or not line.strip():
                continue
            parts = line.split()
            # Format: window_number  start_time(s)  end_time(s)  duration(s)  cc  shift(s)  dlnA
            if len(parts) >= 7:
                win = {
                    'id': int(parts[0]),
                    'start': float(parts[1]),
                    'end': float(parts[2]),
                    'cc': float(parts[4]),
                    'shift': float(parts[5]),
                    'dlnA': float(parts[6])
                }
                windows.append(win)
    return windows

def run_pyflex(obs, syn, min_period, max_period, is_split_phases):
    # Config matching Fortran as closely as possible
    resolution_strategy = 'interval_scheduling' if is_split_phases else 'merge'
    config = pyflex.Config(
        min_period=min_period,
        max_period=max_period,
        min_surface_wave_velocity=2.4, # Match Fortran min_velocity
        stalta_waterlevel=0.08,        # Match Fortran water_level
        tshift_acceptance_level=12.0,  # Match Fortran 0.3 * 40
        dlna_acceptance_level=1.3,     # Match Fortran
        cc_acceptance_level=0.8,       # Match Fortran
        s2n_limit=1.5,                 # Revert to default
        c_0=1.0,                       # Fortran c_0 default is 1.0
        c_1=1.5,                       # Fortran c_1 default is 1.5
        c_2=0.0,                       # Fortran c_2 default is 0.0
        c_3a=4.0,                      # Fortran c_3a default is 4.0
        c_3b=2.5,                      # Fortran c_3b default is 2.5
        c_4a=2.0,                      # Fortran c_4a default is 2.0
        c_4b=10.0,                     # Match Fortran updated c_4b
        resolution_strategy = resolution_strategy,
    )
    
    # Select windows
    windows = pyflex.select_windows(obs, syn, config, plot=False)
    
    print("\nPyflex Windows:")
    print(f"{'Start(s)':>10} {'End(s)':>10} {'Duration':>10} {'CC':>8} {'Shift':>8}")
    print("-" * 50)
    for w in windows:
        # Calculate relative start/end times
        start_time = w.left * w.dt
        end_time = w.right * w.dt
        print(f"{start_time:>10.2f} {end_time:>10.2f} {end_time-start_time:>10.2f} {w.max_cc_value:>8.3f} {w.cc_shift:>8.3f}")
    print("-" * 50)
    print(f"Total Pyflex Windows: {len(windows)}\n")
    
    return windows
    return windows

def main(is_split_phases=True):
    # Load data
    obs_file = "../example_data/obs.II.ABKT.LHZ.sac"
    syn_file = "../example_data/syn.II.ABKT.LHZ.sac"
    
    obs = obspy.read(obs_file)[0]
    syn = obspy.read(syn_file)[0]
    
    # Filter parameters
    min_period = 40.0
    max_period = 120.0
    
    # Run Fortran window selection
    run_fortran_window_selection(obs_file, syn_file, min_period, max_period, is_split_phases=is_split_phases)
    
    # Filter
    obs.filter("bandpass", freqmin=1.0/max_period, freqmax=1.0/min_period, corners=2, zerophase=True)
    syn.filter("bandpass", freqmin=1.0/max_period, freqmax=1.0/min_period, corners=2, zerophase=True)
    
    # Get Pyflex windows
    pyflex_windows = run_pyflex(obs, syn, min_period, max_period, is_split_phases)
    
    # Get Fortran windows
    fortran_windows = read_fortran_windows("window_results.txt")
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    times = obs.times()
    
    # Plot waveforms
    ax1.plot(times, obs.data, 'k', label='Observed', linewidth=1)
    ax1.plot(times, syn.data, 'r', label='Synthetic', linewidth=1)
    ax1.set_ylabel('Amplitude')
    ax1.set_title(f'Waveforms & Pyflex Windows ({len(pyflex_windows)})')
    ax1.legend(loc='upper right')
    
    ax2.plot(times, obs.data, 'k', label='Observed', linewidth=1)
    ax2.plot(times, syn.data, 'r', label='Synthetic', linewidth=1)
    ax2.set_ylabel('Amplitude')
    ax2.set_xlabel('Time (s)')
    ax2.set_title(f'Waveforms & Fortran Windows ({len(fortran_windows)})')
    
    # Overlay Pyflex windows on ax1
    ylim = ax1.get_ylim()
    for win in pyflex_windows:
        rect = patches.Rectangle((win.left, ylim[0]), win.right - win.left, ylim[1] - ylim[0],
                                 linewidth=1, edgecolor='none', facecolor='green', alpha=0.3)
        ax1.add_patch(rect)
        # Add label
        ax1.text((win.left + win.right)/2, ylim[1], f"P", ha='center', va='bottom', fontsize=8, color='green')

    # Overlay Fortran windows on ax2
    ylim = ax2.get_ylim()
    for win in fortran_windows:
        rect = patches.Rectangle((win['start'], ylim[0]), win['end'] - win['start'], ylim[1] - ylim[0],
                                 linewidth=1, edgecolor='none', facecolor='blue', alpha=0.3)
        ax2.add_patch(rect)
        # Add label
        ax2.text((win['start'] + win['end'])/2, ylim[1], f"F{win['id']}", ha='center', va='bottom', fontsize=8, color='blue')

    plt.tight_layout()
    plt.savefig('comparison_plot.png', dpi=150)
    print("Plot saved to comparison_plot.png")

if __name__ == "__main__":
    main()
