#!/usr/bin/python

# FES PLOTTING SCRIPT - INDIVIDUAL COLORBARS WITH PROPER SCALING
# Each system gets its own colorbar with 0-10 kcal/mol range (default)
# Uses inferno colormap for better visualization

import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse

# Constants
KB = 3.2976268E-24  # cal/K
AN = 6.02214179E23  # Avogadro's number

def read_data(infilename):
    """Reads a data file and returns v1 and v2 arrays."""
    v1 = []
    v2 = []
    try:
        with open(infilename, 'r') as ifile:
            for line in ifile:
                if line.strip() == "" or line.startswith(("#", "@", "&")):
                    continue
                newline = line.strip().split()
                if len(newline) >= 2:
                    try:
                        v1.append(float(newline[0]))
                        v2.append(float(newline[1]))
                    except ValueError:
                        continue
    except FileNotFoundError:
        print(f"Error: File not found: {infilename}")
        return None, None
    
    if not v1:
        print(f"Error: No valid data read from {infilename}")
        return None, None
    
    v1_array = np.array(v1)
    v2_array = np.array(v2)
    print(f"  Read {len(v1)} points from {infilename}")
    print(f"  X range: [{v1_array.min():.3f}, {v1_array.max():.3f}]")
    print(f"  Y range: [{v2_array.min():.3f}, {v2_array.max():.3f}]")
    
    return v1_array, v2_array

def calculate_fes(v1_array, v2_array, T, i1, i2, x_min, x_max, y_min, y_max):
    """
    Calculates the Free Energy Surface (DG) for a single system.
    Returns DG array with shape (i1, i2) where DG[ix, iy] corresponds to
    x-coordinate bin ix and y-coordinate bin iy.
    """
    
    # Define the bin edges
    x_bins = np.linspace(x_min, x_max, i1 + 1)
    y_bins = np.linspace(y_min, y_max, i2 + 1)

    # Create the 2D histogram
    H, xedges, yedges = np.histogram2d(v1_array, v2_array, bins=[x_bins, y_bins])
    
    # Count points in range
    points_in_range = np.sum(H)
    print(f"  Points within plot range: {int(points_in_range)} / {len(v1_array)}")

    # Find maximum probability
    Pmax = H.max()

    if Pmax == 0:
        print(f"  Warning: No data found within the specified range!")
        return np.full((i1, i2), 10.0)

    # Calculate Delta G values using barrier cutoff
    # Set a reasonable cutoff for unsampled regions
    barrier_cutoff = 10.0
    
    with np.errstate(divide='ignore'):
        # Where H > 0, calculate -RT ln(P/Pmax)
        # Where H == 0, set to barrier_cutoff
        DG = np.where(H > 0, 
                      -0.001 * AN * KB * T * (np.log(H) - np.log(Pmax)),
                      barrier_cutoff)
    
    # Mask out barrier regions for statistics
    valid_DG = DG[DG < barrier_cutoff]
    if len(valid_DG) > 0:
        print(f"  FES range: [{valid_DG.min():.3f}, {valid_DG.max():.3f}] kcal/mol")
        print(f"  FES mean: {valid_DG.mean():.3f} kcal/mol")
        print(f"  FES std: {valid_DG.std():.3f} kcal/mol")
    print(f"  Bins with data: {np.sum(H > 0)} / {i1 * i2}")
    
    return DG

def main():
    #### Defining flags and help messages ############
    parser = argparse.ArgumentParser(description="Generate a 2x2 FES plot for four systems.")
    # Four input files
    parser.add_argument("-f1", help="Input data file for System 1 (top-left)", required=True)
    parser.add_argument("-f2", help="Input data file for System 2 (top-right)", required=True)
    parser.add_argument("-f3", help="Input data file for System 3 (bottom-left)", required=True)
    parser.add_argument("-f4", help="Input data file for System 4 (bottom-right)", required=True)
    
    parser.add_argument("-o", help="Output file, should be .png", required=True)
    parser.add_argument("-t", help="Temperature in Kelvin", type=float, required=True)
    parser.add_argument("-bx", help="Resolution along x", type=int, required=True)
    parser.add_argument("-by", help="Resolution along y", type=int, required=True)
    parser.add_argument("-lx", help="Label x-axis", default="PC1")
    parser.add_argument("-ly", help="Label y-axis", default="PC2")
    
    # Optional arguments for fixed limits
    parser.add_argument("-xmin", help="Force minimum x-axis limit", type=float, default=None)
    parser.add_argument("-xmax", help="Force maximum x-axis limit", type=float, default=None)
    parser.add_argument("-ymin", help="Force minimum y-axis limit", type=float, default=None)
    parser.add_argument("-ymax", help="Force maximum y-axis limit", type=float, default=None)
    parser.add_argument("-vmin", help="Minimum value for color scale (kcal/mol)", type=float, default=None)
    parser.add_argument("-vmax", help="Maximum value for color scale (kcal/mol)", type=float, default=None)
    parser.add_argument("--autoscale", help="Auto-scale color range to data", action='store_true')
    parser.add_argument("--shared-colorbar", help="Use one shared colorbar instead of individual", action='store_true')
    
    args = parser.parse_args()

    ##### Data Loading ##########
    infiles = [args.f1, args.f2, args.f3, args.f4]
    sub_titles = ["C0", "C1", "C2", "C3"]
    
    all_data = []
    v1_global, v2_global = [], []

    print("\n" + "="*60)
    print("READING DATA FILES")
    print("="*60)
    
    for infile in infiles:
        print(f"\nReading {infile}...")
        v1, v2 = read_data(infile)
        if v1 is None:
            sys.exit(f"Failed to read data from {infile}. Exiting.")
        
        all_data.append((v1, v2))
        v1_global.extend(v1)
        v2_global.extend(v2)

    # Determine plot boundaries
    v1_global = np.array(v1_global)
    v2_global = np.array(v2_global)
    
    plot_xmin = args.xmin if args.xmin is not None else v1_global.min()
    plot_xmax = args.xmax if args.xmax is not None else v1_global.max()
    plot_ymin = args.ymin if args.ymin is not None else -7.0
    plot_ymax = args.ymax if args.ymax is not None else v2_global.max()

    print("\n" + "="*60)
    print("GLOBAL PLOT SETTINGS")
    print("="*60)
    print(f"{args.lx} range: [{plot_xmin:.3f}, {plot_xmax:.3f}]")
    print(f"{args.ly} range: [{plot_ymin:.3f}, {plot_ymax:.3f}]")
    print(f"Resolution: {args.bx} x {args.by} bins")
    print(f"Temperature: {args.t} K")
    
    # Define the extent for all plots [left, right, bottom, top]
    plot_extent = [plot_xmin, plot_xmax, plot_ymin, plot_ymax]

    ##### Calculate all FES first ##########
    print("\n" + "="*60)
    print("CALCULATING FREE ENERGY SURFACES")
    print("="*60)
    
    all_DG = []
    for i in range(4):
        title = sub_titles[i]
        v1_data, v2_data = all_data[i]
        
        print(f"\nProcessing {title} ({infiles[i]})...")
        
        # Calculate FES using the pre-read data and global plot limits
        DG = calculate_fes(v1_data, v2_data, args.t, args.bx, args.by, 
                           plot_xmin, plot_xmax, plot_ymin, plot_ymax)
        all_DG.append(DG)
    
    ##### Determine color scale #####
    # Collect all valid (non-barrier) DG values
    all_valid_DG = []
    barrier_cutoff = 10.0
    for DG in all_DG:
        valid = DG[DG < barrier_cutoff]
        if len(valid) > 0:
            all_valid_DG.extend(valid)
    
    if len(all_valid_DG) == 0:
        print("\nError: No valid FES data calculated!")
        sys.exit(1)
    
    all_valid_DG = np.array(all_valid_DG)
    
    # Determine vmin and vmax
    if args.autoscale or (args.vmin is None and args.vmax is None):
        # Auto-scale: use 0 to 95th percentile or max (whichever is reasonable)
        vmin = 0.0
        max_value = all_valid_DG.max()
        p95 = np.percentile(all_valid_DG, 95)
        
        # If the data is very flat (max < 2), use the actual max
        # Otherwise use 95th percentile to avoid extreme outliers
        if max_value < 2.0:
            vmax = max_value
        else:
            vmax = min(p95 * 1.2, max_value)  # 20% above 95th percentile
        
        # Round to nice number
        if vmax < 1:
            vmax = math.ceil(vmax * 10) / 10  # Round to 0.1
        elif vmax < 5:
            vmax = math.ceil(vmax * 2) / 2  # Round to 0.5
        else:
            vmax = math.ceil(vmax)  # Round to 1
            
    else:
        vmin = args.vmin if args.vmin is not None else 0.0
        vmax = args.vmax if args.vmax is not None else 10.0
    
    print("\n" + "="*60)
    print("COLOR SCALE SETTINGS")
    print("="*60)
    print(f"Global FES range: [{all_valid_DG.min():.3f}, {all_valid_DG.max():.3f}] kcal/mol")
    print(f"Color scale: [{vmin:.3f}, {vmax:.3f}] kcal/mol")
    print(f"Colorbar mode: {'Shared' if args.shared_colorbar else 'Individual per system'}")
    if vmax < 2:
        print("WARNING: Very flat energy landscape detected!")
        print("This may indicate:")
        print("  - Insufficient sampling (need more data)")
        print("  - Poor choice of reaction coordinates")
        print("  - Simulation not equilibrated")
    
    ##### Plotting ##########
    print("\n" + "="*60)
    print("GENERATING PLOTS")
    print("="*60)
    
    # Create 2x2 subplot
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 12))
    ax_list = axes.flatten()
    
    z_label = r'$\Delta G$ [kcal/mol]'
    im = None 
    
    outfilename2 = "FES_all_systems.dat"
    
    # Calculate bin centers for writing to file
    bin_width_x = (plot_xmax - plot_xmin) / args.bx
    bin_width_y = (plot_ymax - plot_ymin) / args.by
    
    x_centers = plot_xmin + (np.arange(args.bx) + 0.5) * bin_width_x
    y_centers = plot_ymin + (np.arange(args.by) + 0.5) * bin_width_y

    with open(outfilename2, 'w') as ofile:
        for i in range(4):
            ax = ax_list[i]
            title = sub_titles[i]
            DG = all_DG[i]
            
            if DG is not None:
                # Plot the FES
                # DG has shape (bx, by) where first index is x, second is y
                # imshow expects (ny, nx) with origin='lower', so we transpose
                im = ax.imshow(DG.T, cmap=cm.jet, 
                               extent=plot_extent, 
                               origin='lower', aspect='auto', 
                               vmin=vmin, vmax=vmax)
                
                # Add contour lines if the range is reasonable
                if vmax - vmin > 0.5:
                    X, Y = np.meshgrid(x_centers, y_centers)
                    # Create contour levels
                    n_contours = min(10, int((vmax - vmin) * 5))  # Adaptive number
                    contour_levels = np.linspace(vmin, vmax, n_contours)
                    ax.contour(X, Y, DG.T, levels=contour_levels, 
                              colors='black', alpha=0.3, linewidths=0.5)
                
                # Add INDIVIDUAL colorbar for each subplot (unless shared mode)
                if not args.shared_colorbar:
                    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                    cbar.set_label(z_label, size=10)
                    cbar.ax.tick_params(labelsize=8)
                
                # Write data to the .dat file
                ofile.write(f"# Data for {title} ({infiles[i]})\n")
                ofile.write(f"# X({args.lx})\tY({args.ly})\tDG(kcal/mol)\n")
                for ix in range(args.bx):
                    for iy in range(args.by):
                        ofile.write(f"{x_centers[ix]:.6f}\t{y_centers[iy]:.6f}\t{DG[ix, iy]:.6f}\n")
                    ofile.write("\n")  # Gnuplot blank line
                ofile.write("\n\n")

            # Set titles and labels
            ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
            if i >= 2:  # Bottom row
                ax.set_xlabel(args.lx, fontsize=11)
            if i % 2 == 0:  # Left column
                ax.set_ylabel(args.ly, fontsize=11)
            
            ax.tick_params(labelsize=9)
            ax.grid(True, alpha=0.2, linestyle='--')

    print(f"\nPlot saved to: {args.o}")
    print(f"Data saved to: {outfilename2}")

    # Add a single shared colorbar ONLY if requested
    if args.shared_colorbar:
        fig.subplots_adjust(right=0.85, hspace=0.25, wspace=0.25)
        cbar_ax = fig.add_axes([0.88, 0.15, 0.03, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label(z_label, size=12)
    else:
        # Adjust spacing for individual colorbars
        fig.subplots_adjust(hspace=0.3, wspace=0.4)

    fig.suptitle("Free Energy Surfaces", fontsize=24, y=0.98,fontweight="bold")
    plt.savefig(args.o, dpi=300, bbox_inches='tight')
    plt.savefig("fes2.svg", format="svg", bbox_inches="tight", pad_inches=0.05)
    plt.close()
    
    print("\nDone!")
    print("="*60 + "\n")
    sys.exit(0)

if __name__ == "__main__":
    main()
