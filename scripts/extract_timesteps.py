#!/usr/bin/env python3
import argparse
import os
import sys
import xml.etree.ElementTree as ET

try:
    import h5py
except ImportError:
    print("Error: 'h5py' is not installed. Run 'pip install h5py' first.", file=sys.stderr)
    sys.exit(1)

def slice_fenics_xdmf(input_file, output_file, stride, start, end):
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing XML from {input_file}...")
    tree = ET.parse(input_file)
    root = tree.getroot()

    # Find the temporal collection grid containing all time steps
    temporal_grid = root.find(".//Grid[@CollectionType='Temporal']")
    if temporal_grid is None:
        temporal_grid = root.find(".//Grid[@GridType='Collection']")
        
    if temporal_grid is None:
        print("Error: Could not find a temporal Grid collection in the XDMF file.", file=sys.stderr)
        sys.exit(1)

    # Extract all individual uniform time steps
    uniform_grids = temporal_grid.findall("./Grid[@GridType='Uniform']")
    total_steps = len(uniform_grids)
    print(f"Found {total_steps} total time steps in the original file.")

    # Determine boundaries and slice indices
    stop_index = total_steps if end is None else min(end, total_steps)
    target_indices = list(range(start, stop_index, stride))

    if not target_indices:
        print("Error: No time steps match your criteria.", file=sys.stderr)
        sys.exit(1)

    print(f"Extracting {len(target_indices)} steps...")

    # Clear out all uniform grids from the temporal collection in the XML
    for grid in uniform_grids:
        temporal_grid.remove(grid)

    # Setup destination paths
    output_dir = os.path.dirname(output_file) or '.'
    output_basename = os.path.splitext(os.path.basename(output_file))[0]
    output_h5_name = f"{output_basename}.h5"
    output_h5_path = os.path.join(output_dir, output_h5_name)
    
    input_dir = os.path.dirname(input_file) or '.'
    src_h5_files = {}

    # Open the new destination HDF5 file
    with h5py.File(output_h5_path, 'w') as dst_h5:
        for idx, k in enumerate(target_indices, start=1):
            grid = uniform_grids[k]
            # Put the selected time step back into our lightened XML structure
            temporal_grid.append(grid)
            
            # Find all HDF5 DataItems for this step and copy their raw datasets
            for data_item in grid.findall(".//DataItem[@Format='HDF']"):
                text = data_item.text.strip()
                src_h5_filename, h5_path = text.split(':')
                
                # Open source HDF5 on demand
                if src_h5_filename not in src_h5_files:
                    full_src_h5_path = os.path.join(input_dir, src_h5_filename)
                    if not os.path.exists(full_src_h5_path):
                        print(f"\nError: Source HDF5 file '{full_src_h5_path}' not found.", file=sys.stderr)
                        sys.exit(1)
                    src_h5_files[src_h5_filename] = h5py.File(full_src_h5_path, 'r')
                
                src_h5 = src_h5_files[src_h5_filename]
                
                if h5_path in src_h5:
                    # Recreate matching nested group structures in the target file
                    parent_path = os.path.dirname(h5_path)
                    basename = os.path.basename(h5_path)
                    
                    parent_group = dst_h5.require_group(parent_path)
                    
                    # Copy dataset if it hasn't been copied already
                    if basename not in parent_group:
                        src_h5.copy(h5_path, parent_group, name=basename)
                
                # Point the XML element to the new .h5 file name
                data_item.text = f"{output_h5_name}:{h5_path}"
            
            # Progress line
            sys.stdout.write(f"\rProgress: [{idx}/{len(target_indices)}] Copied step index {k}")
            sys.stdout.flush()

    # Clean up open source files
    for f in src_h5_files.values():
        f.close()

    # Write out the brand new lightened XDMF metadata file
    tree.write(output_file, encoding='utf-8', xml_declaration=True)
    print(f"\n\nSuccess! Lightened FEniCS dataset saved to:\n  -> {output_file}\n  -> {output_h5_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Directly slice FEniCS/DOLFIN temporal XDMF/HDF5 datasets without corrupting finite element layouts."
    )
    parser.add_argument("input_file", type=str, help="Path to original .xdmf")
    parser.add_argument("output_file", type=str, help="Path to new .xdmf")
    parser.add_argument("-s", "--stride", type=int, default=5, help="Keep every Nth time step")
    parser.add_argument("--start", type=int, default=0, help="Start step index")
    parser.add_argument("--end", type=int, default=None, help="End step index")
    
    args = parser.parse_args()
    slice_fenics_xdmf(args.input_file, args.output_file, args.stride, args.start, args.end)

if __name__ == "__main__":
    main()