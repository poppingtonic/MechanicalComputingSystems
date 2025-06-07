import sys
import os

def format_pdb_line(atom_serial, atom_name, res_name, chain_id, res_seq, x, y, z, occupancy, temp_factor, element_sym):
    """
    Formats an ATOM record line for a PDB file according to specification.
    Ensures proper spacing and justification.
    """
    # PDB format reference: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    # Record Type (1-6) "ATOM  "
    # Atom Serial (7-11) Integer, right-justified
    # Atom Name (13-16) String, left-justified (starting at col 14 for 1-char element)
    # Alt Loc (17) Char
    # Residue Name (18-20) String, right-justified? No, usually left/center.
    # Chain ID (22) Char
    # Residue Sequence (23-26) Integer, right-justified
    # Insertion Code (27) Char
    # X coord (31-38) Float (8.3), right-justified
    # Y coord (39-46) Float (8.3), right-justified
    # Z coord (47-54) Float (8.3), right-justified
    # Occupancy (55-60) Float (6.2), right-justified
    # Temp Factor (61-66) Float (6.2), right-justified
    # Element Symbol (77-78) String, right-justified
    # Charge (79-80) String

    # Adjust atom name based on element length for better alignment sometimes
    # Simple approach: Left justify element in atom name field
    atom_name_pdb = f" {atom_name:<3s}" # e.g., " C  ", " H  ", " FE " - places element starting col 14

    line = (
        f"ATOM  {atom_serial:>5d}"      # Atom serial number (cols 7-11)
        f" {atom_name_pdb:<4s}"         # Atom name (cols 13-16) - Use formatted name
        f"{res_name:<3s} "              # Residue name (cols 18-20) + space for Alt Loc (17)
        f"{chain_id}"                   # Chain identifier (col 22)
        f"{res_seq:>4d}    "            # Residue sequence number (cols 23-26) + space for iCode (27) + 3 spaces
        f"{x:>8.3f}"                    # X coordinate (cols 31-38)
        f"{y:>8.3f}"                    # Y coordinate (cols 39-46)
        f"{z:>8.3f}"                    # Z coordinate (cols 47-54)
        f"{occupancy:>6.2f}"            # Occupancy (cols 55-60)
        f"{temp_factor:>6.2f}      "    # Temperature factor (cols 61-66) + spaces
        f"    {element_sym:>2s}  "      # Element symbol (cols 77-78) + Charge (cols 79-80) blank
    )
    # Ensure line does not exceed 80 characters if needed, though this format fits.
    return line + '\n'


def xyz_to_pdb(input_xyz_file, output_pdb_file):
    """
    Converts a multi-frame XYZ file to a multi-model PDB file.

    Args:
        input_xyz_file (str): Path to the input XYZ file.
        output_pdb_file (str): Path to the output PDB file.
    """
    print(f"Starting conversion: {input_xyz_file} -> {output_pdb_file}")
    model_count = 0
    line_num = 0
    processed_atoms_total = 0

    try:
        with open(input_xyz_file, 'r') as infile, open(output_pdb_file, 'w') as outfile:
            while True:
                # --- Read Header for a Frame ---
                # 1. Number of atoms line
                num_atoms_line = infile.readline()
                line_num += 1
                if not num_atoms_line:
                    break # End of file

                try:
                    num_atoms = int(num_atoms_line.strip())
                except ValueError:
                    print(f"Warning: Skipping invalid number of atoms at line {line_num}: '{num_atoms_line.strip()}'", file=sys.stderr)
                    # Try to find the next potential number of atoms
                    while True:
                        next_line = infile.readline()
                        line_num += 1
                        if not next_line:
                            break # End of file
                        try:
                           num_atoms = int(next_line.strip())
                           line_num += 1 # Count the comment line we are about to read
                           infile.readline() # Consume the assumed comment line
                           break # Found a potential start, continue processing
                        except ValueError:
                           continue # Keep searching
                    if not next_line: # Reached EOF while searching
                         break
                    # If we broke from inner loop, we found a potential num_atoms

                if num_atoms <= 0:
                     print(f"Warning: Skipping frame with non-positive number of atoms ({num_atoms}) at line {line_num}.", file=sys.stderr)
                     # Skip the comment line and atoms for this empty frame
                     infile.readline() # Consume comment line
                     line_num += 1
                     for _ in range(num_atoms): # Should not loop if num_atoms <= 0
                         infile.readline()
                         line_num += 1
                     continue

                # 2. Comment line (often ignored, but read it)
                comment_line = infile.readline() # Can use this for TITLE if needed
                line_num += 1
                if not comment_line:
                    print(f"Warning: Reached end of file unexpectedly after reading atom count at line {line_num-1}.", file=sys.stderr)
                    break

                # --- Write PDB Header for this Model ---
                model_count += 1
                outfile.write(f"MODEL     {model_count:>4d}\n") # Use MODEL record for multi-frame XYZ
                print(f"Processing frame {model_count} with {num_atoms} atoms...")

                # --- Process Atom Lines ---
                atom_serial_pdb = 0 # Reset atom serial for each model in PDB
                for i in range(num_atoms):
                    atom_line = infile.readline()
                    line_num += 1
                    if not atom_line:
                        print(f"Error: Unexpected end of file while reading atoms for frame {model_count}. Expected {num_atoms} atoms, read {i}.", file=sys.stderr)
                        break # Exit the inner loop

                    parts = atom_line.split()
                    if len(parts) < 4:
                        print(f"Warning: Skipping malformed atom line {line_num} in frame {model_count}: '{atom_line.strip()}'", file=sys.stderr)
                        continue # Skip this line and try the next

                    try:
                        element = parts[0]
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Skipping malformed atom line {line_num} in frame {model_count} due to parsing error ({e}): '{atom_line.strip()}'", file=sys.stderr)
                        continue # Skip this line

                    # --- Format and Write PDB ATOM Line ---
                    atom_serial_pdb += 1
                    processed_atoms_total += 1

                    # Basic PDB field defaults - customize if needed
                    atom_name = element # Use element symbol as base for atom name
                    res_name = "MOL"    # Default residue name
                    chain_id = "A"      # Default chain ID
                    res_seq = 1         # Default residue sequence number
                    occupancy = 1.00    # Default occupancy
                    temp_factor = 0.00  # Default temperature factor (B-factor)
                    element_sym = element.upper() # Element symbol for columns 77-78

                    # Handle atom serial number > 99999 (PDB limit for standard field)
                    if atom_serial_pdb > 99999:
                       print(f"Warning: Atom serial number {atom_serial_pdb} exceeds 99999 (PDB limit) in frame {model_count}. Resetting or truncation might occur in some viewers.", file=sys.stderr)
                       # PDB format technically wraps around or uses hex, but simple truncation/modulo is common
                       # Let's just let it write for now, viewer might handle it.
                       # Or cap it: atom_serial_display = 99999
                       atom_serial_display = atom_serial_pdb # Let it write > 5 digits if needed, though format is fixed width
                    else:
                       atom_serial_display = atom_serial_pdb

                    # Format and write the line
                    pdb_line = format_pdb_line(
                        atom_serial=atom_serial_display,
                        atom_name=atom_name,
                        res_name=res_name,
                        chain_id=chain_id,
                        res_seq=res_seq,
                        x=x, y=y, z=z,
                        occupancy=occupancy,
                        temp_factor=temp_factor,
                        element_sym=element_sym[:2] # Ensure element symbol is max 2 chars
                    )
                    outfile.write(pdb_line)

                # --- Write PDB Footer for this Model ---
                outfile.write(f"ENDMDL\n")

            # --- Write Final PDB Record ---
            if model_count > 0:
                outfile.write("END\n") # Optional but good practice
            else:
                 print("Warning: No valid frames found in the input XYZ file.", file=sys.stderr)


    except FileNotFoundError:
        print(f"Error: Input file not found: {input_xyz_file}", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error reading or writing file: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred at line ~{line_num}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"\nConversion complete.")
    print(f"Processed {model_count} frames and {processed_atoms_total} atoms.")
    print(f"Output saved to: {output_pdb_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python xyz2pdb.py <input_xyz_filename> <output_pdb_filename>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Basic check if input file exists
    if not os.path.exists(input_file):
         print(f"Error: Input file '{input_file}' does not exist.", file=sys.stderr)
         sys.exit(1)

    xyz_to_pdb(input_file, output_file)