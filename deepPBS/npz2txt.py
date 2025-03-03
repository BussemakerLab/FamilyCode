import numpy as np
import pandas as pd
import os

def npz_to_4xN_tsv(input_dir):
    """
    Convert all .npz files in a directory into 4xN matrices and save them as .tsv files.

    Parameters:
        input_dir (str): Path to the directory containing .npz files.
    """
    try:
        # List all .npz files in the directory
        for file_name in os.listdir(input_dir):
            if file_name.endswith('.npz'):
                input_file = os.path.join(input_dir, file_name)

                # Load the .npz file
                data = np.load(input_file)

                # Extract the first key's data (assuming one column is required)
                key = list(data.keys())[0]
                array = data[key]

                # Flatten the array and reshape it into a 4xN matrix
                flat_array = array.flatten()
                reshaped_array = flat_array.reshape(-1, 4)

                # Convert to a DataFrame
                df = pd.DataFrame(reshaped_array)

                # Generate the output file name
                output_file = os.path.splitext(input_file)[0] + ".tsv"

                # Save the DataFrame to a TSV file
                df.to_csv(output_file, sep='\t', index=False, header=False)
                print(f"Converted {input_file} to {output_file}.")

    except Exception as e:
        print(f"Error processing files in {input_dir}: {e}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_directory>")
    else:
        input_dir = sys.argv[1]
        npz_to_4xN_tsv(input_dir)