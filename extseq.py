# this is only for .seq files. It basically extracts the sequences from the .seq files and stores them in a list. Then it writes the sequences to multiple files. The number of sequences per file is specified by the max_sequences_per_file variable.
import os

# Set the path to the directory containing the .seq files
directory = 'E:\College\Sem 6\Computational Biophysics\Term Project\Dataset\Database(5sRNA)\Database(5sRNA)\Plastids'

# Initialize a list to store the sequences
sequences = []

# Loop through all files in the directory
for filename in os.listdir(directory):
    # Check if the file is a .seq file
    if filename.endswith('.seq'):
        # Open the file and read the contents
        with open(os.path.join(directory, filename), 'r') as f:
            seq = f.readlines()
        # Extract the sequence from the contents
        name = seq[1].strip()
        sequence = seq[2].strip()
        sequence = sequence[:len(sequence)-1]
        # Add the sequence to the list
        sequences.append(sequence)

# Write the sequences to multiple files
max_sequences_per_file = 100
num_files = (len(sequences) + max_sequences_per_file - 1) // max_sequences_per_file

for file_num in range(num_files):
    # Determine the range of sequences to write to this file
    start_index = file_num * max_sequences_per_file
    end_index = min((file_num + 1) * max_sequences_per_file, len(sequences))
    sequences_to_write = sequences[start_index:end_index]

    # Make directory to store output files
    output_dir = 'seqoutput'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Write the sequences to a new file
    filename = os.path.join('seqoutput', f'all_sequences_{file_num+1}.txt')
    with open(filename, 'w') as f:
        # Write the number of sequences
        f.write(f'{len(sequences_to_write)}\n')
        # Write each sequence to the file
        for i, sequence in enumerate(sequences_to_write):
            f.write(f'{sequence}\n')