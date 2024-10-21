import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_batch(batch_data):
    # Process a batch of data
    sequences = []
    all_indices = []
    
    for record in batch_data:
        sequence = record[0]
        indices = np.array([int(x) for x in record[2].split(',')]) - 1
        count = len(indices)  # Use the length of indices as the count
        
        # Make a list of sequence that is as long as the indices
        sequences.extend([sequence] * count)
        # Make a list of indices
        all_indices.append(indices)
    
    #Put it in a dataframe
    return pd.DataFrame({
        'sequence': sequences,
        'index': np.concatenate(all_indices)
    })

def Barcode_scanner(Combined_Barcode, sampleName, length, distance, ind):
    batch_size = 10000  # Adjust this based on your system's memory
    
    # Step 1: Get the columns from Combined_Barcode
    columns = Combined_Barcode.columns.tolist()
    
    # Step 2: Validate ind
    if ind < 0 or ind >= len(columns):
        raise ValueError(f"Invalid ind value: {ind}. Must be between 0 and {len(columns) - 1}")
    
    # Step 3: Select the column to update
    column = columns[ind]

    # If its not processing 
    if length != "sequence": 
        # Step 4: Read the entire file into memory
        with open(f'{sampleName}_Barcode{length}_d{distance}.txt', 'r') as f:
            data = [line.strip().split('\t') for line in f]
        
        # Step 5: Process data in parallel
        with ThreadPoolExecutor() as executor:
            futures = []
            for i in range(0, len(data), batch_size):
                batch = data[i:i+batch_size]
                futures.append(executor.submit(process_batch, batch))
            
            # Collect results
            results = []
            for future in as_completed(futures):
                results.append(future.result())
    
        # Step 6: Combine all results 
        all_data = pd.concat(results, ignore_index=True)
    
        # Update the Combined_Barcode DataFrame
        # loc = pandas' label-based indexer for selecting data
        # selects the 'index' column from our processed data that corresponds to Combined_Barcode 
        # selects the 'sequence' column from our processed data and converts it to a numpy array with .values and inserting them to Combined_Barcode
        Combined_Barcode.loc[all_data['index'], column] = all_data['sequence'].values
    
    return Combined_Barcode


def process_barcodes(df, ind, columns, length, dist):
    # Ensure ind is a list
    if not isinstance(ind, list):
        ind = [ind]
    
    #Take all columns in the dataset 
    columns = df.columns.tolist()
    
    # Select only the specified columns
    selected_columns = [columns[i] for i in ind]
    
    # Step 1: Separate the data by Sample Name
    grouped = df.groupby('Sample')
    
    for sample, group in grouped:
        # Step 2: Process each specified barcode column
        for column in selected_columns:
            
            # Step 3: Sum counts for same barcode and sort
            processed = group[['Counts', column]].copy()
            # reset_index() turns this Series back into a DataFrame with two columns: 'Barcode' and 'Counts'.
            result = processed.groupby(column)['Counts'].sum().reset_index()
            result = result.sort_values('Counts', ascending=False)
            
            # Step 4: Output to txt file
            output_filename = f"{sample}_final{length}_d{dist}.txt"
            result.to_csv("{}_final{}_{}.txt".format(sample,str(length),str(dist)), sep='\t', index=False, header=False, float_format='%g')
            print(f"Processed data for {sample}, {column} saved to {output_filename}")
