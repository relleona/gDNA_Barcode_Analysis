import numpy as np
import csv
import pandas as pd
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

def update_column(df, indices, sequences, ind):
    df.iloc[indices, ind + 1] = sequences

# This function updates the barcodes in the dataframe based on merging by starcode
def Barcode_scanner(Combined_Barcode, sampleName, length, distance, ind):
    batch_size = 1000  # Adjust this based on your system's memory
    with open(f'{sampleName}_Barcode{length}_d{distance}.txt') as B:
        Barcode = csv.reader(B, delimiter='\t')
        data = list(Barcode)
        num_records = len(data)

        with ThreadPoolExecutor() as executor:
            for i in range(0, num_records, batch_size):
                batch_data = data[i:i + batch_size]
                indices = [np.array([int(x) for x in record[2].split(',')], dtype=object) - 1 for record in batch_data]
                sequences = [record[0] for record in batch_data]

                # Update DataFrame in parallel
                executor.submit(update_column, Combined_Barcode, indices, sequences, ind)

    return Combined_Barcode

# This function creates separate csv files for barcodes of different samples
def Separate_samples(Combined_reads, combined_flag, sampleName, length, dist):
    if(combined_flag):
        sample_list = np.unique(np.array(Combined_reads['Sample']))
        for i in sample_list:
            df = Combined_reads[Combined_reads['Sample'] == i]
            Sequence = np.repeat(df['Sequence'], df['Counts'])
            B30 = np.repeat(df['Barcode_30'], df['Counts'])
            B40 = np.repeat(df['Barcode_40'], df['Counts'])
            B50 = np.repeat(df['Barcode_50'], df['Counts'])
            df = pd.DataFrame({'Sequence': Sequence, 'Barcode_30': B30, 'Barcode_40': B40, 'Barcode_50': B50})
            del B30, B40, B50, Sequence
            final = 'Barcode_' + str(length)
            fin = np.array(count_Barcode(np.array(df[final])), dtype = object)
            df = pd.DataFrame(fin, columns=['Sequence', 'Counts'])
            df['Counts'] = df['Counts'].astype(int)
            df = df.sort_values(by= "Counts", ascending=False)
            df.to_csv("{}_final{}_{}.txt".format(i,str(length),str(dist)),sep='\t', index=False, header=False, float_format='%g')
    else:
        df = Combined_reads
        final = 'Barcode_' + str(length)   
        Sequence = np.repeat(df['Sequencie'], df['Counts'])
        B30 = np.repeat(df['Barcode_30'], df['Counts'])
        B40 = np.repeat(df['Barcode_40'], df['Counts'])
        B50 = np.repeat(df['Barcode_50'], df['Counts'])
        df = pd.DataFrame({'Sequence': Sequence, 'Barcode_30': B30, 'Barcode_40': B40, 'Barcode_50': B50})
        del B30, B40, B50, Sequence
        fin = np.array(count_Barcode(np.array(df[final])), dtype = object)
        df = pd.DataFrame(fin, columns=['Sequence', 'Counts'])
        df['Counts'] = df['Counts'].astype(int)
        df = df.sort_values(by= "Counts", ascending=False)
        df.to_csv("{}_final{}_{}.txt".format(sampleName,str(length),str(dist)), sep='\t', index=False, header=False, float_format='%g')

# This function counts the occurences of a particular barcode to give read counts
def count_Barcode(Barcode_array):
    counter = Counter(Barcode_array)
    result_array = np.array([(element, count) for element, count in counter.items()], dtype=object)
    return result_array


