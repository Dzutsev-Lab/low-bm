import argparse
import pandas as pd
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--host-names", dest="host_names", required=True)
    ap.add_argument("--sep", default="\t")
    ap.add_argument("--id-col", default=None,
                    help="Name of samples ID column in input sequence table. If omitted, will use first column")
    ap.add_argument("--method", required=True)
    ap.add_argument("--offset", default=0.0, type=float)
    args = ap.parse_args()

    # Import sequence table
    seq_table_df = pd.read_csv(args.inp, sep=args.sep, dtype={0: str})
    print('Original Sequence Table')
    print(seq_table_df.info())
    print(seq_table_df.head())
    
    # Identify SampleID column
    if args.id_col is None:
        sample_id_col = seq_table_df.columns[0]
    else:
        sample_id_col = args.id_col
    
    sample_ids = seq_table_df[sample_id_col]
    print("\nSample IDs")
    print(sample_ids)



    # Create count-only data frame for easy calculation
    count_df = seq_table_df.drop(columns=[sample_id_col]).apply(pd.to_numeric)
    print('\nCounts Sequence Table')
    print(count_df.info())
    print(count_df.head())



    #---------------------------
    # Normalization Function
    #---------------------------
    def noNorm(count_df):
        return(count_df)
    
    def rawTSS(count_df):
        # host count df included as parameter to make inclusion in function dict easier
        return count_df.div(count_df.sum(axis='columns').replace(0, np.nan),
                            axis='index').fillna(0.0) #keep all as numerial now that zero divsor averted
    
    def hostTSS(count_df):
        # Retrieve Unmapped ASV names
        with open(args.host_names, 'r') as name_file:
            unmapped_names = [line.strip() for line in name_file]
        print('\nUnmapped ASV IDs')
        print(unmapped_names)

        # Create data frame with only host-mapped ASV columns
        host_count_df = count_df.drop(columns= unmapped_names)
        print('\nHost Mapped Sequence Count Table')
        print(host_count_df.info())
        print(host_count_df.head())

        return count_df.div(host_count_df.sum(axis='columns').replace(0, 1), # replacing with 1 to reward samples of no host alignment
                            axis='index')
    
    def rawTSSpM(count_df):
        return (rawTSS(count_df) * 1000000)
    
    def hostTSSpM(count_df):
        return (hostTSS(count_df) * 1000000)
    
    def log2(count_df):
        return np.log2(count_df + args.offset)

    def rawTSSlog2(count_df):
        return log2(rawTSS(count_df))

    def hostTSSlog2(count_df):
        return log2(hostTSS(count_df))
    
    def rawTSSpMlog2(count_df):
        return log2(rawTSS(count_df) * 1000000)
    
    def hostTSSpMlog2(count_df):
        return log2(hostTSS(count_df) * 1000000)
    
    norm_methods_dict = {
        "noNorm" : noNorm,
        "rawTSS" : rawTSS,
        "hostTSS" : hostTSS,
        "rawTSSpM" : rawTSSpM,
        "hostTSSpM" : hostTSSpM,        
        "log2" : log2,
        "rawTSSlog2" : rawTSSlog2,
        "hostTSSlog2" : hostTSSlog2,
        "rawTSSpMlog2" : rawTSSpMlog2,
        "hostTSSpMlog2" : hostTSSpMlog2,        
    }


    #---------------------------
    # Normalize Sequence Table
    #---------------------------
    normalized_df = norm_methods_dict[args.method](count_df)

    
    # Output Normalized data frame with sample IDs attached
    normalized_df = pd.concat([sample_ids, normalized_df], axis='columns')
    normalized_df.to_csv(args.out, sep=args.sep, index=False)

    print('\nNormalized Sequence Table')
    print(normalized_df.info())
    print(normalized_df.head())

if __name__ == "__main__":
    main()

    
