import os, sys
import pandas as pd



def parse_file(file, sample_name, D):
    with open(file, "r") as f:
        data = f.readlines()
        count = 0
        for line in data:
            L = line.strip().split(" ")
            if count == 0:
                D[sample_name] = L
                percent_align = float(L[1])/float(L[3])
                D[sample_name].append(percent_align)
                
            elif count == 1:
                D[sample_name].extend(L)
            elif count == 2:
                D[sample_name].extend(L)
            else:
                break
            count = count + 1
    return D

def main():
    data_dir = str(sys.argv[1])
    print(data_dir)
    result_dict = {}
    for dir in os.listdir(data_dir):
        if dir.endswith(".stat"):
            sample_name = dir.split(".stat")[0]
            print(sample_name)
            stat_dir = data_dir + dir
            print(stat_dir)
            for file in os.listdir(stat_dir):
                if file.endswith(".cnt"):
                    print(file)
                    cnt_file = stat_dir + "/"+ file
                    parse_file(cnt_file, sample_name, result_dict)
    
    print(result_dict)
    #convert to dataframe and save
    
    df = pd.DataFrame.from_dict(result_dict, orient='index')
    print(df)
    # add colnames
    df.columns = ["Unalignable Reads", "Alignable Reads", "Filtered Reads (too many alignments)", "Total Reads", "Percent Aligned", "Uniquely Mapped Reads", "Multi-Mapped Reads", "Uncertain Reads", "Total Alignments", "Read type"]
    df.to_csv(data_dir + "summarized_rsem_align.csv", sep=',', header=True, index=True)
    
                
    
if __name__ == "__main__":
    main()