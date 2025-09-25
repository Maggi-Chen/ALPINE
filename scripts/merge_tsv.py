import os
import ast
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use('ggplot')


def merge_tsv(list_files, prefix):
    master_list = []
    for n, file in enumerate(list_files):
        with open(file, 'rt') as infile:
            for m, line in enumerate(infile):
                content = line.strip('\n').split('\t')
                col_name, element = content[0].replace(
                    ' ', '_').rstrip(":"), content[1]
                element = element.rstrip(prefix) if m == 0 else element
                if n == 0:
                    master_list.append(col_name + '\t' + element)
                else:
                    master_list[m] += ('\t' + element)

    return master_list


def save_to_tsv(output, list_files, prefix):
    master_list = merge_tsv(list_files, prefix)
    with open(output, 'wt', encoding='UTF-8') as of:
        for row in master_list:
            of.write(row + '\n')

def extract(row, pos):
    return [ast.literal_eval(item)[pos] for item in row]

def divid_two_lists(list1, list2):
    return [x/y for x, y in zip(list1, list2)]

input_dir = 'extracted_samstats'
list_files = [input_dir + "/" +
              file for file in os.listdir(input_dir) if 'txt' in file]

merge_tsv(list_files, '.samstats.txt')
output = input_dir + "/temp_summary_samstats.tsv"
save_to_tsv(output, list_files, ".samstats.txt")

df = pd.read_csv(output, sep='\t')
df.T.to_csv(input_dir + "/summary_stat.tsv",
            sep='\t', header=False, index=True)

input_file = 'multiqc_samtools_idxstats.txt'
df = pd.read_csv(input_file, sep='\t', header=0)
df_new = df[~ df['Sample'].str.contains('filtered_')].copy()
df_new['Sample'] = df_new['Sample'].apply(lambda x: x.rstrip(".idxstat"))

df_new['total_mapped'] = df_new.loc[:, df_new.columns !=
                                    'Sample'].apply(lambda row: sum(extract(row, 0)), axis=1)

cols = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
        'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
        'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
        'chr22', 'chrX', 'chrY', 'JDNA4165_TRAC_KI', 'JDNA4587_B2M_KI']
normalized_size = df_new.loc[:, cols].apply(
                                lambda row: extract(row, 0), axis=1)

df2 = pd.DataFrame(normalized_size.to_list(), columns=cols, index=df_new.index)
df2[['total_mapped']] = df_new[['total_mapped']]
df2.to_csv('total_mapped_reads.csv',
             sep=',',
             header=True,
             index=True)

df2[cols] = df2.loc[:, cols].div(df2['total_mapped'], axis=0)
df2.to_csv('total_mapped_reads2.csv',
             sep=',',
             header=True,
             index=True)

df2_cp = df2[[ 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
              'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
              'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21',
              'chr22', 'chrX', 'chrY', 'JDNA4165_TRAC_KI', 'JDNA4587_B2M_KI']]



samples = df2_cp.index
# Create a single figure with 12 subplots arranged in a 3x4 grid
fig, axes = plt.subplots(3, 4, figsize=(12, 10))
fig.subplots_adjust(wspace=0.4, hspace=0.5)

# Iterate through each sample and create a bar plot in the corresponding subplot

for i, ax in enumerate(axes.flatten(), 0):
    if i < len(samples):
        sample_name = samples[i]
        sample_data = df2_cp.loc[sample_name]
        ax.barh(sample_data.index, sample_data.values, color='salmon')
        ax.set_title(sample_name)

plt.suptitle('Distribution of mapped reads for ALLO samples')
plt.tight_layout()
plt.show()

df_new.loc[:, cols].apply(lambda row: divid_two_lists(
    extract(row, 0), extract(row, 1)), axis=1)

df_new['total_mapped'] = df_new.loc[:, df_new.columns !=
                                    'Sample'].apply(lambda row: sum(extract(row, 0)), axis=1)
df_new.drop(columns='total_mapped', inplace=True)
df_new.loc[:, cols].apply(lambda row: divid_two_lists(
    extract(row, 0), extract(row, 1))/df_new['total_mapped'], axis=1)
