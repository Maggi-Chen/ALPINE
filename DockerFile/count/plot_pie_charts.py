#!/usr/bin/env python3
# coding: utf-8

# import modules
import sys
import math

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.style.use('ggplot')


def plot_varint_count_analysis(tsv_file) -> None:
    """ Plot pie charts to show classification of reads 
    Args:
        tsv_file (str): a tsv output file name
    Return:
        output a pdf figure
    """
    print('Creating a pie chart of variant count result. Reading output table.')
    df = pd.read_csv(tsv_file, sep='\t')
    num_rows, num_cols = 1,1
    
    # Create a figure with subplots and a larger figsize
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(11, 6))
    
    # Define a Seaborn color palette
    colors = sns.color_palette("husl", df.shape[1]-1)
    
    # Initialize empty lists to store legend handles and labels
    legend_handles = []
    legend_labels = []

    print(f'Initializing pie chart of variant count results including {num_rows} rows and {num_cols} columns.')
    
    # Iterate through the rows and create pie charts for each sample
    for index, (_, row) in enumerate(df.iterrows()):
        sample_name = row['#Sample']  # Get the sample name from the 'Sample' column
        values = row.values[1:]  # Exclude the 'Sample' column
        labels = df.columns[1:]  # Use column names as labels
        ax = axes
        
        # Increase the size of the pie chart by adjusting the radius
        pie = ax.pie(values, 
                     labels=['']*len(labels), 
                     autopct='', 
                     startangle=140, 
                     radius=0.9, 
                     colors=colors)
        if index == df.shape[0] - 1:
            legend_handles.extend(pie[0])
            legend_labels.extend(labels)
        
        # Add a title to each subplot
        ax.set_title(f'{sample_name}')
        ax.axis('equal')

    # Remove any unused subplots
    for i in range(num_rows * num_cols, 1):
        fig.delaxes(axes[i])

    plt.legend(legend_handles, legend_labels, loc="best")

    plt.tight_layout()  # Adjust spacing between subplots
    plt.show()
    fig.savefig(tsv_file.split('.')[0] + '.pdf', 
                pad_inches=1, 
                bbox_inches='tight')

# --------------------------------------------------------------------------------------
if __name__ == '__main__':
    tsv_file = sys.argv[0]
    plot_varint_count_analysis(tsv_file)
