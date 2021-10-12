# -*- coding: utf-8 -*-
#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Release date: 2021-01-19
# Last update: 2021-09-29

import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata1", required=True, help="Metadata file 1")
    parser.add_argument("--metadata2", required=True, help="Metadata file 2")
    parser.add_argument("--output", required=True, help="Merged metadata file")
    args = parser.parse_args()

    metadata1 = args.metadata1
    metadata2 = args.metadata2
    output = args.output

    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_impacc/nextstrain/run8_20210402_impacc/'
    # metadata1 = path + 'output_files/metadata/base_metadata.tsv'
    # metadata2 = path + 'output_files/assured_data/sample_metadata.tsv'
    # output = path + 'output_files/metadata/metadata_merged.tsv'

    def load_table(file):
        df = ''
        if str(file).split('.')[-1] == 'tsv':
            separator = '\t'
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] == 'csv':
            separator = ','
            df = pd.read_csv(file, encoding='utf-8', sep=separator, dtype='str')
        elif str(file).split('.')[-1] in ['xls', 'xlsx']:
            df = pd.read_excel(file, index_col=None, header=0, sheet_name=0, dtype='str')
            df.fillna('', inplace=True)
        else:
            print('Wrong file format. Compatible file formats: TSV, CSV, XLS, XLSX')
            exit()
        return df


    # nextstrain metadata
    df1 = load_table(metadata1)
    df1.fillna('', inplace=True)

    # Extra metadata
    df2 = load_table(metadata2)
    df2.fillna('', inplace=True)

    # list_columns = [col for col in dfG.columns.to_list() if col in dfN.columns.to_list()]  # list of columns in common
    # list_columns = list(set(dfG.columns.to_list() + dfN.columns.to_list()))
    # print(list_columns)
    # dfG = dfG[list_columns]

    # merge frames
    frames = [df1, df2]
    result = pd.concat(frames)
    result.fillna('', inplace=True)
    result.to_csv(output, sep='\t', index=False)

    print('\nTSV metadata files successfully merged.\n')
