import pandas as pd
import os


def main(out_path):

    df = pd.read_csv('file_info.csv', dtype={'ID': 'string'}, low_memory=False)
    df = df[pd.notnull(df['ID'])].reset_index(drop=True)

    df['Filename'] = df['Filepath'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
    df['file_usi'] = df.apply(lambda x: f"mzspec:{x['MassIVE']}:{x['Filename']}", axis=1)

    # read tree df
    tree_df = pd.read_csv('tree_df.csv', dtype={'ID': 'string', 'Parent_ID': 'string'}, low_memory=False)

    ### All mzML files should be linked to the lowest leaf nodes
    # remove files with IDs not in the lowest leaf nodes
    parent_ids = list(set(tree_df['Parent_ID']))
    parent_ids = [x for x in parent_ids if pd.notnull(x)]
    df = df[~df['ID'].isin(parent_ids)].reset_index(drop=True)
    
    # dereplicate by file_usi and ID
    df = df.drop_duplicates(subset=['file_usi', 'ID']).reset_index(drop=True)

    df.to_csv(out_path, sep='\t', index=False)


if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    main(out_path='../../data/microbiome_masst_table.tsv')
