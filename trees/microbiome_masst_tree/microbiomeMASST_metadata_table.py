import pandas as pd


def main():

    df = pd.read_csv('file_info.tsv', sep='\t')
    df = df[pd.notnull(df['ID'])].reset_index(drop=True)
    df['ID'] = df['ID'].astype(int)

    df['Filename'] = df['Filename'].apply(lambda x: x.split('.mzML')[0])
    df['file_usi'] = df.apply(lambda x: f"mzspec:{x['MassIVE']}:{x['Filename']}", axis=1)

    # read tree df
    tree_df = pd.read_csv('tree_df.tsv', sep='\t')
    tree_df['ID'] = tree_df['ID'].astype(int)
    tree_df['node_name'] = tree_df.apply(lambda x: x['SampleType'] + ' [' + str(x['Interventions']) + ']' if pd.notnull(x['Interventions']) else x['SampleType'], axis=1)

    ### All mzML files should be linked to the lowest leaf nodes
    # remove files with IDs not in the lowest leaf nodes
    parent_ids = list(set(tree_df['Parent_ID']))
    parent_ids = [int(x) for x in parent_ids if pd.notnull(x)]
    df = df[~df['ID'].isin(parent_ids)].reset_index(drop=True)

    # dict from ID to node_name
    id_to_sample_type = dict(zip(tree_df['ID'], tree_df['node_name']))

    df['node_id'] = df['ID'].map(id_to_sample_type)

    df.to_csv('../../data/microbiome_masst_table.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
