import pandas as pd
import json
from collections import defaultdict
import requests
import io

#make metadata table for masst 
############################################


#get Personal care product MRIs from GNPS
url = "https://datasetcache.gnps2.org/datasette/database/uniquemri.csv"
params = {
    '_filter_column': 'dataset',
    '_filter_op': 'exact',
    '_filter_value': 'MSV000095251',
    '_sort': 'usi',
    '_stream': 'on',
    '_dl': 'on'
}

response = requests.get(url, params=params, stream=True)
content = io.StringIO(response.content.decode('utf-8'))
df_unique_mri_neg = pd.read_csv(content)

params['_filter_value'] = 'MSV000095003'
response = requests.get(url, params=params, stream=True)

content = io.StringIO(response.content.decode('utf-8'))
df_unique_mri_pos = pd.read_csv(content)

#merge positive and negative data
df_unique_mri = pd.concat([df_unique_mri_neg, df_unique_mri_pos], ignore_index=True)

#reformat and filter data
df_unique_mri = df_unique_mri[['usi', 'dataset']]

df_unique_mri = df_unique_mri[df_unique_mri['usi'].str.contains(':peak/')]
df_unique_mri = df_unique_mri[~df_unique_mri['usi'].str.contains('six_mix')]
df_unique_mri = df_unique_mri[~df_unique_mri['usi'].str.contains('low_nitrogen')]
df_unique_mri = df_unique_mri[~df_unique_mri['usi'].str.contains('test')]
df_unique_mri = df_unique_mri[~df_unique_mri['usi'].str.contains('wash')]
df_unique_mri.reset_index(drop=True, inplace=True)

#make id column for this usi dataframe by splitting usi column by / and taking last element
df_unique_mri['Filename'] = df_unique_mri['usi'].str.split('/').str[-1]
df_unique_mri['id'] = df_unique_mri['Filename'].str.split('.').str[0]

#load file associations
df_metadata = pd.read_csv('personalCareProduct_FileAssoc.csv')

#reformat and filter data
df_metadata = df_metadata[['filename', 'Category']]
df_metadata['node_id'] = df_metadata['Category'].str.split(' > ').str[-1]
df_metadata.rename(columns={'filename': 'id'}, inplace=True)


#merge usi dataframe with metadata dataframe on id column and reformat
df_unique_mri = pd.merge(df_unique_mri, df_metadata, on='id', how='inner')
df_unique_mri.rename(columns={'usi': 'file_usi'}, inplace=True)
df_unique_mri.rename(columns={'dataset': 'MassIVE'}, inplace=True)
df_unique_mri.drop(columns=['id'], inplace=True)


df_unique_mri['file_usi'] = df_unique_mri['file_usi'].apply(lambda x: ":".join(x.split(":")[:2]) + ":" + x.split("/")[-1].replace('.mzML', ''))



#save table for special_masst
df_unique_mri.to_csv('../../data/personalCareProduct_masst_table.csv', index=False)
print(df_unique_mri.head())


#make tree for masst
############################################

# Create a dictionary to count the nodes
node_counts = defaultdict(int)

# Split categories and count nodes
for category in df_unique_mri['Category']:
    parts = category.split(' > ')
    for i in range(len(parts)):
        node_counts[parts[i]] += 1

# Create node table with only node names
node_table = pd.DataFrame.from_dict(node_counts, orient='index', columns=['Count']).reset_index()

node_table.columns = ['Node', 'Count']

# Create edge table with full path
edges = []
for category in df_unique_mri['Category']:
    parts = category.split(' > ')
    for i in range(len(parts) - 1):
        edges.append((parts[i], parts[i + 1]))

edge_table = pd.DataFrame(edges, columns=['Source', 'Target']).drop_duplicates()


def build_tree(node_table, edge_table, parent_name=None, parent_id=""):
    if parent_name is None:
        # Find the root node
        root_node = node_table.iloc[0]
        node = {
            "ID": root_node['Node'],
            "duplication": "Y",
            "type": "node",
            "name": root_node['Node'],
            "group_size": int(root_node['Count']),
            "children": build_tree(node_table, edge_table, root_node['Node'], root_node['Node'])
        }
        return [node]
    else:
        children = edge_table[edge_table['Source'] == parent_name]
        if children.empty:
            return []
        result = []
        for _, child in children.iterrows():
            child_node_name = child['Target']
            child_count = node_table[node_table['Node'] == child_node_name]['Count'].values[0]
            child_node = {
                "ID": f"{parent_id}_{child_node_name}",
                "duplication": "Y",
                "type": "node",
                "name": child_node_name,
                "group_size": int(child_count),
                "children": build_tree(node_table, edge_table, child_node_name, f"{parent_id}_{child_node_name}")
            }
            result.append(child_node)
        return result

# Build the tree
tree = build_tree(node_table, edge_table)



with open('../../data/personalCareProduct_masst_tree.json', 'w') as file:
    json.dump(tree[0], file, indent=4)


