import pandas as pd
import json


class TreeNode:
    def __init__(self, id, name, group_size, duplication="Y", node_type="node", children=None):
        self.id = id  # unique identifier
        self.name = name  # name appearing in the tree
        self.group_size = group_size
        self.duplication = duplication
        self.type = node_type
        self.children = children if children is not None else []

    def to_dict(self):
        return {
            "ID": self.id,
            "duplication": self.duplication,
            "type": self.type,
            "name": self.name,
            "group_size": self.group_size,
            "children": [child.to_dict() for child in self.children]
        }

    def add_child(self, child):
        self.children.append(child)


def prepare_tree_df():
    # read tree df
    tree_df = pd.read_csv('tree_df.csv', dtype={'ID': 'string', 'Parent_ID': 'string'}, low_memory=False)
    tree_df['Parent_ID'] = tree_df['Parent_ID'].fillna('-1')

    # read metadata df
    file_info = pd.read_csv('../../data/microbiome_masst_table.tsv', dtype={'ID': 'string'}, sep='\t', low_memory=False)

    # Count files for each ID
    id_counts = file_info['ID'].value_counts()

    ######## add interventions from file_info to tree_df ########
    file_info.drop_duplicates(subset=['ID'], inplace=True)
    # Create a mapping from ID to interventions
    id_to_interventions = dict(zip(file_info['ID'], file_info['Intervention_specific']))
    # Map interventions to tree_df
    tree_df['Interventions'] = tree_df['ID'].map(id_to_interventions)

    # Function to recursively calculate group size
    def calculate_group_size(node_id, id_counts):
        # node = tree_df[tree_df['ID'] == node_id].iloc[0]
        children = tree_df[tree_df['Parent_ID'] == node_id]

        if children.empty:  # Leaf node
            return id_counts.get(node_id, 0)
        else:
            total = sum(calculate_group_size(child_id, id_counts) for child_id in children['ID'])
            return total

    # Fill in Group_Size for all nodes
    for node_id in tree_df['ID']:
        if pd.notnull(node_id):
            size = calculate_group_size(node_id, id_counts)
            tree_df.loc[tree_df['ID'] == node_id, 'Group_Size'] = size

    # Fill NaN values with 0
    tree_df['Group_Size'] = tree_df['Group_Size'].fillna(0).astype(int)

    # save
    tree_df.to_csv('tree_df_group_size.tsv', sep='\t', index=False)

    return tree_df


def build_tree_from_tsv_and_write_json(output_file_path):

    # Read the tree dataframe
    df = pd.read_csv('tree_df_group_size.tsv', sep='\t')

    def build_tree(df, parent_id='-1'):
        children = df[df['Parent_ID'] == parent_id]
        if children.empty:
            return []

        result = []
        for _, child in children.iterrows():
            child_node = {
                "ID": str(child['ID']),
                "NCBI": str(child['NCBITaxonomy']) if pd.notnull(child['NCBITaxonomy']) else None,
                "Interventions": str(child['Interventions']) if pd.notnull(child['Interventions']) else None,
                "Community_composition": str(child['Community_composition']) if pd.notnull(child['Community_composition']) else None,
                "duplication": "Y",
                "type": "node",
                "name": str(child['SampleType']),
                "group_size": int(child['Group_Size']),
                "children": build_tree(df, child['ID'])
            }
            result.append(child_node)
        return result

    # Build the tree structure
    tree = build_tree(df)[0]

    # Write the tree to a JSON file
    with open(output_file_path, 'w') as f:
        json.dump(tree, f, indent=4)


if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    prepare_tree_df()

    build_tree_from_tsv_and_write_json('../../data/microbiome_masst_tree.json')
