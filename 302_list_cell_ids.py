import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv("data/cells_meta.csv", index_col=0)
    indices = df.index.tolist()
    with open("data/cell_ids.txt", "w") as f:
        for cell_id in indices:
            f.write(f"{cell_id}\n")
    print("Saved: cell_list/all_cells.txt")