#Import FISHscale
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*Geopandas not installed.*")
from FISHscale.utils import dataset

#Import other libraries
from pathlib import Path
import pandas as pd
import typer

app = typer.Typer(add_completion=False)
app.pretty_exceptions_enable = False 

@app.command()
def process_data(file_name: Path):

    df_RNA = pd.read_parquet(file_name)
    z_mean, y_mean, x_mean = df_RNA[['global_z', 'global_y', 'global_x']].mean()

    # code from chatGPT to try and put them in the right order
    # desired_gene_order = ["Angpt1", "Crmp1", "Dcx", "Hdac11", "Itgam", "Notch2", "Ptch1", "Serpine1", "Shh","Slc1a3", "Tcf12", "Tek", "Tgfb1", "Tgfbi", "Aif1", "Nnat"]
    # df_RNA['gene_id'] = pd.Categorical(df_RNA['gene_id'],
    #                                 categories=desired_gene_order,
    #                                 ordered=True)
    # df_RNA = df_RNA.sort_values(['gene_id', 'global_z', 'global_y', 'global_x'])

    # # then construct Dataset and visualize as before

    d = dataset.Dataset(
        str(file_name),
        x_label = 'global_x',
        y_label = 'global_y',
        z_label = 'global_z',
        x_offset = -x_mean,
        y_offset = -y_mean,
        z_offset = -z_mean,
        gene_label = 'gene_id',
        pixel_size = '1 micrometer',
        select_valid='cell_id',
        verbose = False,
        reparse = True
    )

    d.visualize()

def main():
    app()

if __name__ == "__main__":
    main()