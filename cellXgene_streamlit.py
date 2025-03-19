import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
from streamlit_option_menu import option_menu
from st_aggrid import AgGrid 

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import plotly.express as px

import pyarrow as pa
import pyarrow.parquet as pq
import sys

import socket
import io

from constants import TOTAL_STREAMLIT_TISSUES, TOTAL_STREAMLIT_CELL_TYPES, TOTAL_STREAMLIT_UNIQ_CELLS



st.set_page_config(layout="wide", page_title="Xencor Single Cell Atlas", page_icon=":dna:", initial_sidebar_state="expanded")
#set dark theme


FONT_SIZE=18
path= Path(__file__).parent
style_css_file= path / "style.css"

def _make_scatter_plot_cellXtissue_subplot(df,  gene_selection, scaling_factor=3.0,):
    """
    
    given a dataframe df with mean expression and percentage of cells expressing a gene,
    make a scatter plot with size of circle proportional to percentage of cells expressing a gene
    and color of circle proportional to mean expression of a gene
    
    Each subplot will correspond to tissue type



    """

    columns = df.columns
   
   
    assert 'gene_symbol' in columns, "check if 'gene_symbol' is in columns"
    assert 'organ' in columns, "check if 'organ' is in columns"
    assert 'pct_of_expressed' in columns, "check if 'pct_of_expressed' is in columns"
    assert 'average_expression_of_expressed' in columns, "check if 'average_expression_of_expressed' is in columns"
    assert 'cell_type' in columns, "check if 'cell_type' is in columns"


    global_max_ln_cppt = df['average_expression_of_expressed'].max()

    
    unique_genes = df['gene_symbol'].unique()
    unique_cell_types = df['cell_type'].unique()
    unique_tissues = df['organ'].unique()

  

    # Determine the number of rows and columns for subplots
    num_rows = len(unique_tissues)
    num_cols = 1
   

    # Create subplots
    factor= len(unique_cell_types)
   
    fig_height = 9 + factor // 3.5
   

    
   
    fig_width = 20
   
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(fig_width,fig_height))# sharex=True)
    
    
    # Flatten the axes if necessary
    if num_rows > 1 and num_cols > 1:
        axes = axes.ravel()

    # Plot each group in a subplot
    for i, tissue in enumerate(unique_tissues):
        group_df = df[df['organ'] == tissue]
        # Set the 'feature_name' column to a Categorical with the desired order
        group_df['gene_symbol'] = pd.Categorical(group_df['gene_symbol'], categories=gene_selection, ordered=True)
        #reverse sort by percentage_cells
        group_df = group_df.sort_values(by=['gene_symbol','pct_of_expressed'], ascending=[True,True])

        
        

         
        unique_genes = group_df['gene_symbol'].unique()
        unique_cell_types = group_df['cell_type'].unique()
        cell_types_df = pd.DataFrame(unique_cell_types, columns=['cell'])
        cell_types_df['organ'] = tissue

        if len(unique_tissues) == 1:
            ax=axes
        else:
            ax = axes[i]
        center_value=unique_genes[0]
       
        ax.set_title(tissue, fontsize=25)
        # set the x-axis angle to 90
        ax.tick_params(axis='x', rotation=90)
        #set font size of x and y axis
        #ax.tick_params(axis='x', labelsize=font_size)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=15)
        ax.grid(True, alpha = 0.3)
        
        
        
        norm = colors.Normalize(vmin=0, vmax=global_max_ln_cppt)

        scatter_plot = ax.scatter(group_df['gene_symbol'], group_df['cell_type'],
                                s=group_df['pct_of_expressed'] * scaling_factor,
                                c=group_df['average_expression_of_expressed'],
                                #c=group_df['percentage_cells'] * scaling_factor,
                                cmap='Reds', norm=norm

                                )
      
        plt.ylim(-0.5, len(unique_cell_types) - 0.5)
        # add color bar to ax

        cbar = plt.colorbar(scatter_plot, ax=ax)
        
        
        
        cbar.ax.tick_params(labelsize=20)
        cbar.set_label('Ln(CPPT+1)', fontsize=20)

        legend_labels = ['10', '25', '50', '75', '100']  # Original values
    
        legend_scatter = ax.scatter([center_value] * 5, [5, 5, 5, 5, 5], sizes=[10 * scaling_factor, 
                                                                                25 * scaling_factor,
                                                                                50 * scaling_factor,
                                                                                75 * scaling_factor,
                                                                                100 * scaling_factor], visible=False)
        legend_elements = legend_scatter.legend_elements(prop="sizes")
        ax.legend(legend_elements[0], legend_labels, loc="center left", title="Percentage", bbox_to_anchor=(1.35,.50),  fontsize=20, title_fontsize=20)
       
     
    plt.tight_layout()
    return fig




def _load_parquet_file(tissue_list:list, gene_list:list, n_cell_types:int) -> pd.DataFrame:  
    """
    Load parquet file with mean expression from Cell x Gene Census

    Note for future implementation, we can move the parquet file to the rsconnect server and load the parquet file from the server.
    This would save time for deployment with rsconnect server.
    """
    
    parquet_path= Path(__file__).parent /  "Parquet" / "20240121"
    
    n_obs_cutoff=50

    filters = [
            ('organ', 'in', tissue_list),
            ('gene_symbol', 'in', gene_list),
            ('pct_of_expressed', '>', 0),
            ('total_num_of_cells', '>', n_obs_cutoff)
            
           ]
    

   
    
    res=pd.read_parquet(parquet_path, filters=filters)
    
    n_cell_types=int(n_cell_types)
    # return the top n cell types by percentage of cells expressing a gene
    top_n_cell_types = res.groupby(['gene_symbol', 'organ']).apply(lambda x: x.nlargest(n_cell_types, 'pct_of_expressed')).reset_index(drop=True)
    return top_n_cell_types

    

def _load_tissue():
    """ 
    Load the parquet path and extract the unique tissue and cell types
    The output will be used to make the dropdown menu for tissue.
     
    """
    path = Path(__file__).parent
    file_path=path / "CELLxGENE_gene_expression_032124_selected.csv"
    #parquet_path =  path / "Parquet" /"20240121"
    res=pd.read_csv(file_path)
    res=res[[ 'organ']].drop_duplicates()
    return res

def _load_gene_names():
    """ 
    Load protein coding gene symbols, these are the genes that users can select from the dropdown menu
    
    """
    path = Path(__file__).parent
    gene_names_file = path / "HGNC_protein_coding_genes.txt"
    gene_names_df = pd.read_csv(gene_names_file, sep="\t")
    #rename column approved symbol to symbol
    gene_names_df.rename(columns={'Approved symbol':'gene_symbol'}, inplace=True)
    gene_names_df=gene_names_df.query('gene_symbol == "MS4A1"')
    return gene_names_df[['gene_symbol']]

tissue_cell_type=_load_tissue()
gene_names=_load_gene_names()

#with open(style_css_file) as f:
#    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)





def app_menu():
    

    with st.sidebar:
        selected = option_menu(
                    menu_title = "Xencor Single Cell Atlas",
                    menu_icon= "map-fill",
                    options=['Main',
                            'Xencor Single Cell Expression',
                             ],
                    icons=['book-fill',
                            'collection-fill',
                            
                        ],
                    default_index=0
                )
    return selected

    


def markdown_introduction():

    """
    Brief introduction to Xencor Single Cell Atlas with helpful links
    """
    
    st.markdown("""
                ### Xencor Single Cell Atlas 
                Xencor Single Cell Atlas ingests single cell expression data of **normal subjects** from [CZ Cell X Gene Discover](https://cellxgene.cziscience.com/) by using [CZ Cell X Gene Census API](https://chanzuckerberg.github.io/cellxgene-census//). 
                Users have the capability to search based on organ type and gene symbol(s) and observe the top N cell types expressing the gene(s).

                The analysis involves the inclusion of cells from assays measuring gene expression without the need for gene-length normalization.
                Read counts are normalized using a log transformation (ln([CPTT](https://www.youtube.com/watch?v=TTUrtCY2k-w)+1)). This normalization mitigates batch effects but does not completely
                eliminate them. See the [CellXGene preprint](https://doi.org/10.1101/2023.10.30.563174) for a detailed discussion.

                After normalization, gene/cell combinations with counts less than or equal to 3 are treated as missing data.
                The [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl), representing the hierarchical relationship between cell types, is used for a rollup operation to aggregate expression values and cell counts for a specific cell type and its descendants.
                This operation accommodates variations in granularity across datasets, offering a more robust measure of average expression for terms like "B cell" or "T cell" that may be represented by multiple cell types in different datasets.


                Please note we are in progress of adding more data to Xencor Single Cell Atlas and will be adding more features to this app. Stay tuned!

                Helpful links:
                - [CZ Cell X Gene Discover](https://cellxgene.cziscience.com/)
                - [CZ Cell X Gene Census API](https://chanzuckerberg.github.io/cellxgene-census//)
                - [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl)
                - [Owlready2](https://owlready2.readthedocs.io/en/latest/intro.html) (ontology manipulation library used in CellXGene Census API)
                - [CellXGene preprint](https://doi.org/10.1101/2023.10.30.563174)
                - [CellXGene Data Processing](https://cellxgene.cziscience.com/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_3__Gene%20Expression%20Data%20Processing)




                

                """, unsafe_allow_html=True)


def cellXGene_mean_expression():

    """
    Display KPIs and make a scatter plot of cellXtissue with mean expression and percentage of cells expressing a gene based on user input collected from the form

    """

    placeholder = st.empty()
    with placeholder.container():
    # create three columns
        kpi1, kpi2, kpi3 = st.columns(3)

        # fill in those three columns with respective metrics or KPIs 
        kpi1.metric(label="Total Cells", value= TOTAL_STREAMLIT_UNIQ_CELLS)
        kpi2.metric(label="Total Tissues", value= TOTAL_STREAMLIT_TISSUES)
        kpi3.metric(label="Total Cell types", value= TOTAL_STREAMLIT_CELL_TYPES )
        
    

    #make a dropdown menu to select organ
    tissue_selections = tissue_cell_type['organ'].unique().tolist()
    tissue_selections.sort()
    tissue_selections=['blood']
        

    #cell_selections = tissue_cell_type['cell_type'].unique().tolist()
    #cell_selections.sort()

    gene_selections = gene_names['gene_symbol'].unique().tolist()

  

    with st.form(key='census_selections'):
        gene_index_default=gene_selections.index('MS4A1')
        
        tissue_index_default=tissue_selections.index('blood')
        #make multi-select dropdown menu to select tissue
        tissue_selection = st.multiselect('Select tissue', (tissue_selections), default=tissue_selections[tissue_index_default])
        # add text box to enter top n cell types
        
        
        # make multi-select dropdown menu to select gene
        gene_selection = st.multiselect('Select genes', (gene_selections), default=gene_selections[gene_index_default])
        n_cell_types = st.text_input('Enter the top N number of cell types (by percentage of cells expressing gene)', '10')

        submitted = st.form_submit_button('Submit')

        if submitted:
           
            
            res=pd.read_csv("CELLxGENE_gene_expression_032124_selected.csv")
            #res=_load_parquet_file(tissue_selection, gene_selection, n_cell_types)
            
            fig=_make_scatter_plot_cellXtissue_subplot(res, gene_selection, scaling_factor=2.5)
            st.pyplot(fig)
            
            AgGrid(res, height=500, width='100%')

  
def main():

    
    selected_app = app_menu()
    if selected_app == 'Main':
        markdown_introduction()
    if selected_app == 'Xencor Single Cell Expression':
        st.markdown("""## Xencor Single Cell  Expression""")
        cellXGene_mean_expression()



if __name__ == '__main__':
    main()