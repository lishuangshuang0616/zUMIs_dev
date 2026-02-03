import plotly.express as px
import pandas as pd
import plotly as py
import plotly.graph_objs as go
from plotly.io import *
from plotly.subplots import make_subplots
import collections,math
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.stats import gaussian_kde
from scipy import stats, spatial as spatial

PLOTLYCONFIG = {
    'modeBarButtonsToRemove': [
        "autoScale2d",
        "hoverClosestCartesian",
        "hoverCompareCartesian",
        "lasso2d",
        "zoomIn2d",
        "zoomOut2d",
        "sendDataToCloud",
        "toggleSpikelines" ,
        "logo"
        ],
    'displaylogo': False,
    }

def draw_and_save_plot(figures, filenames, div_filenames, outdir):
    for fig, filename, div_filename in zip(figures, filenames, div_filenames):
        py.offline.plot(
            fig,
            filename=outdir + '/' + filename,
            auto_open=False,
            config=PLOTLYCONFIG
        )
        figplot = py.offline.plot(
            fig,
            include_plotlyjs=False,
            show_link=False,
            output_type='div',
            config=PLOTLYCONFIG
        )
        with open(outdir + '/' + div_filename, 'w') as fw:
            fw.write(figplot)


# barcoderanks plot density
def segment_log_plot(y_data, x_start, x_end):
    log_max_x = np.log(len(y_data))
    log_max_y = np.log(max(y_data))
    segment_len = 0.0
    segment_idx = [x_start]
    for i in range(x_start, x_end):
        last_i = max(x_start, i-1)
        dx = (np.log(i) - np.log(last_i)) / log_max_x
        dy = (np.log(y_data[i]) - np.log(y_data[last_i])) / log_max_y
        segment_len += np.linalg.norm([dx, dy])
        if segment_len >= 0.02 and i > (segment_idx[-1] + 20):
            segment_idx.append(i+1)
            segment_len = 0.0
    if segment_idx[-1] != x_end:
        segment_idx.append(x_end)
    return segment_idx

# barcoderanks plot density
def plot_cmap(density):
    plot_colors =  [
        "#DDDDDD","#D6D9DC","#CFD6DB","#C8D3DA","#C1D0D9",
        "#BACDD9","#B3C9D8","#ACC6D7","#A5C3D6","#9EC0D6",
        "#97BDD5","#90BAD4","#89B6D3","#82B3D3","#7BB0D2",
        "#74ADD1","#6DAAD0","#66A6CF","#5FA3CF","#58A0CE",
        "#539DCC","#4F99CA","#4C95C8","#4992C6","#458EC3",
        "#428AC1","#3F87BF","#3B83BD","#3880BA","#347CB8",
        "#3178B6","#2E75B4","#2A71B1","#276DAF","#236AAD",
        "#2066AB","#1D62A8","#195FA6","#165BA4","#1358A2"
        ]
    levels = len(plot_colors)
    ind = min(levels - 1, int(math.floor(levels * density)))
    return plot_colors[ind]

def downsample_scatterplot_by_density(points_df, npoints, dim1, dim2):
    if len(points_df) <= npoints:
        return points_df
    dim1_log2 = np.log2(points_df[dim1] + 1)
    dim1_z = stats.zscore(dim1_log2)
    dim2_log2 = np.log2(points_df[dim2] + 1)
    dim2_z = stats.zscore(dim2_log2)
    round_df = pd.DataFrame(
        {"z1": np.round(dim1_z, 2), "z2": np.round(dim2_z, 2)}, index=points_df.index
    )
    np.random.seed(0)
    is_dup = round_df.duplicated()
    ind_unique = round_df.index[is_dup == False]
    ind_dup = round_df.index[is_dup]
    if len(ind_unique) <= npoints:
        samp_dups = np.random.choice(ind_dup, size=npoints - len(ind_unique), replace=False)
        return pd.concat([points_df.loc[ind_unique], points_df.loc[samp_dups]])
    tree = spatial.KDTree(round_df.loc[ind_unique])
    radius = 0.1
    neighbors = tree.query_ball_tree(tree, radius)
    frequency = np.array([len(x) for x in neighbors])
    inv_density = radius ** 2 / frequency

    samp_index = np.random.choice(
        round_df.loc[ind_unique].index,
        size=npoints,
        replace=False,
        p=inv_density / sum(inv_density),
    )
    return points_df.loc[samp_index]

def create_saturation_plot(x, y1, y2, y1_label, y2_label, width, height):
    if len(x) > 2:
        xnew = np.linspace(x.min(), x.max(), 50)
        y1new = make_interp_spline(x, y1)(xnew)
        y2new = make_interp_spline(x, y2)(xnew)
    else:
        xnew, y1new, y2new = x, y1, y2

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(
        go.Scatter(x=xnew, y=y1new, name=y1_label, line=dict(color="blue", width=3)),
        secondary_y=False
    )

    fig.add_trace(
        go.Scatter(x=xnew, y=y2new, name=y2_label, line=dict(color="red", width=3)),
        secondary_y=True
    )

    fig.update_layout(
        autosize=False,
        width=width,
        height=height,
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(gridcolor='lightgray', title="Mean Reads per Cell"),
    )

    fig.update_yaxes(title_text=y1_label, gridcolor='lightgray', zeroline=True, zerolinewidth=1, zerolinecolor='gray', secondary_y=False)
    fig.update_yaxes(title_text=y2_label, gridcolor='lightgray', zeroline=True, zerolinewidth=1, zerolinecolor='gray', secondary_y=True)
    
    fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='gray')
    
    return fig


# plotly all
class plotly_summary:

    @staticmethod
    # df : barcode,UMI,iscell
    def _plot_barcoderanks(dataframe_df, well_df, width, height):
        dataframe_df['is_cell_barcode'] = dataframe_df['CB'].apply(
            lambda x: 1 if (x in well_df['internal_barcodes'].values or x in well_df['umi_barcodes'].values) else 0
        )

        dataframe_df = dataframe_df.sort_values(by="reads" , ascending=False)
        dataframe_df = dataframe_df.reset_index(drop=True)
        dataframe_df['New']=dataframe_df.index
        cell_bc = np.array(dataframe_df[dataframe_df['is_cell_barcode'] == 1].index)
        sorted_bc = np.array(dataframe_df.index)
        sorted_counts = np.array(dataframe_df['reads'])
        total_bc = len(sorted_bc)
        
        if len(cell_bc) == 0:
            ix1 = dataframe_df.drop_duplicates('is_cell_barcode',keep='first').index[0]
            ix2 = dataframe_df.drop_duplicates('is_cell_barcode',keep='last').index[0]
        else:
            ix1 = dataframe_df.drop_duplicates('is_cell_barcode',keep='first').index[1]-1
            ix2 = dataframe_df.drop_duplicates('is_cell_barcode',keep='last').index[0]
        plot_segments = []
        barcodeSegment = collections.namedtuple(
            'barcodeSegment', 
            ['start', 'end', 'density', 'legend']
            )

        plot_segments.append(barcodeSegment(
            start=0, end=ix1, density=1.0, legend=True))
        plot_segments.append(barcodeSegment(
            start=ix2+1, end=total_bc, density=0.0, legend=True))

        mixed_segments = segment_log_plot(
            sorted_counts, ix1, ix2
            )
        for i in range(len(mixed_segments) - 1):
            num_cells = sum(
                [1 for j in range(mixed_segments[i], mixed_segments[i + 1]) if sorted_bc[j] in cell_bc])
            
            density = float(num_cells)/float(mixed_segments[i + 1]-mixed_segments[i])
            plot_segments.append(
                barcodeSegment(
                    start=mixed_segments[i], end = mixed_segments[i + 1], density=density, legend=False
                )
            )

        plot_data = []
        for plot_segment in plot_segments:
            start = max(0, plot_segment.start - 1)
            end = plot_segment.end
            selct_count = dataframe_df[start:end]
            dp_first = set(selct_count[selct_count[["reads"]].duplicated(keep="first")].index)
            dp_last = set(selct_count[selct_count[["reads"]].duplicated(keep="last")].index)
            dp_inter = dp_first & dp_last
            selct_count=selct_count.drop(list(dp_inter),axis=0)
            x = list(selct_count['New'])
            y = list(selct_count['reads'])
            name = 'TRUE' if plot_segment.density > 0 else 'NOISE'
            if plot_segment.density > 0:
                n_barcodes = plot_segment.end - plot_segment.start
                n_cells = int(round(plot_segment.density * n_barcodes))
                hover = "{:.0f}% Cell<br>({}/{})".format(100 * plot_segment.density, n_cells, n_barcodes)
            else:
                hover = "NOISE"

            data_dict = {
                "x": x,"y": y,"name": name, "hoverinfo": "text",
                "text": hover,"type": "scattergl","mode": "lines",
                "line": {
                    "width": 3,
                    "color": plot_cmap(plot_segment.density),
                    },
                "showlegend": plot_segment.legend,
                }
            plot_data.append(data_dict)

        plotly_data = [
            go.Scatter(
                x=dat['x'], y=dat['y'], name=dat['name'], mode=dat['mode'], 
                showlegend=dat['showlegend'],
                marker={
                    'color': dat['line']['color']
                    }, 
                line=dat['line'], text=dat['text']
                ) for dat in plot_data
                ]
        layout = go.Layout(
            xaxis = dict(
                type="log", gridcolor="lightgrey", title="Barcode in Rank-descending Order",
                color="black", showline=True, zeroline=True, linewidth=1, fixedrange= True,
                linecolor="black"),
            yaxis = dict(
                type="log", title="reads counts", gridcolor="lightgrey",
                linewidth=1, fixedrange= True, color="black", linecolor="black"
                ),
            height= height, width= width,
            plot_bgcolor='rgba(0,0,0,0)',hovermode='closest',paper_bgcolor='white',
            legend = dict(
                x=1,y=1,traceorder="normal",
                font = dict(
                    family="Arial",size=12,color="black"
                    ),
                bordercolor="Black",borderwidth=0),
            margin = dict(l=0,r=0,b=0,t=0,pad=1),
            font = dict(size=10))
        fig = go.Figure(
            data=plotly_data, layout=layout
            )

        return fig

    
    @staticmethod
    def _plot_saturation(saturantion_df, width, height, type):
        x = saturantion_df['downsample']
        y1 = saturantion_df['median_gene']
        y2 = saturantion_df['median_umis']
        
        if type == "reads_saturation":
            fig = create_saturation_plot(x, y1, y2, "Genes Number", "Reads Number", width, height)
            return fig
        elif type == "umi_saturation":
            fig = create_saturation_plot(x, y1, y2, "Genes Number", "UMIs Number", width, height)
            return fig
        else:
            raise Exception('Unrecognized.')
        

    @staticmethod
    def _plot_stack(stat, width, height):
        featColors =  {
            "Exon": "#1a5084",
            "Intron+Exon": "#914614",
            "Intron": "#118730",
            "Unmapped": "#545454",
            "Ambiguity": "#ffa54f",
            "Intergenic": "#ffd700",
            "Unused_BC": "#bbbbbb",
            "User": "#cd2626"
        }

        data = {
            'Exon': stat['exon_reads'],
            'Intron': stat['intron_reads'],
            'Intergenic': stat['intergenic_reads'],
            'Ambiguity': stat['ambiguity_reads'],
            'Unmapped': stat['unmapped_reads'],
            'Unused_BC': (stat['raw_reads']-stat['white_reads'])
        }

        # Calculate total reads and percentages
        total_reads = sum(data.values())
        percentages = {category: round((value / total_reads),5) * 100 for category, value in data.items()}

        # Create a stacked bar chart
        fig = go.Figure()

        # Add each category as a segment of the stacked bar
        for category, value in data.items():
            percentage = percentages[category]
            color = featColors.get(category, 'black') 
            fig.add_trace(go.Bar(
                name=category,
                y=[''],  # Keep y as an empty string to hide "Reads" on y-axis
                x=[percentage],
                marker_color=color,
                #text=[f"{percentage:.1f}%"],  # Display percentage as text
                textposition='auto',
                orientation='h',
                hovertemplate=f"{category}: {value} reads<br>{percentage:.1f}% of total"  # Custom hover text
            ))

        # Update layout for the stacked bar chart
        fig.update_layout(
            autosize=False,
            barmode='stack',
            plot_bgcolor='rgba(0,0,0,0)',
            hovermode='closest',
            paper_bgcolor='white',
            margin=dict(l=0, r=0, b=0, t=0, pad=1),
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1
            ),
            height=height,
            width=width,
            xaxis=dict(
                title='Percentage of Reads (%)',
                tickformat='.1f'  # Format x-axis as percentages
            ),
            yaxis=dict(
                showticklabels=False  # Hide y-axis labels (the "Reads" text)
            )
        )

        return fig
    

    @staticmethod
    def _plot_384plate(df, width, height):
        # 提取初始的孔板数据，假设默认显示第一个plate
        df['Row'] = df['wellID'].str.extract('([A-Z])')
        df['Column'] = df['wellID'].str.extract('(\d+)').astype(int)
        plate_names = df['plate'].unique()

        # 默认绘制第一个plate的384孔板图
        initial_plate = plate_names[0]
        filtered_df = df[df['plate'] == initial_plate]

        fig = px.scatter(filtered_df, x='Column', y='Row', color='mappingratio',
                        color_continuous_scale='Viridis', 
                        title=f'384-Well Plate Plot: MappingRatio (Plate: {initial_plate})', 
                        labels={'Column':'Column', 'Row':'Row'},
                        hover_name='wellID', 
                        size_max=10)

        fig.update_traces(
            marker=dict(size=20, line=dict(width=1, color='DarkSlateGrey')),
        )
        fig.update_yaxes(categoryorder='array', categoryarray=[chr(i) for i in range(ord('P'), ord('A')-1, -1)])
        fig.update_xaxes(dtick=1)

        # 更新布局
        fig.update_layout(
            height=height,
            width=width,
            plot_bgcolor='rgba(0,0,0,0)',
            template='plotly_white',
            margin=dict(l=0, r=0, b=0, pad=1),
            coloraxis_colorbar=dict(
                title=''
            )
        )

        # 添加updatemenu，用于切换颜色映射和plate
        fig.update_layout(
            updatemenus=[
                dict(
                    buttons=list([
                        dict(label="MappingRatio",
                            method="update",
                            args=[{"marker.color": [filtered_df['mappingratio']]},
                                {"title": f"384-Well Plate Plot: MappingRatio (Plate: {initial_plate})"}]),
                        dict(label="ExonIntronRatio",
                            method="update",
                            args=[{"marker.color": [filtered_df['exonintronratio']]},
                                {"title": f"384-Well Plate Plot: ExonIntronRatio (Plate: {initial_plate})"}]),
                        dict(label="UMIfrac",
                            method="update",
                            args=[{"marker.color": [filtered_df['UMIfrac']]},
                                {"title": f"384-Well Plate Plot: UMIfrac (Plate: {initial_plate})"}]),
                        dict(label="AllReads",
                            method="update",
                            args=[{"marker.color": [filtered_df['AllReads']]},
                                {"title": f"384-Well Plate Plot: AllReads (Plate: {initial_plate})"}]),
                        dict(label="Genes",
                            method="update",
                            args=[{"marker.color": [filtered_df['Intron_Exon_genes']]},
                                {"title": f"384-Well Plate Plot: Genes (Plate: {initial_plate})"}])
                    ]),
                    direction="down",
                    showactive=True,
                    x=0.8,
                    y=1.15,
                    xanchor="left",
                    yanchor="top"
                ),
                dict(
                    buttons=[
                        dict(label=plate,
                            method="update",
                            args=[
                                {"x": [df[df['plate'] == plate]['Column']],
                                "y": [df[df['plate'] == plate]['Row']],
                                "marker.color": [df[df['plate'] == plate]['mappingratio']]},
                                {"title": f"384-Well Plate Plot: MappingRatio (Plate: {plate})"}
                            ])
                        for plate in plate_names
                    ],
                    direction="down",
                    showactive=True,
                    x=0.2,
                    y=1.15,
                    xanchor="left",
                    yanchor="top"
                )
            ]
        )

        return fig
    


    @staticmethod
    def _plot_384plate_histogram(plate384, width, height):
        df = plate384
        initial_column = 'mappingratio'
        fig = px.histogram(df, x=initial_column, nbins=120,
                   title=f'Distribution of {initial_column}',
                   color_discrete_sequence=['#1a5084'],
                   labels={initial_column: 'Values'})
        
        fig.update_layout(
            height=height,
            width=width,
            plot_bgcolor='rgba(0,0,0,0)',
            template='plotly_white',
            updatemenus=[
                dict(
                    buttons=[
                        dict(label="MappingRatio",
                            method="update",
                            args=[{"x": [df['mappingratio']]},
                                {"title": "Distribution of: MappingRatio"},
                                {"coloraxis.colorbar.title": "Mapping Ratio"}],
                            execute=True),
                        dict(label="ExonIntronRatio",
                            method="update",
                            args=[{"x": [df['exonintronratio']]},
                                {"title": "Distribution of: ExonIntronRatio"},
                                {"coloraxis.colorbar.title": "Exon Intron Ratio"}],
                            execute=True),
                        dict(label="UMIfrac",
                            method="update",
                            args=[{"x": [df['UMIfrac']]},
                                {"title": "Distribution of: UMIfrac"},
                                {"coloraxis.colorbar.title": "UMI Fraction"}],
                            execute=True),
                        dict(label="AllReads",
                            method="update",
                            args=[{"x": [df['AllReads']]},
                                {"title": "Distribution of: AllReads"},
                                {"coloraxis.colorbar.title": "All Reads"}],
                            execute=True),
                        dict(label="Genes",
                            method="update",
                            args=[{"x": [df['Intron_Exon_genes']]},
                                {"title": "Distribution of: Intron_Exon_genes"},
                                {"coloraxis.colorbar.title": "Genes"}],
                            execute=True)
                    ],
                    direction="down",
                    showactive=True,
                    x=0.8,
                    y=1.15,
                    xanchor="left",
                    yanchor="top"
                )
            ]
        )

        return fig

    # @staticmethod
    # def _plot_384plate_histogram(plate384, width, height):
    #     df = plate384
    #     initial_column = 'mappingratio'
    #     plate_names = df['plate'].unique()
    #     initial_plate = plate_names[0]
    #     filtered_df = df[df['plate'] == initial_plate]

    #     fig = px.histogram(filtered_df, x=initial_column, nbins=120,
    #             title=f'Distribution of {initial_column} (Plate: {initial_plate})',
    #             color_discrete_sequence=['#1a5084'],
    #             labels={initial_column: 'Values'})
        
    #     fig.update_layout(
    #         height=height,
    #         width=width,
    #         plot_bgcolor='rgba(0,0,0,0)',
    #         template='plotly_white',
    #         updatemenus=[
    #             dict(
    #                 buttons=[
    #                     dict(label="MappingRatio",
    #                         method="update",
    #                         args=[{"x": [filtered_df['mappingratio']]},
    #                             {"title": f"Distribution of MappingRatio (Plate: {initial_plate})"}],
    #                         execute=True),
    #                     dict(label="ExonIntronRatio",
    #                         method="update",
    #                         args=[{"x": [filtered_df['exonintronratio']]},
    #                             {"title": f"Distribution of ExonIntronRatio (Plate: {initial_plate})"}],
    #                         execute=True),
    #                     dict(label="UMIfrac",
    #                         method="update",
    #                         args=[{"x": [filtered_df['UMIfrac']]},
    #                             {"title": f"Distribution of UMIfrac (Plate: {initial_plate})"}],
    #                         execute=True),
    #                     dict(label="AllReads",
    #                         method="update",
    #                         args=[{"x": [filtered_df['AllReads']]},
    #                             {"title": f"Distribution of AllReads (Plate: {initial_plate})"}],
    #                         execute=True),
    #                     dict(label="Genes",
    #                         method="update",
    #                         args=[{"x": [filtered_df['Intron_Exon_genes']]},
    #                             {"title": f"Distribution of Genes (Plate: {initial_plate})"}],
    #                         execute=True)
    #                 ],
    #                 direction="down",
    #                 showactive=True,
    #                 x=0.8,
    #                 y=1.15,
    #                 xanchor="left",
    #                 yanchor="top"
    #             ),
    #             dict(
    #                 buttons=[
    #                     dict(label=plate,
    #                         method="update",
    #                         args=[
    #                             {"x": [df[df['plate'] == plate][initial_column]]},
    #                             {"title": f"Distribution of {initial_column} (Plate: {plate})"}
    #                         ])
    #                     for plate in plate_names
    #                 ],
    #                 direction="down",
    #                 showactive=True,
    #                 x=0.2,
    #                 y=1.15,
    #                 xanchor="left",
    #                 yanchor="top"
    #             )
    #         ]
    #     )

    #     return fig


    

    @staticmethod
    def _plot_density_trace(plate384, width, height):
        df = plate384

        def create_density_trace(y_column):
            x = df['AllReads']
            y = df[y_column]

            xy = np.vstack([x, y])
            z = gaussian_kde(xy)(xy)

            return go.Scatter(
                x=x, y=y, mode='markers',
                marker=dict(size=5, color=z, colorscale='Viridis', showscale=False),  # Hide color scale
                text=df[y_column],
                hoverinfo='text'
            )

        # Create initial plot with 'mappingratio'
        fig = go.Figure()

        # Add initial trace
        initial_trace = create_density_trace('mappingratio')
        fig.add_trace(initial_trace)

        # Add buttons to update y-axis data
        yaxis_buttons = [
            dict(label="MappingRatio",
                method="update",
                args=[{
                    "x": [df['AllReads']],
                    "y": [df['mappingratio']],
                    "marker": [dict(size=5, color=gaussian_kde(np.vstack([df['AllReads'], df['mappingratio']]))(np.vstack([df['AllReads'], df['mappingratio']])), colorscale='Viridis', showscale=False)],
                    "text": [df['mappingratio']]
                },
                    {"title": "Point Density Plot: MappingRatio"}]),
            dict(label="ExonIntronRatio",
                method="update",
                args=[{
                    "x": [df['AllReads']],
                    "y": [df['exonintronratio']],
                    "marker": [dict(size=5, color=gaussian_kde(np.vstack([df['AllReads'], df['exonintronratio']]))(np.vstack([df['AllReads'], df['exonintronratio']])), colorscale='Viridis', showscale=False)],
                    "text": [df['exonintronratio']]
                },
                    {"title": "Point Density Plot: ExonIntronRatio"}]),
            dict(label="UMIfrac",
                method="update",
                args=[{
                    "x": [df['AllReads']],
                    "y": [df['UMIfrac']],
                    "marker": [dict(size=5, color=gaussian_kde(np.vstack([df['AllReads'], df['UMIfrac']]))(np.vstack([df['AllReads'], df['UMIfrac']])), colorscale='Viridis', showscale=False)],
                    "text": [df['UMIfrac']]
                },
                    {"title": "Point Density Plot: UMIfrac"}]),
            dict(label="UMIs",
                method="update",
                args=[{
                    "x": [df['AllReads']],
                    "y": [df['Intron_Exon_umis']],
                    "marker": [dict(size=5, color=gaussian_kde(np.vstack([df['AllReads'], df['Intron_Exon_umis']]))(np.vstack([df['AllReads'], df['Intron_Exon_umis']])), colorscale='Viridis', showscale=False)],
                    "text": [df['Intron_Exon_umis']]
                },
                    {"title": "Point Density Plot: UMIs"}]),
            dict(label="Genes",
                method="update",
                args=[{
                    "x": [df['AllReads']],
                    "y": [df['Intron_Exon_genes']],
                    "marker": [dict(size=5, color=gaussian_kde(np.vstack([df['AllReads'], df['Intron_Exon_genes']]))(np.vstack([df['AllReads'], df['Intron_Exon_genes']])), colorscale='Viridis', showscale=False)],
                    "text": [df['Intron_Exon_genes']]
                },
                    {"title": "Point Density Plot: Genes"}])
        ]

        # Add buttons to toggle x-axis scale
        xaxis_buttons = [
            dict(label="Log Scale",
                method="relayout",
                args=[{"xaxis.type": "log"}]),
            dict(label="Linear Scale",
                method="relayout",
                args=[{"xaxis.type": "linear"}])
        ]

        fig.update_layout(
            title="Point Density Plot: mappingratio",
            height=height,
            width=width,
            plot_bgcolor='rgba(0,0,0,0)',
            template='plotly_white',
            xaxis_title='All reads',
            xaxis_type='log',
            showlegend=False,  # Remove legend
            updatemenus=[
                dict(
                    buttons=yaxis_buttons,
                    direction="down",
                    showactive=True,
                    x=0.8,
                    y=1.12,
                    xanchor="left",
                    yanchor="top"
                ),
                dict(
                    buttons=xaxis_buttons,
                    direction="down",
                    showactive=True,
                    x=0.3,
                    y=1.12,
                    xanchor="left",
                    yanchor="top"
                )
            ]
        )

        return fig

    # @staticmethod
    # def _plot_density_trace(plate384, width, height):
    #     df = plate384
    #     plate_names = df['plate'].unique()
    #     initial_plate = plate_names[0]
    #     filtered_df = df[df['plate'] == initial_plate]

    #     def create_density_trace(y_column, plate_df):
    #         x = plate_df['AllReads']
    #         y = plate_df[y_column]
    #         xy = np.vstack([x, y])
    #         z = gaussian_kde(xy)(xy)

    #         return go.Scatter(
    #             x=x, y=y, mode='markers',
    #             marker=dict(size=5, color=z, colorscale='Viridis', showscale=False),
    #             text=plate_df[y_column],
    #             hoverinfo='text'
    #         )

    #     fig = go.Figure()

    #     initial_trace = create_density_trace('mappingratio', filtered_df)
    #     fig.add_trace(initial_trace)

    #     yaxis_buttons = [
    #         dict(label="MappingRatio",
    #             method="update",
    #             args=[{
    #                 "x": [filtered_df['AllReads']],
    #                 "y": [filtered_df['mappingratio']],
    #                 "marker": [dict(size=5, color=gaussian_kde(np.vstack([filtered_df['AllReads'], filtered_df['mappingratio']]))(np.vstack([filtered_df['AllReads'], filtered_df['mappingratio']])), colorscale='Viridis', showscale=False)],
    #                 "text": [filtered_df['mappingratio']]
    #             },
    #                 {"title": f"Point Density Plot: MappingRatio (Plate: {initial_plate})"}]),
    #         dict(label="ExonIntronRatio",
    #             method="update",
    #             args=[{
    #                 "x": [filtered_df['AllReads']],
    #                 "y": [filtered_df['exonintronratio']],
    #                 "marker": [dict(size=5, color=gaussian_kde(np.vstack([filtered_df['AllReads'], filtered_df['exonintronratio']]))(np.vstack([filtered_df['AllReads'], filtered_df['exonintronratio']])), colorscale='Viridis', showscale=False)],
    #                 "text": [filtered_df['exonintronratio']]
    #             },
    #                 {"title": f"Point Density Plot: ExonIntronRatio (Plate: {initial_plate})"}]),
    #         dict(label="UMIfrac",
    #             method="update",
    #             args=[{
    #                 "x": [filtered_df['AllReads']],
    #                 "y": [filtered_df['UMIfrac']],
    #                 "marker": [dict(size=5, color=gaussian_kde(np.vstack([filtered_df['AllReads'], filtered_df['UMIfrac']]))(np.vstack([filtered_df['AllReads'], filtered_df['UMIfrac']])), colorscale='Viridis', showscale=False)],
    #                 "text": [filtered_df['UMIfrac']]
    #             },
    #                 {"title": f"Point Density Plot: UMIfrac (Plate: {initial_plate})"}]),
    #         dict(label="Genes",
    #             method="update",
    #             args=[{
    #                 "x": [filtered_df['AllReads']],
    #                 "y": [filtered_df['Intron_Exon_genes']],
    #                 "marker": [dict(size=5, color=gaussian_kde(np.vstack([filtered_df['AllReads'], filtered_df['Intron_Exon_genes']]))(np.vstack([filtered_df['AllReads'], filtered_df['Intron_Exon_genes']])), colorscale='Viridis', showscale=False)],
    #                 "text": [filtered_df['Intron_Exon_genes']]
    #             },
    #                 {"title": f"Point Density Plot: Genes (Plate: {initial_plate})"}])
    #     ]

    #     xaxis_buttons = [
    #         dict(label="Log Scale",
    #             method="relayout",
    #             args=[{"xaxis.type": "log"}]),
    #         dict(label="Linear Scale",
    #             method="relayout",
    #             args=[{"xaxis.type": "linear"}])
    #     ]

    #     plate_buttons = [
    #         dict(label=plate,
    #             method="update",
    #             args=[
    #                 {"x": [df[df['plate'] == plate]['AllReads']],
    #                 "y": [df[df['plate'] == plate]['mappingratio']],
    #                 "marker": [dict(size=5, color=gaussian_kde(np.vstack([df[df['plate'] == plate]['AllReads'], df[df['plate'] == plate]['mappingratio']]))(np.vstack([df[df['plate'] == plate]['AllReads'], df[df['plate'] == plate]['mappingratio']])), colorscale='Viridis', showscale=False)],
    #                 "text": [df[df['plate'] == plate]['mappingratio']]},
    #                 {"title": f"Point Density Plot: MappingRatio (Plate: {plate})"}
    #             ])
    #         for plate in plate_names
    #     ]

    #     fig.update_layout(
    #         title=f"Point Density Plot: MappingRatio (Plate: {initial_plate})",
    #         height=height,
    #         width=width,
    #         plot_bgcolor='rgba(0,0,0,0)',
    #         template='plotly_white',
    #         xaxis_title='All reads',
    #         xaxis_type='log',
    #         showlegend=False,
    #         updatemenus=[
    #             dict(
    #                 buttons=yaxis_buttons,
    #                 direction="down",
    #                 showactive=True,
    #                 x=0.8,
    #                 y=1.12,
    #                 xanchor="left",
    #                 yanchor="top"
    #             ),
    #             dict(
    #                 buttons=xaxis_buttons,
    #                 direction="down",
    #                 showactive=True,
    #                 x=0.3,
    #                 y=1.12,
    #                 xanchor="left",
    #                 yanchor="top"
    #             ),
    #             dict(
    #                 buttons=plate_buttons,
    #                 direction="down",
    #                 showactive=True,
    #                 x=0.5,
    #                 y=1.12,
    #                 xanchor="left",
    #                 yanchor="top"
    #             )
    #         ]
    #     )

    #     return fig

    

    @staticmethod
    def _plotly_genebody(data_inter, data_umi ,width, height):
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=data_inter['Gene Body Percentile'],
            y=data_inter['Coverage'],
            mode='lines',
            line=dict(color='#004c99', width=2),
            name='Inter Coverage'
        ))

        fig.add_trace(go.Scatter(
            x=data_umi['Gene Body Percentile'],
            y=data_umi['Coverage'],
            mode='lines',
            line=dict(color='#FF6347', width=2), 
            name='UMI Coverage'
        ))

        fig.update_layout(
            height=height,
            width=width,
            title='Gene Body Coverage Comparison',
            xaxis_title='Gene Body Percentile (5\'->3\')',
            yaxis_title='Coverage',
            plot_bgcolor='rgba(0,0,0,0)',
            template='plotly_white',
            xaxis=dict(showline=True, showgrid=False, showticklabels=True),
            yaxis=dict(showline=True, showgrid=True, showticklabels=True)
        )

        return fig
    

    def _plot_violin_plot(dataframe, width, height):
        meta = dataframe.copy()
        if 'percent.mt' not in meta.columns:
            meta['percent.mt'] = 0
        fig = make_subplots(rows=1, cols=3, vertical_spacing=0.1, horizontal_spacing=0.1)
        fig.add_trace(
            go.Violin(y=meta['nFeature_RNA'], box_visible=True, line_color='black', line_width=2,
                        box_fillcolor = "white",points=False,spanmode = "hard",
                        meanline_visible=True, fillcolor='#005bac', opacity=0.8,
                        x0='genes'),row=1, col=1
        )

        fig.add_trace(
            go.Violin(y=meta['nCount_RNA'], box_visible=True, line_color='black', line_width=2,
                        box_fillcolor = "white",points=False,spanmode = "hard",
                        meanline_visible=True, fillcolor='#005bac', opacity=0.8,
                        x0='counts'),row=1, col=2
        )

        fig.add_trace(
            go.Violin(y=meta['percent.mt'], box_visible=True, line_color='black', line_width=2,
                        box_fillcolor = "white",points=False,spanmode = "hard",
                        meanline_visible=True, fillcolor='#005bac', opacity=0.8,
                        x0='mito.percent'),row=1, col=3
        )

        fig.update_layout(
            autosize=False,
            width=width,
            height=height,
            margin = dict(l=0,r=0,t=0,b=0),
            plot_bgcolor='#F9F9F9',
            title=dict(
                #text="Violin Summary",
                font=dict(
                    family="Arial",
                    color="black"
                ),
                x=0.15,
                y=0.90,
            ),
            showlegend=False
        )
        fig.update_xaxes(
            zeroline=False, zerolinewidth=1, zerolinecolor='lightgrey', showgrid=False , fixedrange= True, gridcolor = 'lightgrey',
        )
        fig.update_yaxes(
            zeroline=False, zerolinewidth=1, zerolinecolor='lightgrey', showgrid=False, fixedrange= True, gridcolor = 'lightgrey',
        )

        return fig


