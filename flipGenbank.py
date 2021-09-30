from Bio import SeqIO
from Bio.Seq import Seq
import os
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
import subprocess
import random
import string
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash_extensions import Download
import base64

app = dash.Dash(__name__)
server = app.server
def layout():
    layout = dbc.Container([
        html.H1('flipGenbank'),
        html.Div([
            dcc.Upload(
                id='upload-data_fg',
                children=html.Div([
                    'Drag and Drop or ',
                    html.A('select genbank file')
                ]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
                # Allow multiple files to be uploaded
                multiple=False),
                dcc.Store(id='hidden_flipGenbank'),
            html.Div([html.Button("Download flipped genbank", id="btn_txt_fg"), Download(id="download-text_fg")])
            ])
        ],className='mt-4')
    return layout

def readGenbank(tmp_dir):
    filename = os.path.join(tmp_dir, 'user_genbank.gbk')
    records = list()
    with open(filename) as handle:
        print(handle)
        for record in SeqIO.parse(handle, "genbank"):
            records.append(record)
    return records[0] #only accepting a single record genbank

def flipGenbank(record, tmp_dir):
    #flip fna sequence
    record.seq = record.seq.reverse_complement()
    #flip CDS
    total_len = len(record.seq)
    for feature in record.features:
        #feature locations must have be start < end no matter the strand
        old_strand = feature.location.strand
        x = abs(feature.location.start - total_len)
        y = abs(feature.location.end - total_len)
        location = sorted([x,y])
        feature.location = FeatureLocation(location[0], location[1])
        feature.location._strand = -old_strand
    handle = open( os.path.join(tmp_dir, 'flipped.gbk'), 'w')
    SeqIO.write(record, handle, 'genbank')
    handle.close()
def randomString(stringLength=8):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))
def makedir(working_directory):
    f  =  os.path.join(working_directory, randomString(4))
    print(f)
    os.mkdir(f)
    directories = ['tmp']
    for basename in directories:
        os.mkdir( os.path.join(f, basename) )
    return f

if __name__ == '__main__':
    @app.callback(Output('hidden_flipGenbank','data'),
             [Input('upload-data_fg','contents')],
             [State('upload-data_fg','filename'),
              State('upload-data_fg','last_modified')],
              prevent_initial_call=True,)
    def getGenbank(content, list_of_names, list_of_dates):
        try:
            tmp_dir = makedir('/Users/matt/Desktop/flipGenbank_tmp')
            content_type, content_string = content.split(',')
            file_content = base64.b64decode(content_string)
            with open( os.path.join(tmp_dir, 'user_genbank.gbk'), "w+") as handle:
                handle.write(file_content.decode("utf-8"))
            return tmp_dir
        except ValueError:
            raise dash.exceptions.PreventUpdate

    @app.callback(Output("download-text_fg", "data"),
             [Input("btn_txt_fg", "n_clicks")],
             [State('hidden_flipGenbank', 'data')],
             )
    def downloadFlippedGenbank(n_clicks, tmp_dir):
        print('here')
        record = readGenbank(tmp_dir)
        flipGenbank(record, tmp_dir)
        output_filename = os.path.join(tmp_dir, 'flipped.gbk')
        with open(output_filename, 'r') as file:
            data = file.read()
        return dict(content=data, filename='flipped.gbk')
    app.layout = layout()
    app.run_server(host="0.0.0.0", port="8050")
