from glob import glob
import numpy as np
import plotly.graph_objects as go
import streamlit as st
from streamlit.components.v1 import html
import streamlit.components.v1 as components
from io import BytesIO
import json
import pandas as pd

st.set_page_config(layout="centered", page_title="ALeCS Database", page_icon=":milky_way:")
aweights = {'H':1, 'C':12, 'N':14, 'O':16, 'S':32, 'P':31}
latexMols = {}
latexMap = []
molsLatex = {}
with open("molLatex.txt", "r") as lm:
    lines = lm.readlines()
    for (i, li) in enumerate(lines):
        s = li.split()
        latexMols[s[0]] = s[1]
        latexMap.append(s[1])
        molsLatex[i] = s[0]

mols = []
molsLabels = []
mweights = []
xsData = {}
pdbStrings = {}
for file in glob("BEB/*.dat"):
    dat = np.loadtxt(file)
    moli = file.split("/")[-1].split(".")[0]

    mass = -1
    with open("xyz/"+moli+".xyz", "r") as sf:
        mass = 0
        lines = sf.readlines()
        for li in lines:
            s = li.split()
            mass += aweights[s[0]]
    if mass == -1:
        print("ERROR!!!")
        import sys
        sys.exit()

    mols.append(moli)
    molsLabels.append(latexMols[moli])
    mweights.append(mass)
    xsShape = dat.shape
    xsDi = {}
    xsDi['E'] = dat[:,0]
    if xsShape[-1] == 3:
        xsDi['xs'] = dat[:,-1]*0.280028521
    else:
        xsDi['xs'] = dat[:,1]*0.280028521
    xsData[moli] = xsDi
    pdbi = """%s"""%(open("pdbs/"+moli+".pdb", "r").read())
    pdbStrings[moli] = json.dumps(pdbi)

st.markdown("### Select molecules to show the cross sections")
molChoices = st.multiselect("Molecules", [x for _, x in sorted(zip(mweights, mols))], max_selections=5)
htmlStr = """
        <html>
        <body>
        </body>
        </html>
    """
if len(molChoices) > 0:
    htmlStr = """
        <html>
        <body>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha1/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-GLhlTQ8iRABdZLl6O3oVMWSktQOp6b7In1Zl3/Jr59b6EGGoI1aFkw7cmDA6j6gD" crossorigin="anonymous">
        <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/2.2.4/jquery.min.js"></script>
        <script type="text/javascript" src="https://brandt-gaches.space/dev/ChemDoodleWeb.js"></script> 
        <link rel="stylesheet" href="https://brandt-gaches.space/dev/ChemDoodleWeb.css" type="text/css">
        <script>
    """
    divStrngs = """"""
    jsStrings = """"""
    for (i, moli) in enumerate(molChoices):
        jsStrings += """
            var rotator%d = new ChemDoodle.RotatorCanvas3D('canvas3d%d', 125, 125);
            rotator%d.styles.set3DRepresentation('Ball and Stick');
            rotator%d.styles.backgroundColor = 'transparent';
            rotator%d.styles.atoms_useJMOLColors = true;
            var mol = ChemDoodle.readPDB('%s');
            rotator%d.loadMolecule(mol);
            rotator%d.startAnimation();
        """%(i,i,i,i,i,pdbStrings[moli],i,i)
        divStrngs += """
                <div class="col-md-3 col-sm-6 col-xs-12">
                   <canvas id="canvas3d%d" class="ChemDoodleWebComponent"></canvas>
                </div> 
                """%(i)
    htmlStr += jsStrings
    htmlStr += """
        </script>
    """
    htmlStr += """
        <div class="container">
            <center>
            <div class="row">
                %s
            </div>
            </center>
        </div>"""%divStrngs
    htmlStr += """
    </body>
    </html>"""
components.html(htmlStr, height=130)

fig = go.Figure()
for molChoice in molChoices:
    fig.add_trace(go.Scatter(x=xsData[molChoice]['E'], y=xsData[molChoice]['xs'], mode='lines', name=molChoice))
fig.update_layout(title="BEB Cross sections" ,
                xaxis_title="Energy (eV)",
                yaxis_title="Cross section (Å²)",
                font=dict(
                    size=18,
                ))
fig.update_xaxes(type="log", exponentformat="power", showexponent = 'all', showgrid=True, tickfont=dict(
                    size=16,
                ))
fig.update_yaxes(type="log", exponentformat="power", showexponent = 'all', showgrid=True, tickfont=dict(
                    size=16,
                ))
st.plotly_chart(fig, use_container_width=True)

if len(molChoices) > 0:
    df = pd.DataFrame()
    df['Energy'] = xsData[molChoices[0]]['E']
    for mi in molChoices:
        df[mi] = xsData[mi]['xs']
    csvStr = df.to_csv(index=False).encode('utf-8')
    headArr = ['E'] + molChoices
    df.columns = headArr
    st.download_button(
        label="Download data as CSV",
        data=csvStr,
        file_name='cross_sections.csv',
        mime='text/csv',
    )   
