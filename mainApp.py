import os
import glob
import numpy as np
import plotly.graph_objects as go
import streamlit as st
import pandas as pd
import streamlit.components.v1 as components

# --- CONFIG & CONSTANTS ---
st.set_page_config(
    layout="centered", page_title="ALeCS Database", page_icon=":milky_way:"
)

A0 = 0.529177211  # Angstroms
A02 = A0**2


# --- CACHED DATA LOADING ---
@st.cache_data
def load_assets():
    """Loads JS and CSS files into memory to inject into the iframe."""
    try:
        with open("assets/uis/ChemDoodleWeb.js", "r") as f:
            js = f.read()
        return js
    except FileNotFoundError:
        return ""


@st.cache_data
def load_molecule_data():
    """Reads physics data and PDB paths once and caches them."""
    # Load LaTeX mappings
    latex_mols = {}
    if os.path.exists("molLatex.txt"):
        with open("molLatex.txt", "r") as lm:
            for line in lm:
                s = line.split()
                if len(s) >= 2:
                    latex_mols[s[0]] = s[1]

    mols = []
    xs_data = {}
    pdb_paths = {}

    # Scan for data files
    for file_path in glob.glob("BEB/*.xs"):
        # Use basename to prevent path traversal
        filename = os.path.basename(file_path)
        moli = filename.split(".")[0]

        pdb_file = f"pdbs/{moli}.pdb"
        if not os.path.exists(pdb_file):
            continue

        # Load cross-section data
        try:
            dat = np.loadtxt(file_path, skiprows=2)
            if dat.size == 0:
                continue

            # Energy limit logic
            e_vals = dat[:, 0]
            xs_vals = dat[:, 1] * A02
            mask = e_vals <= 5000

            mols.append(moli)
            xs_data[moli] = {"E": e_vals[mask], "xs": xs_vals[mask]}
            # Store the actual PDB content for JS injection
            with open(pdb_file, "r") as f:
                pdb_paths[moli] = f.read().replace("\n", "\\n")
        except Exception as e:
            st.error(f"Error loading {moli}: {e}")

    return sorted(mols), latex_mols, xs_data, pdb_paths


# Initialize data
mols_list, latex_map, xs_data, pdb_contents = load_molecule_data()
chemdoodle_js = load_assets()

# --- UI ---
st.markdown("### Select molecules to show the cross sections")
mol_choices = st.multiselect("Molecules", mols_list, max_selections=5)

if mol_choices:
    # Build JS for multiple canvases
    js_canvases = ""
    div_html = ""

    for i, moli in enumerate(mol_choices):
        pdb_data = pdb_contents.get(moli, "")
        div_html += f'<div class="col-2" id="canvas_container_{i}"><canvas id="canvas3d_{i}"></canvas></div>'
        js_canvases += f"""
            var rotator{i} = new ChemDoodle.RotatorCanvas3D('canvas3d_{i}', 125, 125);
            rotator{i}.styles.set3DRepresentation('Ball and Stick');
            rotator{i}.styles.backgroundColor = 'transparent';
            var mol{i} = ChemDoodle.readPDB('{pdb_data}');
            rotator{i}.loadMolecule(mol{i});
            rotator{i}.startAnimation();
        """

    full_html = f"""
        <html>
        <head>
            <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
            <script>{chemdoodle_js}</script>
            <style>canvas {{ width: 125px; height: 125px; }}</style>
            <style>
                .row {{ display: flex; justify-content: center; align-items: center; flex-wrap: wrap; }}
            </style>
        </head>
        <body>
            <div class="container-fluid"><div class="row">{div_html}</div></div>
            <script>{js_canvases}</script>
        </body>
        </html>
    """
    components.html(full_html, height=150)

    # --- PLOTLY CHART ---
    fig = go.Figure()
    for m in mol_choices:
        fig.add_trace(
            go.Scatter(x=xs_data[m]["E"], y=xs_data[m]["xs"], mode="lines", name=m)
        )

    fig.update_layout(
        xaxis_title="Energy (eV)",
        yaxis_title="Cross section (Å²)",
        xaxis_type="log",
        template="plotly_white",
    )
    st.plotly_chart(fig, use_container_width=True)

    # --- DOWNLOAD DATA ---
    df = pd.DataFrame({"Energy": xs_data[mol_choices[0]]["E"]})
    for m in mol_choices:
        df[m] = xs_data[m]["xs"]

    st.download_button(
        label="Download CSV",
        data=df.to_csv(index=False),
        file_name="data.csv",
        mime="text/csv",
    )
