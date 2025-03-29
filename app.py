import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
from Bio.SeqUtils.ProtParam import ProteinAnalysis

st.set_page_config(layout='wide')
st.sidebar.title('PRo-EDICTOR')
st.sidebar.write("Built By Department Of BIOINFORMATICS")

# Default protein sequence
DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
txt = st.sidebar.text_area('Input sequence', DEFAULT_SEQ, height=275)

# Initialize session state variables
if "pdb_string" not in st.session_state:
    st.session_state.pdb_string = None
if "b_value" not in st.session_state:
    st.session_state.b_value = None
if "show_pdb" not in st.session_state:
    st.session_state.show_pdb = False

# Function to update protein structure prediction
def update(sequence=txt):
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    
    # Store PDB structure in session state
    st.session_state.pdb_string = response.content.decode('utf-8')

    with open('predicted.pdb', 'w') as f:
        f.write(st.session_state.pdb_string)

    struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
    st.session_state.b_value = round(struct.b_factor.mean(), 4)  # Mean B-factor (plDDT confidence score)

# Sidebar Prediction Button
if st.sidebar.button('⏳ Predict Structure'):
    update()

# Only display content if prediction has been run
if st.session_state.pdb_string:
    col1, col2 = st.columns([2, 1])  # Adjust width ratio for structure and table
    
    with col1:
        st.markdown("<h1 style='text-align: center;'>PRo-EDICTOR</h1>", unsafe_allow_html=True)
        st.subheader('🧬 Predicted Protein Structure')
        pdbview = py3Dmol.view()
        pdbview.addModel(st.session_state.pdb_string, 'pdb')
        pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
        pdbview.setBackgroundColor('lightblack')
        pdbview.zoomTo()
        pdbview.zoom(2, 800)
        pdbview.spin(True)
        showmol(pdbview, height=500, width=800)
    
    with col2:
        st.subheader('📈 Protein Properties')
        protein_seq = ProteinAnalysis(txt)
        sequence_length = len(txt)
        mol_weight = round(protein_seq.molecular_weight(), 2)
        hydrophobic_residues = sum(txt.upper().count(res) for res in "AILMFWYV")
        acidic = sum(txt.upper().count(res) for res in "DE")
        basic = sum(txt.upper().count(res) for res in "KRH")
        net_charge = basic - acidic
        
        table_data = {
            "Property": ["Sequence Length", "Molecular Weight (Da)", "Hydrophobic Residues", "Net Charge (pH 7)", "plDDT Confidence Score"],
            "Value": [sequence_length, mol_weight, hydrophobic_residues, net_charge, st.session_state.b_value]
        }
        st.table(table_data)

        # Significance of Colors in Predicted Structure
        st.subheader("📊 Color Interpretation")
        color_table = """
        | **Color**  | **Confidence Score (plDDT)** | **Interpretation** |
        |------------|-----------------------------|---------------------|
        | 🟦 **Blue** | **90-100**                  | Very High Confidence |
        | 🟩 **Green** | **70-90**                   | Good Reliability |
        | 🟨 **Yellow** | **50-70**                   | Low Confidence, Flexible Regions |
        | 🟥 **Red** | **0-50**                     | Very Uncertain, Unstructured |
        """
        st.markdown(color_table, unsafe_allow_html=True)
        

        # Fix: PDB File Content Checkbox (No Auto-refresh)
        # Create empty columns to center the content
    col_left, col_center, col_right = st.columns([2, 4, 2])  # Center column is 3x wider

    with col_center:  # Put the content in the center column
        st.subheader("📂 PDB File Viewer")
    
        show_pdb_file = st.checkbox("Show PDB File Content", value=st.session_state.show_pdb)
    
        if show_pdb_file:
            st.session_state.show_pdb = True
            st.text_area("PDB File", st.session_state.pdb_string, height=700)  # Increased height
        else:
            st.session_state.show_pdb = False

    # Download Button
    st.download_button(
        label="📥 Download PDB File",
        data=st.session_state.pdb_string,
        file_name='predicted.pdb',
        mime='text/plain',
    )
else:
   st.markdown(''' 
    <div style="text-align: center;">
        <h1 style="font-family: 'Dustosmo Roman', 'Times New Roman', serif; color:purple;">
            PRo-EDICTOR
        </h1>
        <p>
            <em>A Streamlit <strong>Component</strong> for creating Speck molecular structures within Streamlit Web app.</em>
        </p>
    </div>
    ''', unsafe_allow_html=True)
   st.warning('⬅️ Enter protein sequence data and click "Predict Structure"!')


