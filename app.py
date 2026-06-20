import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
import pubchempy as pcp
import urllib.parse
import re
import math
import json

st.set_page_config(page_title="Triola", layout="wide")

# ========== ФУНКЦИИ ==========
def safe_smiles_to_mol(smiles):
    original = smiles.strip()
    smiles = original.replace("О", "O").replace("о", "O").replace("С", "C").replace("с", "C")
    smiles = ' '.join(smiles.split())
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            try:
                Chem.SanitizeMol(mol)
                return mol, None
            except:
                pass
    except:
        pass
    
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol:
            try:
                Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_AROMATICITY | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
                return mol, "⚠️ Structure loaded"
            except:
                try:
                    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                    return mol, "⚠️ Structure loaded"
                except:
                    return mol, "⚠️ Structure loaded (basic)"
    except:
        pass
    
    cleaned = re.sub(r'([a-z])(\d+)', r'\1', smiles)
    if cleaned != smiles:
        try:
            mol = Chem.MolFromSmiles(cleaned)
            if mol:
                try:
                    Chem.SanitizeMol(mol)
                    return mol, "⚠️ SMILES fixed"
                except:
                    return mol, "⚠️ SMILES fixed (basic)"
        except:
            pass
    
    return None, f"❌ Invalid SMILES: {original[:60]}"

def get_pubchem_data(smiles):
    try:
        compounds = pcp.get_compounds(smiles, 'smiles')
        if compounds:
            comp = compounds[0]
            return {
                'cid': comp.cid,
                'name': comp.iupac_name,
                'formula': comp.molecular_formula,
            }
    except:
        pass
    return None

def get_cas_from_pubchem(smiles):
    try:
        compounds = pcp.get_compounds(smiles, 'smiles')
        if compounds:
            comp = compounds[0]
            if comp.synonyms:
                for syn in comp.synonyms:
                    if re.match(r'^\d{2,7}-\d{2}-\d$', syn):
                        return syn
            return None
    except:
        pass
    return None

def build_search_urls(name, keywords, cid, smiles):
    urls = {}
    
    if name:
        search_terms = name
        if keywords:
            search_terms = f"({name}) AND ({keywords})"
    elif keywords:
        search_terms = keywords
    else:
        search_terms = None
    
    if search_terms:
        urls['pubmed'] = f"https://pubmed.ncbi.nlm.nih.gov/?term={urllib.parse.quote(search_terms)}"
        urls['scholar'] = f"https://scholar.google.com/scholar?q={urllib.parse.quote(search_terms)}"
        urls['europe_pmc'] = f"https://europepmc.org/search?query={urllib.parse.quote(search_terms)}"
    
    if name:
        urls['patents'] = f"https://patents.google.com/?q={urllib.parse.quote(name)}"
    else:
        urls['patents'] = None
    
    if cid:
        urls['pubchem'] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    else:
        urls['pubchem'] = f"https://pubchem.ncbi.nlm.nih.gov/#query={urllib.parse.quote(smiles[:80])}"
    
    return urls

def find_and_show_analogs(cid, smiles, name, keywords):
    st.divider()
    st.subheader("🔬 Structural Analogs")
    st.caption("Find chemical compounds with similar structure")

    pubchem_similar_url = None
    if smiles:
        encoded_smiles = urllib.parse.quote(smiles)
        pubchem_similar_url = f"https://pubchem.ncbi.nlm.nih.gov/#query={encoded_smiles}&sort=similarity"
    elif cid:
        pubchem_similar_url = f"https://pubchem.ncbi.nlm.nih.gov/#query=CID{cid}&sort=similarity"

    scholar_analog_url = None
    if name:
        scholar_analog_url = f"https://scholar.google.com/scholar?q={urllib.parse.quote(name)}+analog"
    elif keywords:
        scholar_analog_url = f"https://scholar.google.com/scholar?q={urllib.parse.quote(keywords)}+analog"

    st.markdown("**Search for structural analogs manually:**")
    
    col1, col2 = st.columns(2)
    
    with col1:
        if pubchem_similar_url:
            st.link_button("🧬 PubChem Similar", pubchem_similar_url, use_container_width=True)
            st.caption("Similar chemical structures")
        else:
            st.link_button("🧬 PubChem Similar", "", disabled=True, help="SMILES required", use_container_width=True)
    
    with col2:
        if scholar_analog_url:
            st.link_button("🎓 Google Scholar", scholar_analog_url, use_container_width=True)
            st.caption("'name' + analog")
        else:
            st.link_button("🎓 Google Scholar", "", disabled=True, help="Compound name required", use_container_width=True)

# ========== МАГАЗИНЫ ==========
STORES = {
    "Российские": {
        "🧪 ХимМед": "https://chimmed.ru/products/search?name={cas}"
    },
    "Международные": {
        "🏭 Sigma-Aldrich": "https://www.sigmaaldrich.com/RU/en/search/{cas}?focus=products&page=1&perpage=30&sort=relevance&term={cas}&type=cas_number",
        "🏭 Macklin": "https://www.macklin.cn/search/{cas}",
        "🏭 Combi-Blocks": "https://www.combi-blocks.com/cgi-bin/find.cgi?cas={cas}",
        "🏭 Aladdin": "https://www.aladdin-e.com/zh_cn/catalogsearch/result/?q={cas}",
        "🏭 Bide Pharmatech": "https://www.bidepharm.com/products/{cas}.html"
    }
}

# ========== ЭЛЕМЕНТЫ ДЛЯ LEPTORIUM ==========
ELEMENTS = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999,
    'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.065,
    'Cl': 35.45, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996,
    'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.630,
    'As': 74.922, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
    'Nb': 92.906, 'Mo': 95.95, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41,
    'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29, 'Cs': 132.91, 'Ba': 137.33,
    'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24, 'Pm': 145.0, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25,
    'Tb': 158.93, 'Dy': 162.50, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05, 'Lu': 174.97, 'Hf': 178.49,
    'Ta': 180.95, 'W': 183.84, 'Re': 186.21, 'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59,
    'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Po': 209.0, 'At': 210.0, 'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0,
    'Ac': 227.0, 'Th': 232.04, 'Pa': 231.04, 'U': 238.03, 'Np': 237.0, 'Pu': 244.0, 'Am': 243.0, 'Cm': 247.0,
    'Bk': 247.0, 'Cf': 251.0, 'Es': 252.0, 'Fm': 257.0, 'Md': 258.0, 'No': 259.0, 'Lr': 266.0, 'Rf': 267.0,
    'Db': 268.0, 'Sg': 269.0, 'Bh': 270.0, 'Hs': 277.0, 'Mt': 278.0, 'Ds': 281.0, 'Rg': 282.0, 'Cn': 285.0,
    'Nh': 286.0, 'Fl': 289.0, 'Mc': 290.0, 'Lv': 293.0, 'Ts': 294.0, 'Og': 294.0
}

def calculate_molar_mass(formula):
    """Рекурсивный парсинг формулы с поддержкой двухбуквенных элементов и скобок"""
    if not formula:
        return 0.0
    
    total_mass = 0.0
    i = 0
    n = len(formula)
    
    while i < n:
        if formula[i] == '(':
            j = i + 1
            bracket_level = 1
            while j < n and bracket_level > 0:
                if formula[j] == '(':
                    bracket_level += 1
                elif formula[j] == ')':
                    bracket_level -= 1
                j += 1
            
            inside = formula[i+1:j-1]
            
            multiplier = 1
            if j < n and formula[j].isdigit():
                k = j
                while k < n and formula[k].isdigit():
                    k += 1
                multiplier = int(formula[j:k])
                j = k
            
            total_mass += calculate_molar_mass(inside) * multiplier
            i = j
            
        else:
            if formula[i].isupper():
                if i + 1 < n and formula[i+1].islower():
                    element = formula[i:i+2]
                    i += 2
                else:
                    element = formula[i]
                    i += 1
                
                count = 1
                if i < n and formula[i].isdigit():
                    k = i
                    while k < n and formula[k].isdigit():
                        k += 1
                    count = int(formula[i:k])
                    i = k
                
                if element in ELEMENTS:
                    total_mass += ELEMENTS[element] * count
            else:
                i += 1
    
    return total_mass

def show_additional_params(mass, moles, volume, solution_mass, mass_substance=None, volume_substance=None, prefix=""):
    """Показывает дополнительные параметры (концентрация, плотность, массовая доля)"""
    params = []
    
    # Концентрация (из n и V раствора)
    if moles > 0 and volume > 0:
        conc = moles / volume
        params.append(f"Концентрация: C = n / V = {moles:.4f} / {volume:.2f} = {conc:.3f} моль/л")
    
    # Молярная масса (из m и n)
    if mass > 0 and moles > 0:
        molar_mass_calc = mass / moles
        params.append(f"Молярная масса (расчётная): M = m / n = {mass:.2f} / {moles:.4f} = {molar_mass_calc:.3f} г/моль")
    
    # Плотность РАСТВОРА (из массы раствора и объёма раствора)
    if solution_mass > 0 and volume > 0:
        density_solution = solution_mass / volume
        params.append(f"Плотность раствора: ρ = m(р-ра) / V = {solution_mass:.2f} / {volume:.2f} = {density_solution:.3f} г/л ({density_solution/1000:.3f} г/мл)")
    
    # Плотность ВЕЩЕСТВА (из массы вещества и объёма вещества)
    if mass_substance is not None and volume_substance is not None:
        if mass_substance > 0 and volume_substance > 0:
            density_substance = mass_substance / volume_substance
            params.append(f"Плотность вещества: ρ = m(в-ва) / V(в-ва) = {mass_substance:.2f} / {volume_substance:.2f} = {density_substance:.3f} г/мл")
    
    # Массовая доля (из массы вещества и массы раствора)
    if mass > 0 and solution_mass > 0:
        mass_fraction = mass / solution_mass
        params.append(f"Массовая доля: ω = m_в / m_р = {mass:.2f} / {solution_mass:.2f} = {mass_fraction:.3f} ({mass_fraction*100:.1f}%)")
    
    if params:
        st.markdown(f"**📊 Дополнительные параметры {prefix}:**")
        for p in params:
            st.markdown(f'<div class="calc-formula">• {p}</div>', unsafe_allow_html=True)

# ========== ОСНОВНОЙ ИНТЕРФЕЙС ==========
query_params = st.query_params
page = query_params.get("page", "main")

# ========== СТРАНИЦА 1: ГЛАВНАЯ ==========
if page == "main":
    st.markdown("""
    <style>
        .stApp { background: linear-gradient(135deg, #006064 0%, #00BCD4 100%); }
        .main-title { color: white; text-align: center; font-size: 4em; font-weight: bold; margin-top: 60px; text-shadow: 2px 2px 10px rgba(0,0,0,0.3); }
        .main-subtitle { text-align: center; color: #E0F7FA; font-size: 1.2rem; margin-top: -15px; margin-bottom: 40px; }
        .footer { text-align: center; color: #B2EBF2; margin-top: 80px; font-size: 0.9rem; }
        div[data-testid="column"] { display: flex; flex-direction: column; align-items: center; justify-content: flex-start; }
        .stButton { display: flex; justify-content: center; width: 100%; }
        .stButton button {
            background: rgba(255,255,255,0.15) !important; backdrop-filter: blur(10px) !important;
            color: white !important; border: 2px solid rgba(255,255,255,0.3) !important;
            border-radius: 20px !important; padding: 25px 20px !important;
            width: 220px !important; height: 120px !important;
            transition: all 0.3s !important; text-align: center !important;
            white-space: pre-line !important; font-size: 1.1em !important;
            line-height: 1.3 !important; margin: 0 auto !important;
            display: flex !important; flex-direction: column !important;
            justify-content: center !important; align-items: center !important;
        }
        .stButton button:hover { transform: translateY(-10px) !important; border-color: white !important; background: rgba(255,255,255,0.3) !important; box-shadow: 0 10px 30px rgba(0,0,0,0.2) !important; }
        .stButton button p { font-size: 1.6em !important; font-weight: bold !important; color: white !important; margin: 0 !important; }
        .module-status { color: #B2EBF2; font-size: 1.0rem; text-align: center; margin-top: 8px; font-style: italic; }
        .row-widget.stColumns { justify-content: center !important; }
    </style>
    """, unsafe_allow_html=True)
    
    st.markdown("<div class='main-title'>Triola — all begin with a dream</div>", unsafe_allow_html=True)
    st.markdown("<div class='main-subtitle'>Три инструмента для химика</div>", unsafe_allow_html=True)
    
    col1, col2, col3 = st.columns([1, 1, 1], gap="large")
    
    with col1:
        if st.button("📡🧪\nChemRadar", key="btn_chemradar", use_container_width=True):
            st.query_params["page"] = "chemradar"
            st.rerun()
        st.markdown("<div class='module-status'>Investigating!</div>", unsafe_allow_html=True)
    
    with col2:
        if st.button("💎\nLaretc", key="btn_laretc", use_container_width=True):
            st.query_params["page"] = "laretc"
            st.rerun()
        st.markdown("<div class='module-status'>Planning!</div>", unsafe_allow_html=True)
    
    with col3:
        if st.button("⚗️📊\nLeptorium", key="btn_leptorium", use_container_width=True):
            st.query_params["page"] = "leptorium"
            st.rerun()
        st.markdown("<div class='module-status'>Creating!</div>", unsafe_allow_html=True)
    
    st.markdown("<div class='footer'>⚡ Triola — лабораторная платформа для химика</div>", unsafe_allow_html=True)

# ========== СТРАНИЦА 2: CHEMRADAR ==========
elif page == "chemradar":
    st.markdown("""
    <style>
        .stApp { background: #C8E6C9; }
        .stButton button { background: #43A047 !important; color: white !important; border-radius: 10px !important; padding: 10px 20px !important; border: none !important; transition: all 0.3s !important; }
        .stButton button:hover { background: #2E7D32 !important; transform: scale(1.02) !important; }
        .stLinkButton button { background: #66BB6A !important; color: white !important; border-radius: 10px !important; padding: 10px 20px !important; border: none !important; transition: all 0.3s !important; }
        .stLinkButton button:hover { background: #388E3C !important; transform: scale(1.02) !important; }
        .centered-image { display: flex; justify-content: center; align-items: center; width: 100%; }
        .centered-image img { margin: 0 auto; display: block; }
    </style>
    """, unsafe_allow_html=True)
    
    if st.button("← На главную"):
        st.query_params["page"] = "main"
        st.rerun()
    
    APP_VERSION = "v.1.0"
    
    st.markdown(
        f"<div style='text-align: center;'>"
        f"<h1>📡🧪 ChemRadar</h1>"
        f"<p style='font-size: 1.1rem; color: #2E7D32;'>Chemical Literature Scanner — search scientific papers by SMILES and keywords</p>"
        f"<p style='font-size: 0.8rem; color: #558B2F;'>version {APP_VERSION}</p>"
        f"</div>",
        unsafe_allow_html=True
    )
    
    st.markdown("---")
    
    with st.sidebar:
        st.markdown("""
    ## 📡🧪 ChemRadar
    *Chemical Literature Scanner*
    ---
    ## 📖 KEYWORD SEARCH GUIDE
    ### Simple queries (work everywhere)
    | How to enter | Example | What it finds |
    |--------------|---------|----------------|
    | One word | rimonabant | Articles with this word |
    | Multiple words | rimonabant tuberculosis | Articles with ALL words |
    | Word list | rimonabant SR141716 | Articles with ANY word |
    ---
    ### Advanced queries (PubMed, Europe PMC, Google Scholar)
    | Operator | Example | What it finds |
    |----------|---------|----------------|
    | AND | rimonabant AND tuberculosis | Both words |
    | OR | rimonabant OR SR141716 | At least one |
    | NOT | rimonabant NOT review | Excludes reviews |
    | " " (phrase) | "rimonabant derivative" | Exact phrase |
    | ( ) (grouping) | (rimonabant OR SR141716) AND tb | Complex logic |
    | : (range) | rimonabant AND 2015:2025 | By year |
    ---
    ### Query examples
    Rimonabant + tuberculosis: rimonabant AND tuberculosis
    Rimonabant or SR141716: rimonabant OR SR141716
    Rimonabant + tuberculosis + resistance: (rimonabant OR SR141716) AND (tuberculosis OR TB) AND resistance
    Exact phrase for derivatives: "rimonabant derivative" AND antitubercular
    Exclude review articles: rimonabant NOT review NOT editorial
    ---
    ### Important notes
    - Keywords in **English only**
    - Use UPPERCASE for operators (AND, OR, NOT)
    - Space between words means AND
    - Quotes required for exact phrases
    """)
    
    left_col, right_col = st.columns([1, 1], gap="medium")
    
    with left_col:
        st.markdown("### 📝 Enter SMILES formula")
        user_smiles = st.text_area(
            "SMILES:", 
            height=100,
            placeholder="Example: CC1=C(C(=O)NN2CCCCC2)N(N=C1C3=CC=C(C=C3)Cl)C4=C(C=C(C=C4)Cl)Cl",
            label_visibility="collapsed",
            key="chemradar_smiles_input"
        )
        
        st.markdown("### 🔑 Enter keywords (English only)")
        st.caption("💡 Examples: rimonabant tuberculosis | rimonabant AND tuberculosis | rimonabant OR SR141716")
        
        keywords = st.text_input(
            "Keywords:", 
            placeholder="e.g., rimonabant tuberculosis",
            label_visibility="collapsed"
        )
        
        st.caption("💡 Full guide in the left sidebar")
        
        submit_button = st.button("🔍 Search", type="primary", use_container_width=True)
    
    if "processed" not in st.session_state:
        st.session_state.processed = False
    
    should_process = user_smiles and (submit_button or st.session_state.get("processed", False))
    if user_smiles and not st.session_state.get("processed", False) and submit_button:
        st.session_state.processed = True
        should_process = True
    elif not user_smiles:
        st.session_state.processed = False
    
    with right_col:
        if should_process or (user_smiles and st.session_state.get("processed", False)):
            mol, warning = safe_smiles_to_mol(user_smiles)
            
            if mol:
                try:
                    inchikey = Chem.MolToInchiKey(mol)
                except:
                    inchikey = "N/A"
                
                pubchem_data = get_pubchem_data(user_smiles)
                cid = pubchem_data["cid"] if pubchem_data else None
                name = pubchem_data["name"] if pubchem_data else None
                
                st.markdown('<div class="centered-image">', unsafe_allow_html=True)
                st.image(Draw.MolToImage(mol, size=(350, 350)), caption="Molecule Structure")
                st.markdown('</div>', unsafe_allow_html=True)
                
                st.code(f"InChIKey: {inchikey}")
                
                if cid:
                    st.success(f"PubChem CID: {cid}")
                if name:
                    st.info(f"📌 Name: {name[:80]}..." if len(name) > 80 else f"📌 Name: {name}")
                if warning:
                    st.warning(warning)
                if keywords:
                    st.caption(f"🔍 Query: {keywords[:100]}" if len(keywords) > 100 else f"🔍 Query: {keywords}")
                
                st.divider()
                
                urls = build_search_urls(name, keywords, cid, user_smiles)
                
                st.subheader("🔍 Search Literature")
                
                st.markdown("**📚 Academic Databases:**")
                if "pubmed" in urls:
                    st.link_button("📖 PubMed", urls["pubmed"], use_container_width=True)
                if "scholar" in urls:
                    st.link_button("🎓 Google Scholar", urls["scholar"], use_container_width=True)
                if "europe_pmc" in urls:
                    st.link_button("📚 Europe PMC", urls["europe_pmc"], use_container_width=True)
                
                st.markdown("---")
                st.markdown("**📄 Patents (experimental):**")
                if urls.get("patents"):
                    st.link_button("📑 Google Patents", urls["patents"], use_container_width=True)
                    st.caption("Search by IUPAC name (may not always find patents)")
                else:
                    st.link_button("📑 Google Patents", "", disabled=True, help="IUPAC name not available", use_container_width=True)
                
                st.markdown("---")
                st.markdown("**🧬 Compound Info:**")
                st.link_button("🧬 PubChem", urls["pubchem"], use_container_width=True)
                
                find_and_show_analogs(cid, user_smiles, name, keywords)
                
                st.session_state.processed = True
                
            else:
                st.error(warning if warning else "❌ Invalid SMILES format")
                st.info("💡 Try copying SMILES from [PubChem](https://pubchem.ncbi.nlm.nih.gov)")
                st.session_state.processed = False
    
    st.divider()
    st.markdown(
        f"<div style='text-align: center; color: #558B2F; font-size: 0.8rem;'>"
        f"📡🧪 ChemRadar — Chemical Literature Scanner | version {APP_VERSION} | Search by SMILES + keywords"
        f"</div>",
        unsafe_allow_html=True
    )

# ========== СТРАНИЦА 3: LARETC ==========
elif page == "laretc":
    st.markdown("""
    <style>
        .stApp { background: #CFD8DC; }
        .stButton button { 
            background: #546E7A !important; 
            color: white !important; 
            border-radius: 8px !important; 
            padding: 6px 14px !important; 
            border: none !important; 
            transition: all 0.3s !important;
            font-size: 13px !important;
        }
        .stButton button:hover { background: #37474F !important; transform: scale(1.02) !important; }
        .stLinkButton button { background: #78909C !important; color: white !important; border-radius: 10px !important; padding: 10px 20px !important; border: none !important; transition: all 0.3s !important; }
        .stLinkButton button:hover { background: #455A64 !important; transform: scale(1.02) !important; }
        .centered-image { display: flex; justify-content: center; align-items: center; width: 100%; }
        .centered-image img { margin: 0 auto; display: block; }
        .info-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 10px; background: rgba(255,255,255,0.5); border-radius: 10px; padding: 12px; margin: 10px 0; }
        .info-item { padding: 6px 10px; background: white; border-radius: 6px; }
        .info-item strong { color: #37474F; font-size: 13px; }
        .info-item span { color: #546E7A; font-size: 14px; font-weight: 500; }
        .compound-card {
            background: rgba(255,255,255,0.7);
            border-radius: 10px;
            padding: 15px;
            margin: 10px 0;
            border-left: 4px solid #546E7A;
        }
    </style>
    """, unsafe_allow_html=True)
    
    if st.button("← На главную"):
        st.query_params["page"] = "main"
        st.rerun()
    
    st.title("💎 Laretc")
    
    st.markdown("""
    ### 📝 Введите SMILES формулы
    Введите структурные формулы в формате SMILES для поиска информации о соединениях и доступных реактивах.
    """)
    
    num_compounds = st.selectbox(
        "Сколько соединений вы хотите проверить?",
        [1, 2, 3],
        index=0,
        key="num_compounds_laretc"
    )
    
    st.divider()
    
    all_smiles = []
    
    for i in range(num_compounds):
        st.markdown(f"### Соединение {i+1}")
        
        smiles_input = st.text_input(
            f"SMILES {i+1}:",
            placeholder=f"Например: CC(=O)O (уксусная кислота)",
            key=f"laretc_smiles_{i}"
        )
        
        if smiles_input:
            mol, warning = safe_smiles_to_mol(smiles_input)
            if mol:
                all_smiles.append(smiles_input)
                st.success(f"✅ Соединение {i+1} распознано")
            else:
                st.error(f"❌ Ошибка в соединении {i+1}: {warning}")
                st.info("💡 Попробуйте скопировать SMILES с [PubChem](https://pubchem.ncbi.nlm.nih.gov)")
        
        if i < num_compounds - 1:
            st.divider()
    
    if st.button("🔍 Найти информацию", key="search_laretc", use_container_width=True):
        if all_smiles:
            st.session_state.laretc_all_smiles = all_smiles
            st.rerun()
        else:
            st.warning("⚠️ Введите хотя бы одно корректное соединение")
    
    if st.session_state.get("laretc_all_smiles"):
        all_smiles = st.session_state.laretc_all_smiles
        
        for idx, smiles in enumerate(all_smiles):
            st.divider()
            st.markdown(f"### 🔍 Информация о соединении {idx+1}")
            
            pubchem_data = get_pubchem_data(smiles)
            cas_number = get_cas_from_pubchem(smiles)
            mol, _ = safe_smiles_to_mol(smiles)
            
            st.markdown('<div class="compound-card">', unsafe_allow_html=True)
            
            if mol:
                st.markdown('<div class="centered-image">', unsafe_allow_html=True)
                st.image(Draw.MolToImage(mol, size=(250, 250)), caption=f"Структура {idx+1}")
                st.markdown('</div>', unsafe_allow_html=True)
            
            st.markdown('<div class="info-grid">', unsafe_allow_html=True)
            
            col_info1, col_info2 = st.columns(2)
            with col_info1:
                if pubchem_data and pubchem_data.get('name'):
                    st.markdown(f'<div class="info-item"><strong>📌 Название</strong><br><span>{pubchem_data["name"]}</span></div>', unsafe_allow_html=True)
                if pubchem_data and pubchem_data.get('formula'):
                    st.markdown(f'<div class="info-item"><strong>📝 Формула</strong><br><span>{pubchem_data["formula"]}</span></div>', unsafe_allow_html=True)
            
            with col_info2:
                if pubchem_data and pubchem_data.get('cid'):
                    st.markdown(f'<div class="info-item"><strong>🔢 CID</strong><br><span>{pubchem_data["cid"]}</span></div>', unsafe_allow_html=True)
                if cas_number:
                    st.markdown(f'<div class="info-item"><strong>🔢 CAS</strong><br><span>{cas_number}</span></div>', unsafe_allow_html=True)
            
            st.markdown('</div>', unsafe_allow_html=True)
            
            if cas_number:
                st.markdown("### 🛒 Поиск реактивов")
                st.caption(f"Поиск по CAS: {cas_number}")
                
                st.markdown("**🇷🇺 Российские магазины:**")
                cols = st.columns(2)
                for idx_store, (name, url_template) in enumerate(STORES["Российские"].items()):
                    with cols[idx_store % 2]:
                        search_url = url_template.format(cas=cas_number)
                        st.link_button(name, search_url, use_container_width=True)
                
                st.markdown("**🌍 Международные магазины:**")
                cols = st.columns(2)
                for idx_store, (name, url_template) in enumerate(STORES["Международные"].items()):
                    with cols[idx_store % 2]:
                        search_url = url_template.format(cas=cas_number)
                        st.link_button(name, search_url, use_container_width=True)
            else:
                st.info("ℹ️ CAS номер не найден. Поиск по магазинам недоступен.")
            
            st.markdown('</div>', unsafe_allow_html=True)

# ========== СТРАНИЦА 4: LEPTORIUM ==========
elif page == "leptorium":
    st.markdown("""
    <style>
        .stApp { background: #FFE0B2; }
        .stButton button { background: #F57C00 !important; color: white !important; border-radius: 10px !important; padding: 10px 20px !important; border: none !important; transition: all 0.3s !important; }
        .stButton button:hover { background: #E65100 !important; transform: scale(1.02) !important; }
        .stLinkButton button { background: #FFA726 !important; color: white !important; border-radius: 10px !important; padding: 10px 20px !important; border: none !important; transition: all 0.3s !important; }
        .stLinkButton button:hover { background: #EF6C00 !important; transform: scale(1.02) !important; }
        .compound-card { background: rgba(255,255,255,0.85); border-radius: 8px; padding: 12px; margin: 8px 0; border-left: 4px solid #F57C00; }
        .result-panel { background: rgba(255,255,255,0.95); border-radius: 8px; padding: 15px; margin: 12px 0; border-left: 4px solid #E65100; }
        .step { background: rgba(255,255,255,0.7); border-radius: 6px; padding: 10px; margin: 6px 0; border-left: 3px solid #FFA726; }
        .stDivider { margin: 4px 0 !important; }
        .field-label { font-size: 11px; color: #78909C; font-weight: 500; text-align: center; margin-bottom: 2px; }
        .calc-line {
            font-family: 'Courier New', monospace;
            font-size: 14px;
            padding: 2px 0;
            color: #1A237E;
        }
        .calc-formula {
            font-family: 'Courier New', monospace;
            font-size: 14px;
            padding: 4px 0 4px 20px;
            color: #0D47A1;
            background: rgba(255,255,255,0.5);
            border-radius: 4px;
            margin: 2px 0;
        }
        .limiting-highlight {
            background: #FFF3E0;
            border: 2px solid #E65100;
            border-radius: 8px;
            padding: 10px;
            margin: 10px 0;
        }
        .result-title {
            font-size: 18px;
            font-weight: bold;
            color: #BF360C;
            margin: 15px 0 10px 0;
            border-bottom: 2px solid #FFE0B2;
            padding-bottom: 5px;
        }
        .subsection-title {
            font-size: 16px;
            font-weight: bold;
            color: #E65100;
            margin: 10px 0 5px 0;
        }
        .structure-container {
            display: flex;
            justify-content: center;
            margin: 8px 0;
        }
        .structure-container img {
            max-width: 200px;
            height: auto;
        }
    </style>
    """, unsafe_allow_html=True)

    if st.button("← На главную"):
        st.query_params["page"] = "main"
        st.rerun()

    st.title("⚗️📊 Leptorium")

    # ===== ВЫБОР РЕЖИМА РАСЧЕТА =====
    st.markdown("### ⚖️ Выберите режим расчета")
    calc_mode = st.radio(
        "Режим расчета:",
        ["🔬 Расчет по одному реагенту", "⚖️ Расчет по избытку и недостатку", "📦 Расчет по продукту"],
        horizontal=True,
        key="calc_mode"
    )
    st.divider()

    # ===== ВЫБОР КОЛИЧЕСТВА РЕАГЕНТОВ =====
    st.markdown("### 🧪 Реагенты")
    num_reactants = st.selectbox(
        "Количество реагентов:",
        [1, 2, 3, 4, 5],
        index=0,
        key="num_reactants_leptorium"
    )
    st.divider()

    # ===== ПОЛЯ ДЛЯ ВВОДА РЕАГЕНТОВ =====
    reactants_data = []

    for i in range(num_reactants):
        with st.expander(f"Реагент {i+1}", expanded=True):
            col_type, col_input = st.columns([1, 3])
            with col_type:
                r_type = st.selectbox(
                    "Тип:",
                    ["Органический", "Неорганический"],
                    key=f"r_type_{i}"
                )
            with col_input:
                if r_type == "Органический":
                    r_smiles = st.text_input(
                        "SMILES:",
                        placeholder="Например: CC(=O)O (уксусная кислота)",
                        key=f"r_smiles_{i}",
                        label_visibility="collapsed"
                    )
                    r_formula = ""
                    r_name = ""
                    r_molar_mass = 0.0
                    r_mol = None

                    if r_smiles:
                        r_mol, _ = safe_smiles_to_mol(r_smiles)
                        if r_mol:
                            pubchem = get_pubchem_data(r_smiles)
                            if pubchem:
                                r_name = pubchem.get('name', r_smiles)
                                r_formula = pubchem.get('formula', '')
                                r_molar_mass = calculate_molar_mass(r_formula)
                                st.caption(f"📌 {r_name} | ⚖️ {r_molar_mass:.3f} г/моль | 📝 {r_formula}")
                            else:
                                r_name = r_smiles
                                r_formula = ""
                                r_molar_mass = 0.0
                        else:
                            st.warning("⚠️ Некорректный SMILES")

                    if r_mol:
                        st.markdown('<div class="structure-container">', unsafe_allow_html=True)
                        st.image(Draw.MolToImage(r_mol, size=(200, 200)), caption="Структура")
                        st.markdown('</div>', unsafe_allow_html=True)

                    r_input = r_smiles

                else:
                    r_formula = st.text_input(
                        "Формула:",
                        placeholder="Например: Cu(OH)2",
                        key=f"r_formula_{i}",
                        label_visibility="collapsed"
                    )
                    r_smiles = ""
                    r_name = r_formula
                    r_molar_mass = 0.0
                    r_mol = None

                    if r_formula:
                        r_molar_mass = calculate_molar_mass(r_formula)
                        if r_molar_mass > 0:
                            st.caption(f"⚖️ {r_molar_mass:.3f} г/моль")
                        else:
                            st.warning("⚠️ Некорректная формула")

                    r_input = r_formula

            is_solvate = st.checkbox(
                "Это сольват (есть сольватная часть)?",
                key=f"r_is_solvate_{i}",
                value=False
            )

            solvate_formula = ""
            solvate_count = 1
            solvate_molar_mass = 0.0

            if is_solvate:
                col_solv1, col_solv2 = st.columns(2)
                with col_solv1:
                    solvate_formula = st.text_input(
                        "Формула сольвата:",
                        placeholder="Например: H2O, H2CO3, HCl",
                        key=f"r_solvate_formula_{i}"
                    )
                    if solvate_formula:
                        solvate_molar_mass = calculate_molar_mass(solvate_formula)
                        if solvate_molar_mass > 0:
                            st.caption(f"⚖️ M(сольвата) = {solvate_molar_mass:.3f} г/моль")
                        else:
                            st.warning("⚠️ Некорректная формула сольвата")
                with col_solv2:
                    solvate_count = st.number_input(
                        "Количество молекул:",
                        min_value=1,
                        max_value=100,
                        value=1,
                        step=1,
                        key=f"r_solvate_count_{i}"
                    )

            if r_molar_mass > 0:
                if is_solvate and solvate_molar_mass > 0:
                    total_molar_mass = r_molar_mass + solvate_count * solvate_molar_mass
                    st.info(f"📊 M(общая) = {r_molar_mass:.3f} + {solvate_count} × {solvate_molar_mass:.3f} = **{total_molar_mass:.3f} г/моль**")
                else:
                    total_molar_mass = r_molar_mass
            else:
                total_molar_mass = 0.0

            st.markdown("**📐 Параметры реагента:**")

            col1, col2, col3, col4, col5, col6 = st.columns(6)
            with col1:
                st.markdown('<p class="field-label">Масса (г)</p>', unsafe_allow_html=True)
                mass = st.number_input(
                    "Масса",
                    min_value=0.0,
                    step=0.1,
                    key=f"r_mass_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )
            with col2:
                st.markdown('<p class="field-label">Плотность (г/мл)</p>', unsafe_allow_html=True)
                density = st.number_input(
                    "Плотность",
                    min_value=0.0,
                    step=0.01,
                    key=f"r_dens_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )
            with col3:
                st.markdown('<p class="field-label">Объём (л)</p>', unsafe_allow_html=True)
                volume = st.number_input(
                    "Объём",
                    min_value=0.0,
                    step=0.01,
                    key=f"r_vol_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )
            with col4:
                st.markdown('<p class="field-label">Масс. доля</p>', unsafe_allow_html=True)
                mass_fraction = st.number_input(
                    "Массовая доля",
                    min_value=0.0,
                    max_value=1.0,
                    step=0.01,
                    key=f"r_mf_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )
            with col5:
                st.markdown('<p class="field-label">Конц. (моль/л)</p>', unsafe_allow_html=True)
                concentration = st.number_input(
                    "Концентрация",
                    min_value=0.0,
                    step=0.1,
                    key=f"r_conc_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )
            with col6:
                st.markdown('<p class="field-label">Эквивалент</p>', unsafe_allow_html=True)
                equiv = st.number_input(
                    "Эквивалент реагента",
                    min_value=0.0,
                    step=0.1,
                    key=f"r_equiv_{i}",
                    value=1.0,
                    label_visibility="collapsed"
                )

            col7, col8 = st.columns([1, 1])
            with col7:
                st.markdown('<p class="field-label">Масса раствора (г)</p>', unsafe_allow_html=True)
                solution_mass = st.number_input(
                    "Масса раствора",
                    min_value=0.0,
                    step=0.1,
                    key=f"r_sol_mass_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )

            reactants_data.append({
                'type': r_type,
                'smiles': r_smiles,
                'formula': r_formula,
                'name': r_name if r_name else f"Реагент {i+1}",
                'molar_mass': total_molar_mass,
                'molar_mass_main': r_molar_mass,
                'molar_mass_solvate': solvate_molar_mass,
                'solvate_formula': solvate_formula,
                'solvate_count': solvate_count,
                'is_solvate': is_solvate,
                'mass': mass,
                'density': density,
                'volume': volume,
                'mass_fraction': mass_fraction,
                'concentration': concentration,
                'equiv': equiv,
                'solution_mass': solution_mass,
                'mol': r_mol,
                'actual_moles': 0.0,
            })

    st.divider()

    # ===== ВЫБОР КОЛИЧЕСТВА ПРОДУКТОВ =====
    st.markdown("### 🧪 Продукты")
    num_products = st.selectbox(
        "Количество продуктов:",
        [1, 2, 3, 4, 5],
        index=0,
        key="num_products_leptorium"
    )
    st.divider()

    # ===== ПОЛЯ ДЛЯ ВВОДА ПРОДУКТОВ =====
    products_data = []

    for i in range(num_products):
        with st.expander(f"Продукт {i+1}", expanded=True):
            col_type, col_input = st.columns([1, 3])
            with col_type:
                p_type = st.selectbox(
                    "Тип:",
                    ["Органический", "Неорганический"],
                    key=f"p_type_{i}"
                )
            with col_input:
                if p_type == "Органический":
                    p_smiles = st.text_input(
                        "SMILES:",
                        placeholder="Например: CCO (этанол)",
                        key=f"p_smiles_{i}",
                        label_visibility="collapsed"
                    )
                    p_formula = ""
                    p_name = ""
                    p_molar_mass = 0.0
                    p_mol = None

                    if p_smiles:
                        p_mol, _ = safe_smiles_to_mol(p_smiles)
                        if p_mol:
                            pubchem = get_pubchem_data(p_smiles)
                            if pubchem:
                                p_name = pubchem.get('name', p_smiles)
                                p_formula = pubchem.get('formula', '')
                                p_molar_mass = calculate_molar_mass(p_formula)
                                st.caption(f"📌 {p_name} | ⚖️ {p_molar_mass:.3f} г/моль | 📝 {p_formula}")
                            else:
                                p_name = p_smiles
                                p_formula = ""
                                p_molar_mass = 0.0
                        else:
                            st.warning("⚠️ Некорректный SMILES")

                    if p_mol:
                        st.markdown('<div class="structure-container">', unsafe_allow_html=True)
                        st.image(Draw.MolToImage(p_mol, size=(200, 200)), caption="Структура")
                        st.markdown('</div>', unsafe_allow_html=True)

                    p_input = p_smiles

                else:
                    p_formula = st.text_input(
                        "Формула:",
                        placeholder="Например: H2O",
                        key=f"p_formula_{i}",
                        label_visibility="collapsed"
                    )
                    p_smiles = ""
                    p_name = p_formula
                    p_molar_mass = 0.0
                    p_mol = None

                    if p_formula:
                        p_molar_mass = calculate_molar_mass(p_formula)
                        if p_molar_mass > 0:
                            st.caption(f"⚖️ {p_molar_mass:.3f} г/моль")
                        else:
                            st.warning("⚠️ Некорректная формула")

                    p_input = p_formula

            p_is_solvate = st.checkbox(
                "Это сольват (есть сольватная часть)?",
                key=f"p_is_solvate_{i}",
                value=False
            )

            p_solvate_formula = ""
            p_solvate_count = 1
            p_solvate_molar_mass = 0.0

            if p_is_solvate:
                col_solv1, col_solv2 = st.columns(2)
                with col_solv1:
                    p_solvate_formula = st.text_input(
                        "Формула сольвата:",
                        placeholder="Например: H2O",
                        key=f"p_solvate_formula_{i}"
                    )
                    if p_solvate_formula:
                        p_solvate_molar_mass = calculate_molar_mass(p_solvate_formula)
                        if p_solvate_molar_mass > 0:
                            st.caption(f"⚖️ M(сольвата) = {p_solvate_molar_mass:.3f} г/моль")
                with col_solv2:
                    p_solvate_count = st.number_input(
                        "Количество молекул:",
                        min_value=1,
                        max_value=100,
                        value=1,
                        step=1,
                        key=f"p_solvate_count_{i}"
                    )

            if p_molar_mass > 0:
                if p_is_solvate and p_solvate_molar_mass > 0:
                    p_total_molar_mass = p_molar_mass + p_solvate_count * p_solvate_molar_mass
                    st.info(f"📊 M(общая) = {p_molar_mass:.3f} + {p_solvate_count} × {p_solvate_molar_mass:.3f} = **{p_total_molar_mass:.3f} г/моль**")
                else:
                    p_total_molar_mass = p_molar_mass
            else:
                p_total_molar_mass = 0.0

            st.markdown("**📐 Параметры продукта:**")

            col1, col2, col3, col4, col5, col6 = st.columns(6)
            with col1:
                st.markdown('<p class="field-label">Масса (г)</p>', unsafe_allow_html=True)
                p_mass = st.number_input(
                    "Масса",
                    min_value=0.0,
                    step=0.1,
                    key=f"p_mass_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )
            with col2:
                st.markdown('<p class="field-label">Плотность (г/мл)</p>', unsafe_allow_html=True)
                p_density = st.number_input(
                    "Плотность",
                    min_value=0.0,
                    step=0.01,
                    key=f"p_dens_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )
            with col3:
                st.markdown('<p class="field-label">Объём (л)</p>', unsafe_allow_html=True)
                p_volume = st.number_input(
                    "Объём",
                    min_value=0.0,
                    step=0.01,
                    key=f"p_vol_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )
            with col4:
                st.markdown('<p class="field-label">Масс. доля</p>', unsafe_allow_html=True)
                p_mass_fraction = st.number_input(
                    "Массовая доля",
                    min_value=0.0,
                    max_value=1.0,
                    step=0.01,
                    key=f"p_mf_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )
            with col5:
                st.markdown('<p class="field-label">Конц. (моль/л)</p>', unsafe_allow_html=True)
                p_concentration = st.number_input(
                    "Концентрация",
                    min_value=0.0,
                    step=0.1,
                    key=f"p_conc_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )
            with col6:
                st.markdown('<p class="field-label">Эквивалент</p>', unsafe_allow_html=True)
                p_equiv = st.number_input(
                    "Эквивалент продукта",
                    min_value=0.0,
                    step=0.1,
                    key=f"p_equiv_{i}",
                    value=1.0,
                    label_visibility="collapsed"
                )

            col7, col8 = st.columns([1, 1])
            with col7:
                st.markdown('<p class="field-label">Масса раствора (г)</p>', unsafe_allow_html=True)
                p_solution_mass = st.number_input(
                    "Масса раствора",
                    min_value=0.0,
                    step=0.1,
                    key=f"p_sol_mass_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )

            col9, col10 = st.columns([1, 1])
            with col9:
                st.markdown('<p class="field-label">Практическая масса (г)</p>', unsafe_allow_html=True)
                practical_mass = st.number_input(
                    "Практическая масса",
                    min_value=0.0,
                    step=0.1,
                    key=f"p_practical_{i}",
                    value=0.0,
                    label_visibility="collapsed"
                )

            products_data.append({
                'type': p_type,
                'smiles': p_smiles,
                'formula': p_formula,
                'name': p_name if p_name else f"Продукт {i+1}",
                'molar_mass': p_total_molar_mass,
                'molar_mass_main': p_molar_mass,
                'molar_mass_solvate': p_solvate_molar_mass,
                'solvate_formula': p_solvate_formula,
                'solvate_count': p_solvate_count,
                'is_solvate': p_is_solvate,
                'mass': p_mass,
                'density': p_density,
                'volume': p_volume,
                'mass_fraction': p_mass_fraction,
                'concentration': p_concentration,
                'equiv': p_equiv,
                'solution_mass': p_solution_mass,
                'practical_mass': practical_mass,
                'mol': p_mol,
                'actual_moles': 0.0,
            })

    st.divider()

    # ===== КНОПКА РАССЧИТАТЬ =====
    if st.button("🧮 Рассчитать", type="primary", use_container_width=True):

        valid_reactants = [r for r in reactants_data if r['molar_mass'] > 0]
        valid_products = [p for p in products_data if p['molar_mass'] > 0]

        if not valid_reactants:
            st.warning("⚠️ Введите корректные данные для реагентов")
            st.stop()

        st.divider()
        st.markdown('<div class="result-panel">', unsafe_allow_html=True)
        st.markdown('<h3 style="color: #BF360C;">📊 РЕЗУЛЬТАТЫ РАСЧЁТА</h3>', unsafe_allow_html=True)

        # ===== УРАВНЕНИЕ РЕАКЦИИ =====
        st.markdown('<div class="result-title">📌 УРАВНЕНИЕ РЕАКЦИИ</div>', unsafe_allow_html=True)

        r_names = [r['name'] for r in valid_reactants]
        p_names = [p['name'] for p in valid_products] if valid_products else ["продукты"]

        r_str = " + ".join(r_names)
        p_str = " + ".join(p_names)

        st.markdown(f"**{r_str} → {p_str}**")
        st.caption("(Коэффициенты реакции определяются эквивалентами)")

        st.divider()

        # ===== РАСЧЁТ ПАРАМЕТРОВ РЕАГЕНТОВ =====
        st.markdown('<div class="result-title">📌 ПАРАМЕТРЫ РЕАГЕНТОВ</div>', unsafe_allow_html=True)

        reactant_moles = {}  # для режима "по избытку и недостатку" храним n/экв
        reactant_actual_moles = {}  # для режимов "по одному реагенту" и "по продукту" храним фактическое n
        limiting_actual_moles = 0  # фактическое количество лимитирующего реагента для расчёта продуктов

        for idx, r in enumerate(valid_reactants):
            name = r['name']
            M = r['molar_mass']
            M_main = r['molar_mass_main']
            M_solvate = r['molar_mass_solvate']
            solvate_count = r['solvate_count']
            is_solvate = r['is_solvate']

            st.markdown(f"**{idx+1}. {name}**")
            st.markdown("---")

            if r['type'] == 'Органический' and r.get('mol'):
                st.markdown('<div class="structure-container">', unsafe_allow_html=True)
                st.image(Draw.MolToImage(r['mol'], size=(150, 150)), caption="Структура")
                st.markdown('</div>', unsafe_allow_html=True)

            if is_solvate and M_solvate > 0:
                st.markdown(f'<div class="calc-line">M(основная часть) = {M_main:.3f} г/моль</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-line">M(сольватная часть) = {solvate_count} × {M_solvate:.3f} = {solvate_count * M_solvate:.3f} г/моль</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-line">M(общая) = {M_main:.3f} + {solvate_count * M_solvate:.3f} = {M:.3f} г/моль</div>', unsafe_allow_html=True)
            else:
                st.markdown(f'<div class="calc-line">M = {M:.3f} г/моль</div>', unsafe_allow_html=True)

            mass = r['mass']
            solution_mass = r['solution_mass']
            mass_fraction = r['mass_fraction']
            density = r['density']
            volume = r['volume']
            concentration = r['concentration']
            equiv = r['equiv']

            actual_mass = 0.0
            moles = 0.0

            # ===== НОВЫЙ ПРИОРИТЕТ С УЧЁТОМ МАССОВОЙ ДОЛИ =====
            if mass > 0:
                # Приоритет 1: просто масса
                actual_mass = mass
                moles = actual_mass / M
                st.markdown(f'<div class="calc-line">✅ Расчёт по массе:</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-formula">m = {mass:.2f} г</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-formula">n = m / M = {mass:.2f} / {M:.3f} = {moles:.4f} моль</div>', unsafe_allow_html=True)

            elif solution_mass > 0 and mass_fraction > 0:
                # Приоритет 2: масса раствора + массовая доля
                actual_mass = solution_mass * mass_fraction
                moles = actual_mass / M
                st.markdown(f'<div class="calc-line">✅ Расчёт по массе раствора и массовой доле:</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-formula">m = m(раствора) × ω = {solution_mass:.2f} × {mass_fraction:.2f} = {actual_mass:.2f} г</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-formula">n = m / M = {actual_mass:.2f} / {M:.3f} = {moles:.4f} моль</div>', unsafe_allow_html=True)

            elif volume > 0 and density > 0 and mass_fraction > 0:
                # 🆕 Приоритет 3: объём + плотность + массовая доля
                solution_mass_calc = volume * density * 1000
                actual_mass = solution_mass_calc * mass_fraction
                moles = actual_mass / M
                st.markdown(f'<div class="calc-line">✅ Расчёт по объёму, плотности и массовой доле:</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-formula">m(р-ра) = V × ρ × 1000 = {volume:.2f} × {density:.2f} × 1000 = {solution_mass_calc:.2f} г</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-formula">m(в-ва) = m(р-ра) × ω = {solution_mass_calc:.2f} × {mass_fraction:.2f} = {actual_mass:.2f} г</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-formula">n = m / M = {actual_mass:.2f} / {M:.3f} = {moles:.4f} моль</div>', unsafe_allow_html=True)

            elif volume > 0 and density > 0:
                # Приоритет 4: объём + плотность (без массовой доли)
                actual_mass = volume * density * 1000
                moles = actual_mass / M
                st.markdown(f'<div class="calc-line">✅ Расчёт по объёму и плотности:</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-formula">m = V × ρ × 1000 = {volume:.2f} × {density:.2f} × 1000 = {actual_mass:.2f} г</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-formula">n = m / M = {actual_mass:.2f} / {M:.3f} = {moles:.4f} моль</div>', unsafe_allow_html=True)

            elif volume > 0 and concentration > 0:
                # Приоритет 5: объём + концентрация
                moles = concentration * volume
                st.markdown(f'<div class="calc-line">✅ Расчёт по объёму и концентрации:</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-formula">n = C × V = {concentration:.2f} × {volume:.2f} = {moles:.4f} моль</div>', unsafe_allow_html=True)

            else:
                moles = 0.0

            # Сохраняем фактическое количество вещества
            r['actual_moles'] = moles
            reactant_actual_moles[idx] = moles

            # Для режима "по избытку и недостатку" сохраняем n/экв
            if equiv > 0 and moles > 0:
                normalized = moles / equiv
                st.markdown(f'<div class="calc-line">n / экв = {moles:.4f} / {equiv:.1f} = {normalized:.4f} моль</div>', unsafe_allow_html=True)
                reactant_moles[idx] = normalized
            elif moles > 0:
                st.markdown(f'<div class="calc-line">n / экв = {moles:.4f} / 1.0 = {moles:.4f} моль</div>', unsafe_allow_html=True)
                reactant_moles[idx] = moles
            else:
                reactant_moles[idx] = 0.0

            # ===== ДОПОЛНИТЕЛЬНЫЕ ПАРАМЕТРЫ (РАССЧИТАНЫ СИСТЕМОЙ) =====
            show_additional_params(actual_mass, moles, volume, solution_mass, actual_mass, volume, "")

            st.divider()

        # ===== ОПРЕДЕЛЕНИЕ ЛИМИТИРУЮЩЕГО РЕАГЕНТА =====
        st.markdown('<div class="result-title">🔬 ОПРЕДЕЛЕНИЕ ЛИМИТИРУЮЩЕГО РЕАГЕНТА</div>', unsafe_allow_html=True)

        if "⚖️ Расчет по избытку и недостатку" in calc_mode:
            if reactant_moles and any(m > 0 for m in reactant_moles.values()):
                valid_moles = {idx: m for idx, m in reactant_moles.items() if m > 0}
                if valid_moles:
                    limiting_idx = min(valid_moles, key=valid_moles.get)
                    limiting_name = valid_reactants[limiting_idx]['name']
                    limiting_moles = valid_moles[limiting_idx]
                    limiting_actual_moles = valid_reactants[limiting_idx]['actual_moles']

                    st.markdown("**Сравнение реагентов (по моль/экв):**")
                    for idx, r in enumerate(valid_reactants):
                        if idx in reactant_moles:
                            name = r['name']
                            moles = reactant_moles[idx]
                            if moles > 0:
                                if idx == limiting_idx:
                                    st.markdown(f"• **{name}**: {moles:.4f} моль ← **МИНИМУМ**")
                                else:
                                    st.markdown(f"• {name}: {moles:.4f} моль")

                    st.markdown('<div class="limiting-highlight">', unsafe_allow_html=True)
                    st.markdown(f"**⚠️ Лимитирующий реагент: {limiting_name}**")
                    st.markdown(f"**n(экв) = {limiting_moles:.4f} моль**")
                    st.markdown(f"**n(фактическое) = {limiting_actual_moles:.4f} моль**")
                    st.markdown(f"Так как это значение минимальное, {limiting_name} расходуется первым.")
                    st.markdown("</div>", unsafe_allow_html=True)

                    limiting_moles_total = limiting_moles
                    limiting_idx_total = limiting_idx
                    known_moles = limiting_actual_moles
                    known_equiv = valid_reactants[limiting_idx]['equiv']
                    known_name = limiting_name
                else:
                    st.warning("⚠️ Нет данных для определения лимитирующего реагента")
                    limiting_moles_total = 0
                    limiting_idx_total = 0
                    known_moles = 0
                    known_equiv = 1.0
                    known_name = ""
                    limiting_actual_moles = 0
            else:
                st.warning("⚠️ Нет данных для определения лимитирующего реагента")
                limiting_moles_total = 0
                limiting_idx_total = 0
                known_moles = 0
                known_equiv = 1.0
                known_name = ""
                limiting_actual_moles = 0

        elif "🔬 Расчет по одному реагенту" in calc_mode:
            st.info("ℹ️ В режиме 'Расчет по одному реагенту' лимитирующий реагент не определяется")
            
            # Ищем реагент с данными
            found_known = False
            for idx, r in enumerate(valid_reactants):
                has_data = (
                    r['mass'] > 0 or 
                    r['solution_mass'] > 0 or 
                    (r['volume'] > 0 and r['density'] > 0) or 
                    (r['volume'] > 0 and r['concentration'] > 0) or
                    (r['volume'] > 0 and r['density'] > 0 and r['mass_fraction'] > 0)
                )
                if has_data and r['actual_moles'] > 0:
                    limiting_idx_total = idx
                    known_moles = r['actual_moles']
                    known_equiv = r['equiv']
                    known_name = r['name']
                    found_known = True
                    st.success(f"✅ Известный реагент: {r['name']} (n = {known_moles:.4f} моль, экв = {known_equiv:.1f})")
                    break
            
            if not found_known:
                st.warning("⚠️ Нет данных ни для одного реагента")
                known_moles = 0
                known_equiv = 1.0
                known_name = ""
                limiting_idx_total = 0

        else:  # 📦 Расчет по продукту
            st.info("ℹ️ В режиме 'Расчет по продукту' известным считается продукт, реагенты и остальные продукты рассчитываются относительно него")
            
            found_known = False
            known_moles = 0
            known_equiv = 1.0
            known_name = ""
            known_product_idx = 0
            
            for idx, p in enumerate(valid_products):
                has_data = (
                    p['mass'] > 0 or 
                    p['solution_mass'] > 0 or 
                    (p['volume'] > 0 and p['density'] > 0) or 
                    (p['volume'] > 0 and p['concentration'] > 0) or
                    (p['volume'] > 0 and p['density'] > 0 and p['mass_fraction'] > 0)
                )
                if has_data and p['molar_mass'] > 0:
                    M = p['molar_mass']
                    mass = p['mass']
                    solution_mass = p['solution_mass']
                    mass_fraction = p['mass_fraction']
                    density = p['density']
                    volume = p['volume']
                    concentration = p['concentration']
                    equiv = p['equiv']
                    
                    moles = 0.0
                    
                    if mass > 0:
                        moles = mass / M
                    elif solution_mass > 0 and mass_fraction > 0:
                        actual_mass = solution_mass * mass_fraction
                        moles = actual_mass / M
                    elif volume > 0 and density > 0 and mass_fraction > 0:
                        solution_mass_calc = volume * density * 1000
                        actual_mass = solution_mass_calc * mass_fraction
                        moles = actual_mass / M
                    elif volume > 0 and density > 0:
                        actual_mass = volume * density * 1000
                        moles = actual_mass / M
                    elif volume > 0 and concentration > 0:
                        moles = concentration * volume
                    
                    if moles > 0:
                        p['actual_moles'] = moles
                        known_moles = moles
                        known_equiv = equiv
                        known_name = p['name']
                        known_product_idx = idx
                        found_known = True
                        st.success(f"✅ Известный продукт: {p['name']} (n = {known_moles:.4f} моль, экв = {known_equiv:.1f})")
                        break
            
            if not found_known:
                st.warning("⚠️ Нет данных ни для одного продукта")
                known_moles = 0
                known_equiv = 1.0
                known_name = ""
                limiting_idx_total = 0

        st.divider()

        # ===== РАСХОД РЕАГЕНТОВ =====
        st.markdown('<div class="result-title">📊 РАСЧЁТ РАСХОДА РЕАГЕНТОВ</div>', unsafe_allow_html=True)

        if "⚖️ Расчет по избытку и недостатку" in calc_mode:
            if limiting_moles_total > 0:
                st.markdown(f"**Лимитирующий реагент: {known_name}**")
                st.markdown(f"**n(экв) = {limiting_moles_total:.4f} моль**")
                st.markdown(f"**n(фактическое) = {known_moles:.4f} моль**")
                st.markdown("---")

                for idx, r in enumerate(valid_reactants):
                    name = r['name']
                    M = r['molar_mass']
                    equiv = r['equiv']

                    if idx == limiting_idx_total:
                        st.markdown(f"**{name} (лимитирующий):**")
                        st.markdown(f'<div class="calc-line">Расходуется полностью: {known_moles:.4f} моль</div>', unsafe_allow_html=True)
                        st.markdown(f'<div class="calc-line">m = n × M = {known_moles:.4f} × {M:.3f} = {known_moles * M:.2f} г</div>', unsafe_allow_html=True)
                    else:
                        needed_moles = known_moles * (equiv / known_equiv)
                        actual_moles = reactant_actual_moles.get(idx, 0)
                        remaining_moles = actual_moles - needed_moles

                        st.markdown(f"**{name}:**")
                        st.markdown(f'<div class="calc-line">n({name}) = n(лимитирующего_факт) × (экв({name}) / экв(лимитирующего))</div>', unsafe_allow_html=True)
                        st.markdown(f'<div class="calc-line">n({name}) = {known_moles:.4f} × ({equiv:.1f} / {known_equiv:.1f}) = {needed_moles:.4f} моль</div>', unsafe_allow_html=True)
                        st.markdown(f'<div class="calc-line">m = n × M = {needed_moles:.4f} × {M:.3f} = {needed_moles * M:.2f} г</div>', unsafe_allow_html=True)
                        st.markdown(f'<div class="calc-line">Есть: {actual_moles:.4f} моль</div>', unsafe_allow_html=True)
                        if remaining_moles > 0:
                            st.markdown(f'<div class="calc-line">✅ Остаток: {remaining_moles:.4f} моль ({remaining_moles * M:.2f} г)</div>', unsafe_allow_html=True)
                        elif remaining_moles < 0:
                            st.markdown(f'<div class="calc-line">⚠️ Недостаток: {-remaining_moles:.4f} моль ({-remaining_moles * M:.2f} г)</div>', unsafe_allow_html=True)
                        else:
                            st.markdown(f'<div class="calc-line">✅ Точно по стехиометрии</div>', unsafe_allow_html=True)

                        # Дополнительные параметры для продуктов (в режиме по избытку)
                        if idx != limiting_idx_total:
                            mass_product = needed_moles * M
                            show_additional_params(mass_product, needed_moles, r['volume'], r['solution_mass'], mass_product, r['volume'], f"для {name}")

                    st.divider()
            else:
                st.warning("⚠️ Недостаточно данных для расчёта расхода реагентов")

        else:
            if known_moles > 0:
                if "🔬 Расчет по одному реагенту" in calc_mode:
                    st.markdown(f"**Известный реагент: {known_name}**")
                    st.markdown(f"**n = {known_moles:.4f} моль**")
                    st.markdown(f"**экв = {known_equiv:.1f}**")
                    st.markdown("---")

                    for idx, r in enumerate(valid_reactants):
                        name = r['name']
                        M = r['molar_mass']
                        equiv = r['equiv']

                        if idx == limiting_idx_total:
                            st.markdown(f"**{name} (известный):**")
                            st.markdown(f'<div class="calc-line">n = {known_moles:.4f} моль</div>', unsafe_allow_html=True)
                            st.markdown(f'<div class="calc-line">m = n × M = {known_moles:.4f} × {M:.3f} = {known_moles * M:.2f} г</div>', unsafe_allow_html=True)
                        else:
                            needed_moles = known_moles * (equiv / known_equiv)
                            actual_moles = reactant_actual_moles.get(idx, 0)

                            st.markdown(f"**{name}:**")
                            st.markdown(f'<div class="calc-line">n({name}) = n(известного) × (экв({name}) / экв(известного))</div>', unsafe_allow_html=True)
                            st.markdown(f'<div class="calc-line">n({name}) = {known_moles:.4f} × ({equiv:.1f} / {known_equiv:.1f}) = {needed_moles:.4f} моль</div>', unsafe_allow_html=True)
                            st.markdown(f'<div class="calc-line">m = n × M = {needed_moles:.4f} × {M:.3f} = {needed_moles * M:.2f} г</div>', unsafe_allow_html=True)
                            if actual_moles > 0:
                                st.markdown(f'<div class="calc-line">Есть: {actual_moles:.4f} моль</div>', unsafe_allow_html=True)
                                remaining_moles = actual_moles - needed_moles
                                if remaining_moles > 0:
                                    st.markdown(f'<div class="calc-line">✅ Остаток: {remaining_moles:.4f} моль ({remaining_moles * M:.2f} г)</div>', unsafe_allow_html=True)
                                elif remaining_moles < 0:
                                    st.markdown(f'<div class="calc-line">⚠️ Недостаток: {-remaining_moles:.4f} моль ({-remaining_moles * M:.2f} г)</div>', unsafe_allow_html=True)

                            # Дополнительные параметры для продуктов
                            if idx != limiting_idx_total:
                                mass_product = needed_moles * M
                                show_additional_params(mass_product, needed_moles, r['volume'], r['solution_mass'], mass_product, r['volume'], f"для {name}")

                        st.divider()

                else:
                    st.markdown(f"**Известный продукт: {known_name}**")
                    st.markdown(f"**n = {known_moles:.4f} моль**")
                    st.markdown(f"**экв = {known_equiv:.1f}**")
                    st.markdown("---")

                    for idx, r in enumerate(valid_reactants):
                        name = r['name']
                        M = r['molar_mass']
                        equiv = r['equiv']

                        needed_moles = known_moles * (equiv / known_equiv)
                        actual_moles = reactant_actual_moles.get(idx, 0)

                        st.markdown(f"**{name}:**")
                        st.markdown(f'<div class="calc-line">n({name}) = n(продукта) × (экв({name}) / экв(продукта))</div>', unsafe_allow_html=True)
                        st.markdown(f'<div class="calc-line">n({name}) = {known_moles:.4f} × ({equiv:.1f} / {known_equiv:.1f}) = {needed_moles:.4f} моль</div>', unsafe_allow_html=True)
                        st.markdown(f'<div class="calc-line">m = n × M = {needed_moles:.4f} × {M:.3f} = {needed_moles * M:.2f} г</div>', unsafe_allow_html=True)
                        if actual_moles > 0:
                            st.markdown(f'<div class="calc-line">Есть: {actual_moles:.4f} моль</div>', unsafe_allow_html=True)
                            remaining_moles = actual_moles - needed_moles
                            if remaining_moles > 0:
                                st.markdown(f'<div class="calc-line">✅ Остаток: {remaining_moles:.4f} моль ({remaining_moles * M:.2f} г)</div>', unsafe_allow_html=True)
                            elif remaining_moles < 0:
                                st.markdown(f'<div class="calc-line">⚠️ Недостаток: {-remaining_moles:.4f} моль ({-remaining_moles * M:.2f} г)</div>', unsafe_allow_html=True)

                        # Дополнительные параметры для реагентов (в режиме по продукту)
                        mass_product = needed_moles * M
                        show_additional_params(mass_product, needed_moles, r['volume'], r['solution_mass'], mass_product, r['volume'], f"для {name}")

                        st.divider()
            else:
                st.warning("⚠️ Недостаточно данных для расчёта расхода реагентов")

        # ===== ОБРАЗОВАНИЕ ПРОДУКТОВ =====
        st.markdown('<div class="result-title">🧪 ОБРАЗОВАНИЕ ПРОДУКТОВ</div>', unsafe_allow_html=True)

        if known_moles > 0 and valid_products:
            for idx, p in enumerate(valid_products):
                name = p['name']
                M = p['molar_mass']
                equiv = p['equiv']
                M_main = p['molar_mass_main']
                M_solvate = p['molar_mass_solvate']
                solvate_count = p['solvate_count']
                is_solvate = p['is_solvate']
                volume = p['volume']
                solution_mass = p['solution_mass']

                if "📦 Расчет по продукту" in calc_mode and idx == known_product_idx:
                    formed_moles = known_moles
                else:
                    formed_moles = known_moles * (equiv / known_equiv)
                
                formed_mass = formed_moles * M

                st.markdown(f"**{name}:**")
                
                if p['type'] == 'Органический' and p.get('mol'):
                    st.markdown('<div class="structure-container">', unsafe_allow_html=True)
                    st.image(Draw.MolToImage(p['mol'], size=(150, 150)), caption="Структура")
                    st.markdown('</div>', unsafe_allow_html=True)
                
                if "📦 Расчет по продукту" in calc_mode and idx == known_product_idx:
                    st.markdown(f'<div class="calc-line">✅ ИЗВЕСТНЫЙ ПРОДУКТ</div>', unsafe_allow_html=True)
                
                st.markdown(f'<div class="calc-line">n = n(известного) × (экв({name}) / экв(известного))</div>', unsafe_allow_html=True)
                st.markdown(f'<div class="calc-line">n = {known_moles:.4f} × ({equiv:.1f} / {known_equiv:.1f}) = {formed_moles:.4f} моль</div>', unsafe_allow_html=True)

                if is_solvate and M_solvate > 0:
                    st.markdown(f'<div class="calc-line">M(общая) = {M:.3f} г/моль</div>', unsafe_allow_html=True)
                    st.markdown(f'<div class="calc-line">m = n × M = {formed_moles:.4f} × {M:.3f} = {formed_mass:.2f} г</div>', unsafe_allow_html=True)
                    st.markdown(f'<div class="calc-line">Из них: основная часть — {formed_moles * M_main:.2f} г, сольватная часть — {formed_moles * solvate_count * M_solvate:.2f} г</div>', unsafe_allow_html=True)
                else:
                    st.markdown(f'<div class="calc-line">m = n × M = {formed_moles:.4f} × {M:.3f} = {formed_mass:.2f} г</div>', unsafe_allow_html=True)

                # Дополнительные параметры для продуктов
                show_additional_params(formed_mass, formed_moles, volume, solution_mass, formed_mass, volume, f"для {name}")

                st.divider()
        else:
            if not valid_products:
                st.info("ℹ️ Нет продуктов для расчёта")
            else:
                st.warning("⚠️ Недостаточно данных для расчёта образования продуктов")

        # ===== ВЫХОД ПРОДУКТОВ =====
        st.markdown('<div class="result-title">🎯 ВЫХОД ПРОДУКТОВ</div>', unsafe_allow_html=True)

        if known_moles > 0 and valid_products:
            for idx, p in enumerate(valid_products):
                name = p['name']
                M = p['molar_mass']
                equiv = p['equiv']

                if "📦 Расчет по продукту" in calc_mode and idx == known_product_idx:
                    theoretical_mass = known_moles * M
                else:
                    formed_moles = known_moles * (equiv / known_equiv)
                    theoretical_mass = formed_moles * M

                st.markdown(f"**{name}**")
                st.markdown(f'<div class="calc-line">Теоретическая масса: {theoretical_mass:.2f} г</div>', unsafe_allow_html=True)

                practical_mass = p['practical_mass']

                if practical_mass > 0:
                    yield_percent = (practical_mass / theoretical_mass) * 100
                    st.markdown(f'<div class="calc-line">Выход = ({practical_mass:.2f} / {theoretical_mass:.2f}) × 100% = {yield_percent:.1f}%</div>', unsafe_allow_html=True)

                    if yield_percent > 100:
                        st.warning("⚠️ Выход > 100%! Проверьте введённые данные.")
                    elif yield_percent > 90:
                        st.success(f"✅ Отличный выход: {yield_percent:.1f}%")
                    elif yield_percent > 70:
                        st.info(f"✅ Хороший выход: {yield_percent:.1f}%")
                    else:
                        st.warning(f"⚠️ Низкий выход: {yield_percent:.1f}%")
                else:
                    st.info("ℹ️ Введите практическую массу для расчёта выхода")

                st.divider()
        else:
            if not valid_products:
                st.info("ℹ️ Нет продуктов для расчёта выхода")
            else:
                st.warning("⚠️ Недостаточно данных для расчёта выхода")

        st.markdown('</div>', unsafe_allow_html=True)
        st.session_state.calc_done = True
