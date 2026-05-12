import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
import pubchempy as pcp
import urllib.parse
import re

st.set_page_config(page_title="ChemRadar", layout="wide")

col_logo1, col_logo2, col_logo3 = st.columns([1, 10, 1])
with col_logo2:
    st.markdown(
        "<div style='text-align: center;'>"
        "<h1>📡🧪 ChemRadar</h1>"
        "<p style='font-size: 1.1rem; color: #666;'>Chemical Literature Scanner — search scientific papers by SMILES and keywords</p>"
        "</div>", 
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

Rimonabant + tuberculosis:
rimonabant AND tuberculosis

Rimonabant or SR141716:
rimonabant OR SR141716

Rimonabant + tuberculosis + resistance:
(rimonabant OR SR141716) AND (tuberculosis OR TB) AND resistance

Exact phrase for derivatives:
"rimonabant derivative" AND antitubercular

Exclude review articles:
rimonabant NOT review NOT editorial

---

### What works where

| System | Simple | Advanced |
|--------|:------:|:--------:|
| 📖 PubMed | ✅ | ✅ |
| 🎓 Google Scholar | ✅ | ⚠️ |
| 📚 Europe PMC | ✅ | ✅ |
| 📑 Google Patents | ✅ | ❌ |
| 🌐 Google Web | ✅ | ❌ |
| 🧬 PubChem | ✅ | ❌ |

✅ — supports AND, OR, NOT, quotes, brackets, ranges
⚠️ — partial support
❌ — simple words only

---

### Important notes

- Keywords in **English only**
- Use UPPERCASE for operators (AND, OR, NOT)
- Space between words means AND
- Quotes required for exact phrases
""")

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

def get_safe_inchikey(mol):
    try:
        return Chem.MolToInchiKey(mol)
    except:
        try:
            Chem.RemoveStereochemistry(mol)
            return Chem.MolToInchiKey(mol)
        except:
            return "N/A"

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

def build_search_urls(name, keywords, cid, ikey, smiles):
    urls = {}
    
    if keywords:
        search_terms = keywords
    elif name:
        search_terms = name
    else:
        search_terms = None
    
    if search_terms:
        urls['pubmed'] = f"https://pubmed.ncbi.nlm.nih.gov/?term={urllib.parse.quote(search_terms)}"
        urls['scholar'] = f"https://scholar.google.com/scholar?q={urllib.parse.quote(search_terms)}"
        urls['europe_pmc'] = f"https://europepmc.org/search?query={urllib.parse.quote(search_terms)}"
        
        # Исправленная строка google_exclude
        simple_terms = search_terms.split()[0] if search_terms.split() else search_terms
        google_exclude = '-buy -price -pharmacy -shop -store -amazon -ebay -tablet -drugstore -review -rating -\"store\" -price -cost -\"for sale\" -\"online shop\"'
        urls['google'] = f"https://www.google.com/search?q={urllib.parse.quote(simple_terms)}+{urllib.parse.quote(google_exclude)}"
    
    urls['patents'] = f"https://patents.google.com/?q={ikey}"
    
    if cid:
        urls['pubchem'] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    else:
        urls['pubchem'] = f"https://pubchem.ncbi.nlm.nih.gov/#query={urllib.parse.quote(smiles[:80])}"
    
    return urls

def find_and_show_analogs(cid, smiles, name, keywords, ikey):
    st.divider()
    st.subheader("🔬 Structural Analogs")
    st.caption("Compounds with similar chemical structure")
    
    if not cid and not name:
        st.info("No identifier available for analog search")
        return
    
    analogs_list = []
    try:
        with st.spinner("Searching for similar structures in PubChem..."):
            compounds = pcp.get_compounds(smiles, 'similarity', threshold=85)
            analogs_list = [c for c in compounds if c.cid != cid][:4]
    except:
        pass
    
    if analogs_list:
        st.markdown("**Found similar compounds:**")
        cols = st.columns(len(analogs_list))
        
        for idx, comp in enumerate(analogs_list):
            with cols[idx]:
                try:
                    analog_smiles = comp.canonical_smiles
                    analog_mol = Chem.MolFromSmiles(analog_smiles)
                    
                    if analog_mol:
                        st.image(Draw.MolToImage(analog_mol, size=(150, 150)), caption=f"CID: {comp.cid}")
                        
                        if comp.iupac_name:
                            short_name = comp.iupac_name[:40] + "..." if len(comp.iupac_name) > 40 else comp.iupac_name
                            st.caption(short_name)
                        elif comp.synonyms:
                            short_syn = comp.synonyms[0][:40] + "..." if len(comp.synonyms[0]) > 40 else comp.synonyms[0]
                            st.caption(short_syn)
                        
                        st.markdown("**Links:**")
                        
                        pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={comp.cid}[uid]"
                        if keywords:
                            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/?term=({comp.cid}[uid])+AND+({keywords})"
                        st.link_button("📖 PubMed", pubmed_url, use_container_width=True)
                        
                        pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{comp.cid}"
                        st.link_button("🧬 PubChem", pubchem_url, use_container_width=True)
                        
                        scholar_query = comp.iupac_name if comp.iupac_name else str(comp.cid)
                        scholar_url = f"https://scholar.google.com/scholar?q={urllib.parse.quote(scholar_query)}"
                        st.link_button("🎓 Scholar", scholar_url, use_container_width=True)
                        
                except:
                    pass
    else:
        st.info("No similar compounds automatically found. Try searching manually:")
        
        st.markdown("**Search for structural analogs manually:**")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            if cid:
                url_pubchem = f"https://pubchem.ncbi.nlm.nih.gov/#query=CID{cid}&sort=similarity"
                st.link_button("🧬 PubChem Similar", url_pubchem, use_container_width=True)
        
        with col2:
            if name:
                url_scholar = f"https://scholar.google.com/scholar?q={urllib.parse.quote(name)}+analog"
                st.link_button("🎓 Scholar: analog", url_scholar, use_container_width=True)
            elif keywords:
                url_scholar = f"https://scholar.google.com/scholar?q={urllib.parse.quote(keywords)}+analog"
                st.link_button("🎓 Scholar: analog", url_scholar, use_container_width=True)
        
        with col3:
            if smiles:
                url_zinc = f"https://zinc.docking.org/substances/similar/{urllib.parse.quote(smiles[:50])}/"
                st.link_button("⚡ ZINC Similar", url_zinc, use_container_width=True)
        
        with col4:
            if ikey and ikey != "N/A":
                url_chemspider = f"https://www.chemspider.com/Search.aspx?q={ikey}"
                st.link_button("🔬 ChemSpider", url_chemspider, use_container_width=True)
            elif name:
                url_chemspider = f"https://www.chemspider.com/Search.aspx?q={urllib.parse.quote(name)}"
                st.link_button("🔬 ChemSpider", url_chemspider, use_container_width=True)

left_col, right_col = st.columns([1, 1], gap="medium")

with left_col:
    st.markdown("### 📝 Enter SMILES formula")
    user_smiles = st.text_area(
        "SMILES:", 
        height=100,
        placeholder="Example: CC1=C(C(=O)NN2CCCCC2)N(N=C1C3=CC=C(C=C3)Cl)C4=C(C=C(C=C4)Cl)Cl",
        label_visibility="collapsed"
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
            ikey = get_safe_inchikey(mol)
            pubchem_data = get_pubchem_data(user_smiles)
            cid = pubchem_data["cid"] if pubchem_data else None
            name = pubchem_data["name"] if pubchem_data else None
            
            st.image(Draw.MolToImage(mol, size=(300, 300)), caption="Molecule Structure")
            
            st.code(f"InChIKey: {ikey}")
            if cid:
                st.success(f"PubChem CID: {cid}")
            if name:
                st.info(f"📌 Name: {name[:80]}..." if len(name) > 80 else f"📌 Name: {name}")
            if warning:
                st.warning(warning)
            if keywords:
                st.caption(f"🔍 Query: {keywords[:100]}" if len(keywords) > 100 else f"🔍 Query: {keywords}")
            
            st.divider()
            
            urls = build_search_urls(name, keywords, cid, ikey, user_smiles)
            
            st.subheader("🔍 Search Literature")
            
            st.markdown("**📚 Academic Databases:**")
            if "pubmed" in urls:
                st.link_button("📖 PubMed", urls["pubmed"], use_container_width=True)
            if "scholar" in urls:
                st.link_button("🎓 Google Scholar", urls["scholar"], use_container_width=True)
            if "europe_pmc" in urls:
                st.link_button("📚 Europe PMC", urls["europe_pmc"], use_container_width=True)
            
            st.markdown("---")
            st.markdown("**📄 Patents & Compounds:**")
            st.link_button("📑 Google Patents", urls["patents"], use_container_width=True)
            st.link_button("🧬 PubChem", urls["pubchem"], use_container_width=True)
            
            st.markdown("---")
            st.markdown("**🌐 General Web:**")
            if "google" in urls:
                st.link_button("🌐 Google", urls["google"], use_container_width=True)
            
            find_and_show_analogs(cid, user_smiles, name, keywords, ikey)
            
            st.session_state.processed = True
            
        else:
            st.error(warning if warning else "❌ Invalid SMILES format")
            st.info("💡 Try copying SMILES from [PubChem](https://pubchem.ncbi.nlm.nih.gov)")
            st.session_state.processed = False

st.divider()
st.markdown(
    "<div style='text-align: center; color: #888; font-size: 0.8rem;'>"
    "📡🧪 ChemRadar — Chemical Literature Scanner | Search by SMILES + keywords"
    "</div>", 
    unsafe_allow_html=True
)
