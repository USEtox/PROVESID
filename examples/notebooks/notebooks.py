import marimo

__generated_with = "0.20.4"
app = marimo.App()


@app.cell
def _(mo):
    mo.md(r"""
    # PROVESID: Offline Identifier Reconciliation Demo

    This notebook tests the latest offline resolver functions on local data files in this folder:
    - `ids_from_CAS`
    - `ids_from_name`
    - `ids_from_SMILES`

    All lookups use offline sources (`ChebiSDF`, `CompToxID`, `PubChemID`, `ZeroPM`, `CheMBL`).
    """)
    return


@app.cell
def _():
    from pathlib import Path
    import pandas as pd
    from IPython.display import display

    from provesid import PubChemID, CompToxID, ZeroPM, ChebiSDF, CheMBL
    from provesid.tools import ids_from_CAS, ids_from_name, ids_from_SMILES

    return (
        CheMBL,
        ChebiSDF,
        CompToxID,
        Path,
        PubChemID,
        ZeroPM,
        display,
        ids_from_CAS,
        ids_from_SMILES,
        ids_from_name,
        pd,
    )


@app.cell
def _(mo):
    mo.md(r"""
    ## Load local data files

    The notebook searches for the provided files in the current directory, and falls back to `examples/notebooks` if needed.
    """)
    return


@app.cell
def _(Path, display, pd):
    DATA_DIR = Path('.')
    if not (DATA_DIR / 'unique_cas_list.csv').exists():
        DATA_DIR = Path('examples/notebooks')

    cas_path = DATA_DIR / 'unique_cas_list.csv'
    esol_path = DATA_DIR / 'solubility_data_ESOL.csv'

    cas_df = pd.read_csv(cas_path)
    esol_df = pd.read_csv(esol_path)

    print('Using data directory:', DATA_DIR.resolve())
    print('CAS rows:', len(cas_df))
    print('ESOL rows:', len(esol_df))

    display(cas_df.head(3))
    display(esol_df[['chemical_name', 'SMILES']].head(3))
    return DATA_DIR, cas_df, esol_df


@app.cell
def _(mo):
    mo.md(r"""
    ## Initialize offline databases

    On first run, local data artifacts may be downloaded and cached.
    """)
    return


@app.cell
def _(CheMBL, ChebiSDF, CompToxID, PubChemID, ZeroPM):
    pubchem = PubChemID()
    comptox = CompToxID()
    zeropm = ZeroPM()
    chebi = ChebiSDF()
    chembl = CheMBL()

    print('Offline clients initialized.')
    return chebi, chembl, comptox, pubchem, zeropm


@app.cell
def _(mo):
    mo.md(r"""
    ## Test 1: CASRN to identifiers

    This section uses `ids_from_CAS` on a sample from `unique_cas_list.csv`.
    """)
    return


@app.cell
def _(
    cas_df,
    chebi,
    chembl,
    comptox,
    display,
    ids_from_CAS,
    pd,
    pubchem,
    zeropm,
):
    cas_sample = cas_df['input_CAS'].dropna().astype(str).head(25).tolist()

    cas_results = pd.DataFrame([
        ids_from_CAS(
            cas,
            chebi=chebi,
            comptox=comptox,
            pubchem=pubchem,
            zeropm=zeropm,
            chembl=chembl,
        )
        for cas in cas_sample
    ])

    display(cas_results.head(10))
    return (cas_results,)


@app.cell
def _(cas_results, display, pd):
    cas_summary = pd.DataFrame({
        'total': [len(cas_results)],
        'with_smiles': [cas_results['SMILES'].notna().sum()],
        'with_inchikey': [cas_results['InChIKey'].notna().sum()],
        'unique_sources': [cas_results['source'].nunique(dropna=True)]
    })

    display(cas_summary)
    display(cas_results['source'].value_counts(dropna=False).rename_axis('source').to_frame('count'))
    return


@app.cell
def _():
    print('Test 2: name to identifiers with consensus scoring')
    print('Input column: chemical_name from solubility_data_ESOL.csv')
    return


@app.cell
def _(
    chebi,
    chembl,
    comptox,
    display,
    esol_df,
    ids_from_name,
    pd,
    pubchem,
    zeropm,
):
    name_sample = esol_df['chemical_name'].dropna().astype(str).head(30).tolist()

    name_results = pd.DataFrame([
        ids_from_name(
            name,
            chebi=chebi,
            comptox=comptox,
            pubchem=pubchem,
            zeropm=zeropm,
            chembl=chembl,
        )
        for name in name_sample
    ])

    display(name_results[['query', 'CASRN', 'SMILES', 'source', 'consensus_source', 'match_score']].head(15))
    return (name_results,)


@app.cell
def _(display, name_results, pd):
    name_summary = pd.DataFrame({
        'total': [len(name_results)],
        'with_smiles': [name_results['SMILES'].notna().sum()],
        'with_cas': [name_results['CASRN'].notna().sum()],
        'avg_match_score': [name_results['match_score'].mean()]
    })

    display(name_summary)
    display(name_results[['query', 'source', 'consensus_source', 'match_score']].sort_values('match_score').head(10))
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Test 3: SMILES to identifiers with consensus scoring

    This section uses `ids_from_SMILES` on the `SMILES` column from `solubility_data_ESOL.csv`.
    """)
    return


@app.cell
def _(
    chebi,
    chembl,
    comptox,
    display,
    esol_df,
    ids_from_SMILES,
    pd,
    pubchem,
    zeropm,
):
    smiles_sample = esol_df['SMILES'].dropna().astype(str).head(30).tolist()

    smiles_results = pd.DataFrame([
        ids_from_SMILES(
            smi,
            chebi=chebi,
            comptox=comptox,
            pubchem=pubchem,
            zeropm=zeropm,
            chembl=chembl,
        )
        for smi in smiles_sample
    ])

    display(smiles_results[['query', 'CASRN', 'name', 'source', 'consensus_source', 'match_score']].head(15))
    return (smiles_results,)


@app.cell
def _(display, pd, smiles_results):
    smiles_summary = pd.DataFrame({
        'total': [len(smiles_results)],
        'with_cas': [smiles_results['CASRN'].notna().sum()],
        'with_name': [smiles_results['name'].notna().sum()],
        'avg_match_score': [smiles_results['match_score'].mean()]
    })

    display(smiles_summary)
    display(smiles_results[['query', 'source', 'consensus_source', 'match_score']].sort_values('match_score').head(10))
    return


@app.cell
def _(DATA_DIR, cas_results, name_results, smiles_results):
    cas_out = DATA_DIR / 'demo_ids_from_CAS_results.csv'
    name_out = DATA_DIR / 'demo_ids_from_name_results.csv'
    smiles_out = DATA_DIR / 'demo_ids_from_SMILES_results.csv'

    cas_results.to_csv(cas_out, index=False)
    name_results.to_csv(name_out, index=False)
    smiles_results.to_csv(smiles_out, index=False)

    print('Saved:', cas_out)
    print('Saved:', name_out)
    print('Saved:', smiles_out)
    return


if __name__ == "__main__":
    app.run()
