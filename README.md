# PROVESID
Access online services of chemical identifiers from Python. The goal is to have a clean interface to the most important online databases with a clean and simple interface. We offer interfaces to [PubChem](), [NCI chemical identifier resolver](), [CAS Common Chemistry](), and [IUPAC OPSIN]().

TODO: explain [this link](https://pubchem.ncbi.nlm.nih.gov/classification/#hid=72)

# Other tools
  -[PubChemPy]()
  -[CIRpy]()
  -[IUPAC webbook]()
  -more?

# TODO list
We will provide Python interfaces to more online services, including:  
  -[ChEBI](https://www.ebi.ac.uk/chebi/beta/tools) as soon as the REST API documentations become available. The page is empty now!
  - [ZeroPM](https://database.zeropm.eu/) even though there is no web API, the data is available on GitHub. I have written an interface that is not shared here since it can make this codebase too large, and I aim to keep it lean. We will find a way to share it.
  - More? Please open an issue in this repo.

  ## Documentation deployment

  This repository includes a GitHub Actions workflow at `.github/workflows/mkdocs-deploy.yml` which builds the MkDocs documentation and deploys the generated site to the `gh-pages` branch using `peaceiris/actions-gh-pages`.

  - The workflow triggers on pushes to `main` and can be run manually from the Actions tab.
  - Per your request, tests are disabled in this workflow; it only installs documentation dependencies, builds the site, and publishes the `site/` directory.
