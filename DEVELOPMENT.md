# Development Guide

This document provides a comprehensive guide for developers working on PROVESID, including the release workflow, development setup, and contribution guidelines.

## 📋 Table of Contents

1. [Development Setup](#development-setup)
2. [Release Workflow](#release-workflow)
3. [Testing](#testing)
4. [Documentation](#documentation)
5. [Code Quality](#code-quality)
6. [Contributing](#contributing)

## 🛠️ Development Setup

### Prerequisites

- Python 3.8 or higher
- Git
- A GitHub account with access to the USEtox/PROVESID repository

### Initial Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/USEtox/PROVESID.git
   cd PROVESID
   ```

2. **Create a virtual environment:**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install development dependencies:**
   ```bash
   pip install -e ".[dev,test]"
   ```

4. **Verify installation:**
   ```bash
   python -c "import provesid; print('PROVESID installed successfully!')"
   pytest tests/ -v
   ```

## 🚀 Release Workflow

### Step-by-Step Release Process

#### 1. Pre-Release Preparation

**a) Update Version Number**
- Edit `pyproject.toml` and update the version number:
  ```toml
  version = "X.Y.Z"  # e.g., "0.2.0"
  ```

**b) Update Documentation**
- Ensure all new features are documented
- Update `CHANGELOG.md` (if exists) or create release notes
- Verify all examples and tutorials work with the new version

**c) Run Complete Test Suite**
```bash
# Run all tests
pytest tests/ -v

# Run tests with coverage
pytest tests/ --cov=provesid --cov-report=html

# Test package building
python -m build
```

**d) Update Dependencies (if needed)**
- Review and update dependencies in `pyproject.toml`
- Test with updated dependencies

#### 2. Create Release Commit

```bash
# Add all changes
git add .

# Commit with descriptive message
git commit -m "Prepare for vX.Y.Z release: [brief description of changes]"

# Push to main branch
git push origin main
```

#### 3. Create and Push Release Tag

```bash
# Create annotated tag
git tag -a vX.Y.Z -m "Release version X.Y.Z - [brief description]"

# Push the tag to trigger release workflow
git push origin vX.Y.Z
```

#### 4. Monitor Automated Release Process

The GitHub Actions workflow will automatically:
1. **Build the package** (wheel and source distribution)
2. **Run tests** to ensure quality
3. **Create GitHub release** with artifacts
4. **Publish to PyPI** (if configured)

**Check progress at:**
- Actions: https://github.com/USEtox/PROVESID/actions
- Releases: https://github.com/USEtox/PROVESID/releases

#### 5. Post-Release Tasks

**a) Verify Release**
- Check GitHub release page
- Verify PyPI package: https://pypi.org/project/provesid/
- Test installation: `pip install provesid==X.Y.Z`

**b) Update Documentation**
- Ensure documentation reflects new version
- Update installation instructions if needed

**c) Communicate Release**
- Update project README if needed
- Announce in relevant channels/communities

### PyPI Configuration (First-Time Setup)

#### Option 1: Trusted Publishing (Recommended)

1. Go to https://pypi.org/manage/account/publishing/
2. Add a new publisher with:
   - **PyPI Project Name:** `provesid`
   - **Owner:** `USEtox`
   - **Repository:** `PROVESID`
   - **Workflow:** `release.yml`
   - **Environment:** `pypi`

#### Option 2: API Token

1. Create PyPI API token at https://pypi.org/manage/account/token/
2. Add to GitHub repository secrets as `PYPI_API_TOKEN`
3. Go to: Repository Settings → Secrets and variables → Actions

### Emergency Procedures

#### Rollback a Release

```bash
# Delete local tag
git tag -d vX.Y.Z

# Delete remote tag
git push origin :refs/tags/vX.Y.Z

# Create new corrected release
git tag -a vX.Y.Z-fix -m "Fix for version X.Y.Z"
git push origin vX.Y.Z-fix
```

#### Manual PyPI Upload

If automated publishing fails:

```bash
# Install publishing tools
pip install build twine

# Build package
python -m build

# Upload to PyPI
twine upload dist/*
```

## 🧪 Testing

### Running Tests

```bash
# Run all tests
pytest tests/

# Run with verbose output
pytest tests/ -v

# Run specific test file
pytest tests/test_pubchem.py -v

# Run tests with coverage
pytest tests/ --cov=provesid --cov-report=html

# Run tests in parallel
pytest tests/ -n auto
```

### Test Structure

```
tests/
├── conftest.py              # Test configuration and fixtures
├── test_pubchem.py         # PubChem API tests
├── test_pubchemview.py     # PubChem View tests
├── test_cascommonchem.py   # CAS Common Chemistry tests
├── test_classyfire.py      # ClassyFire API tests
├── test_nci_resolver.py    # NCI Resolver tests
├── test_opsin.py          # OPSIN tests
└── test_chebi.py          # ChEBI tests
```

### Adding New Tests

1. Create test file following naming convention: `test_<module>.py`
2. Use pytest fixtures for common setup
3. Include both unit tests and integration tests
4. Test error handling and edge cases
5. Add tests for any new features or bug fixes

## 📚 Documentation

### Building Documentation Locally

```bash
# Install documentation dependencies
pip install -e ".[dev]"

# Serve documentation locally
mkdocs serve

# Build documentation
mkdocs build
```

### Documentation Structure

```
docs/
├── index.md                # Main documentation page
├── installation.md         # Installation instructions
├── quickstart.md          # Quick start guide
├── api/                   # API documentation
│   ├── index.md
│   ├── pubchem.md
│   ├── pubchemview.md
│   └── ...
├── examples/              # Jupyter notebook tutorials
│   ├── pubchem/
│   ├── resolver/
│   └── ...
└── stylesheets/           # Custom CSS
```

### Adding Documentation

1. **API Documentation:** Auto-generated from docstrings using mkdocstrings
2. **Tutorials:** Add Jupyter notebooks to `examples/` and `docs/examples/`
3. **Guides:** Add Markdown files to `docs/`

## 🔍 Code Quality

### Code Style

We use the following tools for code quality:

```bash
# Format code with black
black src/ tests/

# Lint with flake8
flake8 src/ tests/

# Type checking (if using mypy)
mypy src/provesid/
```

### Pre-commit Hooks (Optional)

```bash
# Install pre-commit
pip install pre-commit

# Install hooks
pre-commit install

# Run manually
pre-commit run --all-files
```

## 🤝 Contributing

### Workflow for Contributors

1. **Fork the repository**
2. **Create a feature branch:**
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. **Make changes and add tests**
4. **Run tests and ensure they pass**
5. **Update documentation if needed**
6. **Submit a pull request**

### Pull Request Guidelines

- Include tests for new features
- Update documentation as needed
- Follow existing code style
- Write clear commit messages
- Reference any related issues

### Issue Reporting

When reporting issues, please include:
- Python version
- PROVESID version
- Operating system
- Minimal code example
- Full error traceback

## 📁 Project Structure

```
PROVESID/
├── .github/
│   └── workflows/          # GitHub Actions workflows
├── docs/                   # Documentation source
├── examples/              # Example scripts and notebooks
├── src/
│   └── provesid/          # Main package code
│       ├── __init__.py
│       ├── pubchem.py
│       ├── pubchemview.py
│       ├── resolver.py
│       ├── cascommonchem.py
│       ├── classyfire.py
│       ├── opsin.py
│       ├── chebi.py
│       ├── utils.py
│       └── data/          # Package data files
├── tests/                 # Test suite
├── pyproject.toml        # Package configuration
├── README.md             # Project overview
├── LICENSE               # License file
├── MANIFEST.in           # Package manifest
└── mkdocs.yml           # Documentation configuration
```

## 🔧 Troubleshooting

### Common Issues

**1. Import Errors:**
- Ensure package is installed in development mode: `pip install -e .`
- Check Python path and virtual environment

**2. Test Failures:**
- Check API availability (some tests depend on external services)
- Verify network connectivity
- Check rate limiting

**3. Documentation Build Errors:**
- Ensure all dependencies are installed: `pip install -e ".[dev]"`
- Check for syntax errors in Markdown files

**4. Release Workflow Failures:**
- Verify GitHub Actions permissions
- Check PyPI credentials configuration
- Review workflow logs in GitHub Actions

## 📞 Getting Help

- **Issues:** https://github.com/USEtox/PROVESID/issues
- **Discussions:** https://github.com/USEtox/PROVESID/discussions
- **Documentation:** https://usetox.github.io/PROVESID/

---

For questions not covered in this guide, please open an issue or start a discussion on GitHub.
