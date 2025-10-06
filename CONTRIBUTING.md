# Contributing to EGFR Allosteric Drug Discovery Pipeline

Thank you for your interest in contributing! This document provides guidelines for contributing to this project.

## Code of Conduct

- Be respectful and inclusive
- Focus on constructive feedback
- Help maintain a welcoming environment for all contributors

## How to Contribute

### Reporting Issues

If you find a bug or have a feature request:

1. Check if the issue already exists in the GitHub Issues
2. If not, create a new issue with:
   - Clear, descriptive title
   - Detailed description of the problem/feature
   - Steps to reproduce (for bugs)
   - Expected vs actual behavior
   - Your environment (Python version, OS, etc.)

### Submitting Pull Requests

1. **Fork the repository** and create a new branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following the code style guidelines below

3. **Test your changes**:
   ```bash
   # Run the pipeline on test data
   python scripts/00_prepare_inputs_from_mmcif.py --help
   # Add tests if applicable
   pytest tests/
   ```

4. **Commit your changes** with clear, descriptive messages:
   ```bash
   git commit -m "Add feature: brief description"
   ```

5. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

6. **Create a Pull Request** with:
   - Clear description of changes
   - Reference to related issues
   - Screenshots/examples if applicable

## Code Style Guidelines

### Python Code

- Follow **PEP 8** style guide
- Use **type hints** for function parameters and return values
- Maximum line length: **100 characters**
- Use **docstrings** for all functions, classes, and modules

Example:
```python
def process_molecule(mol_dir: Path, output_base: Path, source: str) -> int:
    """
    Process a single molecule directory.

    Args:
        mol_dir: Path to molecule directory
        output_base: Base output directory
        source: Source type (cmolgpt or enamine_real)

    Returns:
        Number of successfully processed samples
    """
    pass
```

### Documentation

- Use **Markdown** for documentation files
- Keep README.md concise and up-to-date
- Document all parameters and outputs
- Include examples where helpful

### Commit Messages

Follow the conventional commits format:

- `feat:` for new features
- `fix:` for bug fixes
- `docs:` for documentation changes
- `refactor:` for code refactoring
- `test:` for adding tests
- `chore:` for maintenance tasks

Examples:
```
feat: add support for additional file formats
fix: correct IC50 normalization for μM units
docs: update installation instructions
refactor: simplify pose selection logic
```

## Development Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/egfr_allosteric.git
   cd egfr_allosteric
   ```

2. **Create a virtual environment**:
   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies** (including dev dependencies):
   ```bash
   pip install -r requirements.txt
   pip install -r requirements-dev.txt  # If available
   ```

4. **Install pre-commit hooks** (optional but recommended):
   ```bash
   pip install pre-commit
   pre-commit install
   ```

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_prepare_inputs.py

# Run with coverage
pytest --cov=scripts tests/
```

### Writing Tests

- Place tests in the `tests/` directory
- Name test files `test_*.py`
- Use descriptive test function names: `test_convert_mmcif_to_pdb_success()`
- Include both positive and negative test cases

Example:
```python
def test_normalize_ic50_nM():
    """Test IC50 normalization for nanomolar units."""
    result = normalize_ic50(100.0, 'nM')
    assert result == 1e-7

def test_normalize_ic50_invalid():
    """Test IC50 normalization with invalid input."""
    result = normalize_ic50(-10.0, 'nM')
    assert result is None
```

## Project Structure

```
egfr_allosteric/
├── scripts/              # Pipeline scripts (00-05)
├── tests/                # Unit tests
├── data/                 # Raw data (not tracked in git)
├── data_std/             # Standardized data (not tracked)
├── outputs/              # Results (not tracked)
├── analysis/             # Analysis and visualization
├── docs/                 # Documentation
├── requirements.txt      # Python dependencies
├── .gitignore           # Git ignore rules
├── LICENSE              # MIT License
└── README.md            # Main documentation
```

## Areas for Contribution

We welcome contributions in the following areas:

### High Priority
- [ ] Add comprehensive unit tests
- [ ] Improve error handling and validation
- [ ] Add support for additional structure formats
- [ ] Optimize performance for large datasets
- [ ] Improve documentation and examples

### Medium Priority
- [ ] Add visualization tools for results
- [ ] Implement alternative scoring functions
- [ ] Add support for parallel processing
- [ ] Create Docker container for reproducibility
- [ ] Add example datasets

### Low Priority
- [ ] Web interface for result exploration
- [ ] Integration with chemical databases
- [ ] Machine learning-based scoring
- [ ] Interactive Jupyter notebooks

## Questions?

If you have questions about contributing:
- Open a GitHub Discussion
- Email the maintainers
- Check existing documentation in `docs/`

## Attribution

Contributors will be acknowledged in:
- README.md contributors section
- CHANGELOG.md
- GitHub contributors page

Thank you for contributing to scientific research!
