Contributing to PyMEGA
====================

Thank you for your interest in contributing to PyMEGA! This guide will help you get started with contributing code, documentation, or other improvements to the project.

Getting Started
---------------

Development Environment Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Fork the Repository**

   Fork the PyMEGA repository on GitHub and clone your fork:

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/pymega.git
      cd pymega

2. **Create a Virtual Environment**

   .. code-block:: bash

      python -m venv pymega_dev
      source pymega_dev/bin/activate  # On Windows: pymega_dev\Scripts\activate

3. **Install Development Dependencies**

   .. code-block:: bash

      pip install -e .[dev]

   This installs PyMEGA in editable mode with all development dependencies including:
   - Testing tools (pytest, pytest-cov)
   - Code formatting (black, flake8)
   - Type checking (mypy)
   - Documentation tools (sphinx)

4. **Install Pre-commit Hooks**

   .. code-block:: bash

      pre-commit install

   This ensures code quality checks run automatically before each commit.

Code Style and Standards
------------------------

Code Formatting
~~~~~~~~~~~~~~~

PyMEGA uses several tools to maintain code quality:

**Black** for code formatting:

.. code-block:: bash

   black python_scmega/
   
**Flake8** for linting:

.. code-block:: bash

   flake8 python_scmega/
   
**MyPy** for type checking:

.. code-block:: bash

   mypy python_scmega/

**Import Sorting** with isort:

.. code-block:: bash

   isort python_scmega/

All these checks run automatically with pre-commit hooks, but you can run them manually as well.

Code Style Guidelines
~~~~~~~~~~~~~~~~~~~~

1. **Follow PEP 8** for Python code style
2. **Use descriptive variable names** that clearly indicate purpose
3. **Add type hints** for function parameters and return values
4. **Keep functions focused** and do one thing well
5. **Use docstrings** for all public functions and classes

Example of well-formatted code:

.. code-block:: python

   def calculate_correlation_matrix(
       data: np.ndarray,
       method: str = "pearson",
       min_observations: int = 10
   ) -> np.ndarray:
       """
       Calculate correlation matrix for the input data.
       
       Args:
           data: Input data matrix (n_samples, n_features)
           method: Correlation method ("pearson", "spearman")
           min_observations: Minimum number of observations required
           
       Returns:
           Correlation matrix (n_features, n_features)
           
       Raises:
           ValueError: If method is not supported
           
       Example:
           >>> data = np.random.randn(100, 10)
           >>> corr_matrix = calculate_correlation_matrix(data)
           >>> corr_matrix.shape
           (10, 10)
       """
       if method not in ["pearson", "spearman"]:
           raise ValueError(f"Unsupported method: {method}")
       
       # Implementation here
       pass

Documentation Standards
~~~~~~~~~~~~~~~~~~~~~~~

1. **Use Google-style docstrings** for consistency
2. **Include examples** in docstrings when helpful
3. **Document all parameters** and return values
4. **Add type information** in docstrings
5. **Keep docstrings concise** but informative

Testing Guidelines
------------------

Writing Tests
~~~~~~~~~~~~~

PyMEGA uses pytest for testing. All new code should include tests:

1. **Unit tests** for individual functions
2. **Integration tests** for workflows
3. **Property-based tests** for complex algorithms

Test Organization
~~~~~~~~~~~~~~~~~

Tests are organized in the ``tests/`` directory:

.. code-block:: text

   tests/
   ├── test_data_processing.py
   ├── test_cell_annotation.py
   ├── test_network_inference.py
   └── fixtures/
       ├── conftest.py
       └── sample_data.py

Example Test Structure
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pytest
   import numpy as np
   from pymega import MultiomeData, quality_control_multiome
   
   
   class TestQualityControl:
       """Test quality control functions."""
       
       def setup_method(self):
           """Set up test data before each test method."""
           self.test_data = self._create_test_multiome()
       
       def test_basic_qc_filtering(self):
           """Test basic quality control filtering."""
           filtered_data = quality_control_multiome(
               self.test_data,
               min_genes=10,
               min_cells=5
           )
           
           assert filtered_data.n_obs <= self.test_data.n_obs
           assert filtered_data.n_vars <= self.test_data.n_vars
       
       def test_qc_parameter_validation(self):
           """Test parameter validation in QC functions."""
           with pytest.raises(ValueError):
               quality_control_multiome(self.test_data, min_genes=-1)
       
       def _create_test_multiome(self) -> MultiomeData:
           """Create test multiome data."""
           # Implementation here
           pass

Running Tests
~~~~~~~~~~~~~

Run all tests:

.. code-block:: bash

   pytest

Run specific test file:

.. code-block:: bash

   pytest tests/test_data_processing.py

Run with coverage:

.. code-block:: bash

   pytest --cov=python_scmega --cov-report=html

Contribution Workflow
---------------------

Making Changes
~~~~~~~~~~~~~~

1. **Create a Feature Branch**

   .. code-block:: bash

      git checkout -b feature/your-feature-name

2. **Make Your Changes**

   - Write code following the style guidelines
   - Add tests for new functionality
   - Update documentation as needed
   - Ensure all tests pass

3. **Commit Your Changes**

   .. code-block:: bash

      git add .
      git commit -m "Add feature: brief description"

   Use clear, descriptive commit messages following this format:
   
   - ``Add feature: description`` for new features
   - ``Fix bug: description`` for bug fixes
   - ``Update docs: description`` for documentation
   - ``Optimize: description`` for performance improvements

4. **Push to Your Fork**

   .. code-block:: bash

      git push origin feature/your-feature-name

5. **Create a Pull Request**

   - Go to your fork on GitHub
   - Click "New Pull Request"
   - Provide a clear title and description
   - Link to any related issues

Pull Request Guidelines
~~~~~~~~~~~~~~~~~~~~~~

**Before Submitting:**
- [ ] Tests pass locally
- [ ] Code follows style guidelines
- [ ] Documentation is updated
- [ ] Changes are well-tested
- [ ] Performance impact is considered

**PR Description Should Include:**
- Clear description of changes
- Motivation for the changes
- Any breaking changes
- Testing strategy
- Performance impact (if any)

**PR Title Format:**
- ``Add: description`` for new features
- ``Fix: description`` for bug fixes
- ``Docs: description`` for documentation
- ``Perf: description`` for performance improvements

Types of Contributions
----------------------

Code Contributions
~~~~~~~~~~~~~~~~~~

**New Features:**
- Single-cell analysis algorithms
- Visualization functions
- Performance optimizations
- Integration with other tools

**Bug Fixes:**
- Algorithmic corrections
- Performance improvements
- Edge case handling
- Memory optimization

**Code Improvements:**
- Algorithm efficiency
- Memory optimization
- Code quality enhancements

Documentation Contributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**User Documentation:**
- Tutorials and examples
- API documentation
- Installation guides
- Troubleshooting guides

**Developer Documentation:**
- Code architecture
- Algorithm descriptions
- Performance guidelines
- Contributing guides

Testing Contributions
~~~~~~~~~~~~~~~~~~~~

**Test Coverage:**
- Unit tests for uncovered code
- Integration tests for workflows
- Edge case testing

**Test Infrastructure:**
- Testing utilities
- Mock data generation
- Continuous integration

Community Contributions
~~~~~~~~~~~~~~~~~~~~~~~

**Issue Reporting:**
- Bug reports with reproduction steps
- Feature requests with use cases
- Documentation improvements

**Code Review:**
- Review pull requests
- Provide constructive feedback
- Test proposed changes
- Suggest improvements

**Community Support:**
- Answer questions in discussions
- Help with troubleshooting
- Provide examples and tutorials
- Mentor new contributors

Development Practices
---------------------

Version Control
~~~~~~~~~~~~~~~

- Use descriptive commit messages
- Keep commits focused and atomic
- Rebase feature branches before merging
- Use conventional commit format

Code Organization
~~~~~~~~~~~~~~~~

- Keep modules focused and cohesive
- Use clear, descriptive names
- Minimize dependencies between modules
- Follow established project structure

Performance Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~

- Profile code before optimizing
- Consider memory usage in large datasets
- Use appropriate data structures
- Document performance characteristics

Documentation
~~~~~~~~~~~~~

- Keep documentation up to date
- Include examples in docstrings
- Update API documentation for changes
- Test documentation examples

Getting Help
------------

**Development Questions:**
- Open a discussion on GitHub
- Ask in the development channel
- Review existing documentation
- Check similar projects

**Technical Issues:**
- Search existing issues
- Provide minimal reproduction
- Include system information
- Test with latest version

**Design Decisions:**
- Discuss in GitHub issues
- Propose clear alternatives
- Consider backwards compatibility
- Think about user impact

Recognition
-----------

All contributors will be recognized in:

- ``CONTRIBUTORS.md`` file
- Release notes for significant contributions
- Documentation acknowledgments
- Project README

Thank you for contributing to PyMEGA! Your efforts help make single-cell multiome analysis more accessible and powerful for the scientific community.
