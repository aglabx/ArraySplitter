#!/bin/bash
# Setup script for ArraySplitter testing environment

echo "Setting up ArraySplitter testing environment..."

# Check if we're in the correct directory
if [ ! -f "setup.py" ]; then
    echo "Error: setup.py not found. Please run this script from the ArraySplitter root directory."
    exit 1
fi

# Install ArraySplitter in development mode
echo "Installing ArraySplitter in development mode..."
pip install -e .

# Install testing dependencies
echo "Installing testing dependencies..."
if [ -f "requirements-dev.txt" ]; then
    pip install -r requirements-dev.txt
else
    echo "requirements-dev.txt not found, installing minimal test dependencies..."
    pip install pytest pytest-cov
fi

# Verify installation
echo "Verifying installation..."
python -c "import ArraySplitter; print('✓ ArraySplitter imported successfully')" || {
    echo "✗ Failed to import ArraySplitter"
    exit 1
}

pytest --version || {
    echo "✗ pytest not installed properly"
    exit 1
}

echo "✓ Testing environment setup complete!"
echo ""
echo "You can now run tests with:"
echo "  pytest                    # Run all tests"
echo "  pytest -v                 # Run with verbose output"
echo "  pytest tests/test_fs_tree.py  # Run specific test file"
echo ""