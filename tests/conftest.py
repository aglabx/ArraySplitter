#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Pytest configuration and shared fixtures for ArraySplitter tests.
"""

import pytest
import os
import tempfile
from pathlib import Path


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_fasta_file(temp_dir):
    """Create a sample FASTA file for testing."""
    fasta_content = """>array1
ATCGATCGATCGATCGATCG
>array2
CAGCAGCAGCAGCAGCAGCAG
>array3
ATTCCATTCCATTCCATTCC
"""
    fasta_path = temp_dir / "test_arrays.fa"
    fasta_path.write_text(fasta_content)
    return fasta_path


@pytest.fixture
def sample_tandem_repeats():
    """Provide various tandem repeat patterns for testing."""
    return {
        "simple_dinucleotide": "AT" * 20,
        "simple_trinucleotide": "CAG" * 15,
        "alpha_satellite": "ATTCCATTCCATTCCATTCC" * 5,
        "complex_repeat": "ATCGATCG" * 10,
        "with_variation": "CAGCAGCATCAGCAGCAG",  # CAG repeat with one CAT
    }


@pytest.fixture
def sample_monomers():
    """Provide extracted monomers for testing."""
    return [
        "ATCGATCG",
        "ATCGATCG",
        "ATCGATCG",
        "ATCGATCC",  # One with variation
        "ATCGATCG",
    ]


@pytest.fixture
def performance_test_array():
    """Large array for performance testing."""
    # 100kb array
    monomer = "ATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCC"
    return monomer * 1300  # ~100kb


def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "unit: marks tests as unit tests"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add markers."""
    for item in items:
        # Add unit marker to all tests in test files starting with 'test_'
        if item.fspath.basename.startswith("test_"):
            if "integration" not in item.keywords:
                item.add_marker(pytest.mark.unit)
        
        # Mark tests with 'performance' in name as slow
        if "performance" in item.nodeid:
            item.add_marker(pytest.mark.slow)


@pytest.fixture
def mock_progress_bar(monkeypatch):
    """Mock tqdm progress bar for cleaner test output."""
    class MockTqdm:
        def __init__(self, *args, **kwargs):
            self.iterable = args[0] if args else None
        
        def __iter__(self):
            return iter(self.iterable) if self.iterable else iter([])
        
        def __enter__(self):
            return self
        
        def __exit__(self, *args):
            pass
        
        def update(self, n=1):
            pass
        
        def close(self):
            pass
    
    monkeypatch.setattr("tqdm.tqdm", MockTqdm)
    return MockTqdm


# Helper functions for tests

def create_fasta_content(sequences):
    """Helper to create FASTA formatted content."""
    content = []
    for i, seq in enumerate(sequences):
        content.append(f">sequence_{i}")
        content.append(seq)
    return "\n".join(content)


def assert_monomers_similar(monomers, expected_length, tolerance=0.1):
    """Assert that monomers have similar lengths."""
    lengths = [len(m) for m in monomers]
    mean_length = sum(lengths) / len(lengths)
    
    for length in lengths:
        deviation = abs(length - expected_length) / expected_length
        assert deviation <= tolerance, f"Monomer length {length} deviates too much from expected {expected_length}"