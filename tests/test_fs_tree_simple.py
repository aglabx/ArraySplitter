#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simplified tests for the frequency suffix tree module.
Tests only existing functions.
"""

import pytest
from ArraySplitter.core_functions.tools.fs_tree import (
    WeightedValueHeap,
    update,
    iter_fs_tree_from_sequence,
)
from ArraySplitter.decompose import (
    get_fs_tree,
    iterate_hints,
)


class TestWeightedValueHeap:
    """Test the WeightedValueHeap data structure."""
    
    def test_init(self):
        """Test heap initialization."""
        heap = WeightedValueHeap()
        assert heap.heap == []
        assert heap.peek() is None
        assert heap.pop() is None
    
    def test_insert_and_pop(self):
        """Test insertion and popping maintains heap order."""
        heap = WeightedValueHeap()
        
        # Insert values with weights
        heap.insert(5, "five")
        heap.insert(1, "one")
        heap.insert(3, "three")
        heap.insert(2, "two")
        
        # Pop should return in order of weights (min first)
        assert heap.pop() == (1, "one")
        assert heap.pop() == (2, "two")
        assert heap.pop() == (3, "three")
        assert heap.pop() == (5, "five")
        assert heap.pop() is None
    
    def test_peek(self):
        """Test peeking without removing."""
        heap = WeightedValueHeap()
        heap.insert(3, "three")
        heap.insert(1, "one")
        
        # Peek should return min without removing
        assert heap.peek() == (1, "one")
        assert heap.peek() == (1, "one")  # Still there
        assert len(heap.heap) == 2


class TestFSTreeConstruction:
    """Test the main fs_tree construction function."""
    
    def test_simple_repeat_pattern(self):
        """Test fs_tree construction on a simple repeat."""
        sequence = "ATCATCATC"
        
        # Build fs_tree starting from 'A'
        fs_tree = get_fs_tree(sequence, 'A', cutoff=1, depth=3)
        
        # Should find the ATC pattern
        assert len(fs_tree) > 1
        
        # Check that root node exists
        assert 0 in fs_tree
    
    def test_no_pattern(self):
        """Test with sequence without repeating patterns."""
        sequence = "ACGTACGT"
        
        fs_tree = get_fs_tree(sequence, 'A', cutoff=3, depth=3)
        
        # With high cutoff, should only have root
        assert len(fs_tree) == 1
    
    def test_iterate_hints_basic(self):
        """Test hint generation from fs_tree."""
        sequence = "CAGCAGCAGCAGCAG"
        fs_tree = get_fs_tree(sequence, 'C', cutoff=2, depth=5)
        hints = list(iterate_hints(sequence, fs_tree, depth=5))
        
        # Should find hints
        assert len(hints) > 0
        
        # Check hint structure
        for hint in hints:
            assert len(hint) == 3  # (length, sequence, frequency)
            L, seq, N = hint
            assert isinstance(L, int)
            assert isinstance(seq, str)
            assert isinstance(N, int)


class TestUpdate:
    """Test the update function."""
    
    def test_update_function(self):
        """Test the update function for fs_tree construction."""
        fs_tree = {0: (0, [], [0], [0], None, [])}
        queue = []
        
        # Test adding a new node
        cid = update(
            fs_tree=fs_tree,
            queue=queue,
            nucl='A',
            fs_x=[0, 1, 2, 3],  # 4 occurrences
            fs_xp=[1, 2, 3, 4],
            current_cid=0,
            cid=1,
            cutoff=3  # Will pass since len(fs_x) > cutoff
        )
        
        assert cid == 2  # Should increment
        assert 1 in fs_tree
        assert fs_tree[1][1] == ['A']
        assert fs_tree[1][2] == [0, 1, 2, 3]
        assert len(queue) == 1
        assert queue[0] == 1
    
    def test_update_below_cutoff(self):
        """Test update when frequency is below cutoff."""
        fs_tree = {0: (0, [], [0], [0], None, [])}
        queue = []
        
        cid = update(
            fs_tree=fs_tree,
            queue=queue,
            nucl='A',
            fs_x=[0, 1],  # Only 2 occurrences
            fs_xp=[1, 2],
            current_cid=0,
            cid=1,
            cutoff=3  # Won't pass since len(fs_x) <= cutoff
        )
        
        assert cid == 1  # Should not increment
        assert 1 not in fs_tree
        assert len(queue) == 0


@pytest.fixture
def sample_tandem_repeat():
    """Provide a sample tandem repeat sequence."""
    return "CAGCAGCAGCAGCAG"  # CAG repeat


def test_fs_tree_with_tandem(sample_tandem_repeat):
    """Test fs_tree with a known tandem repeat."""
    fs_tree = get_fs_tree(sample_tandem_repeat, 'C', cutoff=2, depth=5)
    
    # Should build a tree
    assert len(fs_tree) >= 1
    
    # Generate hints
    hints = list(iterate_hints(sample_tandem_repeat, fs_tree, depth=5))
    
    # Should find patterns
    assert len(hints) > 0
    
    # Look for the CAG pattern in hints
    hint_sequences = [hint[1] for hint in hints]
    # At least one hint should contain 'CAG' or be 'CAG'
    assert any('CAG' in seq or seq == 'CAG' for seq in hint_sequences)