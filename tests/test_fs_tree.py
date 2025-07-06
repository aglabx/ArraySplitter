#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for the frequency suffix tree module.
"""

import pytest
from ArraySplitter.core_functions.tools.fs_tree import (
    WeightedValueHeap,
    update,
    iter_fs_tree_from_sequence,
    build_fs_tree_from_sequence,
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


class TestFSTreeFunctions:
    """Test frequency suffix tree construction functions."""
    
    def test_get_starting_index_simple(self):
        """Test finding starting indices for a nucleotide."""
        sequence = "ATCGATCG"
        
        # Find all 'A' positions
        indices = get_starting_index(sequence, 'A')
        assert indices == [0, 4]
        
        # Find all 'T' positions
        indices = get_starting_index(sequence, 'T')
        assert indices == [1, 5]
        
        # Find all 'G' positions
        indices = get_starting_index(sequence, 'G')
        assert indices == [3, 7]
        
        # Non-existent nucleotide
        indices = get_starting_index(sequence, 'N')
        assert indices == []
    
    def test_get_starting_index_repeated(self):
        """Test with repeated nucleotides."""
        sequence = "AAATTTCCC"
        
        indices = get_starting_index(sequence, 'A')
        assert indices == [0, 1, 2]
        
        indices = get_starting_index(sequence, 'T')
        assert indices == [3, 4, 5]
    
    @pytest.mark.parametrize("sequence,nucl,expected_count", [
        ("ATCG", "A", 1),
        ("AAAA", "A", 4),
        ("CGCG", "T", 0),
        ("", "A", 0),
    ])
    def test_get_starting_index_parametrized(self, sequence, nucl, expected_count):
        """Parametrized test for various sequences."""
        indices = get_starting_index(sequence, nucl)
        assert len(indices) == expected_count
    
    def test_get_nodes(self):
        """Test node extraction from fs_tree."""
        # Create a simple fs_tree structure
        fs_tree = {
            0: (0, ['A'], [0, 4], [1, 5], None, [1, 2]),
            1: (1, ['A', 'T'], [0], [2], 0, []),
            2: (2, ['A', 'C'], [4], [6], 0, []),
        }
        
        nodes = get_nodes(fs_tree, 0)
        assert len(nodes) == 3
        assert 0 in nodes
        assert 1 in nodes
        assert 2 in nodes
    
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


class TestFSTreeConstruction:
    """Test the main fs_tree construction function."""
    
    def test_simple_repeat_pattern(self):
        """Test fs_tree construction on a simple repeat."""
        sequence = "ATCATCATC"
        
        # Build fs_tree starting from 'A'
        fs_tree = get_fs_tree(sequence, 'A', cutoff=1, depth=3)
        
        # Should find the ATC pattern
        assert len(fs_tree) > 1
        
        # Check that starting positions are correct
        root_node = fs_tree[0]
        assert root_node[2] == [0, 3, 6]  # Positions of 'A'
    
    def test_no_pattern(self):
        """Test with sequence without repeating patterns."""
        sequence = "ACGTACGT"
        
        fs_tree = get_fs_tree(sequence, 'A', cutoff=3, depth=3)
        
        # With high cutoff, should only have root
        assert len(fs_tree) == 1


@pytest.fixture
def sample_tandem_repeat():
    """Provide a sample tandem repeat sequence."""
    return "CAGCAGCAGCAGCAG"  # CAG repeat


@pytest.fixture
def complex_satellite():
    """Provide a more complex satellite array."""
    monomer = "ATCGATCGATCG"
    return monomer * 10  # 10 copies


def test_iterate_hints(sample_tandem_repeat):
    """Test hint generation from fs_tree."""
    fs_tree = get_fs_tree(sample_tandem_repeat, 'C', cutoff=2, depth=5)
    hints = list(iterate_hints(fs_tree, 0))
    
    # Should find CAG pattern
    assert len(hints) > 0
    
    # Check if CAG is among the hints
    hint_sequences = [hint[1] for hint in hints]
    assert any('CAG' in seq for seq in hint_sequences)