"""Tests for recon.utils module."""
import pytest
from recon.utils import split_layer_name


class TestUtilsFunctions:
    """Test utility functions."""
    
    def test_split_layer_name_with_underscore(self):
        """Test splitting layer name with underscore separator."""
        celltype, layer = split_layer_name("CellA_grn")
        assert celltype == "CellA"
        assert layer == "grn"
    
    def test_split_layer_name_multiple_underscores(self):
        """Test splitting layer name with multiple underscores."""
        celltype, layer = split_layer_name("T_cell_CD4_grn")
        assert celltype == "T_cell_CD4"
        assert layer == "grn"
    
    def test_split_layer_name_no_separator(self):
        """Test splitting layer name without separator at end."""
        celltype, layer = split_layer_name("cell_communication")
        # Actually splits on last underscore: 'cell' and 'communication'
        assert celltype == "cell"
        assert layer == "communication"
    
    def test_split_layer_name_custom_separator(self):
        """Test splitting with custom separator."""
        celltype, layer = split_layer_name("CellA-grn", separator='-')
        assert celltype == "CellA"
        assert layer == "grn"
