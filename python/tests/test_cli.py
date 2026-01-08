import pytest
from click.testing import CliRunner
from linkapy.CLI import linkapy, parsing, example


class TestCLI:
    """Test basic CLI functionality."""
    
    def test_linkapy_main_empty(self):
        """Run linkapy without any args"""
        cli = CliRunner()
        clires = cli.invoke(linkapy, [])
        assert clires.exit_code == 0

    def test_linkapy_main(self):
        """Test that the main linkapy group command shows help."""
        cli = CliRunner()
        clires = cli.invoke(linkapy, ['--help'])
        assert clires.exit_code == 0
    
    def test_linkapy_version(self):
        """Test that version flag works."""
        cli = CliRunner()
        clires = cli.invoke(linkapy, ['--version'])
        assert clires.exit_code == 0
        assert 'linkapy, version ' in clires.output

    def test_parsing_nopaths(self):
        """Test parsing command without required paths."""
        cli = CliRunner()
        clires = cli.invoke(parsing, [])
        assert clires.exit_code == 0
        assert "Provide either a methylation path and/or a transcriptome path." in clires.output
    
    def test_parsing_nochroms_noregions(self):
        """Test parsing command with methylation path but no chromsizes or regions."""
        cli = CliRunner()
        clires = cli.invoke(parsing, ['--methylation_path', './'])
        print(clires.output)
        assert clires.exit_code == 0
        assert "Methylation data requires either a chromsizes file or at least one regions file." in clires.output
    
    def test_parsing_nomemode(self, nome_path):
        """Test parsing with --nome flag (sets patterns and names to empty tuples)."""
        cli = CliRunner()
        clires = cli.invoke(parsing, [
            '--methylation_path', str(nome_path),
            '--chromsizes', str(nome_path / 'chromsizes.txt'),
            '--NOMe',
            '--cli_test'
        ])
        assert clires.exit_code == 0
        print(clires)
        print(clires.output)