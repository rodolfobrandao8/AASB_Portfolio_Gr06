import unittest
from bioinf import sequencias

class TestValidacaoDNA(unittest.TestCase):
    def test_dna_valido(self):
        self.assertTrue(sequencias.validar_dna("ACGT"))
        self.assertTrue(sequencias.validar_dna("acgt"))  # minúsculas

    def test_dna_invalido(self):
        self.assertFalse(sequencias.validar_dna("ACGTX"))
        self.assertFalse(sequencias.validar_dna(""))  # vazia

class TestValidacaoRNA(unittest.TestCase):
    def test_rna_valido(self):
        self.assertTrue(sequencias.validar_rna("ACGU"))
        self.assertTrue(sequencias.validar_rna("acgu"))

    def test_rna_invalido(self):
        self.assertFalse(sequencias.validar_rna("ACGT"))
        self.assertFalse(sequencias.validar_rna(""))  # vazia

class TestValidacaoProteina(unittest.TestCase):
    def test_proteina_valida(self):
        self.assertTrue(sequencias.validar_proteina("ACDEFGHIKLMNPQRSTVWY"))
        self.assertTrue(sequencias.validar_proteina("acdefghiklmnpqrstvwy"))

    def test_proteina_invalida(self):
        self.assertFalse(sequencias.validar_proteina("ACDEFX"))
        self.assertFalse(sequencias.validar_proteina(""))  # vazia

class TestTranscricao(unittest.TestCase):
    def test_transcricao_basica(self):
        self.assertEqual(sequencias.transcricao("ATGC"), "AUGC")

    def test_transcricao_invalidas(self):
        self.assertEqual(sequencias.transcricao("ACGTX"), "")
        self.assertEqual(sequencias.transcricao(""), "")

class TestComplemento(unittest.TestCase):
    def test_complemento_basico(self):
        self.assertEqual(sequencias.complemento("ATGC"), "TACG")

    def test_complemento_invalidas(self):
        self.assertEqual(sequencias.complemento("ATGX"), "")
        self.assertEqual(sequencias.complemento(""), "")

class TestReverso(unittest.TestCase):
    def test_reverso_basico(self):
        self.assertEqual(sequencias.reverso("ATGC"), "CGTA")

    def test_reverso_vazio(self):
        self.assertEqual(sequencias.reverso(""), "")

class TestComplementoInverso(unittest.TestCase):
    def test_complemento_inverso_basico(self):
        self.assertEqual(sequencias.complemento_inverso("ATGC"), "GCAT")

    def test_complemento_inverso_invalidas(self):
        self.assertEqual(sequencias.complemento_inverso("ATGX"), "")
        self.assertEqual(sequencias.complemento_inverso(""), "")

class TestSequenciasTamanho1(unittest.TestCase):
    def test_tamanho1_dna(self):
        self.assertTrue(sequencias.validar_dna("A"))
        self.assertEqual(sequencias.transcricao("T"), "U")
        self.assertEqual(sequencias.complemento("A"), "T")
        self.assertEqual(sequencias.reverso("G"), "G")
        self.assertEqual(sequencias.complemento_inverso("C"), "G")

class TestGetORFs(unittest.TestCase): 

    def test_get_orfs(self):
        casos = [
            ("ATGAAATAG", ["ATGAAATAG"]),             
            ("AAACCCGGG", []),                       
            ("ATGAAATAGATGCCCTAA", ["ATGAAATAG", "ATGCCCTAA"])  
        ]

        for dna, esperado in casos:
            self.assertEqual(get_orfs(dna), esperado)

if __name__ == "__main__":
    unittest.main()

