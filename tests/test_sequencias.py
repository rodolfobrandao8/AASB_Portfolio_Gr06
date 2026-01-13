import unittest
from bioinf.sequencias import validar_dna, transcricao, complemento_inverso

class TestSequencias(unittest.TestCase):
    
    def test_dna_valido(self):
        # Testar sequências válidas
        self.assertTrue(validar_dna("ACGT"))
        self.assertTrue(validar_dna("acgt"))

    def test_dna_invalido(self):
        # Testar sequências inválidas
        self.assertFalse(validar_dna("ACGXB"))
        self.assertFalse(validar_dna("123"))
        self.assertFalse(validar_dna("")) # Vazio deve ser falso

    def test_transcricao(self):
        # Testar transcrição (T -> U)
        self.assertEqual(transcricao("ATGC"), "AUGC") 
        self.assertEqual(transcricao("TTTT"), "UUUU")
        self.assertEqual(transcricao(""), "")

    def test_complemento_inverso(self):
        # Testar complemento inverso (Inverte + Complementa)
        # ATCG -> TAGC (Comp) -> CGAT (Inv)
        self.assertEqual(complemento_inverso("ATCG"), "CGAT")
        self.assertEqual(complemento_inverso("aaa"), "TTT")
        self.assertEqual(complemento_inverso("GAATTC"), "GAATTC") 

if __name__ == '__main__':
    unittest.main()